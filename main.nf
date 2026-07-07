// takes VCF
// filename must be: TUMORNAME_vs_NORMALNAME.vcf OR TUMORNAME_vs_NORMALNAME.vcf.gz
// extracts tumor and normal sample names from filename

nextflow.enable.dsl=2

params.vcf         = null
params.ref_fasta   = null
params.chain_file  = null
params.tumor_name  = null
params.normal_name = null

// Consensus filtering options
params.anchor_caller = "M1"
params.min_non_anchor_callers = 2
params.consensus_type = "SM"

// ── process 1: inject ##tumor_sample / ##normal_sample + add missing sample cols ──
process FIX_SAMPLE_HEADERS {
    tag "${sample_id}"
    container "broadinstitute/gatk:4.6.2.0"

    input:
    tuple val(sample_id), val(tumor_name), val(normal_name), path(vcf)

    output:
    tuple val(sample_id), val(tumor_name), val(normal_name), path("${sample_id}.sample_fixed.vcf")
    path "${sample_id}.fix_headers.diagnostics.txt"

    script:
    """
    set -euo pipefail

    DIAG="${sample_id}.fix_headers.diagnostics.txt"

    echo "===== FIX SAMPLE HEADERS =====" > "\$DIAG"
    echo "Sample ID  : ${sample_id}"   >> "\$DIAG"
    echo "Input VCF  : ${vcf}"         >> "\$DIAG"
    echo "Tumor name : ${tumor_name}"  >> "\$DIAG"
    echo "Normal name: ${normal_name}" >> "\$DIAG"
    echo "" >> "\$DIAG"

    python3 - << 'PYEOF' >> "\$DIAG" 2>&1
tumor  = "${tumor_name}"
normal = "${normal_name}"
infile = "${vcf}"
outfile = "${sample_id}.sample_fixed.vcf"

with open(infile) as fh:
    lines = fh.readlines()

chrom_lines = [l for l in lines if l.startswith("#CHROM")]
if len(chrom_lines) != 1:
    raise ValueError(f"Expected exactly one #CHROM line, found {len(chrom_lines)}")

has_tumor_meta  = any(l.startswith("##tumor_sample=")  for l in lines)
has_normal_meta = any(l.startswith("##normal_sample=") for l in lines)

chrom_line = chrom_lines[0]
cols = chrom_line.rstrip("\\n").split("\\t")

has_tumor_col  = tumor in cols
has_normal_col = normal in cols

with open(outfile, "w") as out_fh:
    for line in lines:
        if line.startswith("#CHROM"):
            if not has_tumor_meta:
                out_fh.write(f"##tumor_sample={tumor}\\n")
            if not has_normal_meta:
                out_fh.write(f"##normal_sample={normal}\\n")

            new_cols = cols[:]
            if not has_tumor_col:
                new_cols.append(tumor)
            if not has_normal_col:
                new_cols.append(normal)

            out_fh.write("\\t".join(new_cols) + "\\n")
            continue

        if not line.startswith("#"):
            fields = line.rstrip("\\n").split("\\t")
            if not has_tumor_col:
                fields.append("./.")
            if not has_normal_col:
                fields.append("./.")
            out_fh.write("\\t".join(fields) + "\\n")
            continue

        out_fh.write(line)
PYEOF
    """
}

// ── process 2: consensus filter + caller-specific tumor AF annotation ─────────
process CONSENSUS_FILTER_AND_CALLER_AF {
    tag "${sample_id}"

    container "broadinstitute/gatk:4.5.0.0"
    publishDir "${params.outdir}", mode: "copy"

    input:
    tuple val(sample_id), val(tumor_name), val(normal_name), path(vcf)

    output:
    tuple val(sample_id), val(tumor_name), val(normal_name), path("${sample_id}.consensus.vcf"), emit: consensus_vcf
    path "${sample_id}.consensus_af.tsv",    emit: af_tsv
    path "${sample_id}.consensus_stats.tsv", emit: stats_tsv

    script:
    """
    set -euo pipefail

    python3 - << 'PYEOF'
import re
from collections import Counter

infile = "${vcf}"
out_vcf = "${sample_id}.consensus.vcf"
out_af_tsv = "${sample_id}.consensus_af.tsv"
out_stats_tsv = "${sample_id}.consensus_stats.tsv"
sample_id = "${sample_id}"

anchor_caller = "${params.anchor_caller}".upper()
min_non_anchor_callers = int("${params.min_non_anchor_callers}")
consensus_type = "${params.consensus_type}"

anchor_aliases = {anchor_caller, "M1", "MUTECT1", "MUTECT"}
known_callers = ["HC", "M1", "M2", "VS"]

def parse_info(info_str):
    d = {}
    order = []
    if info_str == ".":
        return d, order
    for item in info_str.split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
        else:
            k, v = item, True
        d[k] = v
        order.append(k)
    return d, order

def rebuild_info(d, order):
    seen = set()
    items = []
    for k in order:
        if k in d and k not in seen:
            v = d[k]
            items.append(k if v is True else f"{k}={v}")
            seen.add(k)
    for k, v in d.items():
        if k not in seen:
            items.append(k if v is True else f"{k}={v}")
    return ";".join(items) if items else "."

def sanitize_info_id(label):
    return re.sub(r"[^A-Za-z0-9_.-]", "_", label)

def get_format_for_caller(format_str, idx, n_callers):
    if "|" in format_str:
        parts = format_str.split("|")
        if len(parts) == n_callers:
            return parts[idx]
    return format_str

def calc_af(fmt, sample_value):
    keys = fmt.split(":")
    vals = sample_value.split(":")
    data = dict(zip(keys, vals))

    ad = data.get("AD")
    dp = data.get("DP")

    if not ad or not dp or ad == "." or dp == ".":
        return None, None, None

    try:
        depths = [int(x) for x in ad.split(",") if x != "."]
        dp_int = int(dp)
        if len(depths) < 2 or dp_int <= 0:
            return None, ad, dp

        alt_depth = sum(depths[1:])
        return alt_depth / dp_int, ad, dp
    except Exception:
        return None, ad, dp

def empty_caller_values():
    vals = {}
    for c in known_callers:
        vals[f"AF_{c}"] = "."
        vals[f"AD_{c}"] = "."
        vals[f"DP_{c}"] = "."
    return vals

n_in = 0
n_type_match = 0
n_anchor_present = 0
n_pass_consensus = 0
n_kept = 0
n_filtered_type = 0
n_filtered_no_anchor = 0
n_filtered_too_few_non_anchor = 0
n_missing_required = 0
n_mismatched_caller_tumor_count = 0
n_missing_af_records = 0
caller_combo_counts = Counter()

base_tsv_cols = ["sample", "CHROM", "POS", "ID", "REF", "ALT", "FILTER", "TYPE", "NB_CALLERS", "CALLERS"]
caller_tsv_cols = []
for c in known_callers:
    caller_tsv_cols.extend([f"AF_{c}", f"AD_{c}", f"DP_{c}"])
tsv_cols = base_tsv_cols + caller_tsv_cols + ["CALLER_AF", "CALLER_AD", "CALLER_DP"]

with open(infile) as fh, open(out_vcf, "w") as vcf_out, open(out_af_tsv, "w") as af_out:
    af_out.write("\\t".join(tsv_cols) + "\\n")

    for line in fh:
        if line.startswith("##"):
            vcf_out.write(line)
            continue

        if line.startswith("#CHROM"):
            vcf_out.write('##INFO=<ID=CALLER_AF,Number=.,Type=String,Description="Tumor allele fraction per caller calculated from NeoDisc INFO/TUMOR as ALT AD / DP; values are caller:AF in CALLERS order">\\n')
            vcf_out.write('##INFO=<ID=CALLER_AD,Number=.,Type=String,Description="Tumor AD per caller from NeoDisc INFO/TUMOR; values are caller:AD in CALLERS order">\\n')
            vcf_out.write('##INFO=<ID=CALLER_DP,Number=.,Type=String,Description="Tumor DP per caller from NeoDisc INFO/TUMOR; values are caller:DP in CALLERS order">\\n')
            for c in known_callers:
                vcf_out.write(f'##INFO=<ID=AF_{c},Number=1,Type=Float,Description="Tumor allele fraction for caller {c} calculated from NeoDisc INFO/TUMOR as ALT AD / DP">\\n')
            vcf_out.write(line)
            continue

        n_in += 1
        fields = line.rstrip("\\n").split("\\t")
        if len(fields) < 8:
            n_missing_required += 1
            continue

        chrom, pos, var_id, ref, alt, qual, filt, info_str = fields[:8]
        info, info_order = parse_info(info_str)

        callers_raw = info.get("CALLERS")
        tumor_raw = info.get("TUMOR")
        format_raw = info.get("FORMAT")

        if not callers_raw or not tumor_raw or not format_raw:
            n_missing_required += 1
            continue

        if consensus_type and consensus_type.lower() != "null":
            if info.get("TYPE") != consensus_type:
                n_filtered_type += 1
                continue
        n_type_match += 1

        callers = callers_raw.split("|")
        tumors = tumor_raw.split("|")
        callers_upper = [c.upper() for c in callers]

        has_anchor = any(c in anchor_aliases for c in callers_upper)
        non_anchor_count = sum(c not in anchor_aliases for c in callers_upper)

        if not has_anchor:
            n_filtered_no_anchor += 1
            continue
        n_anchor_present += 1

        if non_anchor_count < min_non_anchor_callers:
            n_filtered_too_few_non_anchor += 1
            continue

        n_pass_consensus += 1

        if len(tumors) != len(callers):
            n_mismatched_caller_tumor_count += 1
            continue

        caller_values = empty_caller_values()
        caller_af_items = []
        caller_ad_items = []
        caller_dp_items = []

        for i, caller in enumerate(callers):
            fmt_i = get_format_for_caller(format_raw, i, len(callers))
            af, ad, dp = calc_af(fmt_i, tumors[i])
            caller_clean = sanitize_info_id(caller)
            caller_upper = caller.upper()

            if af is None:
                af_str = "."
                n_missing_af_records += 1
            else:
                af_str = f"{af:.6f}"
                info[f"AF_{caller_clean}"] = af_str

            caller_af_items.append(f"{caller}:{af_str}")
            caller_ad_items.append(f"{caller}:{ad if ad is not None else '.'}")
            caller_dp_items.append(f"{caller}:{dp if dp is not None else '.'}")

            if caller_upper in known_callers:
                caller_values[f"AF_{caller_upper}"] = af_str
                caller_values[f"AD_{caller_upper}"] = ad if ad is not None else "."
                caller_values[f"DP_{caller_upper}"] = dp if dp is not None else "."

        caller_af = "|".join(caller_af_items)
        caller_ad = "|".join(caller_ad_items)
        caller_dp = "|".join(caller_dp_items)

        info["CALLER_AF"] = caller_af
        info["CALLER_AD"] = caller_ad
        info["CALLER_DP"] = caller_dp

        fields[7] = rebuild_info(info, info_order)
        vcf_out.write("\\t".join(fields) + "\\n")

        row = {
            "sample": sample_id,
            "CHROM": chrom,
            "POS": pos,
            "ID": var_id,
            "REF": ref,
            "ALT": alt,
            "FILTER": filt,
            "TYPE": info.get("TYPE", "."),
            "NB_CALLERS": info.get("NB_CALLERS", "."),
            "CALLERS": callers_raw,
            "CALLER_AF": caller_af,
            "CALLER_AD": caller_ad,
            "CALLER_DP": caller_dp,
        }
        row.update(caller_values)
        af_out.write("\\t".join(str(row.get(col, ".")) for col in tsv_cols) + "\\n")

        caller_combo_counts[callers_raw] += 1
        n_kept += 1

n_filtered_total = n_in - n_kept

with open(out_stats_tsv, "w") as stats:
    stats.write("metric\\tvalue\\n")
    stats.write(f"input_vcf\\t{infile}\\n")
    stats.write(f"output_vcf\\t{out_vcf}\\n")
    stats.write(f"anchor_caller\\t{anchor_caller}\\n")
    stats.write(f"anchor_aliases\\t{'|'.join(sorted(anchor_aliases))}\\n")
    stats.write(f"min_non_anchor_callers\\t{min_non_anchor_callers}\\n")
    stats.write(f"required_type\\t{consensus_type}\\n")
    stats.write(f"variants_read\\t{n_in}\\n")
    stats.write(f"variants_matching_required_type\\t{n_type_match}\\n")
    stats.write(f"variants_with_anchor\\t{n_anchor_present}\\n")
    stats.write(f"variants_passing_consensus_rule\\t{n_pass_consensus}\\n")
    stats.write(f"variants_kept\\t{n_kept}\\n")
    stats.write(f"variants_filtered_total\\t{n_filtered_total}\\n")
    stats.write(f"filtered_wrong_type\\t{n_filtered_type}\\n")
    stats.write(f"filtered_no_anchor\\t{n_filtered_no_anchor}\\n")
    stats.write(f"filtered_too_few_non_anchor_callers\\t{n_filtered_too_few_non_anchor}\\n")
    stats.write(f"skipped_missing_required_fields\\t{n_missing_required}\\n")
    stats.write(f"skipped_mismatched_caller_tumor_count\\t{n_mismatched_caller_tumor_count}\\n")
    stats.write(f"missing_af_values\\t{n_missing_af_records}\\n")
    for combo, count in caller_combo_counts.most_common():
        stats.write(f"kept_CALLERS_{combo}\\t{count}\\n")
PYEOF
    """
}

// ── process 3: add chr prefixes to ##contig headers and data rows ─────────────
process ADD_CHR_PREFIX {
    tag "${sample_id}"
    container "broadinstitute/gatk:4.5.0.0"

    input:
    tuple val(sample_id), val(tumor_name), val(normal_name), path(vcf)

    output:
    tuple val(sample_id), val(tumor_name), val(normal_name), path("${sample_id}.chr_prefixed.vcf")
    path "${sample_id}.chr_prefix.diagnostics.txt"

    script:
    """
    set -euo pipefail

    DIAG="${sample_id}.chr_prefix.diagnostics.txt"

    echo "===== ADD CHR PREFIX =====" > "\$DIAG"

    awk '
    /^##contig=/ {
        if (\$0 !~ /ID=chr/) sub(/ID=/, "ID=chr")
        print; next
    }
    /^#/ { print; next }
    {
        if (\$1 ~ /^chr/)             { print; next }
        if (\$1 == "MT")              { \$1 = "chrM" }
        else if (\$1 ~ /^[0-9XY]+\$/) { \$1 = "chr" \$1 }
        print
    }
    ' OFS='\\t' ${vcf} > ${sample_id}.chr_prefixed.vcf
    """
}

// ── process 4: download chain file if not provided ────────────────────────────
process DOWNLOAD_CHAIN {
    container "curlimages/curl:8.6.0"

    output:
    path "b37ToHg38.over.chain", emit: chain

    script:
    """
    curl -fsSL \
        https://raw.githubusercontent.com/broadinstitute/gatk/master/scripts/funcotator/data_sources/gnomAD/b37ToHg38.over.chain \
        -o b37ToHg38.over.chain
    """
}

// ── process 5: GATK LiftoverVcf ───────────────────────────────────────────────
process LIFTOVER {
    tag "${sample_id}"
    container "broadinstitute/gatk:4.6.2.0"
    publishDir "${params.outdir}", mode: "copy"

    input:
    tuple val(sample_id), val(tumor_name), val(normal_name), path(vcf)
    path chain
    path ref_fasta
    path ref_dict
    path ref_fai

    output:
    tuple val(sample_id), val(tumor_name), val(normal_name), path("${sample_id}.hg38.vcf"), emit: lifted_vcf
    path "${sample_id}.hg38.rejected.vcf", emit: rejected_vcf
    path "${sample_id}.liftover.diagnostics.txt", emit: diag

    script:
    """
    set -euo pipefail

    DIAG="${sample_id}.liftover.diagnostics.txt"

    echo "PWD: \$(pwd)" > "\$DIAG"
    echo "Files staged:" >> "\$DIAG"
    ls -lh >> "\$DIAG"
    echo "" >> "\$DIAG"

    echo "Checking inputs..." >> "\$DIAG"
    ls -lh ${vcf} ${chain} ${ref_fasta} ${ref_dict} ${ref_fai} >> "\$DIAG" 2>&1 || true
    echo "" >> "\$DIAG"

    echo "VCF header contigs:" >> "\$DIAG"
    grep '^##contig' ${vcf} | head -20 >> "\$DIAG" 2>&1 || true
    echo "" >> "\$DIAG"

    echo "Reference contigs:" >> "\$DIAG"
    cut -f1 ${ref_fai} | head -20 >> "\$DIAG" 2>&1 || true
    echo "" >> "\$DIAG"

    echo "Running LiftoverVcf..." >> "\$DIAG"

    gatk --java-options "-Xmx12g" LiftoverVcf \
        -I ${vcf} \
        -O ${sample_id}.hg38.vcf \
        -R ${ref_fasta} \
        --CHAIN ${chain} \
        --REJECT ${sample_id}.hg38.rejected.vcf \
        --RECOVER_SWAPPED_REF_ALT true \
        --WARN_ON_MISSING_CONTIG true \
        --MAX_RECORDS_IN_RAM 100000 \
        >> "\$DIAG" 2>&1
    """
}

// ── workflow ──────────────────────────────────────────────────────────────────
workflow {
    println "params.vcf         = ${params.vcf}"
    println "params.ref_fasta   = ${params.ref_fasta}"
    println "params.chain_file  = ${params.chain_file}"
    println "params.tumor_name  = ${params.tumor_name}"
    println "params.normal_name = ${params.normal_name}"
    println params.dump()

    if (!params.vcf)         error "Missing required parameter: --vcf"
    if (!params.ref_fasta)   error "Missing required parameter: --ref_fasta"
    if (!params.tumor_name)  error "Missing required parameter: --tumor_name"
    if (!params.normal_name) error "Missing required parameter: --normal_name"

    vcf_ch = Channel.fromPath(params.vcf, checkIfExists: true)
        .map { vcf_file ->

            def tumor_name  = params.tumor_name
            def normal_name = params.normal_name
            def sample_id   = "${tumor_name}_vs_${normal_name}"

            println "Input VCF        = ${vcf_file}"
            println "Tumor name       = ${tumor_name}"
            println "Normal name      = ${normal_name}"
            println "Sample ID        = ${sample_id}"

            tuple(sample_id, tumor_name, normal_name, vcf_file)
        }

    if (params.chain_file) {
        chain_ch = Channel.fromPath(params.chain_file, checkIfExists: true)
    } else {
        chain_ch = DOWNLOAD_CHAIN().chain
    }

    ref_fasta_ch = Channel.fromPath(params.ref_fasta, checkIfExists: true)
    ref_dict_ch  = Channel.fromPath("${params.ref_fasta}".replaceAll(/\.(fa|fasta)$/, ".dict"), checkIfExists: true)
    ref_fai_ch   = Channel.fromPath("${params.ref_fasta}.fai", checkIfExists: true)

    FIX_SAMPLE_HEADERS(vcf_ch)
    CONSENSUS_FILTER_AND_CALLER_AF(FIX_SAMPLE_HEADERS.out[0])
    ADD_CHR_PREFIX(CONSENSUS_FILTER_AND_CALLER_AF.out.consensus_vcf)

    LIFTOVER(
        ADD_CHR_PREFIX.out[0],
        chain_ch,
        ref_fasta_ch,
        ref_dict_ch,
        ref_fai_ch
    )
}// // takes VCF
// // filename must be: TUMORNAME_vs_NORMALNAME.vcf OR TUMORNAME_vs_NORMALNAME.vcf.gz
// // extracts tumor and normal sample names from filename

// nextflow.enable.dsl=2

// // params.outdir     = "$.params.dataset.s3|/data/"
// // params.outdir = "./results"
// params.vcf        = null
// params.ref_fasta   = null
// params.chain_file = null   // optional: auto-downloaded from GATK/Broad if not provided

// // ── process 1: inject ##tumor_sample / ##normal_sample + add missing sample cols ──
// process FIX_SAMPLE_HEADERS {
//     tag "${sample_id}"
//     container "broadinstitute/gatk:4.6.2.0"

//     input:
//     tuple val(sample_id), val(tumor_name), val(normal_name), path(vcf)

//     output:
//     tuple val(sample_id), val(tumor_name), val(normal_name), path("${sample_id}.sample_fixed.vcf")
//     path "${sample_id}.fix_headers.diagnostics.txt"

//     script:
//     """
//     set -euo pipefail

//     DIAG="${sample_id}.fix_headers.diagnostics.txt"

//     echo "===== FIX SAMPLE HEADERS =====" > "\$DIAG"
//     echo "Sample ID  : ${sample_id}"   >> "\$DIAG"
//     echo "Input VCF  : ${vcf}"         >> "\$DIAG"
//     echo "Tumor name : ${tumor_name}"  >> "\$DIAG"
//     echo "Normal name: ${normal_name}" >> "\$DIAG"
//     echo "" >> "\$DIAG"

//     echo "===== VCF #CHROM LINE =====" >> "\$DIAG"
//     grep "^#CHROM" "${vcf}" >> "\$DIAG" || true

//     echo "" >> "\$DIAG"
//     echo "===== EXISTING SAMPLE META LINES =====" >> "\$DIAG"
//     grep "^##tumor_sample\\|^##normal_sample" "${vcf}" >> "\$DIAG" || echo "(none found)" >> "\$DIAG"

//     echo "" >> "\$DIAG"
//     echo "===== RUNNING FIX =====" >> "\$DIAG"

//     python3 - << 'PYEOF' >> "\$DIAG" 2>&1
// tumor  = "${tumor_name}"
// normal = "${normal_name}"
// infile = "${vcf}"
// outfile = "${sample_id}.sample_fixed.vcf"

// with open(infile) as fh:
//     lines = fh.readlines()

// chrom_lines = [l for l in lines if l.startswith("#CHROM")]
// if len(chrom_lines) != 1:
//     raise ValueError(f"Expected exactly one #CHROM line, found {len(chrom_lines)}")

// has_tumor_meta  = any(l.startswith("##tumor_sample=")  for l in lines)
// has_normal_meta = any(l.startswith("##normal_sample=") for l in lines)

// chrom_line = chrom_lines[0]
// cols = chrom_line.rstrip("\\n").split("\\t")

// has_tumor_col  = tumor in cols
// has_normal_col = normal in cols

// print(f"##tumor_sample present : {has_tumor_meta}")
// print(f"##normal_sample present: {has_normal_meta}")
// print(f"tumor col present      : {has_tumor_col}  (looking for '{tumor}')")
// print(f"normal col present     : {has_normal_col}  (looking for '{normal}')")

// with open(outfile, "w") as out_fh:
//     for line in lines:
//         if line.startswith("#CHROM"):
//             if not has_tumor_meta:
//                 out_fh.write(f"##tumor_sample={tumor}\\n")
//             if not has_normal_meta:
//                 out_fh.write(f"##normal_sample={normal}\\n")

//             new_cols = cols[:]

//             if not has_tumor_col:
//                 new_cols.append(tumor)
//             if not has_normal_col:
//                 new_cols.append(normal)

//             out_fh.write("\\t".join(new_cols) + "\\n")
//             continue

//         if not line.startswith("#"):
//             fields = line.rstrip("\\n").split("\\t")

//             if not has_tumor_col:
//                 fields.append("./.")
//             if not has_normal_col:
//                 fields.append("./.")

//             out_fh.write("\\t".join(fields) + "\\n")
//             continue

//         out_fh.write(line)

// print("Done.")
// PYEOF

//     echo "" >> "\$DIAG"
//     echo "===== OUTPUT #CHROM LINE =====" >> "\$DIAG"
//     grep "^#CHROM" "${sample_id}.sample_fixed.vcf" >> "\$DIAG" || true

//     echo "" >> "\$DIAG"
//     echo "===== OUTPUT SAMPLE META LINES =====" >> "\$DIAG"
//     grep "^##tumor_sample\\|^##normal_sample" "${sample_id}.sample_fixed.vcf" >> "\$DIAG" || true
//     """
// }

// // ── process 2: add chr prefixes to ##contig headers and data rows ─────────────
// process ADD_CHR_PREFIX {
//     tag "${sample_id}"

//     container "broadinstitute/gatk:4.5.0.0"

//     input:
//     tuple val(sample_id), val(tumor_name), val(normal_name), path(vcf)

//     output:
//     tuple val(sample_id), val(tumor_name), val(normal_name), path("${sample_id}.chr_prefixed.vcf")
//     path "${sample_id}.chr_prefix.diagnostics.txt"

//     script:
//     """
//     set -euo pipefail

//     DIAG="${sample_id}.chr_prefix.diagnostics.txt"

//     echo "===== ADD CHR PREFIX =====" > "\$DIAG"
//     echo "Sample ID  : ${sample_id}"   >> "\$DIAG"
//     echo "Tumor name : ${tumor_name}"  >> "\$DIAG"
//     echo "Normal name: ${normal_name}" >> "\$DIAG"
//     echo "Input VCF  : ${vcf}"         >> "\$DIAG"
//     echo "" >> "\$DIAG"

//     echo "===== CONTIGS BEFORE =====" >> "\$DIAG"
//     grep -v "^#" "${vcf}" | cut -f1 | sort -u >> "\$DIAG" || true

//     awk '
//     /^##contig=/ {
//         if (\$0 !~ /ID=chr/) sub(/ID=/, "ID=chr")
//         print; next
//     }
//     /^#/ { print; next }
//     {
//         if (\$1 ~ /^chr/)             { print; next }
//         if (\$1 == "MT")              { \$1 = "chrM" }
//         else if (\$1 ~ /^[0-9XY]+\$/) { \$1 = "chr" \$1 }
//         print
//     }
//     ' OFS='\t' ${vcf} > ${sample_id}.chr_prefixed.vcf

//     echo "" >> "\$DIAG"
//     echo "===== CONTIGS AFTER =====" >> "\$DIAG"
//     grep -v "^#" "${sample_id}.chr_prefixed.vcf" | cut -f1 | sort -u >> "\$DIAG" || true

//     echo "" >> "\$DIAG"
//     echo "===== CONTIG HEADER LINES AFTER =====" >> "\$DIAG"
//     grep "^##contig" "${sample_id}.chr_prefixed.vcf" | head -10 >> "\$DIAG" || true
//     """
// }

// // ── process 3: download chain file if not provided ────────────────────────────
// process DOWNLOAD_CHAIN {
//     container "curlimages/curl:8.6.0"

//     output:
//     path "b37ToHg38.over.chain", emit: chain

//     script:
//     """
//     curl -fsSL \
//         https://raw.githubusercontent.com/broadinstitute/gatk/master/scripts/funcotator/data_sources/gnomAD/b37ToHg38.over.chain \
//         -o b37ToHg38.over.chain
//     """
// }

// // ── process 4: GATK LiftoverVcf ───────────────────────────────────────────────
// process LIFTOVER {
//     tag "${sample_id}"

//     container "broadinstitute/gatk:4.6.2.0"

//     publishDir "${params.outdir}", mode: "copy"
//     input:
//     tuple val(sample_id), val(tumor_name), val(normal_name), path(vcf)
//     path chain
//     path ref_fasta
//     path ref_dict
//     path ref_fai

//     output:
//     tuple val(sample_id), val(tumor_name), val(normal_name), path("${sample_id}_vs_${normal_name}.hg38.vcf"), emit: lifted_vcf
//     path "${sample_id}.hg38.rejected.vcf",        emit: rejected_vcf
//     path "${sample_id}.liftover.diagnostics.txt", emit: diag

//     script:
//     """
//     set -euo pipefail

//     DIAG="${sample_id}.liftover.diagnostics.txt"

//     echo "===== LIFTOVER =====" > "\$DIAG"
//     echo "Sample ID  : ${sample_id}"   >> "\$DIAG"
//     echo "Tumor name : ${tumor_name}"  >> "\$DIAG"
//     echo "Normal name: ${normal_name}" >> "\$DIAG"
//     echo "Input VCF  : ${vcf}"         >> "\$DIAG"
//     echo "Chain file : ${chain}"       >> "\$DIAG"
//     echo "Ref FASTA  : ${ref_fasta}"   >> "\$DIAG"
//     echo "" >> "\$DIAG"

//     echo "===== VARIANT COUNT IN =====" >> "\$DIAG"
//     grep -vc "^#" ${vcf} >> "\$DIAG" || true

//     echo "" >> "\$DIAG"
//     echo "===== RUNNING LIFTOVER =====" >> "\$DIAG"

//     gatk LiftoverVcf \
//         -I  ${vcf} \
//         -O  ${sample_id}_vs_${normal_name}.hg38.vcf \
//         -R  ${ref_fasta} \
//         --CHAIN  ${chain} \
//         --REJECT ${sample_id}.hg38.rejected.vcf \
//         --RECOVER_SWAPPED_REF_ALT true \
//         --WARN_ON_MISSING_CONTIG true \
//         --MAX_RECORDS_IN_RAM 500000 \
//         >> "\$DIAG" 2>&1

//     echo "" >> "\$DIAG"
//     echo "===== VARIANT COUNT OUT =====" >> "\$DIAG"
//     grep -vc "^#" ${sample_id}_vs_${normal_name}.hg38.vcf >> "\$DIAG" || true

//     echo "" >> "\$DIAG"
//     echo "===== REJECTED COUNT =====" >> "\$DIAG"
//     grep -vc "^#" ${sample_id}.hg38.rejected.vcf >> "\$DIAG" || true

//     echo "" >> "\$DIAG"
//     echo "===== OUTPUT CONTIG SAMPLE =====" >> "\$DIAG"
//     grep -v "^#" ${sample_id}_vs_${normal_name}.hg38.vcf | head -5 | cut -f1 >> "\$DIAG" || true
//     """
// }

// // ── workflow ──────────────────────────────────────────────────────────────────
// workflow {
//     println "params.vcf        = ${params.vcf}"
//     println "params.ref_fasta   = ${params.ref_fasta}"
//     println "params.chain_file = ${params.chain_file}"
//     println params.dump()

//     if (!params.vcf)      error "Missing required parameter: --vcf"
//     if (!params.ref_fasta) error "Missing required parameter: --ref_fasta"

//     vcf_ch = Channel.fromPath(params.vcf, checkIfExists: true)
//         .map { vcf_file ->

//             def filename = vcf_file.name

//             def m = filename =~ /^(.+)_vs_(.+)\.vcf(?:\.gz)?$/

//             if (!m.matches()) {
//                 error """
// Invalid input VCF filename.

// Expected naming convention:
//     TUMORNAME_vs_NORMALNAME.vcf
// or:
//     TUMORNAME_vs_NORMALNAME.vcf.gz

// Example:
//     GBM1.DFCI4.S1.C4_vs_GBM1.DFCI4.PBMC.vcf

// Got:
//     ${filename}

// The pipeline extracts tumor_name and normal_name from the filename, so --tumor_name and --normal_name should not be supplied.
// """
//             }

//             def tumor_name  = m[0][1]
//             def normal_name = m[0][2]

//             if (!tumor_name?.trim()) {
//                 error "Invalid VCF filename '${filename}': tumor name before '_vs_' is empty."
//             }

//             if (!normal_name?.trim()) {
//                 error "Invalid VCF filename '${filename}': normal name after '_vs_' is empty."
//             }

//             if (tumor_name.contains("_vs_") || normal_name.contains("_vs_")) {
//                 error """
// Invalid VCF filename '${filename}'.

// The delimiter '_vs_' should appear exactly once.
// """
//             }

//             def sample_id = "${tumor_name}_vs_${normal_name}"

//             println "Parsed tumor_name  = ${tumor_name}"
//             println "Parsed normal_name = ${normal_name}"
//             println "Parsed sample_id   = ${sample_id}"

//             tuple(sample_id, tumor_name, normal_name, vcf_file)
//         }

//     if (params.chain_file) {
//         chain_ch = Channel.fromPath(params.chain_file, checkIfExists: true)
//     } else {
//         chain_ch = DOWNLOAD_CHAIN().chain
//     }

//     ref_fasta_ch = Channel.fromPath(params.ref_fasta, checkIfExists: true)
//     ref_dict_ch  = Channel.fromPath("${params.ref_fasta}".replaceAll(/\.(fa|fasta)$/, ".dict"), checkIfExists: true)
//     ref_fai_ch   = Channel.fromPath("${params.ref_fasta}.fai", checkIfExists: true)

//     FIX_SAMPLE_HEADERS(vcf_ch)
//     ADD_CHR_PREFIX(FIX_SAMPLE_HEADERS.out[0])

//     LIFTOVER(
//         ADD_CHR_PREFIX.out[0],
//         chain_ch,
//         ref_fasta_ch,
//         ref_dict_ch,
//         ref_fai_ch
//     )
// }