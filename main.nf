// takes VCF
// filename must be: TUMORNAME_vs_NORMALNAME.vcf OR TUMORNAME_vs_NORMALNAME.vcf.gz
// extracts tumor and normal sample names from filename

nextflow.enable.dsl=2

// params.outdir     = "$.params.dataset.s3|/data/"
// params.outdir = "./results"
params.vcf        = null
params.ref_fasta   = null
params.chain_file = null   // optional: auto-downloaded from GATK/Broad if not provided

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

    echo "===== VCF #CHROM LINE =====" >> "\$DIAG"
    grep "^#CHROM" "${vcf}" >> "\$DIAG" || true

    echo "" >> "\$DIAG"
    echo "===== EXISTING SAMPLE META LINES =====" >> "\$DIAG"
    grep "^##tumor_sample\\|^##normal_sample" "${vcf}" >> "\$DIAG" || echo "(none found)" >> "\$DIAG"

    echo "" >> "\$DIAG"
    echo "===== RUNNING FIX =====" >> "\$DIAG"

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

print(f"##tumor_sample present : {has_tumor_meta}")
print(f"##normal_sample present: {has_normal_meta}")
print(f"tumor col present      : {has_tumor_col}  (looking for '{tumor}')")
print(f"normal col present     : {has_normal_col}  (looking for '{normal}')")

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

print("Done.")
PYEOF

    echo "" >> "\$DIAG"
    echo "===== OUTPUT #CHROM LINE =====" >> "\$DIAG"
    grep "^#CHROM" "${sample_id}.sample_fixed.vcf" >> "\$DIAG" || true

    echo "" >> "\$DIAG"
    echo "===== OUTPUT SAMPLE META LINES =====" >> "\$DIAG"
    grep "^##tumor_sample\\|^##normal_sample" "${sample_id}.sample_fixed.vcf" >> "\$DIAG" || true
    """
}

// ── process 2: add chr prefixes to ##contig headers and data rows ─────────────
process ADD_CHR_PREFIX {
    tag "${sample_id}"

    container "biocontainers/samtools:v1.9-4-deb_cv1"

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
    echo "Sample ID  : ${sample_id}"   >> "\$DIAG"
    echo "Tumor name : ${tumor_name}"  >> "\$DIAG"
    echo "Normal name: ${normal_name}" >> "\$DIAG"
    echo "Input VCF  : ${vcf}"         >> "\$DIAG"
    echo "" >> "\$DIAG"

    echo "===== CONTIGS BEFORE =====" >> "\$DIAG"
    grep -v "^#" "${vcf}" | cut -f1 | sort -u >> "\$DIAG" || true

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
    ' OFS='\t' ${vcf} > ${sample_id}.chr_prefixed.vcf

    echo "" >> "\$DIAG"
    echo "===== CONTIGS AFTER =====" >> "\$DIAG"
    grep -v "^#" "${sample_id}.chr_prefixed.vcf" | cut -f1 | sort -u >> "\$DIAG" || true

    echo "" >> "\$DIAG"
    echo "===== CONTIG HEADER LINES AFTER =====" >> "\$DIAG"
    grep "^##contig" "${sample_id}.chr_prefixed.vcf" | head -10 >> "\$DIAG" || true
    """
}

// ── process 3: download chain file if not provided ────────────────────────────
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

// ── process 4: GATK LiftoverVcf ───────────────────────────────────────────────
process LIFTOVER {
    tag "${sample_id}"

    container "broadinstitute/gatk:4.6.2.0"


    input:
    tuple val(sample_id), val(tumor_name), val(normal_name), path(vcf)
    path chain
    path ref_fasta
    path ref_dict
    path ref_fai

    output:
    tuple val(sample_id), val(tumor_name), val(normal_name), path("${sample_id}.hg38.vcf"), emit: lifted_vcf
    path "${sample_id}.hg38.rejected.vcf",        emit: rejected_vcf
    path "${sample_id}.liftover.diagnostics.txt", emit: diag

    script:
    """
    set -euo pipefail

    DIAG="${sample_id}.liftover.diagnostics.txt"

    echo "===== LIFTOVER =====" > "\$DIAG"
    echo "Sample ID  : ${sample_id}"   >> "\$DIAG"
    echo "Tumor name : ${tumor_name}"  >> "\$DIAG"
    echo "Normal name: ${normal_name}" >> "\$DIAG"
    echo "Input VCF  : ${vcf}"         >> "\$DIAG"
    echo "Chain file : ${chain}"       >> "\$DIAG"
    echo "Ref FASTA  : ${ref_fasta}"   >> "\$DIAG"
    echo "" >> "\$DIAG"

    echo "===== VARIANT COUNT IN =====" >> "\$DIAG"
    grep -vc "^#" ${vcf} >> "\$DIAG" || true

    echo "" >> "\$DIAG"
    echo "===== RUNNING LIFTOVER =====" >> "\$DIAG"

    gatk LiftoverVcf \
        -I  ${vcf} \
        -O  ${sample_id}.hg38.vcf \
        -R  ${ref_fasta} \
        --CHAIN  ${chain} \
        --REJECT ${sample_id}.hg38.rejected.vcf \
        --RECOVER_SWAPPED_REF_ALT true \
        --WARN_ON_MISSING_CONTIG true \
        --MAX_RECORDS_IN_RAM 500000 \
        >> "\$DIAG" 2>&1

    echo "" >> "\$DIAG"
    echo "===== VARIANT COUNT OUT =====" >> "\$DIAG"
    grep -vc "^#" ${sample_id}.hg38.vcf >> "\$DIAG" || true

    echo "" >> "\$DIAG"
    echo "===== REJECTED COUNT =====" >> "\$DIAG"
    grep -vc "^#" ${sample_id}.hg38.rejected.vcf >> "\$DIAG" || true

    echo "" >> "\$DIAG"
    echo "===== OUTPUT CONTIG SAMPLE =====" >> "\$DIAG"
    grep -v "^#" ${sample_id}.hg38.vcf | head -5 | cut -f1 >> "\$DIAG" || true
    """
}

// ── workflow ──────────────────────────────────────────────────────────────────
workflow {
    println "params.vcf        = ${params.vcf}"
    println "params.ref_fasta   = ${params.ref_fasta}"
    println "params.chain_file = ${params.chain_file}"
    println params.dump()

    if (!params.vcf)      error "Missing required parameter: --vcf"
    if (!params.ref_fasta) error "Missing required parameter: --ref_fasta"

    vcf_ch = Channel.fromPath(params.vcf, checkIfExists: true)
        .map { vcf_file ->

            def filename = vcf_file.name

            def m = filename =~ /^(.+)_vs_(.+)\.vcf(?:\.gz)?$/

            if (!m.matches()) {
                error """
Invalid input VCF filename.

Expected naming convention:
    TUMORNAME_vs_NORMALNAME.vcf
or:
    TUMORNAME_vs_NORMALNAME.vcf.gz

Example:
    GBM1.DFCI4.S1.C4_vs_GBM1.DFCI4.PBMC.vcf

Got:
    ${filename}

The pipeline extracts tumor_name and normal_name from the filename, so --tumor_name and --normal_name should not be supplied.
"""
            }

            def tumor_name  = m[0][1]
            def normal_name = m[0][2]

            if (!tumor_name?.trim()) {
                error "Invalid VCF filename '${filename}': tumor name before '_vs_' is empty."
            }

            if (!normal_name?.trim()) {
                error "Invalid VCF filename '${filename}': normal name after '_vs_' is empty."
            }

            if (tumor_name.contains("_vs_") || normal_name.contains("_vs_")) {
                error """
Invalid VCF filename '${filename}'.

The delimiter '_vs_' should appear exactly once.
"""
            }

            def sample_id = "${tumor_name}_vs_${normal_name}"

            println "Parsed tumor_name  = ${tumor_name}"
            println "Parsed normal_name = ${normal_name}"
            println "Parsed sample_id   = ${sample_id}"

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
    ADD_CHR_PREFIX(FIX_SAMPLE_HEADERS.out[0])

    LIFTOVER(
        ADD_CHR_PREFIX.out[0],
        chain_ch,
        ref_fasta_ch,
        ref_dict_ch,
        ref_fai_ch
    )
}