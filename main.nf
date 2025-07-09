#!/usr/local/bin nextflow

process prepareInputData {
    publishDir "${params.outdir}/intermediate/extracted_input"
    label 'low'

    input:
    path input_data_path

    output:
    path "**/*.bam", emit: extracted_bams

    script:
    """
    # improves error handling by catching errors early
    set -euo pipefail

    tar -xf ${input_data_path}
    """
}

process convertToHG38BED {
    publishDir "${params.outdir}/intermediate/converted_beds", pattern: "*_hg38.bed"
    publishDir "${params.outdir}/intermediate/unmapped", pattern: "*unMapped"
    tag "${sample_id}"
    label 'low'

    input:
    tuple val(sample_id), path(bam_file)

    output:
    tuple val(sample_id), path("*_hg38.bed"), emit: converted_hg38_bed
    path "*unMapped", emit: unmapped_bed, optional: true

    script:
    """
    set -euo pipefail

    bedtools bamtobed -i ${bam_file} > ${sample_id}.bed

    # selected UCSC LiftOver chain file because the files were initially aligned to a UCSC reference
    # genome that is prepared for next-gen sequencing alignment pipelines
    liftOver ${sample_id}.bed \$LIFTOVER_CHAIN_PATH ${sample_id}_hg38.bed ${sample_id}_unMapped

    if [[ ! -s "${sample_id}_hg38.bed" ]]; then
        echo "ERROR: Output file is missing or empty" >&2
        exit 1
    fi
    """
}

process convertToWIG {
    publishDir "${params.outdir}/intermediate/tumor_wigs", pattern: "*.wig"
    tag "${sample_id}"
    label 'low'

    input:
    tuple val(sample_id), path(bed_file)

    output:
    tuple val(sample_id), path("*_tumor.wig"), emit: tumor_wig

    script:
    """
    set -euo pipefail

    # filter out any rows with alternative contigs from the hg38.chrom.sizes because ichorCNA expects canonical chromosomes
    awk '\$1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)\$/' \$CHROM_SIZES_PATH > hg38_filtered.chrom.sizes
    bedtools makewindows -g hg38_filtered.chrom.sizes -w 1000000 > genome_1Mb_windows.bed
    bedtools intersect -a genome_1Mb_windows.bed -b ${bed_file} -c > ${sample_id}_sample_counts.bed
    
    awk '
        BEGIN { step=1000000; span=1000000 }
        /^chr/ {
            if (\$1 != prev_chr) {
                print "fixedStep chrom=" \$1 " start=1 step=" step " span=" span
                prev_chr = \$1
            }
            print \$4
        }' ${sample_id}_sample_counts.bed > ${sample_id}_tumor.wig

    if [[ ! -s "${sample_id}_tumor.wig" ]]; then
        echo "ERROR: Output file is missing or empty" >&2
        exit 1
    fi
    """
}

process predictTumorFraction {
    publishDir "${params.outdir}/output"
    tag "${sample_id}"
    label 'high'

    input:
    tuple val(sample_id), path(tumor_wig)

    output:
    path('*.params.txt'), emit: ichorCNA_params
    path('*')

    script:
    """
    set -euo pipefail

    ichorCNA_extdata="/ichorCNA/inst/extdata"
    Rscript /ichorCNA/scripts/runIchorCNA.R \\
        --id ${sample_id} \\
        --WIG ${tumor_wig} \\
        --gcWig \${ichorCNA_extdata}/gc_hg38_1000kb.wig \\
        --mapWig \${ichorCNA_extdata}/map_hg38_1000kb.wig \\
        --centromere \${ichorCNA_extdata}/GRCh38.GCA_000001405.2_centromere_acen.txt \\
        --normalPanel \${ichorCNA_extdata}/HD_ULP_PoN_hg38_1Mb_median_normAutosome_median.rds \\
        --includeHOMD FALSE \\
        --chrs "c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX')" \
        --chrTrain "c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22')" \
        --genomeBuild hg38 \\
        --genomeStyle UCSC \\
        --outDir .
    
    if [[ ! -s "${sample_id}.params.txt" ]]; then
        echo "ERROR: Output file is missing or empty" >&2
        exit 1
    fi
    """
}

process generateOutputTable {
    publishDir "${params.outdir}/output"
    label 'high'

    input:
    path sample_params

    output:
    path "*.tsv", emit: output_table

    script:
    """
    set -euo pipefail

    echo -e "sample identifier\ttumor fraction" > sample_tumor_fractions.tsv
    for file in ${sample_params}; do
        # Extract sample identifier and tumor fraction from the second row
        awk 'NR==2 { print \$1 "\\t" \$2 }' \$file >> sample_tumor_fractions.tsv
    done
    """
}

workflow {
    // Prepare input data and extract BAM files
    extracted = prepareInputData(params.input_data_path)

    // Collect extracted BAM files and create named tuples
    named_bams = extracted.extracted_bams.flatten().map { bam -> tuple(bam.baseName, bam) } 

    // Convert extracted HG19 BAM files to HG38 BED format
    converted_hg38_beds = convertToHG38BED(named_bams)

    // Convert BED files to WIG format
    tumor_wigs = convertToWIG(converted_hg38_beds.converted_hg38_bed)

    // Predict tumor fraction using the WIG files
    ichorCNA_output = predictTumorFraction(tumor_wigs.tumor_wig)

    // Generate output table from all of the .params.txt files
    generateOutputTable(ichorCNA_output.ichorCNA_params.collect())
}