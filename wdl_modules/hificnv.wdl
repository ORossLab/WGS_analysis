version 1.0

task hificnv {
    input {
        String sample_id
        String? sex

        File bam
        File bam_index

        File reference
        File reference_index
        String reference_id

        File exclude_bed
        File exclude_bed_index 
        File expected_bed_male
        File expected_bed_female
    }

    Boolean sex_defined = defined(sex)
    File expected_bed = if select_first([sex, "FEMALE"]) == "MALE" then expected_bed_male else expected_bed_female

    Int threads = 20
    Int mem_gb = 32

    command <<<
        set -euo pipefail

        echo ~{if sex_defined then "" else "Sex is not defined for ~{sample_id}. Defaulting to karyotype XX for HiFiCNV."}

        hificnv --version

        hificnv \
            --threads ~{threads} \
            --bam ~{bam} \
            --ref ~{reference} \
            --exclude ~{exclude_bed} \
            --expected-cn ~{expected_bed} \
            --output-prefix ~{sample_id}.~{reference_id}.hificnv

        bcftools index --tbi ~{sample_id}.~{reference_id}.hificnv.~{sample_id}.vcf.gz

    >>>

    output {
        File cnv_vcf = "~{sample_id}.~{reference_id}.hificnv.~{sample_id}.vcf.gz"
        File cnv_vcf_index = "~{sample_id}.~{reference_id}.hificnv.~{sample_id}.vcf.gz.tbi"
        File copynum_bedgraph = "~{sample_id}.~{reference_id}.hificnv.~{sample_id}.copynum.bedgraph"
        File depth_bw = "~{sample_id}.~{reference_id}.hificnv.~{sample_id}.depth.bw"
    }

    runtime {
        cpu: threads
        memory: mem_gb + " GB"
    }
}
