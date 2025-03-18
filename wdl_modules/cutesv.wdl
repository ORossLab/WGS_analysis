version 1.0

task cutesv {
    input {
        String sample_id
        File aligned_bam
        File aligned_bam_index

        File reference
        String reference_id

        String aligner
    }

    Int threads = 8
    Int mem_gb = 16

    command <<<
        set -euo pipefail

        cuteSV --version
        cuteSV \
            --max_size -1 --threads ~{threads} \
            --max_cluster_bias_INS 1000 \
            --diff_ratio_merging_INS 0.9 \
            --max_cluster_bias_DEL 1000 \
            --diff_ratio_merging_DEL 0.5 \
            --genotype \
            ~{aligned_bam} ~{reference} "temp.vcf" "./"

        bcftools --version
        echo "~{sample_id}.~{reference_id}.~{aligner}.cutesv" > new_name.txt
        bcftools reheader -s new_name.txt -o "~{sample_id}.~{reference_id}.~{aligner}.cutesv.vcf"  "temp.vcf"

        bgzip --version
        bgzip "~{sample_id}.~{reference_id}.~{aligner}.cutesv.vcf"
        tabix --version
        tabix -p vcf "~{sample_id}.~{reference_id}.~{aligner}.cutesv.vcf.gz"
    >>>

    output {
        File cutesv_vcf = "~{sample_id}.~{reference_id}.~{aligner}.cutesv.vcf.gz"
        File cutesv_vcf_index =  "~{sample_id}.~{reference_id}.~{aligner}.cutesv.vcf.gz.tbi"
    }

    runtime {
        cpu: threads
        memory: mem_gb + " GB"
    }
}
