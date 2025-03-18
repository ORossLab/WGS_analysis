version 1.0

task pbsv {
    input {
        String sample_id
        File aligned_bam
        File aligned_bam_index

        File reference
        String reference_id

        String aligner
    }

    Int threads = 32
    Int mem_gb = 98

    command <<<
        set -euo pipefail

        pbsv --version

        pbsv discover \
            --log-level INFO \
            --hifi \
            ~{aligned_bam} \
            ~{sample_id}.svsig.gz

        pbsv call \
            --hifi \
            -j ~{threads} \
            --log-level INFO \
            --min-sv-length 20 \
            --max-dup-length 20M \
            --max-ins-length 1M \
            ~{reference} \
            ~{sample_id}.svsig.gz \
            "temp.vcf"

        bcftools --version
        echo "~{sample_id}.~{reference_id}.~{aligner}.pbsv" > new_name.txt
        bcftools reheader -s new_name.txt -o "~{sample_id}.~{reference_id}.~{aligner}.pbsv.vcf"  "temp.vcf"


        bgzip --version
        bgzip "~{sample_id}.~{reference_id}.~{aligner}.pbsv.vcf"
        tabix --version
        tabix -p vcf "~{sample_id}.~{reference_id}.~{aligner}.pbsv.vcf.gz"
    >>>

    output {
        File pbsv_vcf = "~{sample_id}.~{reference_id}.~{aligner}.pbsv.vcf.gz"
        File pbsv_vcf_index =  "~{sample_id}.~{reference_id}.~{aligner}.pbsv.vcf.gz.tbi"
    }

    runtime {
        cpu: threads
        memory: mem_gb + " GB"
    }
}
