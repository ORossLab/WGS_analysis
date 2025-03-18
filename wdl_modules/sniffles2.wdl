version 1.0

task sniffles2 {
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

        sniffles --version 
        sniffles \
            -i ~{aligned_bam} \
            --reference ~{reference} \
            --no-qc \
            --threads ~{threads} \
            --output-rnames \
            -v "temp.vcf"

        bcftools --version
        echo "~{sample_id}.~{reference_id}.~{aligner}.sniffles2" > new_name.txt
        bcftools reheader -s new_name.txt -o "~{sample_id}.~{reference_id}.~{aligner}.sniffles2.vcf"  "temp.vcf"

        bgzip --version
        bgzip "~{sample_id}.~{reference_id}.~{aligner}.sniffles2.vcf"
        tabix --version
        tabix -p vcf "~{sample_id}.~{reference_id}.~{aligner}.sniffles2.vcf.gz"
    >>>

    output {
        File sniffles2_vcf = "~{sample_id}.~{reference_id}.~{aligner}.sniffles2.vcf.gz"
        File sniffles2_vcf_index =  "~{sample_id}.~{reference_id}.~{aligner}.sniffles2.vcf.gz.tbi"
    }
    runtime {
        cpu: threads
        memory: mem_gb + " GB"
    }
}
