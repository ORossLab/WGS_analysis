version 1.0

task sniffles2_vacmap_filter {
    input {
        String sample_id
        File aligned_bam
        File aligned_bam_index

        File reference
        String reference_id

        String script_path
    }

    Int threads = 8
    Int mem_gb = 16

    command <<<
        set -euo pipefail

        sniffles --version
        sniffles \
            -i ~{aligned_bam} \
            --reference ~{reference} \
            --threads ~{threads} \
            --output-rnames \
            -v "temp.vcf"

        bcftools --version
        echo "~{sample_id}.~{reference_id}.vacmap.sniffles2" > new_name.txt
        bcftools reheader -s new_name.txt -o "~{sample_id}.~{reference_id}.vacmap.sniffles2.vcf"  "temp.vcf"



        bcftools view "~{sample_id}.~{reference_id}.vacmap.sniffles2.vcf" | grep -E "SVTYPE=DUP|SVTYPE=INV" > "~{sample_id}.~{reference_id}.vacmap.sniffles2.INV.DUP.vcf"


        bash ~{script_path}/vacmap_filter.sh --sample_id "~{sample_id}" --input_vcf "~{sample_id}.~{reference_id}.vacmap.sniffles2.INV.DUP.vcf" --ref "~{reference}" --ref_name "~{reference_id}"
        python ~{script_path}/vacmap_analysis.py -f "~{sample_id}.~{reference_id}.vacmap.sniffles2.INV.DUP.collapsed.csv"

        bgzip --version
        bgzip "~{sample_id}.~{reference_id}.vacmap.sniffles2.INV.DUP.collapsed.vcf"
        tabix --version
        tabix -p vcf "~{sample_id}.~{reference_id}.vacmap.sniffles2.INV.DUP.collapsed.vcf.gz"

    >>>

    output {
        File vacmap_sniffles2_vcf = "~{sample_id}.~{reference_id}.vacmap.sniffles2.INV.DUP.collapsed.vcf.gz"
        File vacmap_sniffles2_vcf_index =  "~{sample_id}.~{reference_id}.vacmap.sniffles2.INV.DUP.collapsed.vcf.gz.tbi"
    }

    runtime {
        cpu: threads
        memory: mem_gb + " GB"
    }
}
