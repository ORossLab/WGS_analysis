version 1.0

task trgt {
    input {
        String sample_id
        String? sex

        File bam
        File bam_index

        File reference
        File reference_index
        String reference_id
    
        File tandem_repeat_bed
    }

    Boolean sex_defined = defined(sex)
    String karyotype = if select_first([sex, "FEMALE"]) == "MALE" then "XY" else "XX"

    Int threads = 20
    Int mem_gb = 32

    command <<<
        set -euo pipefail

        echo ~{if sex_defined then "" else "Sex is not defined for ~{sample_id}. Defaulting to karyotype XX for TRGT."}

        trgt --version

        trgt genotype \
            --threads ~{threads} \
            --karyotype ~{karyotype} \
            --genome ~{reference} \
            --repeats ~{tandem_repeat_bed} \
            --reads ~{bam} \
            --output-prefix ~{sample_id}.~{reference_id}.trgt


        bcftools --version
        bcftools sort -Oz -o ~{sample_id}.~{reference_id}.trgt.sorted.vcf.gz ~{sample_id}.~{reference_id}.trgt.vcf.gz
        bcftools index --tbi ~{sample_id}.~{reference_id}.trgt.sorted.vcf.gz

        samtools --version
        samtools sort -@10 -o ~{sample_id}.~{reference_id}.trgt.spanning.sorted.bam ~{sample_id}.~{reference_id}.trgt.spanning.bam
        samtools index -@10 ~{sample_id}.~{reference_id}.trgt.spanning.sorted.bam
    >>>

    output {
        File repeat_vcf = "~{sample_id}.~{reference_id}.trgt.sorted.vcf.gz"
        File repeat_vcf_index = "~{sample_id}.~{reference_id}.trgt.sorted.vcf.gz.tbi"
        File spanning_reads = "~{sample_id}.~{reference_id}.trgt.spanning.sorted.bam"
        File spanning_reads_index = "~{sample_id}.~{reference_id}.trgt.spanning.sorted.bam.bai"
    }
    
    runtime {
        cpu: threads
        memory: mem_gb + " GB"
    }
}
