version 1.0

task vacmap_alignment {
    input {
        String sample_id
        File fastq
        File reference
        File reference_index
        String reference_id
    }

    Int threads = 68
    Int mem_gb = 90

    command <<<
        set -euo pipefail

        # Run alignment
        vacmap \
            -ref ~{reference} -read ~{fastq} -mode L -t ~{threads-8} --eqx --MD --cs --L | samtools sort -@8 > ~{sample_id}.~{reference_id}.vacmap.bam

        # Index alignment
        samtools index -@20 ~{sample_id}.~{reference_id}.vacmap.bam
    >>>

    output {
        File aligned_bam = "~{sample_id}.~{reference_id}.vacmap.bam"
        File aligned_bam_index = "~{sample_id}.~{reference_id}.vacmap.bam.bai"
    }

    runtime {
        cpu: threads
        memory: mem_gb + "GB"
    }
}
