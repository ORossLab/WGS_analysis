version 1.0

task mosdepth {
    input {
        File aligned_bam
        File aligned_bam_index
        String sample_id
    }

    Int threads = 32
    Int mem_gb = 64

    command <<<
        set -euo pipefail

        mosdepth --version

        mosdepth \
            --threads ~{threads - 1} \
            --by 500 \
            --no-per-base \
            --use-median \
            ~{sample_id} \
            ~{aligned_bam}
    >>>

    output {
        File summary = "~{sample_id}.mosdepth.summary.txt"
        File region_bed = "~{sample_id}.regions.bed.gz"
    }

    runtime {
        cpu: threads
        memory: mem_gb + "GB"
    }
}
