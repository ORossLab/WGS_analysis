version 1.0

task pbmm2_alignment {
    input {
        String sample_id
        File ubam
        File reference
        File reference_index
        String reference_id
        String REF_PATH_DIR
    }

    Int threads = 60
    Int mem_gb = 90

    command <<<
        set -euo pipefail

        pbmm2 --version

        # Create index if necessary 
        if [ ! -f ~{reference}.mmi ]; then
            echo "Reference index does not exist. Creating file index ~{reference}.mmi."
            pbmm2 index \
                ~{reference} \
                ~{reference}.mmi 
        else
            echo "Reference index exists: ~{reference}.mmi"
        fi

        # Run alignment
        pbmm2 align \
            -j ~{threads} \
            --sort \
            --preset HIFI \
            --log-level INFO \
            ~{reference}.mmi \
            ~{ubam} \
            --sample ~{sample_id} \
            ~{sample_id}.~{reference_id}.pbmm2.bam    >>>
    
    output {
        File aligned_bam = "~{sample_id}.~{reference_id}.pbmm2.bam"
        File aligned_bam_index = "~{sample_id}.~{reference_id}.pbmm2.bam.bai"
    }

    runtime {
        cpu: threads
        memory: mem_gb + "GB"
    }
}
