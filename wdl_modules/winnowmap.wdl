version 1.0

task winnowmap_alignment {
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
        if [ ! -f "repetitive_k15_GRCh38.txt" ]; then
            echo "Pre-computing high frequency k-mers index does not exist. Creating them."
            meryl count threads=16 k=15 output meryl_GRCh38_DB $grch38
            meryl print greater-than distinct=0.9998 meryl_GRCh38_DB > repetitive_k15_GRCh38.txt
        else
            echo "Pre-computing high frequency k-mers exists"
        fi

        winnowmap \
            -t ~{threads} \
            -W repetitive_k15_GRCh38.txt \
            -ax map-pb -Y -L --cs \
            ~{reference} \
            ~{fastq} > unsorted.sam

        samtools view -hb -@~{threads} unsorted.sam > unsorted.bam
        samtools sort -@~{threads} unsorted.bam -o ~{sample_id}.~{reference_id}.winnowmap.bam 

        samtools index -@16 ~{sample_id}.~{reference_id}.winnowmap.bam
    >>>

    output {
        File aligned_bam = "~{sample_id}.~{reference_id}.winnowmap.bam"
        File aligned_bam_index = "~{sample_id}.~{reference_id}.winnowmap.bam.bai"
    }

    runtime {
        cpu: threads
        memory: mem_gb + "GB"
    }
}
