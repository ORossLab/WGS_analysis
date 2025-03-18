version 1.0

task deepvariant {
    input {
        String sample_id
        File aligned_bam
        File aligned_bam_index
        File reference
        File reference_index
        String REF_PATH_DIR
    }

    Int threads = 94
    Int mem_gb = 128

    command <<<
        set -euo pipefail

        WD_PATH=$PWD
        DATA_PATH=$PWD
        INPUTS_PATH_ref=$(dirname ~{reference})
        INPUTS_PATH_bam=$(dirname ~{aligned_bam})

        # Create output dir
        # mkdir -p ${WD_PATH}/~{sample_id}

        GVCF_FILE=${WD_PATH}/~{sample_id}.gvcf
        VCF_FILE=${WD_PATH}/~{sample_id}.vcf
        LOG_FILE=${WD_PATH}/LOG_OUT
        TEMP_FILE=${WD_PATH}/TEMP_OUT

        echo ${WD_PATH}
        echo ${INPUTS_PATH_ref}
        echo ${INPUTS_PATH_bam}
        echo ${DATA_PATH}

        # Run DeepVariant
        apptainer exec --unsquash --bind ${WD_PATH},${DATA_PATH},${INPUTS_PATH_ref},${INPUTS_PATH_bam} \
            /research/bsi/tools/biotools/deepvariant/1.6.0/deepvariant_1.6.0.sif \
            /opt/deepvariant/bin/run_deepvariant \
            --intermediate_results_dir ${TEMP_FILE} \
            --logging_dir ${LOG_FILE} \
            --model_type PACBIO \
            --num_shards=32 \
            --output_gvcf ${GVCF_FILE} \
            --output_vcf ${VCF_FILE} \
            --reads ~{aligned_bam} \
            --ref ~{reference} \
            --sample_name ~{sample_id}

        # Remove temp files
        rm -rf ${TEMP_FILE}
        rm -rf ${LOG_FILE}
    >>>

    output {
        File gvcf_file = "~{sample_id}.gvcf"
        File vcf_file = "~{sample_id}.vcf"
    }

    runtime {
        cpu: threads
        memory: mem_gb + "GB"
    }
}

task bcftools {
    input {
        File vcf
        String stats_params
        File reference
    }

    String vcf_basename = basename(vcf, ".vcf.gz")
    Int threads = 16
    Int mem_gb = 32
    Int reference_size = if (defined(reference)) then ceil(size(reference, "GB")) else 0

    command <<<
        set -euo pipefail

        bcftools --version

        bcftools stats \
            --threads ~{threads - 1} \
            ~{stats_params} \
            ~{"--fasta-ref " + reference} \
            ~{vcf} \
        > ~{vcf_basename}.vcf.stats.txt

        bcftools roh \
            --threads ~{threads - 1} \
            --AF-dflt 0.4 \
            ~{vcf} \
        > ~{vcf_basename}.bcftools_roh.out

        echo -e "#chr\\tstart\\tend\\tqual" > ~{vcf_basename}.roh.bed
        awk -v OFS='\t' '$1=="RG" {{ print $3, $4, $5, $8 }}' \
            ~{vcf_basename}.bcftools_roh.out \
        >> ~{vcf_basename}.roh.bed
    >>>

    output {
        File stats = "~{vcf_basename}.vcf.stats.txt"
        File roh_out = "~{vcf_basename}.bcftools_roh.out"
        File roh_bed = "~{vcf_basename}.roh.bed"
    }

    runtime {
        cpu: threads
        memory: mem_gb + " GB"
    }
}
