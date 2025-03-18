version 1.0

# mamba activate wgs_workflow_env
# mamba activate --stack vacmap_env
# mamba activate --stack Winnowmap
# module load java apptainer

import "input_structs.wdl"
import "/home/ext_gavrielatos_marios_mayo_edu/tools/WGS_analysis_pipeline/wdl_modules/pbmm2.wdl" as Pbmm2
import "/home/ext_gavrielatos_marios_mayo_edu/tools/WGS_analysis_pipeline/wdl_modules/vacmap.wdl" as VACmap
import "/home/ext_gavrielatos_marios_mayo_edu/tools/WGS_analysis_pipeline/wdl_modules/winnowmap.wdl" as Winnowmap
import "/home/ext_gavrielatos_marios_mayo_edu/tools/WGS_analysis_pipeline/wdl_modules/pbsv.wdl" as Pbsv
import "/home/ext_gavrielatos_marios_mayo_edu/tools/WGS_analysis_pipeline/wdl_modules/deepvariant.wdl" as DeepVariant
import "/home/ext_gavrielatos_marios_mayo_edu/tools/WGS_analysis_pipeline/wdl_modules/cutesv.wdl" as Cutesv
import "/home/ext_gavrielatos_marios_mayo_edu/tools/WGS_analysis_pipeline/wdl_modules/sniffles2.wdl" as Sniffles2
import "/home/ext_gavrielatos_marios_mayo_edu/tools/WGS_analysis_pipeline/wdl_modules/sniffles2_vacmap_filter.wdl" as Sniffles2_Vacmap_filter
import "/home/ext_gavrielatos_marios_mayo_edu/tools/WGS_analysis_pipeline/wdl_modules/hificnv.wdl" as Hificnv
import "/home/ext_gavrielatos_marios_mayo_edu/tools/WGS_analysis_pipeline/wdl_modules/trgt.wdl" as Trgt
import "/home/ext_gavrielatos_marios_mayo_edu/tools/WGS_analysis_pipeline/wdl_modules/mosdepth.wdl" as Mosdepth

workflow individual_sample_analysis {
    input {
        SampleInfo sample_info
        ReferenceData reference_data
        Paths paths
    }

    # Group 1: pbmm2_alignment and vacmap_alignment run in parallel
    scatter (alignment in [1]) {
        # fixed
        call Pbmm2.pbmm2_alignment {
            input:
                sample_id = sample_info.sample_id,
                ubam = sample_info.ubam,
                reference = reference_data.reference,
                reference_index = reference_data.reference_index,
                reference_id = reference_data.reference_id,
                REF_PATH_DIR = reference_data.REF_PATH_DIR
            }

        # fixed
        call VACmap.vacmap_alignment {
            input:
                sample_id = sample_info.sample_id,
                fastq = sample_info.fastq,
                reference = reference_data.reference,
                reference_index = reference_data.reference_index,
                reference_id = reference_data.reference_id,
            }

        call Winnowmap.winnowmap_alignment {
            input:
                sample_id = sample_info.sample_id,
                fastq = sample_info.fastq,
                reference = reference_data.reference,
                reference_index = reference_data.reference_index,
                reference_id = reference_data.reference_id,
            }
    }

    # Extract the first result from each scatter block array
    File pbmm2_aligned_bam = pbmm2_alignment.aligned_bam[0]
    File pbmm2_aligned_bam_index = pbmm2_alignment.aligned_bam_index[0]

    File vacmap_aligned_bam = vacmap_alignment.aligned_bam[0]
    File vacmap_aligned_bam_index = vacmap_alignment.aligned_bam_index[0]

    File winnowmap_aligned_bam = winnowmap_alignment.aligned_bam[0]
    File winnowmap_aligned_bam_index = winnowmap_alignment.aligned_bam_index[0]

    # Group 2: Run deepvariant along with pbsv, sniffles2, cutesv, and sniffles2_vacmap_filter in parallel
    scatter (variant_calling in [1]) {

        # fixed
        call DeepVariant.deepvariant {
            input:
                sample_id = sample_info.sample_id,
                aligned_bam = pbmm2_aligned_bam,
                aligned_bam_index = pbmm2_aligned_bam_index,
                reference = reference_data.reference,
                reference_index = reference_data.reference_index,
                REF_PATH_DIR = reference_data.REF_PATH_DIR,
        }

        # fixed
        call Pbsv.pbsv {
            input:
                sample_id = sample_info.sample_id,
                aligned_bam = pbmm2_aligned_bam,
                aligned_bam_index = pbmm2_aligned_bam_index,
                reference = reference_data.reference,
                reference_id = reference_data.reference_id,
                aligner = 'pbmm2'
        }

        call Sniffles2.sniffles2 as sniffles2_pbmm2 {
            input:
                sample_id = sample_info.sample_id,
                aligned_bam = pbmm2_aligned_bam,
                aligned_bam_index = pbmm2_aligned_bam_index,
                reference = reference_data.reference,
                reference_id = reference_data.reference_id,
                aligner = 'pbmm2'
        }
        call Sniffles2.sniffles2 as sniffles2_winnowmap {
            input:
                sample_id = sample_info.sample_id,
                aligned_bam = winnowmap_aligned_bam,
                aligned_bam_index = winnowmap_aligned_bam_index,
                reference = reference_data.reference,
                reference_id = reference_data.reference_id,
                aligner = 'winnowmap'
        }

        call Cutesv.cutesv as cutesv_pbmm2 {
            input:
                sample_id = sample_info.sample_id,
                aligned_bam = pbmm2_aligned_bam,
                aligned_bam_index = pbmm2_aligned_bam_index,
                reference = reference_data.reference,
                reference_id = reference_data.reference_id,
                aligner = 'pbmm2'
        }
        call Cutesv.cutesv as cutesv_winnowmap {
            input:
                sample_id = sample_info.sample_id,
                aligned_bam = winnowmap_aligned_bam,
                aligned_bam_index = winnowmap_aligned_bam_index,
                reference = reference_data.reference,
                reference_id = reference_data.reference_id,
                aligner = 'winnowmap'
        }


        call Sniffles2_Vacmap_filter.sniffles2_vacmap_filter {
            input :
                sample_id = sample_info.sample_id,
                aligned_bam = vacmap_aligned_bam,
                aligned_bam_index = vacmap_aligned_bam_index,
                reference = reference_data.reference,
                reference_id = reference_data.reference_id,
                script_path = paths.script_path
        }
    }

    File snp_vcf = deepvariant.vcf_file[0]
    File snp_gvcf = deepvariant.gvcf_file[0]

    call DeepVariant.bcftools {
            input:
                vcf = snp_vcf,
                stats_params = "--apply-filters PASS --samples ~{sample_info.sample_id}",
                reference = reference_data.reference,
        }

    call combine_vcfs {
        input:
            pbmm2_pbsv_vcf = pbsv.pbsv_vcf,
            pbmm2_cutesv_vcf = cutesv_pbmm2.cutesv_vcf,
            pbmm2_sniffles2_vcf = sniffles2_pbmm2.sniffles2_vcf,
            winnowmap_cutesv_vcf = cutesv_winnowmap.cutesv_vcf,
            winnowmap_sniffles2_vcf = sniffles2_winnowmap.sniffles2_vcf,
            sniffles2_vacmap_vcf = sniffles2_vacmap_filter.vacmap_sniffles2_vcf,
            sample_id = sample_info.sample_id,
            reference_id = reference_data.reference_id
    }

    call Mosdepth.mosdepth {
        input:
            aligned_bam = pbmm2_aligned_bam,
            aligned_bam_index = pbmm2_aligned_bam_index,
            sample_id = sample_info.sample_id,
    }



    call Trgt.trgt {
        input:
            sample_id = sample_info.sample_id,
            sex = sample_info.sex,

            bam = pbmm2_aligned_bam,
            bam_index = pbmm2_aligned_bam_index,

            reference = reference_data.reference,
            reference_index = reference_data.reference_index,
            reference_id = reference_data.reference_id,

            tandem_repeat_bed = reference_data.tandem_repeat_bed,
    }



    call Hificnv.hificnv {
        input:
            sample_id = sample_info.sample_id,
            sex = sample_info.sex,

            bam = pbmm2_aligned_bam,
            bam_index = pbmm2_aligned_bam_index,

            # phased_vcf = hiphase.phased_vcfs[0].data,
            # phased_vcf_index = hiphase.phased_vcfs[0].data_index,

            reference = reference_data.reference,
            reference_index = reference_data.reference_index,
            reference_id = reference_data.reference_id,

            exclude_bed = reference_data.exclude_bed,
            exclude_bed_index = reference_data.exclude_bed_index,
            expected_bed_male = reference_data.expected_bed_male,
            expected_bed_female = reference_data.expected_bed_female,
    }




    call organize_outputs {
        input:
            sample_id = sample_info.sample_id,
            reference_id = reference_data.reference_id,
            
            # Alignments
            pbmm2_bam = pbmm2_aligned_bam,
            pbmm2_bam_index = pbmm2_aligned_bam_index,
            vacmap_bam = vacmap_aligned_bam,
            vacmap_bam_index = vacmap_aligned_bam_index,
            winnowmap_bam = winnowmap_aligned_bam,
            winnowmap_bam_index = winnowmap_aligned_bam_index,
            
            # Coverage
            mosdepth_summary = mosdepth.summary,
            mosdepth_region_bed = mosdepth.region_bed,
            
            # Small variants
            deepvariant_vcf = snp_vcf,
            deepvariant_gvcf = snp_gvcf,
            bcftools_stats = bcftools.stats,
            bcftools_roh_out = bcftools.roh_out,
            bcftools_roh_bed = bcftools.roh_bed,
            
            # Structural variants
            pbsv_vcf = pbsv.pbsv_vcf[0],
            pbsv_vcf_index = pbsv.pbsv_vcf_index[0],
            pbmm2_sniffles2_vcf = sniffles2_pbmm2.sniffles2_vcf[0],
            pbmm2_sniffles2_vcf_index = sniffles2_pbmm2.sniffles2_vcf_index[0],
            winnowmap_sniffles2_vcf = sniffles2_winnowmap.sniffles2_vcf[0],
            winnowmap_sniffles2_vcf_index = sniffles2_winnowmap.sniffles2_vcf_index[0],
            pbmm2_cutesv_vcf = cutesv_pbmm2.cutesv_vcf[0],
            pbmm2_cutesv_vcf_index = cutesv_pbmm2.cutesv_vcf_index[0],
            winnowmap_cutesv_vcf = cutesv_winnowmap.cutesv_vcf[0],
            winnowmap_cutesv_vcf_index = cutesv_winnowmap.cutesv_vcf_index[0],
            vacmap_sniffles2_vcf = sniffles2_vacmap_filter.vacmap_sniffles2_vcf[0],
            vacmap_sniffles2_vcf_index = sniffles2_vacmap_filter.vacmap_sniffles2_vcf_index[0],
            merged_sv_vcf = combine_vcfs.merged_vcf,
            merged_sv_vcf_index = combine_vcfs.merged_vcf_index,
            
            # Repeats
            trgt_vcf = trgt.repeat_vcf,
            trgt_vcf_index = trgt.repeat_vcf_index,
            trgt_spanning_reads = trgt.spanning_reads,
            trgt_spanning_reads_index = trgt.spanning_reads_index,
            
            # CNV
            cnv_vcf = hificnv.cnv_vcf,
            cnv_vcf_index = hificnv.cnv_vcf_index,
            cnv_bedgraph = hificnv.copynum_bedgraph,
            cnv_depth = hificnv.depth_bw
    }

    output {
        File organized_outputs = organize_outputs.organized_outputs
        File output_md5sums = organize_outputs.md5sums
    }
}


task combine_vcfs {
        input {
        Array[File] pbmm2_pbsv_vcf
        Array[File] pbmm2_cutesv_vcf
        Array[File] pbmm2_sniffles2_vcf
        Array[File] winnowmap_cutesv_vcf
        Array[File] winnowmap_sniffles2_vcf
        Array[File] sniffles2_vacmap_vcf
        String sample_id
        String reference_id
    }

    Int threads = 8
    Int mem_gb = 16

    command <<<
        set -euo pipefail

        # Merge VCFs from arrays
        svdb --merge --vcf ~{sep=" " pbmm2_pbsv_vcf} \
                          ~{sep=" " pbmm2_cutesv_vcf} \
                          ~{sep=" " pbmm2_sniffles2_vcf} \
                          ~{sep=" " winnowmap_cutesv_vcf} \
                          ~{sep=" " winnowmap_sniffles2_vcf} \
                          ~{sep=" " sniffles2_vacmap_vcf} \
             --no_intra --overlap 1 > ~{sample_id}.~{reference_id}.merged.vcf

        bgzip --version
        bgzip "~{sample_id}.~{reference_id}.merged.vcf"
        tabix --version
        tabix -p vcf "~{sample_id}.~{reference_id}.merged.vcf.gz"
    >>>

    output {
        File merged_vcf = "~{sample_id}.~{reference_id}.merged.vcf.gz"
        File merged_vcf_index = "~{sample_id}.~{reference_id}.merged.vcf.gz.tbi"
    }

    runtime {
        cpu: threads
        memory: mem_gb + " GB"
    }
}

task organize_outputs {
    input {
        String sample_id
        String reference_id
        
        # Alignments
        File pbmm2_bam
        File pbmm2_bam_index
        File vacmap_bam
        File vacmap_bam_index
        File winnowmap_bam
        File winnowmap_bam_index
        
        # Coverage
        File mosdepth_summary
        File mosdepth_region_bed
        
        # Small variants
        File deepvariant_vcf
        File deepvariant_gvcf
        File bcftools_stats
        File bcftools_roh_out
        File bcftools_roh_bed
        
        # Structural variants
        File pbsv_vcf
        File pbsv_vcf_index
        File pbmm2_sniffles2_vcf
        File pbmm2_sniffles2_vcf_index
        File winnowmap_sniffles2_vcf
        File winnowmap_sniffles2_vcf_index
        File pbmm2_cutesv_vcf
        File pbmm2_cutesv_vcf_index
        File winnowmap_cutesv_vcf
        File winnowmap_cutesv_vcf_index
        File vacmap_sniffles2_vcf
        File vacmap_sniffles2_vcf_index
        File merged_sv_vcf
        File merged_sv_vcf_index
        
        # Repeats
        File trgt_vcf
        File trgt_vcf_index
        File trgt_spanning_reads
        File trgt_spanning_reads_index
        
        # CNV
        File cnv_vcf
        File cnv_vcf_index
        File cnv_bedgraph
        File cnv_depth
    }
    
    command <<<
        # Create directory structure
        mkdir -p ~{sample_id}/{alignments/{pbmm2,vacmap,winnowmap},coverage,variants/{small_variants,structural_variants/{individual_callers,merged}},repeats,cnv}
        
        # Copy alignments
        cp ~{pbmm2_bam} ~{pbmm2_bam_index} ~{sample_id}/alignments/pbmm2/
        cp ~{vacmap_bam} ~{vacmap_bam_index} ~{sample_id}/alignments/vacmap/
        cp ~{winnowmap_bam} ~{winnowmap_bam_index} ~{sample_id}/alignments/winnowmap/
        
        # Copy coverage
        cp ~{mosdepth_summary} ~{mosdepth_region_bed} ~{sample_id}/coverage/
        
        # Copy small variants
        cp ~{deepvariant_vcf} ~{deepvariant_gvcf} ~{bcftools_stats} ~{bcftools_roh_out} ~{bcftools_roh_bed} ~{sample_id}/variants/small_variants/
        
        # Copy structural variants
        cp ~{pbsv_vcf} ~{pbsv_vcf_index} \
           ~{pbmm2_sniffles2_vcf} ~{pbmm2_sniffles2_vcf_index} \
           ~{winnowmap_sniffles2_vcf} ~{winnowmap_sniffles2_vcf_index} \
           ~{pbmm2_cutesv_vcf} ~{pbmm2_cutesv_vcf_index} \
           ~{winnowmap_cutesv_vcf} ~{winnowmap_cutesv_vcf_index} \
           ~{vacmap_sniffles2_vcf} ~{vacmap_sniffles2_vcf_index} \
           ~{sample_id}/variants/structural_variants/individual_callers/
           
        cp ~{merged_sv_vcf} ~{merged_sv_vcf_index} ~{sample_id}/variants/structural_variants/merged/
        
        # Copy repeats
        cp ~{trgt_vcf} ~{trgt_vcf_index} ~{trgt_spanning_reads} ~{trgt_spanning_reads_index} ~{sample_id}/repeats/
        
        # Copy CNV
        cp ~{cnv_vcf} ~{cnv_vcf_index} ~{cnv_bedgraph} ~{cnv_depth} ~{sample_id}/cnv/
        
        # Create MD5sums
        cd ~{sample_id}
        find . -type f -exec md5sum {} \; > md5sums.txt
        
        # Create tarball of organized outputs
        cd ..
        tar -czf ~{sample_id}_organized_outputs.tar.gz ~{sample_id}
    >>>
    
    output {
        File organized_outputs = "~{sample_id}_organized_outputs.tar.gz"
        File md5sums = "~{sample_id}/md5sums.txt"
    }
    
    runtime {
        cpu: 1
        memory: "4 GB"
    }
}