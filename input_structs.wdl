version 1.0

# Struct for SampleInfo
struct SampleInfo {
    String sample_id
    File ubam
    File fastq
    String sex
}

# Struct for ReferenceData
struct ReferenceData {
    File reference
    File reference_index
    String REF_PATH_DIR
    String reference_id
    File exclude_bed
    File exclude_bed_index
    File expected_bed_male
    File expected_bed_female
    File tandem_repeat_bed
}

# Struct for Paths
struct Paths {
    String script_path
}
