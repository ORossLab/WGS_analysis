#!/bin/bash

# Function to display usage instructions
usage() {
    echo "Usage: $0 --sample_id <sample_id> --input_vcf <input_vcf> --ref <reference> --ref_name <ref_name>"
    exit 1
}

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --sample_id) sample_id="$2"; shift ;;
        --input_vcf) input_vcf="$2"; shift ;;
        --ref) ref="$2"; shift ;;
        --ref_name) ref_name="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; usage ;;
    esac
    shift
done

# Check if all required arguments are provided
if [[ -z "$sample_id" || -z "$input_vcf" || -z "$ref" || -z "$ref_name" ]]; then
    echo "Error: Missing required arguments."
    usage
fi

pctsize=0.8
pctovl=0.2

merge_vcf="${sample_id}.${ref_name}.vacmap.sniffles2.INV.DUP.merge.vcf"
collapsed_vcf="${sample_id}.${ref_name}.vacmap.sniffles2.INV.DUP.collapsed.vcf"
merge_jl="${sample_id}.${ref_name}.vacmap.sniffles2.INV.DUP.merge.jl"
collapsed_jl="${sample_id}.${ref_name}.vacmap.sniffles2.INV.DUP.collapsed.jl"

truvari collapse -i $input_vcf -o $merge_vcf \
-c $collapsed_vcf -f ${ref} --pctsize $pctsize --pctovl $pctovl \
--sizemax 3000000 --passonly

truvari vcf2df -i -f $merge_vcf $merge_jl
truvari vcf2df -i -f $collapsed_vcf $collapsed_jl


python -c "
import pandas as pd
import joblib
import glob

# Load all .jl files
merge_jl = glob.glob('*merge.jl')
merge_df = joblib.load(merge_jl[0])
merge_df.to_csv(f'{merge_jl[0][:-3]}.csv', sep='#', index=False)

collapsed_jl = glob.glob('*collapsed.jl')
collapsed_df = joblib.load(collapsed_jl[0])
collapsed_df.to_csv(f'{collapsed_jl[0][:-3]}.csv', sep='#', index=False)
"
