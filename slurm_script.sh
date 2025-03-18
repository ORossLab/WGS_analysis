#!/bin/bash

#SBATCH --mail-user=gavrielatos.marios@mayo.edu
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_50          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --job-name=wgs_wdl                                 # Job name
#SBATCH --partition=cpu-short                       ####
#SBATCH --nodes=1                                    # Number of nodes
#SBATCH --ntasks=10                                   # Useful if you have commands that you want to run in parallel
# #SBATCH --cpus-per-task=10                           # Number of CPUs required per task
# #SBATCH --ntasks-per-node=4                         # How many tasks on each node
# #SBATCH --ntasks-per-socket=2                       # How many tasks on each CPU or socke

#SBATCH --time=80:00:00                                    # Time requested
#SBATCH --mem=12G                                    # Memory/node (or mem-per-cpu)
# #SBATCH --mem-per-cpu=8G
#SBATCH --chdir /home/mayo/m301955/my_data/my_projects/Genome_Assembly/WGS_workflow_runs/WGS_analysis_pipeline
#SBATCH --output logs/%x.%N.%j.stdout
#SBATCH --error logs/%x.%N.%j.stderr

module load java
module load apptainer
module load conda
source activate /home/mayo/m301955/mambaforge/envs/wdl_workflow_env
source activate --stack /home/mayo/m301955/mambaforge/envs/Winnowmap
source activate --stack /home/mayo/m301955/mambaforge/envs/vacmap_env

java -Dconfig.file=/home/mayo/m301955/.config/cromwell_permCache.conf -jar \
/home/mayo/m301955/my_data/tools/cromwell-87.jar \
run /home/mayo/m301955/my_data/my_projects/Genome_Assembly/WGS_workflow_runs/WGS_analysis_pipeline/individual_sample_analysis.wdl \
-i /home/mayo/m301955/my_data/my_projects/Genome_Assembly/WGS_workflow_runs/WGS_analysis_pipeline/sample_test_input.json