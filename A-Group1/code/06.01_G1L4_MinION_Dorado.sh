#!/bin/bash
## Job Name
#SBATCH --job-name=G1L4_MinION_Dorado
## Allocation Definition
#SBATCH --account=srlab-ckpt
#SBATCH --partition=ckpt
## Resources
## GPU
#SBATCH --gres=gpu:a40:1
## Nodes
#SBATCH --nodes=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=0-12:00:00
## Memory per node
#SBATCH --mem=120G
##turn on e-mail notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kdurkin1@uw.edu
## Specify the working directory for this job
#SBATCH --chdir=/gscratch/srlab/kdurkin1/SIFP-nanopore/A-Group1/output/06.01-G1-Library4-MinION-Dorado-recall-GPU/

## Script for running ONT Dorado to perform
## basecalling (i.e. convert raw ONT pod5 to FastQ) of NanaPore data generated
## Summer 2025, as psrt of K.Durkin SIFP project

## This script utilizes a GPU node. These nodes are only available as part of the checkpoint
## partition/account. Since we don't own a GPU node, our GPU jobs are lowest priority and
## can be interrupted at any time if the node owner submits a new job.

###################################################################################
# These variables need to be set by user

wd=$(pwd)

# Programs array
# declare -A programs_array
# programs_array=(
# [dorado]="apptainer exec --nv --bind /gscratch /gscratch/srlab/containers/srlab-R4.4-bioinformatics-container-3886a1c.sif dorado"
# )


# Establish variables for more readable code

# Input files directory
raw_pod5_dir=/gscratch/srlab/kdurkin1/SIFP-nanopore/A-Group1/data/06.01-G1-Library4-MinION-Dorado-recall-GPU/
output_dir=/gscratch/srlab/kdurkin1/SIFP-nanopore/A-Group1/output/06.01-G1-Library4-MinION-Dorado-recall-GPU/
genome_file=/gscratch/srlab/kdurkin1/SIFP-nanopore/data/GCA_965233905.1_jaEunKnig1.1/GCA_965233905.1_jaEunKnig1.1_genomic.fna

# Output directory
out_dir=${wd}

# CPU threads
threads=28

# Sequencing kit used
kit="SQK-NBD114-96"

# Flow Cell ID
flow_cell_id="FBD09922"

# GPU devices setting
GPU_devices=auto

###################################################################################

# Exit script if any command fails
set -e

# Load CUDA GPU module
module load cuda/12.9.1

apptainer exec \
--nv \
--home "$PWD" \
--bind /mmfs1/home/ \
--bind /gscratch \
/gscratch/srlab/containers/srlab-R4.4-bioinformatics-container-3886a1c.sif \
dorado basecaller \
hac \
-r ${raw_pod5_dir}/ \
--kit-name SQK-NBD114-96 \
--trim 'all' \
--reference ${genome_file} \
--modified-bases 5mCG_5hmCG 6mA \
--device ${GPU_devices} \
> ${output_dir}/FBD09922_pass_recalled.bam


###################################################################################

# Document programs in PATH (primarily for program version ID)
{
date
echo ""
echo "System PATH for $SLURM_JOB_ID"
echo ""
printf "%0.s-" {1..10}
echo "${PATH}" | tr : n
} >> system_path.log


# Capture program options
for program in "${!programs_array[@]}"
do
	{
  echo "Program options for ${program}: "
	echo ""
	${programs_array[$program]} --help
	echo ""
	echo ""
	echo "----------------------------------------------"
	echo ""
	echo ""
} &>> program_options.log || true
done
