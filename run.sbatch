#!/bin/bash --login

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --partition=work
#SBATCH --time=1-00:00:00
#SBATCH --account=y95
#SBATCH --mail-type=ALL
#SBATCH --export=NONE

module load nextflow/22.10.0
module load singularity/3.11.4-slurm


NXF_ANSI_LOG=false nextflow run ./main.nf \
  -profile singularity,pawsey_setonix \
  -resume \
  --genomes "test/*.fasta" \
  --outdir "test/results"

#  --dfam_hmm "/scratch/y95/kgagalova/pante2/dfam_db/dfam.hmm.gz"
# --repbase downloads/RepBaseRepeatMaskerEdition-20181026.tar.gz
#  --rm_meta downloads/RepeatMaskerMetaData-20181026.tar.gz
# --species "Ascomycota"
# --rnammer
