NXF_ANSI_LOG=false /home/ubuntu/src/nextflow-22.10.7/nextflow run ./main.nf \
  -profile singularity \
  -resume \
  --genomes "test/*.fasta" \
  --outdir "test/results"

# --repbase downloads/RepBaseRepeatMaskerEdition-20181026.tar.gz
#  --rm_meta downloads/RepeatMaskerMetaData-20181026.tar.gz
# --species "Ascomycota"
# --rnammer
