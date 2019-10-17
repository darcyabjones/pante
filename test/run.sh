NXF_ANSI_LOG=false nextflow run ./main.nf \
  -profile singularity,standard \
  -resume \
  --genomes "test/*.fasta" \
  --outdir "test/results"

# --repbase downloads/RepBaseRepeatMaskerEdition-20181026.tar.gz
#  --rm_meta downloads/RepeatMaskerMetaData-20181026.tar.gz
# --species "Ascomycota"
# --rnammer
