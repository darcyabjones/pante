NXF_ANSI_LOG=false nextflow run ./main.nf \
  -profile singularity,standard \
  -resume \
  --genomes "test/*.fasta" \
  --outdir "test/results"

# --rnammer \

#  --repbase containers/downloads/RepBaseRepeatMaskerEdition-20181026.tar.gz \

#  --rm_meta containers/downloads/RepeatMaskerMetaData-20181026.tar.gz \

