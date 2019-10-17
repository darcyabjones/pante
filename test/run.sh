NXF_ANSI_LOG=false nextflow run ./main.nf \
  -profile dev,standard \
  -resume \
  --genomes "test/*.fasta" \
  --repbase containers/downloads/RepBaseRepeatMaskerEdition-20181026.tar.gz \
  --rm_meta containers/downloads/RepeatMaskerMetaData-20181026.tar.gz \
  --pfam ./pfam.stk \
  --noinfernal \
  --rnammer \
  --outdir "test/results"




