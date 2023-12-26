
NXF_ANSI_LOG=false nextflow-22.10.3-all run ./main.nf \
  -with-docker "kristinagagalova/pante:pante2-v1.0.0" \
  -profile standard \
  -resume \
  --genomes "test/*.fasta" \
  --outdir "test/results" #\

#--dfam_hmm "/scratch/y95/kgagalova/pante2/dfam_db/dfam.hmm.gz"
# --repbase downloads/RepBaseRepeatMaskerEdition-20181026.tar.gz
#  --rm_meta downloads/RepeatMaskerMetaData-20181026.tar.gz
# --species "Ascomycota"
# --rnammer
