
NXF_ANSI_LOG=false nextflow-22.10.3-all run ./main.nf \
  -profile docker \
  -resume \
  --genomes "test/*.fasta" \
  --outdir "test/results" #\

#--dfam_hmm "/scratch/y95/kgagalova/pante2/dfam_db/dfam.hmm.gz"
# --repbase downloads/RepBaseRepeatMaskerEdition-20181026.tar.gz
#  --rm_meta downloads/RepeatMaskerMetaData-20181026.tar.gz
# --species "Ascomycota"
# --rnammer
