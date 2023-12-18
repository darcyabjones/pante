
NXF_ANSI_LOG=false /home/ubuntu/src/nextflow-22.10.7/build/releases/nextflow-22.10.7-all run KristinaGagalova\pante2 -r master \
  -profile singularity,standard \
  -resume \
  --genomes "test/*.fasta" \
  --outdir "test/results" \
  --dfam_hmm "/scratch/y95/kgagalova/pante2/dfam_db/dfam.hmm.gz"

# --repbase downloads/RepBaseRepeatMaskerEdition-20181026.tar.gz
#  --rm_meta downloads/RepeatMaskerMetaData-20181026.tar.gz
# --species "Ascomycota"
# --rnammer
