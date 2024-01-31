


NXF_ANSI_LOG=false /home/kgagalova/src/nextflow-22.10.7/nextflow run ./main.nf \
  -profile standard \
  -resume \
  --rnammer true \
  --genomes "test/*.fasta" \
  --outdir "test/results" \
  --dfam_h5 "/home/kgagalova/pante-debug/dfam38_full.0.h5.gz" \
  --repbase "/home/kgagalova/pante-debug/RepBaseRepeatMaskerEdition-20181026.tar.gz" \
  --rm_meta "/home/kgagalova/pante-debug/RepeatMaskerMetaData-20181026.tar.gz" \
  -with-singularity "containers/singularity/pante2-rnammer.sif" 

#   --rm_species "fungi" \
#--dfam_hmm "/scratch/y95/kgagalova/pante2/dfam_db/dfam.hmm.gz"
# --repbase downloads/RepBaseRepeatMaskerEdition-20181026.tar.gz
#  --rm_meta downloads/RepeatMaskerMetaData-20181026.tar.gz
# --species "Ascomycota"
# --rnammer
