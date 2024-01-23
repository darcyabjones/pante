


NXF_ANSI_LOG=false /home/kgagalova/src/nextflow-22.10.7/nextflow run ./main.nf \
  -profile standard \
  -resume \
  --genomes "test/*.fasta" \
  --outdir "test/results" \
  --dfam_h5 "/home/kgagalova/pante-debug/dfam38_full.0.h5.gz" \
  -with-singularity "containers/singularity/pante2.sif" 

#--dfam_hmm "/scratch/y95/kgagalova/pante2/dfam_db/dfam.hmm.gz"
# --repbase downloads/RepBaseRepeatMaskerEdition-20181026.tar.gz
#  --rm_meta downloads/RepeatMaskerMetaData-20181026.tar.gz
# --species "Ascomycota"
# --rnammer
