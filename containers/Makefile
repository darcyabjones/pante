PWD=$(shell pwd)
DEBIAN_VERSION=buster-20190708-slim
TAG=pante-v0.0.1
DOCKER_USER=darcyabjones
BASE_IMAGE=$(DOCKER_USER)/base:$(TAG)

INSTALL_BASE=/opt

define make_image_name
$(DOCKER_USER)/$(notdir $1):$(TAG)
endef

# SOFTWARE_VERSIONs
GENOMETOOLS_VERSION=1.5.10
GENOMETOOLS_URL=http://genometools.org/pub/genometools-$(GENOMETOOLS_VERSION).tar.gz
GENOMETOOLS_PREFIX_ARG=$(INSTALL_BASE)/genometools/$(GENOMETOOLS_VERSION)

GFFPAL_TAG=v0.0.2
GFFPAL_REPO=https://github.com/darcyabjones/gffpal.git
GFFPAL_PREFIX_ARG=$(INSTALL_BASE)/gffpal/$(GFFPAL_TAG)

SEQRENAMER_TAG=v0.0.1
SEQRENAMER_REPO=https://github.com/darcyabjones/seqrenamer.git
SEQRENAMER_PREFIX_ARG=$(INSTALL_BASE)/seqrenamer/$(SEQRENAMER_TAG)

RMBLAST_VERSION=2.9.0+
RMBLAST_URL=ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/ncbi-blast-2.9.0+-src.tar.gz
RMBLAST_PATCH_URL=http://www.repeatmasker.org/isb-2.9.0+-rmblast.patch.gz
RMBLAST_PREFIX_ARG=$(INSTALL_BASE)/rmblast/$(RMBLAST_VERSION)

TRF_VERSION=409
TRF_URL=http://tandem.bu.edu/trf/downloads/trf409.linux64
TRF_PREFIX_ARG=$(INSTALL_BASE)/trf/$(TRF_VERSION)

NSEG_URL=ftp://ftp.ncbi.nih.gov/pub/seg/nseg/
NSEG_PREFIX_ARG=$(INSTALL_BASE)/nseg

REPEATSCOUT_VERSION=1.0.5
REPEATSCOUT_URL=http://www.repeatmasker.org/RepeatScout-1.0.5.tar.gz
REPEATSCOUT_PREFIX_ARG=$(INSTALL_BASE)/repeatscout/$(REPEATSCOUT_VERSION)

RECON_VERSION=1.08
RECON_URL=http://www.repeatmasker.org/RepeatModeler/RECON-1.08.tar.gz
RECON_PREFIX_ARG=$(INSTALL_BASE)/recon/$(RECON_VERSION)

RMASK_VERSION=4-0-9
RMASK_URL=http://www.repeatmasker.org/RepeatMasker-open-$(RMASK_VERSION).tar.gz
RMASK_PREFIX_ARG=$(INSTALL_BASE)/repeatmasker/$(RMASK_VERSION)

RMOD_VERSION=1.0.11
RMOD_URL=http://www.repeatmasker.org/RepeatModeler/RepeatModeler-open-$(RMOD_VERSION).tar.gz
RMOD_PREFIX_ARG=$(INSTALL_BASE)/repeatmodeller/$(RMOD_VERSION)


HMMER3_VERSION=3.2.1
HMMER3_URL=http://eddylab.org/software/hmmer/hmmer-3.2.1.tar.gz
HMMER3_PREFIX_ARG=$(INSTALL_BASE)/hmmer3/$(HMMER3_VERSION)

MMSEQS_TAG=9-d36de
MMSEQS_REPO=https://github.com/soedinglab/MMseqs2.git
MMSEQS_PREFIX_ARG=$(INSTALL_BASE)/mmseqs/$(MMSEQS_TAG)

INFERNAL_VERSION=1.1.2
INFERNAL_URL=http://eddylab.org/infernal/infernal-1.1.2.tar.gz
INFERNAL_PREFIX_ARG=$(INSTALL_BASE)/infernal/$(INFERNAL_VERSION)

ARAGORN_VERSION=1.2.38
ARAGORN_URL=http://130.235.244.92/ARAGORN/Downloads/aragorn1.2.38.tgz
ARAGORN_PREFIX_ARG=$(INSTALL_BASE)/aragorn/$(ARAGORN_VERSION)

VIENNA_VERSION=2.4.13
VIENNA_URL=https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.13.tar.gz
VIENNA_PREFIX_ARG=$(INSTALL_BASE)/vienna/$(VIENNA_VERSION)

STRUCTRNAFINDER_COMMIT=3d5479b
STRUCTRNAFINDER_REPO=https://github.com/viniciusmaracaja/structRNAfinder.git
STRUCTRNAFINDER_PREFIX_ARG=$(INSTALL_BASE)/structrnafinder/$(STRUCTRNAFINDER_COMMIT)

TRNASCAN_VERSION=2.0.3
TRNASCAN_URL=http://trna.ucsc.edu/software/trnascan-se-2.0.3.tar.gz
TRNASCAN_PREFIX_ARG=$(INSTALL_BASE)/trnascan/$(TRNASCAN_VERSION)

HMMER2_VERSION=2.3.2
HMMER2_URL=http://eddylab.org/software/hmmer/hmmer-2.3.2.tar.gz
HMMER2_PREFIX_ARG=$(INSTALL_BASE)/hmmer2/$(HMMER2_VERSION)

RNAMMER_VERSION=1.2
RNAMMER_PREFIX_ARG=$(INSTALL_BASE)/rnammer/$(RNAMMER_VERSION)
RNAMMER_TAR=sources/rnammer-1.2.src.tar.Z

OCCULTERCUT_VERSION=1.1
OCCULTERCUT_URL=https://downloads.sourceforge.net/project/occultercut/OcculterCut_v1.1.tar.gz
OCCULTERCUT_PREFIX_ARG=$(INSTALL_BASE)/occultercut/$(OCCULTERCUT_VERSION)


MITEFINDER_COMMIT=833754b
MITEFINDER_REPO=https://github.com/screamer/miteFinder.git
MITEFINDER_PREFIX_ARG=$(INSTALL_BASE)/mitefinder/$(MITEFINDER_COMMIT)


EAHELITRON_COMMIT=c4c3dca
EAHELITRON_REPO=https://github.com/dontkme/EAHelitron.git
EAHELITRON_PREFIX_ARG=$(INSTALL_BASE)/eahelitron/$(EAHELITRON_COMMIT)


VSEARCH_VERSION=2.14.1
VSEARCH_URL=https://github.com/torognes/vsearch/archive/v2.14.1.tar.gz
VSEARCH_PREFIX_ARG=$(INSTALL_BASE)/vsearch/$(VSEARCH_VERSION)

DECIPHER_VERSION=2.10.0
DECIPHER_URL=http://www2.decipher.codes/DECIPHER_2.10.0.tar.gz
DECIPHER_PREFIX_ARG=$(INSTALL_BASE)/decipher/$(DECIPHER_VERSION)


# Downloads

RM_LIB_ARG=/data/rmlib

RMASK_META_TAR=RepeatMaskerMetaData-20181026.tar.gz
RMASK_META_URL=http://www.repeatmasker.org/libraries/RepeatMaskerMetaData-20181026.tar.gz

REPBASE_TAR=RepBaseRepeatMaskerEdition-20181026.tar.gz
REPBASE_URL=http://www.girinst.org/server/RepBase/protected/repeatmaskerlibraries/RepBaseRepeatMaskerEdition-20181026.tar.gz
REPBASE_USERNAME:=
REPBASE_PASSWORD:=

GENERAL=genometools python3 posix gffpal seqrenamer
NONCODING=aragorn structrnafinder trnascan
ALIGNMENT=rmblast vsearch decipher mmseqs
ANNOTATION=repeatmasker eahelitron mitefinder
PROPRIETARY=rnammer

OPEN=$(GENERAL) $(NONCODING) $(ALIGNMENT) $(ANNOTATION)

all: docker/all singularity/all

# Docker targets
docker: docker/pante
docker/all: docker/general docker/alignment docker/annotation docker/noncoding
	@echo "WARNING: proprietary software will not be built by default."
	@echo "WARNING: call `make docker/proprietary` explicitly to create them."

docker/general: $(addprefix docker/,$(GENERAL))
docker/alignment: $(addprefix docker/,$(ALIGNMENT))
docker/annotation: $(addprefix docker/,$(ANNOTATION))
docker/noncoding: $(addprefix docker/,$(NONCODING))

# Proprietary software won't be built by default.
docker/proprietary: $(addprefix docker/,$(PROPRIETARY))

# Singularity targets
singularity: singularity/pante.sif
singularity/all: singularity/general singularity/alignment singularity/annotation singularity/noncoding
	@echo "WARNING: proprietary software will not be built by default."
	@echo "WARNING: call `make singularity/proprietary` explicitly to create them."

singularity/general: $(addprefix singularity/,$(addsuffix .sif,$(GENERAL)))
singularity/alignment: $(addprefix singularity/,$(addsuffix .sif,$(ALIGNMENT)))
singularity/annotation: $(addprefix singularity/,$(addsuffix .sif,$(ANNOTATION)))
singularity/noncoding: $(addprefix singularity/,$(addsuffix .sif,$(NONCODING)))

# Proprietary software won't be built by default.
singularity/proprietary: $(addprefix singularity/,$(addsuffix .sif,$(PROPRIETARY)))

# Pull containers for dockerhub
dockerpull: dockerpull/pante
dockerpull/all: dockerpull/general dockerpull/alignment dockerpull/annotation dockerpull/noncoding
dockerpull/general: $(addprefix dockerpull/,$(GENERAL))
dockerpull/alignment: $(addprefix dockerpull/,$(ALIGNMENT))
dockerpull/annotation: $(addprefix dockerpull/,$(ANNOTATION))
dockerpull/noncoding: $(addprefix dockerpull/,$(NONCODING))


clean:
	docker stop $(shell docker ps -aq) || true
	docker rm $(shell docker ps -a -q) || true
	docker rmi $(shell docker images -q) || true
	singularity cache clean -a
	rm -rf -- tmp
	rm -rf -- docker

tidy:
	# Just removes untagged images to free up space
	rm -rf -- tmp
	docker image prune -f

# Generic rules
singularity/%.sif: docker/%
	@mkdir -p "$(dir $@)"
	@mkdir -p tmp tmp/cache
	SINGULARITY_TMPDIR="${PWD}/tmp" SINGULARITY_CACHEDIR="${PWD}/tmp/cache" singularity build \
	  --force \
	  --tmpdir "${PWD}/tmp" \
	  "$@" \
	  "docker-daemon://$(call make_image_name,$(notdir $(basename $@)))"

docker/%.tar.gz: docker/%
	@mkdir -p "$(dir $@)"
	docker save "$(call make_image_name,$(basename $(basename $@)))" | gzip > $@

dockerpull/%:
	docker pull "$(call make_image_name,$(notdir $@))"
	@mkdir -p "$(dir $@)"
	@touch docker/$(notdir $@)


docker/base: base.Dockerfile
	docker build \
	  --build-arg DEBIAN_VERSION=$(DEBIAN_VERSION) \
	  --build-arg MAINTAINER="$(MAINTAINER)" \
	  --build-arg PIPELINE_VERSION="$(TAG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/python3: python3.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/posix: posix.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/genometools: genometools.Dockerfile docker/base docker/hmmer3
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg GENOMETOOLS_VERSION="$(GENOMETOOLS_VERSION)" \
	  --build-arg GENOMETOOLS_URL="$(GENOMETOOLS_URL)" \
	  --build-arg GENOMETOOLS_PREFIX_ARG="$(GENOMETOOLS_PREFIX_ARG)" \
	  --build-arg HMMER3_IMAGE="$(call make_image_name,hmmer3)" \
	  --build-arg HMMER3_VERSION="$(HMMER3_VERSION)" \
	  --build-arg HMMER3_URL="$(HMMER3_URL)" \
	  --build-arg HMMER3_PREFIX_ARG="$(HMMER3_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/gffpal: gffpal.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg GFFPAL_TAG="$(GFFPAL_TAG)" \
	  --build-arg GFFPAL_REPO="$(GFFPAL_REPO)" \
	  --build-arg GFFPAL_PREFIX_ARG="$(GFFPAL_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/seqrenamer: seqrenamer.Dockerfile docker/base docker/gffpal
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg GFFPAL_IMAGE="$(call make_image_name,gffpal)" \
	  --build-arg SEQRENAMER_TAG="$(SEQRENAMER_TAG)" \
	  --build-arg SEQRENAMER_REPO="$(SEQRENAMER_REPO)" \
	  --build-arg SEQRENAMER_PREFIX_ARG="$(SEQRENAMER_PREFIX_ARG)" \
	  --build-arg GFFPAL_TAG="$(GFFPAL_TAG)" \
	  --build-arg GFFPAL_PREFIX_ARG="$(GFFPAL_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


# Repeat masking stuff


docker/mitefinder: mitefinder.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg MITEFINDER_COMMIT="$(MITEFINDER_COMMIT)" \
	  --build-arg MITEFINDER_REPO="$(MITEFINDER_REPO)" \
	  --build-arg MITEFINDER_PREFIX_ARG="$(MITEFINDER_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/eahelitron: eahelitron.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg EAHELITRON_COMMIT="$(EAHELITRON_COMMIT)" \
	  --build-arg EAHELITRON_REPO="$(EAHELITRON_REPO)" \
	  --build-arg EAHELITRON_PREFIX_ARG="$(EAHELITRON_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/vsearch: vsearch.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg VSEARCH_VERSION="$(VSEARCH_VERSION)" \
	  --build-arg VSEARCH_URL="$(VSEARCH_URL)" \
	  --build-arg VSEARCH_PREFIX_ARG="$(VSEARCH_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/decipher: decipher.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg DECIPHER_VERSION="$(DECIPHER_VERSION)" \
	  --build-arg DECIPHER_URL="$(DECIPHER_URL)" \
	  --build-arg DECIPHER_PREFIX_ARG="$(DECIPHER_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/rmblast: rmblast.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg RMBLAST_VERSION="$(RMBLAST_VERSION)" \
	  --build-arg RMBLAST_URL="$(RMBLAST_URL)" \
	  --build-arg RMBLAST_PATCH_URL="$(RMBLAST_PATCH_URL)" \
	  --build-arg RMBLAST_PREFIX_ARG="$(RMBLAST_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/repeatmasker: repeatmasker.Dockerfile docker/base docker/rmblast
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg RMBLAST_IMAGE="$(call make_image_name,rmblast)" \
	  --build-arg RMBLAST_VERSION="$(RMBLAST_VERSION)" \
	  --build-arg RMBLAST_PREFIX_ARG="$(RMBLAST_PREFIX_ARG)" \
	  --build-arg TRF_VERSION="$(TRF_VERSION)" \
	  --build-arg TRF_URL="$(TRF_URL)" \
	  --build-arg TRF_PREFIX_ARG="$(TRF_PREFIX_ARG)" \
	  --build-arg NSEG_URL="$(NSEG_URL)" \
	  --build-arg NSEG_PREFIX_ARG="$(NSEG_PREFIX_ARG)" \
	  --build-arg REPEATSCOUT_VERSION="$(REPEATSCOUT_VERSION)" \
	  --build-arg REPEATSCOUT_URL="$(REPEATSCOUT_URL)" \
	  --build-arg REPEATSCOUT_PREFIX_ARG="$(REPEATSCOUT_PREFIX_ARG)" \
	  --build-arg RECON_VERSION="$(RECON_VERSION)" \
	  --build-arg RECON_URL="$(RECON_URL)" \
	  --build-arg RECON_PREFIX_ARG="$(RECON_PREFIX_ARG)" \
	  --build-arg RMASK_VERSION="$(RMASK_VERSION)" \
	  --build-arg RMASK_URL="$(RMASK_URL)" \
	  --build-arg RMOD_VERSION="$(RMOD_VERSION)" \
	  --build-arg RMOD_URL="$(RMOD_URL)" \
	  --build-arg RMOD_PREFIX_ARG="$(RMOD_PREFIX_ARG)" \
	  --build-arg RM_LIB_ARG="$(RM_LIB_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/aragorn: aragorn.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg ARAGORN_VERSION="$(ARAGORN_VERSION)" \
	  --build-arg ARAGORN_URL="$(ARAGORN_URL)" \
	  --build-arg ARAGORN_PREFIX_ARG="$(ARAGORN_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/infernal: infernal.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg INFERNAL_VERSION="$(INFERNAL_VERSION)" \
	  --build-arg INFERNAL_URL="$(INFERNAL_URL)" \
	  --build-arg INFERNAL_PREFIX_ARG="$(INFERNAL_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/structrnafinder: structrnafinder.Dockerfile docker/base docker/infernal
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg STRUCTRNAFINDER_COMMIT="$(STRUCTRNAFINDER_COMMIT)" \
	  --build-arg STRUCTRNAFINDER_REPO="$(STRUCTRNAFINDER_REPO)" \
	  --build-arg STRUCTRNAFINDER_PREFIX_ARG="$(STRUCTRNAFINDER_PREFIX_ARG)" \
	  --build-arg VIENNA_VERSION="$(VIENNA_VERSION)" \
	  --build-arg VIENNA_URL="$(VIENNA_URL)" \
	  --build-arg VIENNA_PREFIX_ARG="$(VIENNA_PREFIX_ARG)" \
	  --build-arg INFERNAL_IMAGE="$(call make_image_name,infernal)" \
	  --build-arg INFERNAL_VERSION="$(INFERNAL_VERSION)" \
	  --build-arg INFERNAL_PREFIX_ARG="$(INFERNAL_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/trnascan: trnascan.Dockerfile docker/base docker/infernal
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg TRNASCAN_VERSION="$(TRNASCAN_VERSION)" \
	  --build-arg TRNASCAN_URL="$(TRNASCAN_URL)" \
	  --build-arg TRNASCAN_PREFIX_ARG="$(TRNASCAN_PREFIX_ARG)" \
	  --build-arg INFERNAL_IMAGE="$(call make_image_name,infernal)" \
	  --build-arg INFERNAL_VERSION="$(INFERNAL_VERSION)" \
	  --build-arg INFERNAL_PREFIX_ARG="$(INFERNAL_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/occultercut: occultercut.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg OCCULTERCUT_VERSION="$(OCCULTERCUT_VERSION)" \
	  --build-arg OCCULTERCUT_URL="$(OCCULTERCUT_URL)" \
	  --build-arg OCCULTERCUT_PREFIX_ARG="$(OCCULTERCUT_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/hmmer3: hmmer3.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg HMMER3_VERSION="$(HMMER3_VERSION)" \
	  --build-arg HMMER3_URL="$(HMMER3_URL)" \
	  --build-arg HMMER3_PREFIX_ARG="$(HMMER3_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/mmseqs: mmseqs.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg MMSEQS_TAG="$(MMSEQS_TAG)" \
	  --build-arg MMSEQS_REPO="$(MMSEQS_REPO)" \
	  --build-arg MMSEQS_PREFIX_ARG="$(MMSEQS_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/hmmer2: hmmer2.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg HMMER2_VERSION="$(HMMER2_VERSION)" \
	  --build-arg HMMER2_URL="$(HMMER2_URL)" \
	  --build-arg HMMER2_PREFIX_ARG="$(HMMER2_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@

docker/rnammer: rnammer.Dockerfile docker/base docker/hmmer2
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg RNAMMER_VERSION="$(RNAMMER_VERSION)" \
	  --build-arg RNAMMER_TAR="$(RNAMMER_TAR)" \
	  --build-arg RNAMMER_PREFIX_ARG="$(RNAMMER_PREFIX_ARG)" \
	  --build-arg HMMER2_IMAGE="$(call make_image_name,hmmer2)" \
	  --build-arg HMMER2_VERSION="$(HMMER2_VERSION)" \
	  --build-arg HMMER2_PREFIX_ARG="$(HMMER2_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@

docker/pante: pante.Dockerfile docker/python3 docker/mmseqs docker/decipher
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg PYTHON3_IMAGE="$(call make_image_name,python3)" \
	  --build-arg MMSEQS_IMAGE="$(call make_image_name,mmseqs)" \
	  --build-arg MMSEQS_TAG="$(MMSEQS_TAG)" \
	  --build-arg MMSEQS_PREFIX_ARG="$(MMSEQS_PREFIX_ARG)" \
	  --build-arg DECIPHER_IMAGE="$(call make_image_name,decipher)" \
	  --build-arg DECIPHER_VERSION="$(DECIPHER_VERSION)" \
	  --build-arg DECIPHER_PREFIX_ARG="$(DECIPHER_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@



downloads/Dfam.hmm.gz:
	mkdir -p downloads
	cd downloads \
	&& wget http://dfam.org/releases/Dfam_3.0/families/Dfam.hmm.gz


downloads/Dfam.embl.gz:
	mkdir -p downloads
	cd downloads \
	&& wget http://dfam.org/releases/Dfam_3.0/families/Dfam.embl.gz


downloads/$(RMASK_META_TAR):
	mkdir -p downloads
	cd downloads \
	&& wget -O $@ $(RMASK_META_URL)


downloads/$(REPBASE_TAR):
	mkdir -p downloads
	cd downloads \
	&& wget \
	     --user $(REPBASE_USERNAME) \
	     --password $(REPBASE_PASSWORD) \
	     -O $@ \
	     $(REPBASE_URL)


rmlib: downloads/Dfam.hmm.gz downloads/Dfam.embl.gz
rmlib: downloads/$(RMASK_META_TAR) downloads/$(REPBASE_TAR)
rmlib: docker/repeatmasker
	docker run \
	  --rm \
	  -v "$(PWD)/downloads:/downloads:rw" \
	  -v "$(PWD)/rmlib:$(RM_LIB_ARG):rw" \
	  -w /downloads \
	  $(call make_image_name,repeatmasker) \
	  prep_repeatmasker_lib.sh \
	    $(REPBASE_TAR) \
	    $(RMASK_META_TAR) \
	    Dfam.hmm.gz \
	    Dfam.embl.gz \
	    $(RM_LIB_ARG)
