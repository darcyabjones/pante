PWD=$(shell pwd)
DEBIAN_VERSION=stretch-20190228-slim
TAG=pante-v0.0.1
DOCKER_USER=darcyabjones
BASE_IMAGE=$(DOCKER_USER)/base:$(TAG)

INSTALL_BASE=/opt

define make_image_name
$(DOCKER_USER)/$(notdir $1):$(TAG)
endef

# SOFTWARE_VERSIONs
RED_TAG=v2.0
RED_REPO=https://github.com/TulsaBioinformaticsToolsmith/Red.git
RED_PREFIX_ARG=$(INSTALL_BASE)/red/$(RED_TAG)

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

RM_LIB_ARG=/data/rmlib

RMASK_META_TAR=RepeatMaskerMetaData-20181026.tar.gz
RMASK_META_URL=http://www.repeatmasker.org/libraries/RepeatMaskerMetaData-20181026.tar.gz

REPBASE_TAR=RepBaseRepeatMaskerEdition-20181026.tar.gz
REPBASE_URL=http://www.girinst.org/server/RepBase/protected/repeatmaskerlibraries/RepBaseRepeatMaskerEdition-20181026.tar.gz
REPBASE_USERNAME:=
REPBASE_PASSWORD:=


all: docker

docker: docker/repeatmasker

# Singularity is distributed as a monolith
# Target not implemented yet
# singularity: singularity/pante.sif


clean:
	docker stop $(shell docker ps -aq) || true
	docker rm $(shell docker ps -a -q) || true
	docker rmi $(shell docker images -q) || true
	rm -rf -- docker

tidy:
	# Just removes untagged images to free up space
	#docker rmi -f $(shell docker images -f "dangling=true" -q) || true
	docker image prune -f

# Generic rules
singularity/%.sif: docker/%
	@mkdir -p "$(dir $@)"
	singularity build "$@" "docker-daemon://$(call make_image_name,$(notdir $(basename $@)))"

docker/%.tar.gz: docker/%
	@mkdir -p "$(dir $@)"
	docker save "$(call make_image_name,$(basename $(basename $@)))" | gzip > $@


docker/base: base.Dockerfile
	docker build \
	  --build-arg DEBIAN_VERSION=$(DEBIAN_VERSION) \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/python3: python3.Dockerfile
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/posix: posix.Dockerfile
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


# Repeat masking stuff

docker/edge_base: base.Dockerfile
	docker build \
	  --build-arg DEBIAN_VERSION="buster-20190228-slim" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/red: red.Dockerfile docker/edge_base docker/base
	docker build \
	  --build-arg BUILD_IMAGE="$(call make_image_name,edge_base)" \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg RED_TAG="$(RED_TAG)" \
	  --build-arg RED_REPO="$(RED_REPO)" \
	  --build-arg RED_PREFIX_ARG="$(RED_PREFIX_ARG)" \
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
