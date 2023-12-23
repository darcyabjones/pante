PWD=$(shell pwd)
TAG=pante2-v1.0.0
DOCKER_USER=kristinagagalova
BASE_IMAGE=$(DOCKER_USER)/base:$(TAG)

INSTALL_BASE=/opt

define make_image_name
$(DOCKER_USER)/$(notdir $1):$(TAG)
endef

# Other tools
# -----------------------------------------------------------
OCCULTERCUT_VERSION=1.1
OCCULTERCUT_URL=https://downloads.sourceforge.net/project/occultercut/OcculterCut_v1.1.tar.gz
OCCULTERCUT_PREFIX_ARG=$(INSTALL_BASE)/occultercut/$(OCCULTERCUT_VERSION)

# Targets
# -----------------------------------------------------------
all: docker/all singularity/all

docker: docker/pante
docker/all: docker/general docker/alignment docker/annotation docker/noncoding
	@echo "WARNING: proprietary software will not be built by default."
	@echo "WARNING: call 'make docker/proprietary' explicitly to create them."

# Proprietary software won't be built by default.
docker/proprietary: $(addprefix docker/,$(PROPRIETARY))

# Singularity targets
singularity: singularity/pante.sif
singularity/all: singularity/general singularity/alignment singularity/annotation singularity/noncoding
	@echo "WARNING: proprietary software will not be built by default."
	@echo "WARNING: call 'make singularity/proprietary' explicitly to create them."

# Proprietary software won't be built by default.
singularity/proprietary: $(addprefix singularity/,$(addsuffix .sif,$(PROPRIETARY)))

# Pull containers for dockerhub
dockerpull: dockerpull/pante

clean:
	docker stop $(shell docker ps -aq) || true
	docker rm $(shell docker ps -a -q) || true
	docker rmi $(shell docker images -q) || true
	singularity cache clean -f
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


# Basic tools


# THE MONOLITH(s)

docker/pante: pante_env.Dockerfile env_yml/pante_env.yml
	docker build \
	  --build-arg OCCULTERCUT_VERSION="$(OCCULTERCUT_VERSION)" \
	  --build-arg OCCULTERCUT_URL="$(OCCULTERCUT_URL)" \
	  --build-arg OCCULTERCUT_PREFIX_ARG="$(OCCULTERCUT_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/pante-plus: pante-plus.Dockerfile docker/pante docker/hmmer2 docker/rnammer
	docker build \
	  --build-arg PANTE_IMAGE="$(call make_image_name,pante)" \
	  --build-arg HMMER2_IMAGE="$(call make_image_name,hmmer2)" \
	  --build-arg HMMER2_VERSION="$(HMMER2_VERSION)" \
	  --build-arg HMMER2_PREFIX_ARG="$(HMMER2_PREFIX_ARG)" \
	  --build-arg RNAMMER_IMAGE="$(call make_image_name,rnammer)" \
	  --build-arg RNAMMER_VERSION="$(RNAMMER_VERSION)" \
	  --build-arg RNAMMER_PREFIX_ARG="$(RNAMMER_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@
