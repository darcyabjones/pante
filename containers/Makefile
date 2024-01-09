PWD=$(shell pwd)
TAG=v1.0.0
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

MITEFINDER_COMMIT=833754b
MITEFINDER_REPO=https://github.com/screamer/miteFinder.git
MITEFINDER_PREFIX_ARG=$(INSTALL_BASE)/mitefinder/$(MITEFINDER_COMMIT)

EAHELITRON_COMMIT=769ae07
EAHELITRON_REPO=https://github.com/dontkme/EAHelitron.git
EAHELITRON_PREFIX_ARG=$(INSTALL_BASE)/eahelitron/$(EAHELITRON_COMMIT)

# rnammer + hmm2
RNAMMER_VERSION=1.2
RNAMMER_PREFIX=$(INSTALL_BASE)/rnammer/$(RNAMMER_VERSION)

HMMER2_VERSION=2.3.2
HMMER2_URL=http://eddylab.org/software/hmmer/hmmer-2.3.2.tar.gz
HMMER2_PREFIX=$(INSTALL_BASE)/hmmer2/$(HMMER2_VERSION)

# Targets
# -----------------------------------------------------------
docker: docker/pante2

# Singularity targets
singularity: singularity/pante2.sif
singularity: singularity/pante2-rnammer.sif

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
	singularity build \
	  --force \
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
docker/pante2: pante2_env.Dockerfile env_yml/pante2_env.yml
	docker build \
	  --build-arg OCCULTERCUT_VERSION="$(OCCULTERCUT_VERSION)" \
	  --build-arg OCCULTERCUT_URL="$(OCCULTERCUT_URL)" \
	  --build-arg OCCULTERCUT_PREFIX_ARG="$(OCCULTERCUT_PREFIX_ARG)" \
	  --build-arg MITEFINDER_COMMIT="$(MITEFINDER_COMMIT)" \
	  --build-arg MITEFINDER_REPO="$(MITEFINDER_REPO)" \
	  --build-arg MITEFINDER_PREFIX_ARG="$(MITEFINDER_PREFIX_ARG)" \
	  --build-arg EAHELITRON_COMMIT="$(EAHELITRON_COMMIT)" \
	  --build-arg EAHELITRON_REPO="$(EAHELITRON_REPO)" \
	  --build-arg EAHELITRON_PREFIX_ARG="$(EAHELITRON_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/pante2-rnammer: pante2-rnammer_env.Dockerfile env_yml/pante2-rnammer_env.yml
	if [ -z $(RNAMMER_TAR) ]; then \
		echo "error RNAMMER_TAR is not set for pante2-rnammer. Please provide a value for RNAMMER_TAR. ex "sudo make docker/pante2-rnammer RNAMMER_TAR=/path/to/rnammertar")"; \
	else \
		docker build \
                    --build-arg OCCULTERCUT_VERSION="$(OCCULTERCUT_VERSION)" \
                    --build-arg OCCULTERCUT_URL="$(OCCULTERCUT_URL)" \
                    --build-arg OCCULTERCUT_PREFIX_ARG="$(OCCULTERCUT_PREFIX_ARG)" \
                    --build-arg MITEFINDER_COMMIT="$(MITEFINDER_COMMIT)" \
                    --build-arg MITEFINDER_REPO="$(MITEFINDER_REPO)" \
                    --build-arg MITEFINDER_PREFIX_ARG="$(MITEFINDER_PREFIX_ARG)" \
                    --build-arg EAHELITRON_COMMIT="$(EAHELITRON_COMMIT)" \
                    --build-arg EAHELITRON_REPO="$(EAHELITRON_REPO)" \
                    --build-arg EAHELITRON_PREFIX_ARG="$(EAHELITRON_PREFIX_ARG)" \
                    --build-arg RNAMMER_TAR="$(RNAMMER_TAR)"\
                    --build-arg RNAMMER_VERSION="$(RNAMMER_VERSION)" \
                    --build-arg RNAMMER_PREFIX="$(RNAMMER_PREFIX)" \
		    --build-arg HMMER2_VERSION="$(HMMER2_VERSION)" \
		    --build-arg HMMER2_URL="$(HMMER2_URL)" \
		    --build-arg HMMER2_PREFIX=$(INSTALL_BASE)/hmmer2/$(HMMER2_VERSION) \
                    --file $< \
                    -t "$(call make_image_name,$@)" \
                    . ; \
		mkdir -p "$(dir $@)" ; \
		touch $@; \
	fi
