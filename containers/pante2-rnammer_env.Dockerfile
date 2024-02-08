FROM mambaorg/micromamba:1.5.3

LABEL image.author.name "Kristina K. Gagalova"
LABEL image.author.email "kristina.gagalova@curtin.edu.au"

COPY base.sh /build/base.sh
USER 0

# install most of dep
COPY --chown=$MAMBA_USER:$MAMBA_USER env_yml/pante2-rnammer_env.yml /tmp/pante2-rnammer_env.yml

RUN micromamba create -n pante2-rnammer

RUN micromamba install -y -n pante2-rnammer -f /tmp/pante2-rnammer_env.yml && \
    micromamba clean --all --yes

# install new repeatmodeler
COPY --chown=$MAMBA_USER:$MAMBA_USER env_yml/repeatmodeler_env.yml /tmp/repeatmodeler_env.yml

RUN micromamba create -n repeatmodeler

RUN micromamba install -y -n repeatmodeler -f /tmp/repeatmodeler_env.yml && \
    micromamba clean --all --yes

ENV PATH="/opt/conda/envs/pante2-rnammer/bin:/opt/conda/envs/repeatmodeler/bin:${PATH}"

ENV PATH="/opt/conda/envs/pante2-rnammer/lib/python3.10/bin:${PATH}"
ENV PYTHONPATH="/opt/conda/envs/pante2-rnammer/lib/python3.10/site-packages/:${PYTHONPATH}"

ENV RMASK_PREFIX="/opt/conda/envs/pante2-rnammer/share/RepeatMasker/"
ENV PATH="${RMASK_PREFIX}:${RMASK_PREFIX}/util:${PATH}"

ARG MITEFINDER_COMMIT
ARG MITEFINDER_REPO
ARG MITEFINDER_PREFIX_ARG
ENV MITEFINDER_PREFIX="${MITEFINDER_PREFIX_ARG}"

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y \
	ca-certificates \
	build-essential \
	procps \
	wget \
	libgl1 \
	gzip \
	zip \
	unzip \
  && git clone "${MITEFINDER_REPO}" \
  && cd miteFinder \
  && git checkout "${MITEFINDER_COMMIT}" \
  && make \
  && mkdir -p "${MITEFINDER_PREFIX}" \
  && mkdir -p "${MITEFINDER_PREFIX}/bin" \
  && mv miteFinder "${MITEFINDER_PREFIX}/bin" \
  && mv profile "${MITEFINDER_PREFIX}"


ENV MITEFINDER_PROFILE="${MITEFINDER_PREFIX}/profile/pattern_scoring.txt"
ENV PATH="${MITEFINDER_PREFIX}/bin:${PATH}"

ARG EAHELITRON_COMMIT
ARG EAHELITRON_REPO
ARG EAHELITRON_PREFIX_ARG
ENV EAHELITRON_PREFIX="${EAHELITRON_PREFIX_ARG}"

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
       ca-certificates \
  && update-ca-certificates \
  && mkdir -p "${EAHELITRON_PREFIX}/bin" \
  && cd "${EAHELITRON_PREFIX}/bin" \
  && git clone "${EAHELITRON_REPO}" \
  && cd EAHelitron \
  && git checkout "${EAHELITRON_COMMIT}" \
  && rm -rf -- .git README.md ChromInfo.txt \
  && mv LICENSE "${EAHELITRON_PREFIX}" \
  && chmod a+x * \
  && sed -i '1d' EAHelitron \
  && sed -i '1i #!/usr/bin/perl' EAHelitron

ENV PATH="${EAHELITRON_PREFIX}/bin/EAHelitron:${PATH}"

# rnammer + hmm2
ARG HMMER2_VERSION
ARG HMMER2_URL
ARG HMMER2_PREFIX
ENV HMMER_NCPU=1

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && rm -rf /var/lib/apt/lists/* \
  && wget -O hmmer.tar.gz "${HMMER2_URL}" \
  && tar xf hmmer.tar.gz \
  && rm hmmer.tar.gz \
  && cd hmmer-*/ \
  && CFLAGS="-g -O3" ./configure --prefix="${HMMER2_PREFIX}" --enable-threads \
  && make \
  && make install

ENV PATH="${HMMER2_PREFIX}/bin:${PATH}"

ARG RNAMMER_TAR
ARG RNAMMER_VERSION
ARG RNAMMER_PREFIX

COPY ${RNAMMER_TAR} /tmp/rnammer-${RNAMMER_VERSION}.tar.gz

WORKDIR ${RNAMMER_PREFIX}
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && tar xf /tmp/rnammer-${RNAMMER_VERSION}.tar.gz \
  && rm /tmp/rnammer-${RNAMMER_VERSION}.tar.gz \
  && rm -rf -- man example \
  && sed -i "s~/usr/cbs/bio/src/rnammer-${RNAMMER_VERSION}~${RNAMMER_PREFIX}~" rnammer \
  && sed -i "s~/usr/cbs/bio/bin/linux64/hmmsearch~${HMMER2_PREFIX}/bin/hmmsearch~" rnammer

ENV PATH="${RNAMMER_PREFIX}:${PATH}"

#---- new version of RepeatMasker 4.1.6 which has a new db structure (lighter!). Will wait until this gets on conda), in the meanwhile compiling manually

ARG RMASK_PREFIX
ARG RMASK_URL

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && update-ca-certificates \
  && mkdir -p "${RMASK_PREFIX}" \
  && cd "${RMASK_PREFIX}" \
  && wget -O repmasker.tar.gz "${RMASK_URL}" \
  && tar -zxf repmasker.tar.gz \
  && chmod a+w RepeatMasker/Libraries \
  && chmod a+w RepeatMasker/Libraries/famdb \
  && mv RepeatMasker RepeatMasker2 \
  && mv RepeatMasker2/* ./ && rm -r RepeatMasker2 \
  && cp RepeatMaskerConfig.pm RepeatMaskerConfig.pm.bak \
  && gawk -i inplace '/TRF_PRGM/{ n=NR+3 } NR==n{ sub(/'\''value'\'' => '\'''\''/, "'\''value'\'' => '\''/opt/conda/envs/pante2-rnammer/bin/trf'\''") }1' RepeatMaskerConfig.pm \
  && gawk -i inplace '/DEFAULT_SEARCH_ENGINE/{ n=NR+5 } NR==n{ sub(/'\''value'\'' => '\'''\''/, "'\''value'\'' => '\''rmblast'\''") }1' RepeatMaskerConfig.pm \
  && gawk -i inplace '/HMMER_DIR/{ n=NR+7 } NR==n{ sub(/'\''value'\'' => '\'''\''/, "'\''value'\'' => '\''/opt/conda/envs/pante2-rnammer/bin'\''") }1' RepeatMaskerConfig.pm \
  && gawk -i inplace '/RMBLAST_DIR/{ n=NR+12 } NR==n{ sub(/'\''value'\'' => '\'''\''/, "'\''value'\'' => '\''/opt/conda/envs/pante2-rnammer/bin'\''") }1' RepeatMaskerConfig.pm \
  && gawk -i inplace -v rmdir="${RMASK_PREFIX}" '/LIBDIR/{ n=NR+9 } NR==n{ sub(/'\''value'\'' => '\'''\''/, "'\''value'\'' => '\''" rmdir "/Libraries'\''") }1' RepeatMaskerConfig.pm \
  && rm repmasker.tar.gz

ENV RMASK_PREFIX="${RMASK_PREFIX}"
ENV PATH="${RMASK_PREFIX}:${RMASK_PREFIX}/util:${PATH}"
