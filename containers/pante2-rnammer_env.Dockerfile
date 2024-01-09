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

# independent tools

ARG OCCULTERCUT_VERSION
ARG OCCULTERCUT_URL
ARG OCCULTERCUT_PREFIX_ARG
ENV OCCULTERCUT_PREFIX="${OCCULTERCUT_PREFIX_ARG}"

WORKDIR /tmp
RUN  set -eu \
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
  && update-ca-certificates \
  && wget -O occultercut.tar.gz "${OCCULTERCUT_URL}" \
  && tar -zxf occultercut.tar.gz \
  && cd OcculterCut* \
  && make -f makefile \
  && mkdir -p "${OCCULTERCUT_PREFIX}/bin" \
  && cp -r OcculterCut "${OCCULTERCUT_PREFIX}/bin" \
  && cd .. \
  && rm -rf -- OcculterCut* occultercut*

ENV PATH="${OCCULTERCUT_PREFIX}/bin:${PATH}"

ARG MITEFINDER_COMMIT
ARG MITEFINDER_REPO
ARG MITEFINDER_PREFIX_ARG
ENV MITEFINDER_PREFIX="${MITEFINDER_PREFIX_ARG}"

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
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
