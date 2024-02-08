FROM mambaorg/micromamba:1.5.3

LABEL image.author.name "Kristina K. Gagalova"
LABEL image.author.email "kristina.gagalova@curtin.edu.au"

COPY base.sh /build/base.sh
USER 0

# install most of dep
COPY --chown=$MAMBA_USER:$MAMBA_USER env_yml/pante2_env.yml /tmp/pante2_env.yml

RUN micromamba create -n pante2

RUN micromamba install -y -n pante2 -f /tmp/pante2_env.yml && \
    micromamba clean --all --yes

# install new repeatmodeler
COPY --chown=$MAMBA_USER:$MAMBA_USER env_yml/repeatmodeler_env.yml /tmp/repeatmodeler_env.yml

RUN micromamba create -n repeatmodeler

RUN micromamba install -y -n repeatmodeler -f /tmp/repeatmodeler_env.yml && \
    micromamba clean --all --yes

ENV PATH="/opt/conda/envs/pante2/bin:/opt/conda/envs/repeatmodeler/bin:${PATH}"

ENV PATH="/opt/conda/envs/pante2/lib/python3.10/bin:${PATH}"
ENV PYTHONPATH="/opt/conda/envs/pante2/lib/python3.10/site-packages/:${PYTHONPATH}"

# independent tools

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
  && gawk -i inplace '/TRF_PRGM/{ n=NR+3 } NR==n{ sub(/'\''value'\'' => '\'''\''/, "'\''value'\'' => '\''/opt/conda/envs/pante2/bin/trf'\''") }1' RepeatMaskerConfig.pm \  
  && gawk -i inplace '/DEFAULT_SEARCH_ENGINE/{ n=NR+5 } NR==n{ sub(/'\''value'\'' => '\'''\''/, "'\''value'\'' => '\''rmblast'\''") }1' RepeatMaskerConfig.pm \
  && gawk -i inplace '/HMMER_DIR/{ n=NR+7 } NR==n{ sub(/'\''value'\'' => '\'''\''/, "'\''value'\'' => '\''/opt/conda/envs/pante2/bin'\''") }1' RepeatMaskerConfig.pm \
  && gawk -i inplace '/RMBLAST_DIR/{ n=NR+12 } NR==n{ sub(/'\''value'\'' => '\'''\''/, "'\''value'\'' => '\''/opt/conda/envs/pante2/bin'\''") }1' RepeatMaskerConfig.pm \ 
  && gawk -i inplace -v rmdir="${RMASK_PREFIX}" '/LIBDIR/{ n=NR+9 } NR==n{ sub(/'\''value'\'' => '\'''\''/, "'\''value'\'' => '\''" rmdir "/Libraries'\''") }1' RepeatMaskerConfig.pm \
  && rm repmasker.tar.gz

ENV RMASK_PREFIX="${RMASK_PREFIX}"
ENV PATH="${RMASK_PREFIX}:${RMASK_PREFIX}/util:${PATH}"
