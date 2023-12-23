FROM mambaorg/micromamba:1.5.3

LABEL image.author.name "Kristina K. Gagalova"
LABEL image.author.email "kristina.gagalova@curtin.edu.au"

ARG OCCULTERCUT_VERSION
ARG OCCULTERCUT_URL
ARG OCCULTERCUT_PREFIX_ARG
ENV OCCULTERCUT_PREFIX="${OCCULTERCUT_PREFIX_ARG}"

COPY base.sh /build/base.sh
USER 0

COPY --chown=$MAMBA_USER:$MAMBA_USER env_yml/pante_env.yml /tmp/pante_env.yml

RUN micromamba create -n pante

RUN micromamba install -y -n pante -f /tmp/pante_env.yml && \
    micromamba clean --all --yes

ENV PATH /opt/conda/envs/pante/bin:$PATH

WORKDIR /tmp
RUN  set -eu \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y \
       ca-certificates \
       build-essential \
       wget \
  && update-ca-certificates \
  && wget -O occultercut.tar.gz "${OCCULTERCUT_URL}" \
  && tar -zxf occultercut.tar.gz \
  && cd OcculterCut* \
  && make -f makefile \
  && mkdir -p "${OCCULTERCUT_PREFIX}/bin" \
  && cp -r OcculterCut "${OCCULTERCUT_PREFIX}/bin" \
  && cd .. \
  && rm -rf -- OcculterCut* occultercut*
