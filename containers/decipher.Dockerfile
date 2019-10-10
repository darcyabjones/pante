ARG IMAGE="darcyabjones/base"

FROM "${IMAGE}" as decipher_builder

ARG DECIPHER_VERSION
ARG DECIPHER_URL
ARG DECIPHER_PREFIX_ARG
ENV DECIPHER_PREFIX="${DECIPHER_PREFIX_ARG}"


WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
       build-essential \
       ca-certificates \
       wget \
  && add_runtime_dep \
       r-base \
       r-cran-ape \
       r-cran-phangorn \
       r-cran-rsqlite \
       r-bioc-biostrings \
       r-cran-optparse \
       r-cran-magrittr \
  && apt_install_from_file "${APT_REQUIREMENTS_FILE}" \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates \
  && wget -O "decipher.tar.gz" "${DECIPHER_URL}" \
  && tar -zxf decipher.tar.gz \
  && cd DECIPHER*/ \
  && R CMD build --no-build-vignettes --no-manual . \
  && mkdir -p "${DECIPHER_PREFIX}" \
  && R CMD INSTALL --build --library="${DECIPHER_PREFIX}" ./DECIPHER_*.tar.gz


FROM "${IMAGE}"

ARG DECIPHER_VERSION
ARG DECIPHER_PREFIX_ARG
ENV DECIPHER_PREFIX="${DECIPHER_PREFIX_ARG}"
LABEL decipher.version="${DECIPHER_VERSION}"

ENV R_LIBS_USER="${DECIPHER_PREFIX}:${R_LIBS_USER:-}"

COPY --from=decipher_builder "${DECIPHER_PREFIX}" "${DECIPHER_PREFIX}"
COPY --from=decipher_builder "${APT_REQUIREMENTS_FILE}" /build/apt/decipher.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
