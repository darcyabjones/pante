ARG IMAGE

FROM "${IMAGE}" as bedtools_builder

ARG BEDTOOLS_VERSION
ARG BEDTOOLS_URL
ARG BEDTOOLS_PREFIX_ARG
ENV BEDTOOLS_PREFIX="${BEDTOOLS_PREFIX_ARG}"

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
       build-essential \
       ca-certificates \
       libbz2-dev \
       liblzma-dev \
       python \
       wget \
       zlib1g-dev \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates \
  && wget -O bedtools.tar.gz "${BEDTOOLS_URL}" \
  && tar zxf bedtools.tar.gz \
  && cd bedtools*/ \
  && make \
  && make install prefix="${BEDTOOLS_PREFIX}" \
  && add_runtime_dep zlib1g libbz2-1.0 lzma
# CA cert stuff sometime required for git clone https


FROM "${IMAGE}"

ARG BEDTOOLS_VERSION
ARG BEDTOOLS_PREFIX_ARG
ENV BEDTOOLS_PREFIX="${BEDTOOLS_PREFIX_ARG}"
LABEL bedtools.version="${BEDTOOLS_VERSION}"

ENV PATH "${BEDTOOLS_PREFIX}/bin:${PATH}"

COPY --from=bedtools_builder "${BEDTOOLS_PREFIX}" "${BEDTOOLS_PREFIX}"
COPY --from=bedtools_builder "${APT_REQUIREMENTS_FILE}" /build/apt/bedtools.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
