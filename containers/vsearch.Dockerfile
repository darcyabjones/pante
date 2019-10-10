ARG IMAGE

FROM "${IMAGE}" as vsearch_builder

ARG VSEARCH_VERSION
ARG VSEARCH_URL
ARG VSEARCH_PREFIX_ARG
ENV VSEARCH_PREFIX="${VSEARCH_PREFIX_ARG}"

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
       autoconf \
       automake \
       build-essential \
       ca-certificates \
       libbz2-dev \
       wget \
       zlib1g-dev \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates \
  && wget -O vsearch.tar.gz "${VSEARCH_URL}" \
  && tar -zxf vsearch.tar.gz \
  && cd vsearch*/ \
  && ./autogen.sh \
  && ./configure --prefix="${VSEARCH_PREFIX}" --disable-pdfman \
  && make \
  && make install \
  && add_runtime_dep \
       libbz2-1.0 \
       zlib1g


FROM "${IMAGE}"

ARG VSEARCH_VERSION
ARG VSEARCH_PREFIX_ARG
ENV VSEARCH_PREFIX="${VSEARCH_PREFIX_ARG}"

ENV PATH="${VSEARCH_PREFIX}/bin:${PATH}"

LABEL vsearch.version="${VSEARCH_VERSION}"

COPY --from=vsearch_builder "${VSEARCH_PREFIX}" "${VSEARCH_PREFIX}"
COPY --from=vsearch_builder "${APT_REQUIREMENTS_FILE}" /build/apt/vsearch.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
