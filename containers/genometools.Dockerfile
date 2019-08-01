ARG IMAGE

FROM "${IMAGE}" as builder

ARG GENOMETOOLS_VERSION
ARG GENOMETOOLS_URL="http://genometools.org/pub/genometools-${GENOMETOOLS_VERSION}.tar.gz"
ARG GENOMETOOLS_PREFIX_ARG="/opt/genometools/${GENOMETOOLS_VERSION}"
ENV GENOMETOOLS_PREFIX="${GENOMETOOLS_PREFIX_ARG}"

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y \
       build-essential \
       libcairo2-dev \
       libpango1.0-dev \
       wget \
  && rm -rf /var/lib/apt/lists/* \
  && wget -O genometools.tar.gz "${GENOMETOOLS_URL}" \
  && tar zxf genometools.tar.gz \
  && rm genometools.tar.gz \
  && cd genometools*/ \
  && sed -i 's/-Wall//g' Makefile \
  && make errorcheck=no \
  && make errorcheck=no prefix="${GENOMETOOLS_PREFIX}" install \
  && add_runtime_dep libcairo2 libpango-1.0-0 libpangocairo-1.0-0


FROM "${IMAGE}"

ARG GENOMETOOLS_VERSION
ARG GENOMETOOLS_PREFIX_ARG="/opt/genometools/${GENOMETOOLS_VERSION}"
ENV GENOMETOOLS_PREFIX="${GENOMETOOLS_PREFIX_ARG}"
LABEL genometools.version="${GENOMETOOLS_VERSION}"

ENV PATH="${GENOMETOOLS_PREFIX}/bin:${PATH}"
ENV INCLUDE="${GENOMETOOLS_PREFIX}/include:${INCLUDE}"
ENV CPATH="${GENOMETOOLS_PREFIX}/include:${CPATH}"
ENV LIBRARY_PATH="${GENOMETOOLS_PREFIX}/lib:${LIBRARY_PATH}"
ENV LD_LIBRARY_PATH="${GENOMETOOLS_PREFIX}/lib:${LD_LIBRARY_PATH}"
ENV LD_RUN_PATH="${GENOMETOOLS_PREFIX}/lib:${LD_RUN_PATH}"

COPY --from=builder "${GENOMETOOLS_PREFIX}" "${GENOMETOOLS_PREFIX}"
COPY --from=builder "${APT_REQUIREMENTS_FILE}" /build/apt/genometools.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
