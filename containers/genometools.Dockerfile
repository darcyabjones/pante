ARG IMAGE
ARG HMMER3_IMAGE

FROM "${HMMER3_IMAGE}" as hmmer3_builder

FROM "${IMAGE}" as genometools_builder

ARG GENOMETOOLS_VERSION
ARG GENOMETOOLS_URL="http://genometools.org/pub/genometools-${GENOMETOOLS_VERSION}.tar.gz"
ARG GENOMETOOLS_PREFIX_ARG="/opt/genometools/${GENOMETOOLS_VERSION}"
ENV GENOMETOOLS_PREFIX="${GENOMETOOLS_PREFIX_ARG}"

ARG HMMER3_URL

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
  && wget -O genometools.tar.gz "${GENOMETOOLS_URL}"

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && tar -zxf genometools.tar.gz \
  && rm genometools.tar.gz \
  && cd genometools*/ \
  && sed -i 's/-Wall//g' Makefile \
  && cd src/external \
  && wget -O hmmer.tar.gz "${HMMER3_URL}" \
  && tar -zxf hmmer.tar.gz \
  && cd ../.. \
  && make errorcheck=no threads=yes with-hmmer=yes \
  && make errorcheck=no threads=yes with-hmmer=yes prefix="${GENOMETOOLS_PREFIX}" install \
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

COPY --from=genometools_builder "${GENOMETOOLS_PREFIX}" "${GENOMETOOLS_PREFIX}"
COPY --from=genometools_builder "${APT_REQUIREMENTS_FILE}" /build/apt/genometools.txt

ARG HMMER3_VERSION
ARG HMMER3_PREFIX_ARG
ENV HMMER3_PREFIX="${HMMER3_PREFIX_ARG}"

ENV PATH="${HMMER3_PREFIX}/bin:${PATH}"

LABEL hmmer3.version="${HMMER3_VERSION}"

COPY --from=hmmer3_builder "${HMMER3_PREFIX}" "${HMMER3_PREFIX}"
COPY --from=hmmer3_builder "${APT_REQUIREMENTS_FILE}" /build/apt/hmmer3.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
