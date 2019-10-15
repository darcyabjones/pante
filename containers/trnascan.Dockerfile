ARG IMAGE
ARG INFERNAL_IMAGE

FROM "${INFERNAL_IMAGE}" as infernal_builder

FROM "${IMAGE}" as trnascan_builder

ARG TRNASCAN_VERSION
ARG TRNASCAN_URL
ARG TRNASCAN_PREFIX_ARG
ENV TRNASCAN_PREFIX="${TRNASCAN_PREFIX_ARG}"

ARG INFERNAL_VERSION
ARG INFERNAL_PREFIX_ARG="/opt/infernal/${INFERNAL_VERSION}"
ENV INFERNAL_PREFIX="${INFERNAL_PREFIX_ARG}"

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
       build-essential \
       ca-certificates \
       wget \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates \
  && wget -O trnascan.tar.gz "${TRNASCAN_URL}" \
  && tar zxf trnascan.tar.gz \
  && cd tRNAscan-*/ \
  && ./configure --prefix="${TRNASCAN_PREFIX}" \
  && make \
  && sed -i "s~infernal_dir: {bin_dir}~infernal_dir: ${INFERNAL_PREFIX}/bin~" tRNAscan-SE.conf \
  && make install \
  && add_runtime_dep perl


FROM "${IMAGE}"

ARG TRNASCAN_VERSION
ARG TRNASCAN_PREFIX_ARG
ENV TRNASCAN_PREFIX="${TRNASCAN_PREFIX_ARG}"
LABEL trnascan.version="${TRNASCAN_VERSION}"

ENV PATH="${TRNASCAN_PREFIX}/bin:${PATH}"
ENV INCLUDE="${TRNASCAN_PREFIX}/include:${INCLUDE}"
ENV CPATH="${TRNASCAN_PREFIX}/include:${CPATH}"
ENV LIBRARY_PATH="${TRNASCAN_PREFIX}/lib:${LIBRARY_PATH}"
ENV LD_LIBRARY_PATH="${TRNASCAN_PREFIX}/lib:${LD_LIBRARY_PATH}"
ENV LD_RUN_PATH="${TRNASCAN_PREFIX}/lib:${LD_RUN_PATH}"

COPY --from=trnascan_builder "${TRNASCAN_PREFIX}" "${TRNASCAN_PREFIX}"
COPY --from=trnascan_builder "${APT_REQUIREMENTS_FILE}" /build/apt/trnascan.txt

ARG INFERNAL_VERSION
ARG INFERNAL_PREFIX_ARG="/opt/infernal/${INFERNAL_VERSION}"
ENV INFERNAL_PREFIX="${INFERNAL_PREFIX_ARG}"
LABEL infernal.version="${INFERNAL_VERSION}"

ENV PATH="${INFERNAL_PREFIX}/bin:${PATH}"

COPY --from=infernal_builder "${INFERNAL_PREFIX}" "${INFERNAL_PREFIX}"
COPY --from=infernal_builder "${APT_REQUIREMENTS_FILE}" /build/apt/infernal.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
