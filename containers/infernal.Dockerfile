ARG IMAGE

FROM "${IMAGE}" as infernal_builder

ARG INFERNAL_VERSION
ARG INFERNAL_URL="http://eddylab.org/infernal/infernal-1.1.2.tar.gz"
ARG INFERNAL_PREFIX_ARG="/opt/infernal/${INFERNAL_VERSION}"
ENV INFERNAL_PREFIX="${INFERNAL_PREFIX_ARG}"

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
       build-essential \
       wget \
  && rm -rf /var/lib/apt/lists/* \
  && wget -O infernal.tar.gz "${INFERNAL_URL}" \
  && tar -xf infernal.tar.gz \
  && cd infernal*/ \
  && ./configure --prefix="${INFERNAL_PREFIX}" \
  && make \
  && make install


FROM "${IMAGE}"

ARG INFERNAL_VERSION
ARG INFERNAL_PREFIX_ARG
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
