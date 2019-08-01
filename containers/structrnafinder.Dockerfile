ARG IMAGE
ARG INFERNAL_IMAGE

FROM "${INFERNAL_IMAGE}" as infernal_builder

FROM "${IMAGE}" as vienna_builder

ARG VIENNA_VERSION
ARG VIENNA_URL="https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.13.tar.gz"
ARG VIENNA_PREFIX_ARG="/opt/vienna/${VIENNA_VERSION}"
ENV VIENNA_PREFIX="${VIENNA_PREFIX_ARG}"

# Have to disable simd because of dynamic libraries, harder to automatically select.

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
       build-essential \
       ca-certificates \
       libgsl-dev \
       libmpfr-dev \
       wget \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates \
  && wget -O vienna.tar.gz "${VIENNA_URL}" \
  && tar -xf vienna.tar.gz \
  && cd ViennaRNA*/ \
  && ./configure --disable-simd --prefix="${VIENNA_PREFIX}" \
  && make \
  && make install \
  && add_runtime_dep libgomp1 libgsl23 libmpfr6


FROM "${IMAGE}" as structrnafinder_builder

ARG STRUCTRNAFINDER_COMMIT
ARG STRUCTRNAFINDER_REPO="https://github.com/viniciusmaracaja/structRNAfinder.git"
ARG STRUCTRNAFINDER_PREFIX_ARG="/opt/structrnafinder/${STRUCTRNAFINDER_COMMIT}"
ENV STRUCTRNAFINDER_PREFIX="${STRUCTRNAFINDER_PREFIX_ARG}"

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
       build-essential \
       ca-certificates \
       git \
       libgd-dev \
       libbio-graphics-perl \
       zlib1g-dev \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates \
  && git clone "${STRUCTRNAFINDER_REPO}" . \
  && git checkout "${STRUCTRNAFINDER_COMMIT}" \
  && sed -i "s~^open (IN2,\"/usr/local~open (IN2,\"${STRUCTRNAFINDER_PREFIX}~" bin/SRF_infernal2table \
  && mkdir -p "${STRUCTRNAFINDER_PREFIX}" \
  && cp -r bin "${STRUCTRNAFINDER_PREFIX}" \
  && cp -r share "${STRUCTRNAFINDER_PREFIX}" \
  && add_runtime_dep libgd-perl libbio-graphics-perl perl

# CA cert stuff sometime required for git clone https


FROM "${IMAGE}"

ARG STRUCTRNAFINDER_COMMIT
ARG STRUCTRNAFINDER_PREFIX_ARG="/opt/structrnafinder/${STRUCTRNAFINDER_COMMIT}"
ENV STRUCTRNAFINDER_PREFIX="${STRUCTRNAFINDER_PREFIX_ARG}"
LABEL structrnafinder.version="${STRUCTRNAFINDER_COMMIT}"

ENV PATH="${STRUCTRNAFINDER_PREFIX}/bin:${PATH}"

COPY --from=structrnafinder_builder "${STRUCTRNAFINDER_PREFIX}" "${STRUCTRNAFINDER_PREFIX}"
COPY --from=structrnafinder_builder "${APT_REQUIREMENTS_FILE}" /build/apt/structrnafinder.txt

ARG VIENNA_VERSION
ARG VIENNA_PREFIX_ARG="/opt/vienna/${VIENNA_VERSION}"
ENV VIENNA_PREFIX="${VIENNA_PREFIX_ARG}"
LABEL vienna.version="${VIENNA_VERSION}"

ENV PATH="${VIENNA_PREFIX}/bin:${PATH}"
ENV INCLUDE="${VIENNA_PREFIX}/include:${INCLUDE}"
ENV CPATH="${VIENNA_PREFIX}/include:${CPATH}"
ENV LIBRARY_PATH="${VIENNA_PREFIX}/lib:${LIBRARY_PATH}"
ENV LD_LIBRARY_PATH="${VIENNA_PREFIX}/lib:${LD_LIBRARY_PATH}"
ENV LD_RUN_PATH="${VIENNA_PREFIX}/lib:${LD_RUN_PATH}"

COPY --from=vienna_builder "${VIENNA_PREFIX}" "${VIENNA_PREFIX}"
COPY --from=vienna_builder "${APT_REQUIREMENTS_FILE}" /build/apt/vienna.txt

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
