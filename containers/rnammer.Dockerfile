ARG IMAGE="darcyabjones/base"
ARG HMMER2_IMAGE

FROM "${HMMER2_IMAGE}" as hmmer2_builder

FROM "${IMAGE}"

ARG RNAMMER_VERSION
ARG RNAMMER_PREFIX_ARG="/opt/rnammer/${RNAMMER_VERSION}"
ARG RNAMMER_TAR="sources/rnammer-1.2.src.tar.Z"
ENV RNAMMER_PREFIX="${RNAMMER_PREFIX_ARG}"
ENV PATH="${RNAMMER_PREFIX}:${PATH}"
LABEL rnammer.version="${RNAMMER_VERSION}"

ARG HMMER2_VERSION
ARG HMMER2_PREFIX_ARG="/opt/hmmer2/${HMMER2_VERSION}"
ENV HMMER2_PREFIX="${HMMER2_PREFIX_ARG}"
ENV HMMER_NCPU=1
ENV PATH="${HMMER2_PREFIX}/bin:${PATH}"
LABEL hmmer2.version="${HMMER2_VERSION}"

COPY --from=hmmer2_builder "${HMMER2_PREFIX}" "${HMMER2_PREFIX}"
COPY --from=hmmer2_builder "${APT_REQUIREMENTS_FILE}" /build/apt/hmmer2.txt

COPY "${RNAMMER_TAR}" "${RNAMMER_PREFIX}/rnammer.tar.Z"
WORKDIR ${RNAMMER_PREFIX}
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && add_runtime_dep \
       libxml-simple-perl\
       perl \
  && apt-get update \
  && apt_install_from_file "${APT_REQUIREMENTS_FILE}" \
  && tar xf rnammer.tar.Z \
  && rm rnammer.tar.Z \
  && rm -rf -- man example \
  && sed -i "s~/usr/cbs/bio/src/rnammer-${RNAMMER_VERSION}~${RNAMMER_PREFIX}~" rnammer \
  && sed -i "s~/usr/cbs/bio/bin/linux64/hmmsearch~${HMMER2_PREFIX}/bin/hmmsearch~" rnammer

WORKDIR /
