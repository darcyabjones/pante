ARG PANTE_IMAGE
ARG HMMER2_IMAGE
ARG RNAMMER_IMAGE

FROM "${HMMER2_IMAGE}" as hmmer2_builder
FROM "${RNAMMER_IMAGE}" as rnammer_builder


FROM "${PANTE_IMAGE}"

ARG HMMER2_VERSION
ARG HMMER2_PREFIX_ARG
ENV HMMER2_PREFIX="${HMMER2_PREFIX_ARG}"
ENV HMMER_NCPU=1
ENV PATH="${HMMER2_PREFIX}/bin:${PATH}"

LABEL hmmer2.version="${HMMER2_VERSION}"

COPY --from=hmmer2_builder "${HMMER2_PREFIX}" "${HMMER2_PREFIX}"
COPY --from=hmmer2_builder "${APT_REQUIREMENTS_FILE}" /build/apt/hmmer2.txt


ARG RNAMMER_VERSION
ARG RNAMMER_PREFIX_ARG
ENV RNAMMER_PREFIX="${RNAMMER_PREFIX_ARG}"
ENV PATH="${RNAMMER_PREFIX}:${PATH}"

LABEL rnammer.version="${RNAMMER_VERSION}"

COPY --from=rnammer_builder "${RNAMMER_PREFIX}" "${RNAMMER_PREFIX}"
COPY --from=rnammer_builder "${APT_REQUIREMENTS_FILE}" /build/apt/rnammer.txt


RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"

WORKDIR /
