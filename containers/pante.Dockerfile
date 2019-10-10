ARG IMAGE
ARG PYTHON3_IMAGE
ARG MMSEQS_IMAGE
ARG DECIPHER_IMAGE

FROM "${PYTHON3_IMAGE}" as python3_builder
FROM "${MMSEQS_IMAGE}" as mmseqs_builder
FROM "${DECIPHER_IMAGE}" as decipher_builder

FROM "${IMAGE}"


LABEL maintainer="darcy.ab.jones@gmail.com"

ARG MMSEQS_TAG
ARG MMSEQS_PREFIX_ARG
ENV MMSEQS_PREFIX="${MMSEQS_PREFIX_ARG}"
LABEL mmseqs.version="${MMSEQS_TAG}"

ENV PATH="${MMSEQS_PREFIX}/bin:${PATH}"

COPY --from=mmseqs_builder "${MMSEQS_PREFIX}" "${MMSEQS_PREFIX}"
COPY --from=mmseqs_builder "${APT_REQUIREMENTS_FILE}" /build/apt/mmseqs.txt


ARG DECIPHER_VERSION
ARG DECIPHER_PREFIX_ARG
ENV DECIPHER_PREFIX="${DECIPHER_PREFIX_ARG}"
LABEL decipher.version="${DECIPHER_VERSION}"

ENV R_LIBS_USER="${DECIPHER_PREFIX}:${R_LIBS_USER:-}"

COPY --from=decipher_builder "${DECIPHER_PREFIX}" "${DECIPHER_PREFIX}"
COPY --from=decipher_builder "${APT_REQUIREMENTS_FILE}" /build/apt/decipher.txt

COPY --from=python3_builder "${APT_REQUIREMENTS_FILE}" /build/apt/python3.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
