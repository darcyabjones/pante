ARG IMAGE

FROM "${IMAGE}" as occultercut_builder

ARG OCCULTERCUT_VERSION
ARG OCCULTERCUT_URL
ARG OCCULTERCUT_PREFIX_ARG
ENV OCCULTERCUT_PREFIX="${OCCULTERCUT_PREFIX_ARG}"

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y \
      build-essential \
      wget \
  && rm -rf /var/lib/apt/lists/* \
  && wget -O occultercut.tar.gz "${OCCULTERCUT_URL}" \
  && tar -zxf occultercut.tar.gz \
  && cd OcculterCut* \
  && make -f makefile \
  && mkdir -p "${OCCULTERCUT_PREFIX}/bin" \
  && cp -r OcculterCut "${OCCULTERCUT_PREFIX}/bin" \
  && cd .. \
  && rm -rf -- OcculterCut* occultercut*

#add_runtime_dep libgomp1 python python-biopython gawk


FROM "${IMAGE}"
ARG OCCULTERCUT_VERSION
ARG OCCULTERCUT_PREFIX_ARG
ENV OCCULTERCUT_PREFIX="${OCCULTERCUT_PREFIX_ARG}"
LABEL occultercut.version="${OCCULTERCUT_VERSION}"

ENV PATH="${OCCULTERCUT_PREFIX}/bin:${PATH}"

COPY --from=occultercut_builder "${OCCULTERCUT_PREFIX}" "${OCCULTERCUT_PREFIX}"
COPY --from=occultercut_builder "${APT_REQUIREMENTS_FILE}" /build/apt/occultercut.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
