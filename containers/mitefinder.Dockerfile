ARG IMAGE

FROM "${IMAGE}" as mitefinder_builder

ARG MITEFINDER_COMMIT="833754b0ff1899e8cb0f260d6d5011d3583b3012"
ARG MITEFINDER_REPO="https://github.com/screamer/miteFinder.git"
ARG MITEFINDER_PREFIX_ARG
ENV MITEFINDER_PREFIX="${MITEFINDER_PREFIX_ARG}"

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y \
       build-essential \
       git \
  && rm -rf /var/lib/apt/lists/* \
  && git clone "${MITEFINDER_REPO}" . \
  && git checkout "${MITEFINDER_COMMIT}" \
  && make \
  && mkdir -p "${MITEFINDER_PREFIX}" \
  && mkdir -p "${MITEFINDER_PREFIX}/bin" \
  && mv miteFinder "${MITEFINDER_PREFIX}/bin" \
  && mv profile "${MITEFINDER_PREFIX}"


FROM "${IMAGE}"

ARG MITEFINDER_COMMIT="833754b0ff1899e8cb0f260d6d5011d3583b3012"
ARG MITEFINDER_PREFIX_ARG
ENV MITEFINDER_PREFIX="${MITEFINDER_PREFIX_ARG}"

ENV MITEFINDER_PROFILE="${MITEFINDER_PREFIX}/profile/pattern_scoring.txt"
ENV PATH="${MITEFINDER_PREFIX}/bin:${PATH}"

LABEL mitefinder.version = "${MITEFINDER_COMMIT}"

COPY --from=mitefinder_builder "${MITEFINDER_PREFIX}" "${MITEFINDER_PREFIX}"
COPY --from=mitefinder_builder "${APT_REQUIREMENTS_FILE}" /build/apt/mitefinder.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
