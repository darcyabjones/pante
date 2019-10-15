ARG IMAGE

FROM "${IMAGE}" as eahelitron_builder

ARG EAHELITRON_COMMIT
ARG EAHELITRON_REPO
ARG EAHELITRON_PREFIX_ARG
ENV EAHELITRON_PREFIX="${EAHELITRON_PREFIX_ARG}"
LABEL eahelitron.version="${EAHELITRON_COMMIT}"


WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
       build-essential \
       ca-certificates \
       git \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates \
  && mkdir -p "${EAHELITRON_PREFIX}/bin" \
  && cd "${EAHELITRON_PREFIX}/bin" \
  && git clone "${EAHELITRON_REPO}" . \
  && git checkout "${EAHELITRON_COMMIT}" \
  && rm -rf -- .git README.md ChromInfo.txt \
  && mv LICENSE "${EAHELITRON_PREFIX}" \
  && chmod a+x * \
  && sed -i '1d' EAHelitron \
  && sed -i '1i #!/usr/bin/perl' EAHelitron \
  && add_runtime_dep \
       libparallel-forkmanager-perl \
       perl

ENV PATH="${EAHELITRON_PREFIX}/bin:${PATH}"

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && mkdir -p /build/apt \
  && cp "${APT_REQUIREMENTS_FILE}" /build/apt/eahelitron.txt \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
