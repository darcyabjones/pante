ARG IMAGE="darcyabjones/base"

FROM "${IMAGE}" as builder

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && add_runtime_dep bash curl gawk sed wget \
  && apt-get update \
  && apt_install_from_file "${APT_REQUIREMENTS_FILE}" \
  && rm -rf /var/lib/apt/lists/*
