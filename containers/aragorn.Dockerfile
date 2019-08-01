ARG IMAGE

FROM "${IMAGE}" as builder

ARG ARAGORN_VERSION
ARG ARAGORN_URL="http://130.235.244.92/ARAGORN/Downloads/aragorn1.2.38.tgz"
ARG ARAGORN_PREFIX_ARG="/opt/aragorn/${ARAGORN_VERSION}"
ENV ARAGORN_PREFIX="${ARAGORN_PREFIX_ARG}"

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y \
       build-essential \
       wget \
  && rm -rf /var/lib/apt/lists/* \
  && wget -O aragorn.tar.gz "${ARAGORN_URL}" \
  && tar -xf aragorn.tar.gz \
  && cd aragorn*/ \
  && gcc -O3 -ffast-math -finline-functions -o aragorn aragorn*.c \
  && mkdir -p "${ARAGORN_PREFIX}/bin" \
  && mv aragorn "${ARAGORN_PREFIX}/bin"


FROM "${IMAGE}"

ARG ARAGORN_VERSION
ARG ARAGORN_PREFIX_ARG="/opt/aragorn/${ARAGORN_VERSION}"
ENV ARAGORN_PREFIX="${ARAGORN_PREFIX_ARG}"
LABEL aragorn.version="${ARAGORN_VERSION}"

ENV PATH="${ARAGORN_PREFIX}/bin:${PATH}"

COPY --from=builder "${ARAGORN_PREFIX}" "${ARAGORN_PREFIX}"
COPY --from=builder "${APT_REQUIREMENTS_FILE}" /build/apt/aragorn.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
