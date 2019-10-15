ARG IMAGE

FROM "${IMAGE}" as hmmer2_builder

ARG HMMER2_VERSION
ARG HMMER2_URL
ARG HMMER2_PREFIX_ARG
ENV HMMER2_PREFIX="${HMMER2_PREFIX_ARG}"
ENV HMMER_NCPU=1

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
       build-essential \
       wget \
  && rm -rf /var/lib/apt/lists/* \
  && wget -O hmmer.tar.gz "${HMMER2_URL}" \
  && tar xf hmmer.tar.gz \
  && rm hmmer.tar.gz \
  && cd hmmer-*/ \
  && CFLAGS="-g -O3" ./configure --prefix="${HMMER2_PREFIX}" --enable-threads \
  && make \
  && make install


FROM "${IMAGE}"

ARG HMMER2_VERSION
ARG HMMER2_PREFIX_ARG
ENV HMMER2_PREFIX="${HMMER2_PREFIX_ARG}"
ENV HMMER_NCPU=1
ENV PATH="${HMMER2_PREFIX}/bin:${PATH}"

LABEL hmmer2.version="${HMMER2_VERSION}"

COPY --from=hmmer2_builder "${HMMER2_PREFIX}" "${HMMER2_PREFIX}"
COPY --from=hmmer2_builder "${APT_REQUIREMENTS_FILE}" /build/apt/hmmer2.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
