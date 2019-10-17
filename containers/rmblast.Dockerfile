ARG IMAGE

FROM "${IMAGE}" as rmblast_builder

ARG RMBLAST_VERSION
ARG RMBLAST_URL=
ARG RMBLAST_PATCH_URL
ARG RMBLAST_PREFIX_ARG="/opt/rmblast/${RMBLAST_VERSION}"
ENV RMBLAST_PREFIX="${RMBLAST_PREFIX_ARG}"


WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
       build-essential \
       cpio \
       libblas-dev \
       libboost-dev \
       libboost-test-dev \
       libboost-thread-dev \
       libboost-regex-dev \
       libboost-iostreams-dev \
       libboost-program-options-dev \
       libboost-system-dev \
       libboost-test-dev \
       libboost-filesystem-dev \
       libbz2-dev \
       libcurl4-gnutls-dev \
       libgmp-dev \
       libgnutls28-dev \
       liblapack-dev \
       liblzma-dev \
       libperl-dev \
       libxml2-dev \
       libxslt1-dev \
       zlib1g-dev \
       wget \
  && rm -rf /var/lib/apt/lists/* \
  && wget "${RMBLAST_URL}" \
  && wget "${RMBLAST_PATCH_URL}" \
  && tar zxf "ncbi-blast-${RMBLAST_VERSION}-src.tar.gz" \
  && rm "ncbi-blast-${RMBLAST_VERSION}-src.tar.gz" \
  && gunzip "isb-${RMBLAST_VERSION}-rmblast.patch.gz" \
  && cd "ncbi-blast-${RMBLAST_VERSION}-src" \
  && patch -p1 <  "../isb-${RMBLAST_VERSION}-rmblast.patch" \
  && cd "c++" \
  && ./configure \
       --with-mt \
       --with-openmp \
       --without-debug \
       --without-krb5 \
       --without-openssl \
       --with-projects=scripts/projects/rmblastn/project.lst \
       --prefix="${RMBLAST_PREFIX}" \
  && make -j $(grep -c ^processor /proc/cpuinfo) \
  && make install \
  && cd .. \
  && rm -rf -- "c++" \
  && rm ${RMBLAST_PREFIX}/bin/*_unit_test \
  && rm ${RMBLAST_PREFIX}/bin/test_* \
  && rm ${RMBLAST_PREFIX}/bin/lmdb* \
  && rm ${RMBLAST_PREFIX}/bin/seqdb* \
  && add_runtime_dep \
       libblas3 \
       libbz2-1.0 \
       libcurl3-gnutls \
       libgmp10 \
       libgnutlsxx28 \
       libgomp1 \
       liblapack3 \
       libstdc++6 \
       libxml2 \
       libxslt1.1 \
       lzma \
       perl \
       zlib1g


FROM "${IMAGE}"

ARG RMBLAST_VERSION
ARG RMBLAST_PREFIX_ARG
ENV RMBLAST_PREFIX="${RMBLAST_PREFIX_ARG}"

ENV PATH="${RMBLAST_PREFIX}/bin:${PATH}"
ENV CPATH="${RMBLAST_PREFIX}/include"
ENV LIBRARY_PATH="${RMBLAST_PREFIX}/lib:${LIBRARY_PATH}"
ENV LD_LIBRARY_PATH="${RMBLAST_PREFIX}/lib:${LD_LIBRARY_PATH}"

LABEL rmblast.version="${RMBLAST_VERSION}"

COPY --from=rmblast_builder "${RMBLAST_PREFIX}" "${RMBLAST_PREFIX}"
COPY --from=rmblast_builder "${APT_REQUIREMENTS_FILE}" /build/apt/rmblast.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
