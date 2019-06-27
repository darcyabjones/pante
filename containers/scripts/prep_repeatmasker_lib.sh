#!/usr/bin/env sh

set -eux

REPBASE_ARCHIVE="$1"
META_ARCHIVE="$2"
DFAM_HMM_GZ="$3"
DFAM_CONSENSUS_GZ="$4"
OUT_DIR="$5"

# These tar archives unpack to a folder "Libraries"
tar zxf "${REPBASE_ARCHIVE}"
tar zxf "${META_ARCHIVE}"

mv Libraries/* "${OUT_DIR}"
rmdir Libraries

# Repeat peps gets distributed with the repeat masker executables.
# Not available elsewhere.
cp ${RMASK_PREFIX}/RepeatPeps.lib ${OUT_DIR}
cp "${DFAM_HMM_GZ}" "${DFAM_CONSENSUS_GZ}" ${OUT_DIR}

cd ${OUT_DIR}

gunzip "${DFAM_HMM_GZ}"
# DFAM consensus filenames hard-coded into the script, so we move it.
mv "${DFAM_CONSENSUS_GZ}" "DfamConsensus.embl.gz"
gunzip "DfamConsensus.embl.gz"
gunzip "taxonomy.dat.gz"

# This basically concatenates the different blastable (i.e. not hmm) databases together.
perl -I ${RMASK_PREFIX} -e "use LibraryUtils; LibraryUtils::rebuildMainLibrary( \"${OUT_DIR}\" );"

buildRMLibFromEMBL.pl RepeatMaskerLib.embl > RepeatMasker.lib

makeblastdb -dbtype nucl -in RepeatMasker.lib
makeblastdb -dbtype prot -in RepeatPeps.lib

# As far as I can tell, nothing needs to be done with dfam hmms?
# Why not hmmpress?
exit 0
