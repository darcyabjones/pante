#!/usr/bin/env bash

set -euo pipefail

if [ $# -eq 1 ]
then
  INFILE="/dev/stdin"
else
  INFILE="$2"
fi

awk -v name="$1" '
  /^>/ {
    seqid=gensub(/^>.*\[seqid ([^:]+):.*$/, "\\1", "g", $0);
    start=gensub(/^>.*\[seqid [^:]+:([0-9]+).*$/, "\\1", "g", $0);
    end=gensub(/^>.*\[seqid [^:]+:[0-9]+-([0-9]+).*$/, "\\1", "g", $0);
    strand=gensub(/^>.*(.)\]$/, "\\1", "g", $0);
    if ( strand == "-" ) sid=name ":" seqid ":" end "-" start;
    else sid=name ":" seqid ":" start "-" end;
    $0 = ">" sid;
  }
  { print }
' < "${INFILE}"
