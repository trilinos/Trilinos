#!/bin/bash -e

# This script is used to generate the the index.html file.

ARGS=$@

echo
echo "Generating HTML ..."
echo
../python/generate-docutils-output.py \
  --file-base=index \
  --generate-latex= \
  $ARGS
