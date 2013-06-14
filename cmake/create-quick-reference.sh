#!/bin/bash -e

# This simple script regenerates the TrilinosQuickstart.* files from the
# template file TrilinosQuickstartTemplate.rst and the standard
# TribitsQuickstartBody.rst file.
#
# This script gets run from the base Trilinos source directory as:
#
#   ./cmake/generate-quickstart.sh
#
# If the programs for rst2html and rst2latex need to be specified differently
# (such as rst2html.py and rst2latex.py) these can be set as:
#
#   env RST2HTML=rst2html.py RST2LATEX=rst2latex.py \
#        ./cmake/generate-quickstart.sh
#  

ARGS=""

if [ "${RST2HTML}" != "" ] ; then
  ARGS="$ARGS --generate-html=${RST2HTML}"
fi

if [ "${RST2LATEX}" != "" ] ; then
  ARGS="$ARGS --generate-latex=${RST2LATEX}"
fi
#echo "ARGS = $ARGS"

./cmake/tribits/doc/create-project-user-quickstart.py $ARGS
