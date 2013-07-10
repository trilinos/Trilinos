#!/bin/bash -e

# This simple script regenerates the TrilinosBuildQuickReference.* files from
# the template file TrilinosBuildQuickReferenceTemplate.rst and the standard
# TribitsBuildQuickReferenceBody.rst file.
#
# This script gets run from the base Trilinos source directory as:
#
#   ./cmake/create-build-quickref.sh
#
# if any options need to be overridden (for example to point to different
# versions of rst2thml or rst2latex), they can be passed into this scirpt, for
# example, as:
#
#   ./cmake/create-build-quickref.sh --generate-html=rst2html \
#     --generate-latex=rst2latex
#  

ARGS=$@

./cmake/tribits/doc/create-project-build-quickref.py $ARGS
