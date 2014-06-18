#!/bin/bash -e

# This simple script regenerates the (read-only) TrilinosBuildQuickReference.*
# files from the Trilinos-specific template file:
#
#   cmake/TrilinosBuildQuickReferenceTemplate.rst
#
# and the standard TriBITS template file:
#
#   cmake/tribits/doc/TribitsBuildQuickReferenceBody.rst file.
#
# This generated files TrilinosBuildQuickReference.* (rst, html, and pdf) are
# read-only to avoid direct modification.  To update this output, one must
# modify the input template files described above.
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

./cmake/tribits/doc/build_quick_ref/create-project-build-quickref.py $ARGS
