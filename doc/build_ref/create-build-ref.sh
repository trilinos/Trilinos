#!/bin/bash -e

# This simple script regenerates the (read-only) TrilinosBuildReference.*
# files from the Trilinos-specific template file:
#
#   TrilinosBuildQuickReferenceTemplate.rst
#
# and the standard TriBITS template file:
#
#   TribitsBuildQuickReferenceBody.rst
#
# The generated files TrilinosBuildReference.* (rst, html, and pdf) are
# read-only to avoid direct modification.  To update this output, one must
# modify the input template files described above.
#
# This script gets run from the base Trilinos source directory as:
#
#   cd Trilinos/doc/build_ref/
#   ./create-build-ref.sh
#  

echo
echo "Generating TrilinosBuildReference HTML and PDF files ..."
echo

make
