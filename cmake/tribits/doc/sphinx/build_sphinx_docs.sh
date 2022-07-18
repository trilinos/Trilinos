#!/bin/bash -e
#
# Build all of the TriBITS-related Sphinx documentation
#
# To build this documentation, from any directory, run:
#
#   <this-dir>/build_sphinx_docs.sh
#
# WARNING: This has to remove any prior-created Sphinx-build documentation at
# the beginning using 'git clean -xdf -- .'  Therefore, before running this
# script make sure and commit your work (or modify previously tracked files)
# to make sure that you don't loose it!
#

ABS_FILE_PATH=`readlink -f $0`
SPHINX_DOC_BASE_DIR=`dirname $ABS_FILE_PATH`

cd $SPHINX_DOC_BASE_DIR

git clean -xdf -- .

python3 sphinx_rst_generator.py --copy-base-dir=../../..
