#!/bin/bash

# This script is used to generate the the TribitsDevelopersGuide.(html,pdf)
# files using a new script in TriBITS for that purpose.  You just run it from
# this directory as:
#
#   cd <this-dir>
#   ./generate-dev-guide.sh
#
# If you want to override what utils are used to generate the files, you can
# pass in, for exmaple:
#
#   ./generate-dev-guide.sh \
#      --generate-html=rst2html.py --generate-latex=rst2latex.py
#
# Note that if you are debugging the parsing, you can disable the generation
# of the latex and pdf by setting:
#
#    --generate-latex=
#
# Enjoy!

ARGS=$@

../../python/generate-docutils-output.py \
  --file-base=TribitsDevelopersGuide \
  --generate-html=rst2html.py --generate-latex=rst2latex.py \
  $ARGS

# NOTE: This above invocation by default overrides the DocUtils html and latex
# generators to use the versions without the *.py prefix.  When you install
# docutils from source, you get the *.py prefix.  When you install using the
# yum package on a Linux machine, you don't get the *.py prefix.  What a pain!
