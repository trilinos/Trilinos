#!/bin/bash -e

#  Build a generic version of <Project>BuildQuickRef for general
#  documentation.

ARGS=$@

./create-project-build-quickref.py \
--project-name="<Project>" \
--project-template-file=TribitsBuildQuickRefTemplate.rst \
--file-base=TribitsBuildQuickRef \
$ARGS
