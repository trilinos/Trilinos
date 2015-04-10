#!/bin/sh
find . -type f -name *.hpp -exec scripts/autoheader {} copyright.txt \;
find . -type f -name *.cpp -exec scripts/autoheader {} copyright.txt \;
find . -type f -name \*.am -exec scripts/autoheader {} copyright.txt \;
find . -type f -name \*.ac -exec scripts/autoheader {} copyright.txt \;
scripts/autoheader Makefile.export.phalanx.in copyright.txt
scripts/autoheader README copyright.txt
