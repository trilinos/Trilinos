#!/bin/sh
find . -type f -name *.hpp -exec maintenance/autoheader {} copyright.txt \;
find . -type f -name *.cpp -exec maintenance/autoheader {} copyright.txt \;
find . -type f -name \*.am -exec maintenance/autoheader {} copyright.txt \;
find . -type f -name \*.ac -exec maintenance/autoheader {} copyright.txt \;
maintenance/autoheader Makefile.export.phalanx.in copyright.txt
maintenance/autoheader README copyright.txt
