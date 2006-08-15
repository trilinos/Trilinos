#!/bin/sh
find . -type f -name *.h -exec maintenance/autoheader {} maintenance/header.txt \;
find . -type f -name *.cpp -exec maintenance/autoheader {} maintenance/header.txt \;
find . -type f -name *.H -exec maintenance/autoheader {} maintenance/header.txt \;
find . -type f -name *.C -exec maintenance/autoheader {} maintenance/header.txt \;
find . -type f -name \*.am -exec maintenance/autoheader {} maintenance/header.txt \;
find . -type f -name \*.ac -exec maintenance/autoheader {} maintenance/header.txt \;
maintenance/autoheader Makefile.export.nox.in maintenance/header.txt
