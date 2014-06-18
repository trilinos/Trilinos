#!/bin/sh
../checkin-test.py --no-eg-git-version-check --st-extra-builds=MPI_ss,SERIAL_ss --do-all --disable-packages=TrilinosCouplings -j32

