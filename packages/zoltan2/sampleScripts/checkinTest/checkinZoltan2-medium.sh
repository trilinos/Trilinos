#!/bin/sh
../checkin-test.py --no-eg-git-version-check --ss-extra-builds=MPI_ss,SERIAL_ss --do-all --disable-packages=TrilinosCouplings -j32

