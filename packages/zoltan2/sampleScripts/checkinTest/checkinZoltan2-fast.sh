#!/bin/sh
../checkin-test.py --no-eg-git-version-check --st-extra-builds=MPI_ss,SERIAL_ss --do-all --no-enable-fwd-packages --enable-packages=Zoltan2 -j32

