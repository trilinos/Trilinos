#!/bin/bash

for d in */ ; do
  cd $d;
  rm samples.txt obj_samples.txt;
  find . -type f -name 'samples_*.txt' -exec cat {} + >> samples.txt;
  find . -type f -name 'obj_samples_*.txt' -exec cat {} + >> obj_samples.txt;
  cd ..;
done
