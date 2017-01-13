#!/bin/bash

find . -type f -name 'samples_*' -exec cat {} + >> samples.txt;

find . -type f -name 'obj_samples_*' -exec cat {} + >> obj_samples.txt;
