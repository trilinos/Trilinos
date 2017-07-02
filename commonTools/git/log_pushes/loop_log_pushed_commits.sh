#!/bin/bash -e

# See README.md file in this directory for usage

_ABS_FILE_PATH=`readlink -f $0`
_SCRIPT_DIR=`dirname $_ABS_FILE_PATH`
_TRILINOS_BASE_DIR=$_SCRIPT_DIR/../../../..

cd $_TRILINOS_BASE_DIR

echo "Running in Trilinos base directory: $PWD"

./Trilinos/cmake/tribits/python_utils/generic-looping-demon.py \
--command=$_SCRIPT_DIR/log_pushed_commits.sh \
--loop-interval=60s --today-run-till=23:59:00
