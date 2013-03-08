#!/bin/sh
#
# bash script to compare Epetra and Tpetra timings from MueLu_ScalingTest.exe side-by-side.
#
# syntax:
#    mueprof.sh [-h] file
#
# prerequisites:
#
# 1) MueLu must be configured with -D MueLu_ENABLE_PROFILING:STRING=ON.
# 2) Timer parent/child relations must be written out by putting
#
#       MueLu::MutuallyExclusiveTime<MueLu::BaseClass>::PrintParentChildPairs();
#
#    at the end of your program.
#
# optional arguments:
#   -h           help
#
# required arguments:
#
#   file  epetra or tpetra screen dumps
#
# example:
#    mueprof.sh screen.tpetra
#
# Borrowed the command line parsing from
# http://blog.mafr.de/2007/08/05/cmdline-options-in-shell-scripts

USAGE="Usage: `basename $0` [-h] file"
OPTDESCR="\n  -h  help\n"

numReqd=1;
# Parse command line options.
while getopts ahsl OPT; do
    case "$OPT" in
        h)
            echo $USAGE
            echo -e $OPTDESCR
            exit 0
            ;;
        \?)
            # getopts issues an error message
            echo $USAGE >&2
            echo "2"
            exit 1
            ;;
    esac
done

# Remove the switches we parsed above.
shift `expr $OPTIND - 1`

# We want at least $numReqd non-option arguments. 
# Remove this block if you don't need it.
if (( $# < $numReqd )); then
    echo $USAGE >&2
    exit 1
fi

file=$1;

# Access additional arguments as usual through 
# variables $@, $*, $1, $2, etc. or using this loop:
#for PARAM; do
#    echo $PARAM
#done

#if not defined, set path so that we can find mueprof.awk
set ${AWKPATH:=$PATH}
export AWKPATH

#ttt=`awk --version`
#echo "awk info: $ttt"

awk -f mueprof.awk $file
