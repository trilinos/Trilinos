#!/bin/sh
#
# bash script to compare Epetra and Tpetra timings from MueLu_ScalingTest.exe side-by-side.
#
# syntax:
#    compareTE.sh [-hasl] file1 file2
#
# example:
#    compareTE.sh -a screen.epetra screen.tpetra
#
# optional arguments:
#   -h           help
#   -a           all timings (no by default)
#   -s           timings summed over all levels, i.e., no breakdown by level (no by default)
#   -l           timings level-by-level (no by default)
#
# required arguments
#
#   file1,file2  epetra/tpetra screen dumps  (order doesn't matter)
#
# You should see something like the following:
#                                                                                          Tpetra       Epetra    T/E ratio
#                                                                                          --------------------------------
#                          MueLu: IfpackSmoother: Setup Smoother (total, level=0)  ==>     0.0325       0.0119       2.72
#                MueLu: AmalgamationFactory: AmalgamationFactory (total, level=2)  ==>     0.0173       0.0130       1.34
#                              MueLu: TransPFactory: Transpose P (total, level=4)  ==>     0.0021       0.0006       3.70
#                         MueLu: AmalgamationFactory: AmalgamationFactory (total)  ==>     0.3636       0.3431       1.06
#                                   MueLu: AmesosSmoother: Setup Smoother (total)  ==>     0.0001       0.0008       0.10
#
# The first two columns are times in seconds, the last column is the ratio of the first two columns.

# Borrowed the command line parsing from
# http://blog.mafr.de/2007/08/05/cmdline-options-in-shell-scripts

# Find the script location
SCRIPT=$(readlink -f $0)
SCRIPTPATH=$(dirname $SCRIPT)


USAGE="Usage: `basename $0` [-hasl] file1 file2"
OPTDESCR="\n  -h  help\n   -a  all timings [off]\n  -s  timings summed over all levels [off]\n  -l  timings level-by-level[off]\n"

numReqd=2;
printAll=1;
printSummedStats=0;
printLevelStats=0;
# Parse command line options.
while getopts ahsl OPT; do
    case "$OPT" in
        h)
            echo $USAGE
            echo -e $OPTDESCR
            exit 0
            ;;
        a)
            printAll=1;
            ;;
        s)
            printSummedStats=1;
            ;;
        l)
            printLevelStats=1;
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

# We want at least two non-option arguments.
# Remove this block if you don't need it.
if (( $# < $numReqd )); then
    echo $USAGE >&2
    exit 1
fi

file1=$1;
file2=$2;

# Access additional arguments as usual through
# variables $@, $*, $1, $2, etc. or using this loop:
#for PARAM; do
#    echo $PARAM
#done

if (( $printAll == 1 )); then
  printSummedStats=1;
  printLevelStats=1;
fi

#if not defined, set path so that we can find compareTE.awk
set ${AWKPATH:=$PATH}
export AWKPATH

#ttt=`awk --version`
#echo "awk info: $ttt"

awk -f $SCRIPTPATH/compareTE.awk -v printSummedStats=${printSummedStats} -v printLevelStats=${printLevelStats} $file1 $file2
