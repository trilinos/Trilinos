#!/bin/sh
#
# bash script to compare Epetra and Tpetra timings side-by-side.
#
# syntax:
#    mueprof.sh [-h] [ARGS] file1 file2
#
# example:
#    mueprof.sh screen.tpetra screen.epetra
#
# Borrowed the command line parsing from
# http://blog.mafr.de/2007/08/05/cmdline-options-in-shell-scripts

# Find the script location
SCRIPT=$(readlink -f $0)
SCRIPTPATH=$(dirname $SCRIPT)

USAGE="Usage: `basename $0` [-h] [ARGS] file1 file2"
OPTDESCR="\n\t-h \thelp\n\t-n \twhich solver block to analyze [1]\n
          \t-s \twhich linear algebra lib to sort by: \"Epetra\",[\"Tpetra\"]\n
          \t-t \thow to display delta in times: \"ratio\",[\"diff\"]\n
          \t-a \tallows comparison of two Epetra or two Tpetra files\n
          \t-d \tdebug mode\n"

numReqd=1;
blockNumber=1;
sortByLib="Tpetra"
deltaDisplay="diff"
agnostic=0
debug=0
# Parse command line options.
while getopts hn:s:t:ad OPT; do
    case "$OPT" in
        d)
            debug=1
            ;;
        a)
            agnostic=1
            ;;
        t)
            deltaDisplay=$OPTARG
            ;;
        s)
            sortByLib=$OPTARG
            ;;
        n)
            blockNumber=$OPTARG
            ;;
        h)
            echo $USAGE
            echo -e $OPTDESCR
            exit 0
            ;;
        \?)
            # getopts issues an error message
            echo $USAGE >&2
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

file1=$1;
file2=$2;

# Access additional arguments as usual through
# variables $@, $*, $1, $2, etc. or using this loop:
#for PARAM; do
#    echo $PARAM
#done

#if not defined, set path so that we can find mueprof.awk
set ${AWKPATH:=./}
export AWKPATH

#ttt=`awk --version`
#echo "awk info: $ttt"

awk -v "blockNumber=$blockNumber" -v "sortByLib=$sortByLib" -v "etDelta=$deltaDisplay" -v "agnostic=$agnostic" -v "debug=$debug" -f $SCRIPTPATH/mueprof.awk $file1 $file2
