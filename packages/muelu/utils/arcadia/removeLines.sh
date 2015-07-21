#/bin/sh
# summary:   delete all instances of a pattern and some number of lines following
# usage:     removeLines.sh [-h] -p <pattern> [-n #lines] <file>
#
# example:
#    removeLines.sh "product AB .Pmat_" -n 4 run_1/screen_1.ml
#       will find all instances of "produce AB (Pmat_"  and delete it and the following 4 lines (5 lines total)
#
#    removeLines.sh "product AB .Pmat_|AB .Rebalancing" -n 4 run_1/screen_1.ml
#       will find all instances of "produce AB (Pmat_" *or* "AB (Rebalancing", delete them and
#       the following 4 lines (5 lines total)
#
# [Borrowed the command line parsing from
# http://blog.mafr.de/2007/08/05/cmdline-options-in-shell-scripts]

USAGE="Usage: `basename $0` [-h] -p pattern [-n numberOfLinesToCut] file"
OPTDESCR="\n  -h        help\n  -p <pattern>   pattern to match\n  -n       <n>   number of lines to excise (default=4)\n"

COUNT=4
numReqd=1;
# Parse command line options.
while getopts hp:n: OPT; do
    case "$OPT" in
        p)
            PATTERN=$OPTARG
            ;;

        n)
            COUNT=$OPTARG
            ;;

        h)
            echo
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

FILE=$1;

# Access additional arguments as usual through
# variables $@, $*, $1, $2, etc. or using this loop:
#for PARAM; do
#    echo "##  $PARAM"
#done

ttt=`echo ${PATTERN} xxx | sed "s/ //g"`
if test `expr ${ttt}` = 'xxx'
then
  echo "you must provide a pattern"
  exit
fi

awk -v pat="$PATTERN" -v N="$COUNT" '
BEGIN{}
{
  #match found, stop parsing and get the next line
  if ($0 ~ pat)
  {
    n=N
    next
  }

  #nothing was found previously, print this line
  if (n==0)
  {
    print
  }

  #pattern was found previously, skip n lines
  if (n>0) {n--;}
}
END{}
' $FILE
