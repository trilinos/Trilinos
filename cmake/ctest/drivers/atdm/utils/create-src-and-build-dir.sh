#
# Helper script that creates the SRC_AND_BUILD subdir and moves in
# 

if [ "${WORKSPACE}" == ""  ] ; then
  echo "Error, must set WORKSPACE var before calling!"
  exit 1
fi

if [ "${JOB_NAME}" == ""  ] ; then
  echo "Error, must set JOB_NAME var before calling!"
  exit 1
fi

# Shiller/hansen settings for automated job (rest are set in ./environment.sh)
ulimit -c 0

export SUBDIR=SRC_AND_BUILD
if [ ! -e $SUBDIR ] ; then
  echo "Making $SUBDIR"
  mkdir $SUBDIR
fi

cd $SUBDIR/
echo "Current dir: $PWD"

