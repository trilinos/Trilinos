#!/bin/bash

# To use:
# 1) Make sure you've compiled with the CreateOperator tests and that they run.
# 2) Run ctest
# 3) cd $BUILDDIR/packages/muelu/test/interface/Output
# 4) Run: $SOURCEDIR/packages/muelu/test/interface/default/Output/rebase.sh
# 5) Don't forget the heavy tests, which do not run during checkin. (See comments at the end.)

echo "********************************************"
echo "This is the rebase script for SCALAR=COMPLEX."
echo "********************************************"

RESULTSDIR=`pwd`
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Sanity
if [ "$RESULTSDIR" == $SCRIPTDIR ] || ! [[ "$RESULTSDIR" =~ "/Output" ]]; then
    echo "$0: Please run this from the \$BUILDDIR/packages/muelu/test/interface/Output directory"
    exit -1;
fi

for file in *.out; do
    GOLDFILE=${file%%.out}.gold

    # Diff test the filtered files
    diff -u -w -I"^\\s*$" ${file}_filtered ${GOLDFILE}_filtered >& /dev/null
    returncode=$?

    # Only rebase diffing files"
    if [ "$returncode" -eq 1 ]; then
	echo "$file diffs, rebasing"
	cp ${file}_filtered $SCRIPTDIR/$GOLDFILE
    fi
done

echo ""
for GOLDFILE in *.gold; do
    OUTFILE=${GOLDFILE/.gold/.out}
    if [ ! -f $OUTFILE ]; then
        echo "The test corresponding to $GOLDFILE has not been run."
    fi
done

echo "This is the rebase script for SCALAR=COMPLEX."
echo "If you are trying to rebase a double build, use the script from test/interface/default/Output/"
echo ""
echo ""
echo "Did you remember to rebase the \"heavy\" interface tests?"
echo "To do so, run the following:"
echo ""
echo "  MueLu_ParameterListInterpreter.exe --linAlgebra=Epetra --heavytests"
echo "  MueLu_ParameterListInterpreter.exe --linAlgebra=Tpetra --heavytests"
echo "  mpirun -np 4 MueLu_ParameterListInterpreter.exe --linAlgebra=Epetra --heavytests"
echo "  mpirun -np 4 MueLu_ParameterListInterpreter.exe --linAlgebra=Tpetra --heavytests"
echo ""
echo " and then rerun this script."
