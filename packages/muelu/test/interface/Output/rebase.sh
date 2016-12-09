#!/bin/sh

# To use: 
# 1) Make sure you've compiled with the CreateOperator tests and that they run.
# 2) Run ctest 
# 3) cd $BUILDDIR/packages/muelu/test/interface/Output
# 4) Run: $SOURCEDIR/packages/muelu/test/interface/Output/rebase.sh



RESULTSDIR=`pwd`
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Sanity
if [ "$RESULTSDIR" == $SCRIPTDIR ]; then
    echo "$0: Please run this from the packages/muelu/test/interface/Output directory"
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
	
	# Remove trailing whitespace 
	sed -i 's/[ \t]*$//' $file

	# Copy
	cp $file $SCRIPTDIR/$GOLDFILE
    fi
done

