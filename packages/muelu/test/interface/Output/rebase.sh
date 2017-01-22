#!/bin/sh

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
	cp $file $SCRIPTDIR/$GOLDFILE
    fi
done

