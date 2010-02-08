#!/usr/bin/env sh
#
# usage: remove-double-templating.sh file_name
#
# This script removes the last remnants of double templating
# <RangeScalar,DomainScalar> and replaces it with just <Scalar>.
# 
# Running this scrit on the same file multiple times is fine.
#

FILENAME=$1
#echo $FILENAME

perl -i -pe 's/template *< *class +RangeScalar *, *class +DomainScalar *= *RangeScalar *>/template<class Scalar>/g' $FILENAME
perl -i -pe 's/template *< *class +RangeScalar *, *class +DomainScalar *>/template<class Scalar>/g' $FILENAME
perl -i -pe 's/< *RangeScalar *, *DomainScalar *>/<Scalar>/g' $FILENAME
perl -i -pe 's/< *Scalar *, *Scalar *>/<Scalar>/g' $FILENAME
perl -i -pe 's/RangeScalar/Scalar/g' $FILENAME
perl -i -pe 's/DomainScalar/Scalar/g' $FILENAME
