#!/bin/bash

TESTCASES="Hund-1000 Metis Scotch Hund-100 Hund-10000 HUND-1000s"

for filename in $@ ; do
    myname=${filename/.mtx.gz/}
    myname=${myname/.mtx/}
    for NUMP in 1 2 4 6 8 ; do
    cat >> $myname.pbs <<EOF
echo "MATRIX $myname with MMD on $NUMP threads"
 ./bin/pdlinsol -p $NUMP matrices/$myname.mtx.gz
EOF
    for mycase in $TESTCASES ; do
	cat >> $myname.pbs <<EOF
echo "MATRIX $myname with $mycase on $NUMP threads "
./bin/pdlinsol -p $NUMP -c ordering/$myname-$mycase.ord matrices/$myname.mtx.gz
EOF
    done
    done
done
