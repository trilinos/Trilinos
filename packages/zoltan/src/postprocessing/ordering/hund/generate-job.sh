#!/bin/bash

NC=2
NR=2
NUMP=$(($NC*$NR))
TESTCASES="Hund-1000 Metis Scotch Hund-100 Hund-10000"

for filename in $@ ; do
    myname=${filename/.mtx/}
    sed "s/NUMP/$NUMP/" superlu.pbs > $myname.$NUMP.pbs ;
    cat >> $myname.$NUMP.pbs <<EOF
echo "MATRIX $myname with MMD "
aprun -n $NUMP ./bin/pddrive_mmd -c $NC -r $NR matrices/$myname.mtx
EOF
    for mycase in $TESTCASES ; do
	cat >> $myname.$NUMP.pbs <<EOF
echo "MATRIX $myname with $mycase "
aprun -n $NUMP ./bin/pddrive -c $NC -r $NR matrices/$myname.mtx -p ordering/$myname-$mycase.ord
EOF
    done
done
