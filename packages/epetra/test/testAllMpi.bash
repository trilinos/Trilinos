#!/usr/bin/bash
file=log.mpi.`eval date +%d%b%Y_%H%M%S`
rm -f $file
echo $file
echo "Date: " `eval date` >> $file
echo `uname -a` >> $file
DIRS="BlockMap Comm CrsGraph CrsMatrix CrsRectMatrix FECrsGraph FECrsMatrix FEVbrMatrix FEVector ImportExport Map MapColoring MultiVector Object VbrMatrix Vector"
for f in $DIRS
do
	cd "$f"
	make -f classicMakefile TRILINOS_ID= TRILINOS_COMM=MPI clobber; make -f classicMakefile TRILINOS_ID= TRILINOS_COMM=MPI
	for g in *.exe
     do
		echo "############" "$g" "##############" >> ../$file
		mpirun -np 4 "$g" -v >> ../$file
	done
	make -f classicMakefile TRILINOS_ID= TRILINOS_COMM=MPI clobber
	make -f classicMakefile TRILINOS_ID= TRILINOS_COMM=MPI clobber
	cd ..
done
