//
//  This code populates an empty EpetraCrsMatrix with a tridiagonal with 
//  -1 on the off-diagonals and 2 on the diagonal.  
//
//  This code was plaguerized from epetra/example/petra_power_method/cxx_main.cpp
//  presumably written by Mike Heroux.
//
//  Adapted by Ken Stanley - Aug 2003 
//
#include "Epetra_CrsMatrix.h"

int CreateTridi(Epetra_CrsMatrix& A) ;
