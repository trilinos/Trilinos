//
//  CreateTridi populates an empty EpetraCrsMatrix with a tridiagonal with 
//  -1 on the off-diagonals and 2 on the diagonal.  
//
//  CreateTridiPlus creates the same matrix as CreateTridi except that it adds
//  -1 in the two off diagonal corners.
//
//  This code was plaguerized from epetra/example/petra_power_method/cxx_main.cpp
//  presumably written by Mike Heroux.
//
//  Adapted by Ken Stanley - Aug 2003 
//
#include "Epetra_CrsMatrix.h"

int CreateTridi(Epetra_CrsMatrix& A) ;

int CreateTridiPlus(Epetra_CrsMatrix& A) ;
