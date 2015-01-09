//@HEADER
// ************************************************************************
// 
//         Claps: A Collection of Domain Decomposition Preconditioners
//                and Solvers
//         Copyright (2006) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Clark R. Dohrmann (crdohrm@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#ifndef SPARSE_LU_HPP
#define SPARSE_LU_HPP
#include <stdio.h>
#include <fstream>
#include <iomanip>
#include <assert.h>
//#include "my_feti_sparse_solver.hpp"
#include "Epetra_LAPACK.h"

class CLAPS_sparse_lu {
 public:
  CLAPS_sparse_lu();
  ~CLAPS_sparse_lu();
  void cleanup();
  int factor(int N_, int NNZ_, int COLPTR[], int ROWIDX[], double ANZ[],
	     int scale_flag_ = 0);
  int sol(int NRHS, double RHS[], double SOL[], double TEMP[]);
 private: //
  void getnrm(int n, int colptr[], int rowidx[], 
	      double values[], double &anorm);
  void inpnv(int &n , int colptr[], int rowidx[], double values[], 
	     int perm[], int invp [], int &nsuper, int xsuper[], 
	     int xlindx[], int lindx[], int xlnz[], double lnz[],
	     int offset[]);
  int small_factor(int rowbeg[], int colidx[], double vals[]);
  int small_solve(int NRHS, double RHS[], double SOL[]);
  Epetra_LAPACK EL;
  int N, DEFBLK, NSUPER, NDEF, LBDEF, max_small, scale_flag;
  int *XSUPER;
  int *XLINDX;
  int *LINDX;
  int *XLNZ;
  int *PERM;
  int *INVP;
  int *IPROW;
  int *IPCOL;
  int *DEF;
  double* LNZ;
  double* NS;
  double *SCALE;
};
#endif // SPARSE_LU_HPP
  
