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
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
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
  
