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

#ifndef CLOP_SUB_HPP
#define CLOP_SUB_HPP
#include <stdio.h>
#include <fstream>
#include <iomanip>
#include <math.h>
#include "sparse_lu.hpp"
#include "myzero.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_IntVector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_LAPACK.h"
#include "Epetra_BLAS.h"
#include "Epetra_Comm.h"

class CLOP_sub {
 public:
  CLOP_sub();
  ~CLOP_sub();
  void getmatrix_nnz(int subdofs_[], int ndof_sub, const Epetra_CrsMatrix *A, 
		     int imap[], const Epetra_Comm* Comm, int &nnz,
		     int scale_option_);
  void factormatrix(const Epetra_CrsMatrix *A, int imap[], int rowbeg[], 
		    int colidx[], double K[]);
  void genpu(const Epetra_IntVector *LD, const Epetra_MultiVector *Coords, 
	     double rhs[], double sol[], double temp[], int atype_sub, 
	     int ndim_sub, double WORK[], int LWORK, double Edof_sub[],
	     int & nneg);
  void normalpu(double Edof_sub[], unsigned char nsubdof[]);
  void construct_coarse1(const Epetra_CrsMatrix *A, double rhs[], 
	     double sol[], double temp[], int rowbeg[], int colidx[], 
             double K[], int imap[], unsigned char nsubdof[], int & csdimP, 
	     int & ndof_rot_);
  void get_cdof(int cs_local[], int & csdimP, double & xcent_,
		double & ycent_, double & zcent_);
  void sum_scalar_multiply(double Edof_sub[], int rbm, double alpha);
  void correct(double Edof_sub[], int rbm, unsigned char nsubdof[], int ldof[],
	       int lbound, int ubound);
  void statcond(unsigned char nsubdof[], unsigned char on_sub_bound[], 
	  int imap[], int rowbeg[], int colidx[], double K[], double rhs[], 
	  double sol[], double temp[], const Epetra_CrsMatrix *A);
  void get_Phi_ptr(double* & Phi_ptr);
  void subpre(double rr[], double zz[], double rhs_work[], 
	      double sol_work[], double tmp_work[]);
 private: // variables
  const Epetra_Comm *Comm;
  int ndof, csdim, csdim_max, ndim, atype, INFO, MyPID, scale_option;
  double xcent, ycent, zcent, *x, *y, *z, *Edof, *Phi;
  CLAPS_sparse_lu *A_sub;
  int *jpvt, *subdofs, *locdof;
};
#endif // CLOP_SUB_HPP
  
