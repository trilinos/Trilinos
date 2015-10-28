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

#ifndef CRD_UTILS_HPP
#define CRD_UTILS_HPP
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include "Epetra_CrsMatrix.h"
#include "Epetra_MpiComm.h"
#include "Epetra_BLAS.h"
#include "Epetra_LAPACK.h"
#include "sparse_lu.hpp"
#include "myzero.hpp"

namespace CRD_utils {
  int find_index(int a[], int n, int gdof);
  void sort_and_cull(int a[], int n, int & m);
  void Epetra_datfile(const Epetra_CrsMatrix* A, char fname[]);
  void Epetra_datfile(int* A, int N, char fname[]);
  void spmat_datfile(int nrow, int rowbegp [], int colidxp [],
		     double val[], char fname[]);
  void scale_columns(Epetra_CrsMatrix* A, 
		     const int norm_opt, 
		     const int blocksize);
  void get_column_norm(Epetra_CrsMatrix* A,
		       const int norm_opt,
		       const int blocksize,
		       double *col_norm);
  void tie_down_coarse(int n, 
		       int rowbeg[], 
		       int colidx[], 
		       double vals[], 
		       int ne, 
		       int scale_option,
		       int & num_tied_down,
		       int* & tied_down, 
		       CLAPS_sparse_lu* & AA, 
		       double* & Xvecs);
  void subspace_iteration(int n, 
			  int rowbeg[], 
			  int colidx[], 
			  double vals[], 
			  bool bound_flag[],
			  int scale_option,
			  int & nextra, 
			  int* & extra_corner,
			  CLAPS_sparse_lu* & A, 
			  double* & Xvecs, 
			  int ne = 0);
  class Graph_class 
  {
  public:
    Graph_class(int N_, int A1_[], int A2_[]);
    ~Graph_class();
    void Components(int component_[], int & comp_num);
  private:
    void DFS(const int v, const int comp_num);
    int N;
    int *A1, *A2;
    int *component;
  };
}
#endif // CRD_UTILS_HPP
