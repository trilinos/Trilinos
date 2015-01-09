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
