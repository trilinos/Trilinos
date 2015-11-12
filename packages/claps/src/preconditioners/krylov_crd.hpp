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

#ifndef KRYLOV_CRD_HPP
#define KRYLOV_CRD_HPP

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <math.h>
#include <assert.h>
#include "myzero.hpp"
#include "preconditioner_crd.hpp"

using namespace std;

class krylov_crd 
{
public: // functions
  krylov_crd(preconditioner_crd* pptr_,
	     int krylov_method_,
	     int max_iter_,
	     double solver_tol_,
	     int num_search_dirs_,
	     int vector_length_,
	     int print_flag_,
	     int MyPID_=0);
  ~krylov_crd();
  int solve(double u [], const double f []);
private:
  void zero_pointers();
  void search_dir_correction(double rhs [], double sol [] , int nn);
  void calculate_condition(int miter);
  void store_search_dir(const double pAp, double pvec [], double Apvec []);
  void factor_matrix(const int num_search_orig);
  int solve_pcg(double u [], const double f []);
  
  preconditioner_crd *pptr;
  int krylov_method, max_iter, max_search_dirs, vector_length, print_flag;
  const double solver_tol;
  int num_search_used, num_iter, MyPID;
  double *rcurra, *rhoa, *betaa, *pApa, *Dtri, *Etri, *econa;
  double *AP_matrix, *P_matrix, *p, *r, *z, *rwork, *PAP_matrix;
  double *work_search, *work_search_sum, *PAP_matrix_sum;
protected:
  
};
#endif // KRYLOV_CRD_HPP
