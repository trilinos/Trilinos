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
