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

#include "krylov_crd.hpp"
#include "Claps_ConfigDefs.hpp"  // for definition of F77_FUNC
#define DSTEV_F77 F77_FUNC(dstev, DSTEV)
#define DGEMM_F77 F77_FUNC(dgemm, DGEMM)
#define DGEMV_F77 F77_FUNC(dgemv, DGEMV)
#define DPOTRF_F77 F77_FUNC(dpotrf, DPOTRF)
#define DPOTRS_F77 F77_FUNC(dpotrs, DPOTRS) 

extern "C"{
  void DGEMM_F77(char* TRANSA, char* TRANSB, int* M, int* N, int* K,
		 double* ALPHA, double A[], int* LDA, double B[],
		 int* LDB, double* BETA, double C[], int*LDC);
  void DGEMV_F77(char* TRANS, int* M, int* N, double* ALPHA,
		 double A[], int* LDA, double X[], int* INCX,
		 double* BETA, double Y[], int* INCY);
  void DSTEV_F77(char* JOBZ, int* N, double D[], double E[],
		 double Z[], int* LDZ, double WORK[], int* INFO);
  void DPOTRF_F77(char* UPLO, int* N, double A[], int* LDA, 
		  int* INFO);
  void DPOTRS_F77(char* UPLO, int* N, int* NRHS, double A[],
		  int* LDA, double B[], int* LDB, int* INFO);
}

krylov_crd::krylov_crd(preconditioner_crd* pptr_,
		       int krylov_method_,
		       int max_iter_,
		       double solver_tol_,
		       int max_search_dirs_,
		       int vector_length_,
		       int print_flag_,
		       int MyPID_)
  : pptr(pptr_), krylov_method(krylov_method_), max_iter(max_iter_), 
    solver_tol(solver_tol_), max_search_dirs(max_search_dirs_),
    vector_length(vector_length_), print_flag(print_flag_)
{
  MyPID = MyPID_;
  if (MyPID != 0) print_flag = 0;
  zero_pointers();
  num_search_used = 0;
  //
  // allocate memory
  //
  if (krylov_method == 0) {
    rcurra          = new double[max_iter+1];
    rhoa            = new double[max_iter+1];
    betaa           = new double[max_iter+1];
    pApa            = new double[max_iter+1];
    Dtri            = new double[max_iter+1];
    Etri            = new double[max_iter+1];
    econa           = new double[max_iter+1];
    AP_matrix       = new double[vector_length*max_search_dirs];
    P_matrix        = new double[vector_length*max_search_dirs];
    p               = new double[vector_length];
    r               = new double[vector_length];
    z               = new double[vector_length];
    rwork           = new double[vector_length];
    work_search     = new double[max_search_dirs];
    work_search_sum = new double[max_search_dirs];
    PAP_matrix      = new double[max_search_dirs*max_search_dirs];
    PAP_matrix_sum  = new double[max_search_dirs*max_search_dirs];
  }
}

krylov_crd::~krylov_crd()
{
  delete [] rcurra; delete [] rhoa; delete [] betaa; delete [] pApa;
  delete [] Dtri; delete [] Etri; delete [] AP_matrix; delete [] P_matrix;
  delete [] p; delete [] r; delete [] z; delete [] work_search;
  delete [] work_search_sum; delete [] rwork; delete [] econa;
  delete [] PAP_matrix; delete [] PAP_matrix_sum;
}

void krylov_crd::zero_pointers()
{
  rcurra = 0; rhoa = 0; betaa = 0; pApa = 0; Dtri = 0; Etri = 0;
  AP_matrix = 0; P_matrix = 0; p = 0; r = 0; z = 0; work_search = 0;
  work_search_sum = 0; rwork = 0; econa = 0; PAP_matrix = 0;
  PAP_matrix_sum = 0;
}

int krylov_crd::solve(double u [], const double f[])
{
  int solver_status;
  if (krylov_method == 0) solver_status = solve_pcg(u, f);
  return solver_status;
}

int krylov_crd::solve_pcg(double u [], const double f [])
{
  int i, sod, num_search_orig, init_val;
  double rorig, rcurr, dprod, beta, roldzold, alpha, pAp, ractual;

  sod = sizeof(double);
  memcpy(r, f, vector_length*sod);
  myzero(u, vector_length);
  rorig = pptr->norm2(r, vector_length);
  rcurr = rorig;
  if (print_flag != 0) cout << "initial residual = " << rorig << endl;
  init_val = pptr->initialize_solve(u, r);
  if (init_val == 1) {
    pptr->A_times_x(u, z);
    for (i=0; i<vector_length; i++) r[i] = r[i] - z[i];
    rcurr = pptr->norm2(r, vector_length);
  }
  //
  // initial correction
  //
  num_search_orig = num_search_used;
  if (num_search_orig > 0) {
    search_dir_correction(r, u, num_search_orig);
    pptr->A_times_x(u, z);
    for (i=0; i<vector_length; i++) r[i] -= z[i];
    rcurr = pptr->norm2(r, vector_length);
    if (print_flag != 0) {
      cout << "number of search dirs used       = " << num_search_orig << endl;
      cout << "residual after using search dirs = " << rcurr << endl;
    }
  }
  if (rcurr/rorig < solver_tol) return 0;
  for (int iter=0; iter<max_iter; iter++) {
    if (print_flag != 0)
      cout << "iteration " << iter+1 << " of maxiter = " << max_iter << endl;
    //
    // calculate preconditioned residual
    //
    pptr->apply_preconditioner(r, z);
    if (num_search_orig > 0) {
      pptr->A_times_x(z, rwork);
      for (i=0; i<vector_length; i++) rwork[i] = r[i] - rwork[i];
      search_dir_correction(rwork, z, num_search_orig);
    }
    //
    // standard pcg stuff
    //
    dprod = pptr->dotprod(r, z, vector_length);
    //    assert(dprod >= 0);
    rhoa[iter] = sqrt(fabs(dprod));
    if (iter == 0) {
      beta = 0;
        memcpy(p, z, vector_length*sod);
    }
    else {
      beta = dprod/roldzold;
      for (i=0; i<vector_length; i++) p[i] = beta*p[i] + z[i];
    }
    betaa[iter] = beta;
    roldzold = dprod;
    pptr->A_times_x(p, z);
    pAp = pptr->dotprod(p, z, vector_length);
    pApa[iter] = pAp;
    alpha = dprod/pAp;
    if (print_flag == -1) {
      cout << "dprod, pAp, alpha, beta = " << dprod << " " << pAp 
	   << " " << alpha << " " << beta << endl;
    }
    for (i=0; i<vector_length; i++) {
      u[i] += alpha*p[i];
      r[i] -= alpha*z[i];
    }
    rcurr = pptr->norm2(r, vector_length);
    rcurra[iter+1] = rcurr;
    num_iter = iter+1;
    //
    // store search direction
    //
    store_search_dir(pAp, p, z);
    if (rcurr/rorig <= solver_tol) break;
  }
  //
  // factor matrix for stored vectors
  //
  factor_matrix(num_search_orig);
  pptr->A_times_x(u, z);
  memcpy(r, f, vector_length*sod);
  for (i=0; i<vector_length; i++) r[i] -= z[i];
  ractual = pptr->norm2(r, vector_length);
  if (print_flag != 0) {
    cout << "num_iter = " << num_iter << endl;
    if (num_iter > 0) {
      cout << "rcurr(recursive)      = " << rcurr << endl;
      cout << "rcurr(actual)         = " << ractual << endl;
      cout << "number of iterations  = " << num_iter << endl;
      cout << "solver tolerance      = " << solver_tol << endl;
      cout << "condition # estimate      relative residual" 
	   << "   iteration" << endl;
      calculate_condition(num_iter);
      cout << setiosflags(ios::scientific | ios::uppercase);
      for (i=0; i<num_iter; i++) {
	double ee = 0;
	ee = econa[i]; 
	cout << " " 
	     << setw(17) << setprecision(10) << ee 
	     << "       " 
	     << setw(17) << setprecision(10) << rcurra[i+1]/rorig
	     << "        " 
	     << i+1 << endl;
      }
    }
    cout << resetiosflags(ios::scientific);
    cout << resetiosflags(ios::uppercase);
    cout << setprecision(6);
  }
  if (rcurr/rorig <= solver_tol) return 0;
  else return 1;
}

void krylov_crd::search_dir_correction(double rhs [], double sol [] , int nn)
{
  char TRANS('T'), UPLO('U');
  int INCX(1), INCY(1), NRHS(1), INFO;
  double ALPHA(1), BETA(0);
  DGEMV_F77(&TRANS, &vector_length, &nn, &ALPHA,
	    P_matrix, &vector_length, rhs, &INCX, &BETA,
	    work_search, &INCY);
  pptr->sum_vectors(work_search, nn, work_search_sum);
  DPOTRS_F77(&UPLO, &nn, &NRHS, PAP_matrix_sum, &nn,
	     work_search_sum, &nn, &INFO);
  assert(INFO == 0);
  TRANS = 'N'; BETA = 1;
  DGEMV_F77(&TRANS, &vector_length, &nn, &ALPHA,
	    P_matrix, &vector_length, work_search_sum, &INCX,
	    &BETA, sol, &INCY);
}

void krylov_crd::calculate_condition(int miter)
{
  int i, j, INFO, ip1, one(1);
  double Z, WORK;
  char N = 'N';
  if (miter == 1) {
    econa[0] = 1;
    return;
  }
  for (i=0; i<miter; i++) {
    ip1 = i + 1;
    Dtri[0] = pApa[0]/rhoa[0]/rhoa[0];
    for (j=1; j<ip1; j++) {
      Dtri[j]   = (pApa[j-1]*betaa[j]*betaa[j]+pApa[j])/rhoa[j]/rhoa[j];
      Etri[j-1] = -pApa[j-1]*betaa[j]/rhoa[j-1]/rhoa[j];
    }
    DSTEV_F77(&N, &ip1, Dtri, Etri, &Z, &one, &WORK, &INFO); 
    if (INFO != 0) {
      cout << "error in call to dstev in amge_solver::calculate_condition" 
	   << endl;
      cout << "INFO = " << INFO << endl;
    }
    econa[i] = Dtri[i]/Dtri[0];
  }
  if (print_flag == 2) {
    cout << "eigenvalue estimates = " << endl;
    for (i=0; i<miter; i++) cout << Dtri[i] << endl;
  }
}

void krylov_crd::store_search_dir(const double pAp, double pvec [],
				  double Apvec [])
{
  int i, ibeg;
  double sfac;
  if (num_search_used < max_search_dirs) {
    ibeg = vector_length*num_search_used;
    sfac = 1/sqrt(pAp);
    for (i=0; i<vector_length; i++) {
      P_matrix[ ibeg+i] = sfac*pvec[i];
      AP_matrix[ibeg+i] = sfac*Apvec[i];
    }
    num_search_used++;
  }
}

void krylov_crd::factor_matrix(const int num_search_orig)
{
  char TRANSA('T'), TRANSB('N'), UPLO('U');
  double ALPHA(1), BETA(0);
  int M, N, K, INFO;
  if ((num_search_used > 0) && (num_search_used != num_search_orig)) {
    //
    // calculate P^T * A * P
    //
    M = N = num_search_used;
    K = vector_length;
    DGEMM_F77(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, P_matrix,
	      &K, AP_matrix, &K, &BETA, PAP_matrix, &M);
    pptr->sum_vectors(PAP_matrix, M*N, PAP_matrix_sum);
    //
    // factor P^T * A * P
    //
    DPOTRF_F77(&UPLO, &N, PAP_matrix_sum, &N, &INFO);
    assert(INFO == 0);
  }
}
