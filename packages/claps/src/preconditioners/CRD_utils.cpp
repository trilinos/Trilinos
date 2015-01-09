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

#include "CRD_utils.hpp"

namespace CRD_utils {

int find_index(int a[], int n, int gdof)
{
  int *loc, loc_col;
  loc = std::lower_bound(&a[0], &a[n], gdof);
  loc_col = (int) (loc - &a[0]);
  if ((loc_col < 0) || (loc_col > n-1)) return -1;
  if (a[loc_col] != gdof) return -1;
  return loc_col;
}
  
void sort_and_cull(int a[], int n, int & m)
{
  if (n == 0) {
    m = 0;
    return;
  }
  std::sort(&a[0], &a[n]);
  int previd = a[0];
  m = 1;
  for (int i=1; i<n; i++) {
    if (a[i] != previd) {
      previd = a[i];
      a[m] = a[i];
      m++;
    }
  }
}

void Epetra_datfile(const Epetra_CrsMatrix* A, char fname[])
{
  int i, j, NumEntries, *Indices;
  double *Values;
  std::ofstream fout;
  sprintf(fname, "%s.dat", fname);
  fout.open(fname);
  for (i=0; i<A->NumMyRows(); i++) {
    A->ExtractMyRowView(i, NumEntries, Values, Indices);
    for (j=0; j<NumEntries; j++) {
      fout << A->GRID(i) << " " << A->GCID(Indices[j]) << std::setw(22) << std::setprecision(15)
	   << Values[j] << std::endl;
    }
  }
  fout.close();
}

void spmat_datfile(int nrow, int rowbegp [], int colidxp [],
				 double val[], char fname[])
{
  int i, j;
  std::ofstream fout;
  sprintf(fname, "%s.dat", fname);
  fout.open(fname);
  for (i=0; i<nrow; i++) {
    for (j=rowbegp[i]; j<rowbegp[i+1]; j++) {
      fout << i << " " << colidxp[j] << std::setw(22) << std::setprecision(15)
	   << val[j] << std::endl;
    }
  }
  fout.close();
}

void Epetra_datfile(int* A, int N, char fname[])
{
  int i;
  std::ofstream fout;
  sprintf(fname, "%s.dat", fname);
  fout.open(fname);
  for (i=0; i<N; i++) fout << A[i] << std::endl;
  fout.close();
}

Graph_class::Graph_class(int N_, int A1_[], int A2_[])
{
  N = N_;
  A1 = A1_;
  A2 = A2_;
  component = 0;
  component = new int[N];
  for (int i=0; i<N; i++) component[i] = -1;
}

Graph_class::~Graph_class()
{
  delete [] component; component = 0;
}

void Graph_class::Components(int component_[], int & comp_num)
{
  comp_num = 0;
  for (int i=0; i<N; i++) {
    if (component[i] == -1) {
      DFS(i, comp_num);
      comp_num++;
    }
  }
  for (int i=0; i<N; i++) component_[i] = component[i];
}

void Graph_class::DFS(const int v, const int comp_num)
{
  component[v] = comp_num;
  for (int i=A2[v];i<A2[v+1];i++) {
    int adj_vertex = A1[i];
    if (component[adj_vertex] == -1) DFS(adj_vertex, comp_num);
  }
}

void scale_columns(Epetra_CrsMatrix* A, 
		   const int norm_opt, 
		   const int blocksize)
{
  //
  // scales columns of matrix A according to
  //  norm_opt = 1 (max absolute value of all entries in column is 1)
  //
  int i, j, NumEntries, *Indices;
  double *Values;
  int nrow = A->NumMyRows();
  int ncol = A->NumMyCols();
  double *col_norm = new double[ncol];
  get_column_norm(A, norm_opt, blocksize, col_norm);
  for (i=0; i<nrow; i++) {
    A->ExtractMyRowView(i, NumEntries, Values, Indices);
    for (j=0; j<NumEntries; j++) {
      int col = Indices[j];
      if (col_norm[col] > 0) Values[j] /= col_norm[col];
    }
  }
  delete [] col_norm;
}

void get_column_norm(Epetra_CrsMatrix* A,
		     const int norm_opt,
		     const int blocksize,
		     double *col_norm)
{
  int i, j, NumEntries, *Indices, jpos;
  double *Values;
  int nrow = A->NumMyRows();
  int ncol = A->NumMyCols();
  int ncol_global = A->NumGlobalCols();
  double *vals_block_all = new double[blocksize];
  double *vals_block_loc = new double[blocksize];
  double *vals_col = new double[ncol]; myzero(vals_col, ncol);
  myzero(vals_block_loc, blocksize);
  myzero(vals_block_all, blocksize);
  for (i=0; i<nrow; i++) {
    A->ExtractMyRowView(i, NumEntries, Values, Indices);
    for (j=0; j<NumEntries; j++) {
      if (norm_opt == 1) {
	if (fabs(Values[j]) > vals_col[Indices[j]])
	  vals_col[Indices[j]] = fabs(Values[j]);
      }
    }
  }
  int nmpi = ncol_global/blocksize;
  int nremain = ncol_global - nmpi*blocksize;
  if (nremain > 0) nmpi++;
  for (i=0; i<nmpi; i++) {
    int ibeg = i*blocksize;
    int nsend = blocksize;
    if ((i == (nmpi-1)) && (nremain > 0)) nsend = nremain;
    for (j=0; j<ncol; j++) {
      jpos = A->GCID(j) - ibeg;
      if ((jpos >= 0) && (jpos < nsend)) vals_block_loc[jpos] = vals_col[j];
    }
    if (norm_opt == 1) {
      A->Comm().MaxAll(vals_block_loc, vals_block_all, nsend);
    }
    for (j=0; j<ncol; j++) {
      jpos = A->GCID(j) - ibeg;
      if ((jpos >= 0) && (jpos < nsend)) {
	col_norm[j] = vals_block_all[jpos];
      }
    }
  }
  delete [] vals_block_all; delete [] vals_block_loc; delete [] vals_col;
}

void tie_down_coarse(int n, 
		     int rowbeg[], 
		     int colidx[], 
		     double vals[], 
		     int ne, 
		     int scale_option,
                     int & num_tied_down, 
		     int* & tied_down, 
		     CLAPS_sparse_lu* & AA, 
		     double* & Xvecs)
{
  int i, j, nextra;
  bool *bound_flag = new bool[n];
  for (i=0; i<n; i++) bound_flag[i] = true;
  subspace_iteration(n, rowbeg, colidx, vals, bound_flag, scale_option,
                     num_tied_down, tied_down, AA, Xvecs, ne);
  if (num_tied_down > 0) {
    if (num_tied_down != ne) {
      std::cout << "Error: number of actual rigid body modes less than ";
      std::cout << " num_rigid_mode specified in CLAPS block" << std::endl;
      assert (num_tied_down == ne);
    }
  }
  for (i=0; i<n; i++) bound_flag[i] = false;
  for (i=0; i<num_tied_down; i++) bound_flag[tied_down[i]] = true;
  double min_diag(1e40);
  for (i=0; i<n; i++)
    for (j=rowbeg[i]; j<rowbeg[i+1]; j++)
      if ((colidx[j] == i) && (fabs(vals[j]) < min_diag))
	min_diag = fabs(vals[i]);
  for (i=0; i<n; i++) {
    for (j=rowbeg[i]; j<rowbeg[i+1]; j++) {
      if (colidx[j] == i) {
	if (bound_flag[i] == true) vals[j] = min_diag;
      }
      else
	if ((bound_flag[i] == true) || (bound_flag[colidx[j]] == true))
	  vals[j] = 0;
    }
  }
  delete [] bound_flag; delete AA;
  /*
  subspace_iteration(n, rowbeg, colidx, vals, bound_flag, scale_option, 
                     nextra, tied_down, AKc, Xvecs, ne);
  std::ofstream ffout;
  ffout.open("coarse_mat.dat");
  for (i=0; i<n; i++) 
    for (j=rowbeg[i]; j<rowbeg[i+1]; j++)
      ffout << i+1 << " " << colidx[j]+1 << " " << vals[j] << std::endl;
  ffout.close();
  */
}

void subspace_iteration(int n, int rowbeg[], int colidx[], 
			double vals[], bool bound_flag[], int scale_option,
			int & nextra, int* & extra_corner,
			CLAPS_sparse_lu* & A, double* & Xvecs, int ne)
{
  int i, j, k, nnz, p, q, kbeg;
  double dmin, dmax(-1), dtol(1e-4), sitol(1e-4), etol(1e-4), sum;

  if (n == 0) {
    nextra = 0;
    extra_corner = new int[nextra];
    return;
  }
  //
  // add small multiple of identity to diagonal to handle positive
  // semidefinite matrices
  //
  for (i=0; i<n; i++) 
    for (j=rowbeg[i]; j<rowbeg[i+1]; j++)
      if ((colidx[j] == i) && (fabs(vals[j]) > dmax)) dmax = fabs(vals[j]);
  dmin = dmax;
  for (i=0; i<n; i++) 
    for (j=rowbeg[i]; j<rowbeg[i+1]; j++)
      if ((colidx[j] == i) && (fabs(vals[j]) < dmin)) dmin = fabs(vals[j]);
  for (i=0; i<n; i++)
    for (j=rowbeg[i]; j<rowbeg[i+1]; j++) 
      if (colidx[j] == i) vals[j] += dtol*dmin;
  /*
  if (MyPID == 1) {
    std::ofstream fout;
    fout.open("A.dat");
    fout << setiosflags(ios::scientific | ios::uppercase);
    for (i=0; i<n; i++) {
      for (j=rowbeg[i]; j<rowbeg[i+1]; j++) {
	fout << i+1 << " " << colidx[j]+1 << " " << std::setw(22) <<
	  std::setprecision(15) << vals[j] << std::endl;
      }
    }
    fout.close();
  }
  */
  //
  // factor matrix
  //
  A = new CLAPS_sparse_lu();
  nnz = rowbeg[n];
  A->factor(n, nnz, rowbeg, colidx, vals, scale_option);
  //
  // adjust vals array back to original values
  //
  for (i=0; i<n; i++)
    for (j=rowbeg[i]; j<rowbeg[i+1]; j++) 
      if (colidx[j] == i) vals[j] -= dtol*dmin;
  //
  // subspace iteration stuff
  //
  int maxiter_si(20);
  p = 20;
  if (ne > 0) p = ne + 1;
  if (p > n) p = n;
  q = 2*p; if (2*p > (p+8)) q = p+8; if (q > n) q = n;
  double *X = new double[n*q]; myzero(X, n*q);
  double *SOL = new double[n*q];
  double *TEMP = new double[n*q];
  double *LAMBDA = new double[q];
  Epetra_BLAS EB;
  Epetra_LAPACK EL;
  char TRANSA = 'T'; char TRANSB = 'N'; double ALPHA(1); double BETA(0);
  char JOBZ = 'V'; char UPLO = 'U';
  int ITYPE(1), LWORK, INFO, nj;
  LWORK = 3*q;
  double *WORK = new double[LWORK];
  //
  // initialize first column of X to 1 and remaining columns to random values
  //
  for (i=0; i<n; i++) X[i] = 1;
  for (j=1; j<q; j++) {
    nj = n*j;
    for (i=0; i<n; i++) X[nj+i] = 0.997*rand()/RAND_MAX;
  }
  //
  // subspace iterations
  //
  for (int iter=0; iter<maxiter_si; iter++) {
    A->sol(q, X, SOL, TEMP);
    //
    // SOL^T * A * SOL and SOL^T * SOL
    //
    EB.GEMM(TRANSA, TRANSB, q, q, n, ALPHA, SOL, n,   X, n, BETA, TEMP, q);
    EB.GEMM(TRANSA, TRANSB, q, q, n, ALPHA, SOL, n, SOL, n, BETA,    X, q);   
    /*
    if (MyPID == 1) {
      std::ofstream fout;
      fout.open("AA.m");
      fout << "TEMP = zeros(" << q << "," << q << ");" << std::endl;
      fout << "X    = zeros(" << q << "," << q << ");" << std::endl;
      fout << setiosflags(ios::scientific | ios::uppercase);
      for (i=0; i<q; i++) {
	for (j=0; j<q; j++) {
	  fout << "TEMP(" << i+1 << "," << j+1 << ") = " << std::setw(22) <<
	    std::setprecision(15) <<  TEMP[i+q*j] << ";" << std::endl;
	  fout << "X(   " << i+1 << "," << j+1 << ") = " << std::setw(22) <<
	    std::setprecision(15) <<     X[i+q*j] << ";" << std::endl;
	}
      }
      fout.close();
    }
    */
    EL.SYGV(ITYPE, JOBZ, UPLO, q, TEMP, q, X, q, LAMBDA, WORK, LWORK, &INFO);
    //    assert (INFO == 0);
    if (INFO != 0) {
      delete [] SOL; delete [] TEMP; delete [] WORK; delete [] X;
      delete [] LAMBDA;
      nextra = 0;
      extra_corner = new int[nextra];
      return;
    }
    EB.GEMM(TRANSB, TRANSB, n, q, q, ALPHA, SOL, n, TEMP, q, BETA, X, n);
    sum = 0;
    for (i=0; i<q; i++) sum += TEMP[q*(p-1)+i]*TEMP[q*(p-1)+i];
    sum = sqrt(1-LAMBDA[p-1]*LAMBDA[p-1]/sum);
    /*
    if (MyPID == 0) {
      std::cout << "LAMBDA = " << std::endl;
      for (i=0; i<q; i++) std::cout << LAMBDA[i] << std::endl;
      std::cout << "eigenerror = " << sum << std::endl;
    }
    */
    if (sum <= sitol) break;
  }
  delete [] SOL; delete [] TEMP; delete [] WORK;
  for (i=0; i<q; i++) LAMBDA[i] -= dtol*dmin;
  //
  // determine additional corners needed to remove singularities using
  // Gaussian elimination of X^T with column pivoting
  //
  /*
  if (MyPID == 0) {
    std::cout << "LAMBDA = " << std::endl;
    for (i=0; i<q; i++) std::cout << LAMBDA[i] << std::endl;
  }
  */
  nextra = 0;
  if (ne == 0) {
    for (i=0; i<q; i++) if (fabs(LAMBDA[i]) <= dtol*dmin) nextra++;
  }
  if (ne > 0) {
    double aaa = dtol*fabs(LAMBDA[ne]);
    for (i=0; i<ne; i++) if (fabs(LAMBDA[i]) <= aaa) nextra++;
    Xvecs = new double[n*nextra];
    for (i=0; i<(n*nextra); i++) Xvecs[i] = X[i];
  }
  /*  
  if (ne > 0) {
    for (i=0; i<q; i++) std::cout << "lambda[" << i << "]= " << LAMBDA[i] << std::endl;
    std::ofstream ffout;
    ffout.open("coarse_vec.dat");
    for (i=0; i<n; i++) {
      for (j=0; j<ne; j++) ffout << X[i+j*n] << " ";
      ffout << std::endl;
    }
    ffout.close();
  }
  */
  extra_corner = new int[nextra];
  double maxval, sfac;
  int col, ibeg;
  for (i=0; i<nextra; i++) {
    ibeg = i*n;
    maxval = 0;
    for (j=0; j<n; j++) {
      if ((fabs(X[ibeg+j]) > maxval) && (bound_flag[j] == true)) {
	maxval = fabs(X[ibeg+j]);
	col = j;
      }
    }
    assert (maxval > 0);
    extra_corner[i] = col;
    sfac = 1/X[ibeg+col];
    for (j=0; j<n; j++) X[ibeg+j] *= sfac;
    for (k=i+1; k<nextra; k++) {
      kbeg = k*n;
      for (j=0; j<n; j++) X[kbeg+j] -= X[kbeg+col]*X[ibeg+j];
    }
  }
  /*
  if (MyPID == 0) {
    std::cout << "nextra = " << nextra << std::endl;
    std::cout << "extra constained dofs = :";
    for (i=0; i<nextra; i++) std::cout << extra_corner[i] << " ";
    std::cout << std::endl;
    for (i=0; i<n; i++) {
      for (j=0; j<nextra; j++) std::cout << X[j*n+i] << " ";
      std::cout << std::endl;
    }
  }
  */
  delete [] X; delete [] LAMBDA;
}


} // end of namespace CRD_utils
