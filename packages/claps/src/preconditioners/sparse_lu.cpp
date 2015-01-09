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

#include "sparse_lu.hpp"
#include "my_feti_sparse_solver.hpp"
#include <metis.h>
#include <vector>

CLAPS_sparse_lu::CLAPS_sparse_lu()
{
  XSUPER = 0; XLINDX = 0; LINDX = 0; XLNZ = 0; PERM = 0; INVP = 0;
  IPROW = 0;  IPCOL = 0; DEF = 0; LNZ = 0; NS = 0; SCALE = 0;
}

CLAPS_sparse_lu::~CLAPS_sparse_lu()
{
  delete [] XSUPER; delete [] XLINDX; delete [] LINDX; delete [] XLNZ; 
  delete [] PERM; delete [] INVP; delete [] IPROW; delete [] IPCOL;
  delete [] DEF; delete [] LNZ; delete [] NS; delete [] SCALE;
}

void CLAPS_sparse_lu::cleanup()
{
  delete [] XSUPER; delete [] XLINDX; delete [] LINDX; delete [] XLNZ;
  delete [] PERM; delete [] INVP; delete [] IPROW; delete [] IPCOL;
  delete [] DEF; delete [] LNZ; delete [] NS;
}

int CLAPS_sparse_lu::factor(int N_, int NNZ, int COLPTR[], int ROWIDX[], 
		      double ANZ[], int scale_flag_)
{
  //
  // Input:
  //   N_ = number of equations
  //   NNZ = number of nonzeros in full matrix
  //   ROWIDX[COLPTR[i]:COLPTR[i+1]-1]  = nonzero row numbers for column i
  //   ANZ[   COLPTR[i]:COLPTR[i+1]-1]  = nonzero values in column i
  //   Note: C-style numbering of inputs is assumed (i.e. row and column
  //         numbers are all between 0 and N_-1
  //
  // solver parameters
  //
  DEFBLK = 1;
  int order_opt(2);
  max_small = 4;
  //
  N = N_;
  if (N == 0) return 0;
  if (N <= max_small) {
    int INFO;
    INFO = small_factor(COLPTR, ROWIDX, ANZ);
    return INFO;
  }
  scale_flag = scale_flag_;
  //
  // scale matrix entries if requested
  //
  double *ANZ_SCALED;
  if (scale_flag == 1) {
    int i, j;
    SCALE = new double[N];
    for (i=0; i<N; i++) SCALE[i] = 0;
    for (i=0; i<N; i++) {
      for (j=COLPTR[i]; j<COLPTR[i+1]; j++)
	if (ROWIDX[j] == i) SCALE[i] = sqrt(fabs(ANZ[j]));
      assert (SCALE[i] != 0);
    }
    ANZ_SCALED = new double[NNZ];
    for (i=0; i<N; i++) SCALE[i] = 1/SCALE[i];
    for (i=0; i<N; i++) {
      for (j=COLPTR[i]; j<COLPTR[i+1]; j++) {
	ANZ_SCALED[j] = SCALE[i]*SCALE[ROWIDX[j]]*ANZ[j];
      }
    }
  }
  else {
    ANZ_SCALED = ANZ;
  }
  int NNZA, NADJ, IWMAX, IWSIZE, IFLAG, MAXSUP, NTOT, RWSIZE, LDNS;
  int NNZL, NSUB, NLNZ, TMPSIZ, MAXDEF, ASDEF;
  double ANORM, EPS, TOL;
  int OPTIONS[8];
  //
  OPTIONS[0]=0; // use default values for options
  OPTIONS[1]=0;
  MAXSUP=150;   // maximum supernode size
  NTOT=NNZ;     // total number of nonzeros in matrix
  NADJ=NTOT-N;  // number of nonzeros in full matrix minus those on diagonal
  NNZA=NADJ/2+N;// number of nonzeros in lower triangle including diagonal
  IWMAX=7*N+3;  // dimension for integer working array
  if ((MAXSUP+2*N+1) > IWMAX) IWMAX=MAXSUP+2*N+1;
  if ((3*N+2*MAXSUP) > IWMAX) IWMAX=3*N+2*MAXSUP;
  
  MAXDEF=10;
  ASDEF=0;
  //  std::cout << "N    = " << N << std::endl;
  //  std::cout << "NTOT = " << NTOT << std::endl;
  //  std::cout << "NADJ = " << NADJ << std::endl;
  //  std::cout << "NNZA = " << NNZA << std::endl;
  int* ADJ = new int[NTOT];
  int* XADJ = new int[N+1];
  LINDX = new int[NTOT];
  PERM = new int[N];
  INVP = new int[N];
  int* IWORK = new int[IWMAX];
  int* COLCNT = new int[N];
  int* SNODE = new int[N];
  XSUPER = new int[N+1];
  XLINDX = new int[N+1];
  XLNZ = new int[N+1];
  DEF = new int[N];
  IPROW = new int[N];
  IPCOL = new int[N];
  //
  // determine adjacency structure of matrix
  // ADJ[XADJ[i]:XADJ[i+1]-1] = dofs adjacent to dof i (but not including 
  //                            dof i itself)
  //
  XADJ[0]=0;
  int nnzADJ=0;
  for (int i=0;i<N;i++) {
    for (int j=COLPTR[i];j<COLPTR[i+1];j++) {
      int row=ROWIDX[j];
      if (row != i) {
	ADJ[nnzADJ]=row;
	nnzADJ++;
      }
    }
    XADJ[i+1]=nnzADJ;
  }
  //
  // compute L-infinity norm of matrix
  //
  getnrm(N, COLPTR, ROWIDX, ANZ_SCALED, ANORM);
  //
  // convert COLPTR, ROWIDX, XADJ, and ADJ to Fortran numbering
  //
  for (int i=0;i<=N;i++) COLPTR[i]++;
  for (int i=0;i<NNZ;i++) ROWIDX[i]++;
  for (int i=0;i<=N;i++) XADJ[i]++;
  for (int i=0;i<nnzADJ;i++) ADJ[i]++;
  //
  NADJ=XADJ[N]-1;
  for (int i=0;i<=N;i++) XLINDX[i]=XADJ[i];
  for (int i=0;i<nnzADJ;i++) LINDX[i]=ADJ[i];
  //
  // multiple minimum degree ordering
  //
  IWSIZE=4*N;
  if (order_opt == 1) {
    ORDMMD2_F77(N,XLINDX,LINDX,INVP,PERM,IWSIZE,IWORK,NSUB,IFLAG);
    if (IFLAG != 0) {
      std::cout << "error in call to ordmmd2 in CLAPS_sparse_lu::factor" << std::endl;
      std::cout << "ORDMMD2 IFLAG=" << IFLAG << std::endl;
      return -1;
    }
  }
  //
  // Metis ordering
  //
  if (order_opt == 2) {
    int numflag=1;
    std::vector<idx_t> options(METIS_NOPTIONS,0);
    METIS_SetDefaultOptions(&options[0]);

    options[METIS_OPTION_NUMBERING] = numflag;

    METIS_NodeND(&N,XLINDX,LINDX,0,&options[0],PERM,INVP);
  }
  //
  // symbolic factorization initialization
  //
  IWSIZE=7*N+3;
  SFINIT_F77(N,NADJ,XADJ,ADJ,PERM,INVP,MAXSUP,DEFBLK,COLCNT,
	     NNZL,NSUB,NSUPER,XSUPER,SNODE,IWSIZE,IWORK,IFLAG);
  if (IFLAG != 0) {
    std::cout << "error in call to sfinit in CLAPS_sparse_lu::factor" << std::endl;
    std::cout << "SFINIT IFLAG=" << IFLAG << std::endl;
    return -1;
  }
  //
  // supernodal symbolic factorization
  //
  IWSIZE=NSUPER+2*N+1;
  if (NSUB>NTOT) {
    delete [] LINDX;
    LINDX = new int[NSUB];
  }
  SYMFCT_F77(N,NADJ,XADJ,ADJ,PERM,INVP,COLCNT,NSUPER,XSUPER,
	     SNODE,NSUB,XLINDX,LINDX,XLNZ,IWSIZE,IWORK,IFLAG);
  if (IFLAG != 0) {
    std::cout << "error in call to symfct in CLAPS_sparse_lu::factor" << std::endl;
    std::cout << "SYMFCT IFLAG=" << IFLAG << std::endl;
    return -1;
  }
  //
  // input numerical values into data structures
  //
  NLNZ=XLNZ[N];
  //  std::cout << "number of nonzeros in LU factorization = " << NLNZ << std::endl;
  //  std::cout << "NLNZ = " << NLNZ << std::endl;

  // later, sparsepak will call dgemm on some of the "panels" in this data. We need to have
  // memory on the end of the array to ensure we don't overrun. Without the extra "N", we
  // may not have all of the last column allocated.
  LNZ = new double[NLNZ+N];
  for (int i=0;i<NLNZ;i++) LNZ[i]=0;
  for (int i=0;i<N;i++) DEF[i]=0;
  NDEF=0;
  LBDEF=0;
  inpnv(N,COLPTR,ROWIDX,ANZ_SCALED,PERM,INVP,NSUPER,XSUPER,XLINDX,LINDX,XLNZ,
	LNZ,IWORK);
  if (scale_flag == 1) delete [] ANZ_SCALED;
  //
  // numerical factorization
  //
  BFINIT_F77(NSUPER,XSUPER,SNODE,XLINDX,LINDX,TMPSIZ,RWSIZE);
  if (TMPSIZ < 1) TMPSIZ=1;
  TMPSIZ=2*TMPSIZ;
  double* TMPVEC = new double[TMPSIZ];
  int RWORKdim=N;
  if (RWSIZE > N) RWORKdim=RWSIZE;
  double* RWORK = new double[RWORKdim];
  //  EPS=1e-10;
  EPS = 1e-12;
  //  EPS=DLAMCH('EPS');
  TOL=EPS*ANORM;
  IWSIZE = 3*N + 2*NSUPER;
  BLKLDL_F77(NSUPER,XSUPER,SNODE,XLINDX,LINDX,XLNZ,LNZ,DEFBLK,ASDEF,NDEF,
	     LBDEF,DEF,TOL,IPROW,IPCOL,TMPSIZ,TMPVEC,IWSIZE,IWORK,
	     RWSIZE,RWORK,IFLAG);
  if (IFLAG != 0) {
    std::cout << "error in call to blkldl in CLAPS_sparse_lu::factor" << std::endl;
    std::cout << "BLKLDL IFLAG=" << IFLAG << std::endl;
    return -1;
  }
  //  if (DEFBLK == 0) {
  //    LBDEF=0;
  //    NDEF=0;
  //  }
  if ((NDEF != 0) && (NDEF <= MAXDEF)) {
    //
    // compute null space
    //
    NS = new double [N*NDEF];
    LDNS=N;
    BLKNS_F77(NSUPER,XSUPER,XLINDX,LINDX,XLNZ,LNZ,DEFBLK,NDEF,LBDEF,
	      DEF,IPCOL,INVP,NS,LDNS,RWORK);
    std::cout << "null space dimension = " << NDEF << std::endl;
    /*
    std::cout << "NS = " << std::endl;
    for (int i=0;i<NDEF;i++) {
      for (int j=0;j<N;j++) std::cout << NS[j+N*i] << " ";
      std::cout << std::endl;
    }
    */
  }
  delete [] ADJ;ADJ=0;
  delete [] XADJ;XADJ=0;
  delete [] IWORK;IWORK=0;
  delete [] COLCNT;COLCNT=0;
  delete [] SNODE;SNODE=0;
  delete [] TMPVEC;TMPVEC=0;
  delete [] RWORK;RWORK=0;
  //
  // convert COLPTR and ROWIDX back to C numbering
  //
  for (int i=0;i<=N;i++) COLPTR[i]--;
  for (int i=0;i<NNZ;i++) ROWIDX[i]--;
  return 0;
}

int CLAPS_sparse_lu::sol(int NRHS, double RHS[], double SOL[], double TEMP[])
{
  //
  // numerical solution
  //
  if (N == 0 ) return 0;
  if (N <= max_small) {
    int INFO;
    INFO = small_solve(NRHS, RHS, SOL);
    return INFO;
  }
  if (scale_flag == 1) {
    for (int j=0; j<NRHS; j++) {
      int jbeg = j*N;
      for (int i=0; i<N; i++) RHS[jbeg+i] *= SCALE[i];
    }
  }
  int LRHS=N;
  int LSOL=N;
  //  std::cout << "LBDEF = " << LBDEF << std::endl;
  //  std::cout << "NDEF  = " << NDEF << std::endl;
  BLKSLVN_F77(NSUPER,XSUPER,XLINDX,LINDX,XLNZ,LNZ,DEFBLK,NDEF,LBDEF,
	      DEF,IPROW,IPCOL,PERM,INVP,LRHS,NRHS,RHS,LSOL,SOL,N,TEMP);
  if (scale_flag == 1) {
    for (int j=0; j<NRHS; j++) {
      int jbeg = j*N;
      for (int i=0; i<N; i++) {
	SOL[jbeg+i] *= SCALE[i];
	RHS[jbeg+i] /= SCALE[i];
      }
    }
  }
  return 0;
}

void CLAPS_sparse_lu::getnrm(int n, int colptr[], int rowidx[], 
		       double values[], double &anorm)
{
  int i,j;
  double t;
  anorm  = 0;
  for (i=0;i<n;i++) {
    t=0;
    for (j=colptr[i];j<colptr[i+1];j++) t += fabs(values[j]);
    if (t > anorm) anorm=t;
  }
}

void CLAPS_sparse_lu::inpnv(int &n , int colptr[], int rowidx[], 
                 double values[], int perm[], int invp [], int &nsuper, 
		 int xsuper[], int xlindx[], int lindx[], int xlnz[], 
                 double lnz[], int offset[])
{
  //
  // input numerical values for cholesky factorization
  //
  int ii,lxbeg,fstcol,jsuper,lxend,jlen,irow,lstcol,jcol,oldj,lastl;
  for (ii=1;ii<xlnz[n];ii++) lnz[ii-1]=0;
  lxbeg  = xlindx[0];
  fstcol = xsuper[0];
  for (jsuper=1;jsuper<=nsuper;jsuper++) {
    lxend = xlindx[jsuper];
    jlen  = lxend - lxbeg;
    for (ii=lxbeg;ii<lxend;ii++) {
      irow = lindx[ii-1];
      jlen--;
      offset[irow-1]=jlen;
    }
    lstcol=xsuper[jsuper];
    for (jcol=fstcol;jcol<lstcol;jcol++) {
      oldj = perm[jcol-1];
      lastl = xlnz[jcol]-1;
      for (ii=colptr[oldj-1];ii<colptr[oldj];ii++) {
	irow=invp[rowidx[ii-1]-1];
	if (irow >= fstcol) {
	  lnz[lastl-offset[irow-1]-1]=values[ii-1];
	}
      }
    }
    lxbeg  = lxend;
    fstcol = lstcol;
  }
}

int CLAPS_sparse_lu::small_factor(int rowbeg[], int colidx[], double vals[])
{
  int i, j, INFO, col;
  LNZ = new double[N*N]; for (i=0; i<N*N; i++) LNZ[i] = 0;
  XSUPER = new int[N];
  for (i=0; i<N; i++) {
    for (j=rowbeg[i]; j<rowbeg[i+1]; j++) {
      col = colidx[j];
      LNZ[N*col+i] = vals[j];
    }
  }
  EL.GETRF(N, N, LNZ, N, XSUPER, &INFO);
  return INFO;
}

int CLAPS_sparse_lu::small_solve(int NRHS, double RHS[], double SOL[])
{
  int INFO;
  char TRANS = 'N';
  memcpy(SOL, RHS, N*NRHS*sizeof(double));
  EL.GETRS(TRANS, N, NRHS, LNZ, N, XSUPER, SOL, N, &INFO);
  return INFO;
}
