#include "CLOP_sub.hpp"
#include "CRD_utils.hpp"
//#include "../include/sort_prototypes.h"
#include <assert.h>
extern "C"{
  void dgelss_(int* M, int* N, int* NRHS, double A[], int* LDA, double B[],
	       int* LDB, double S[], double* RCOND, int* RANK, double WORK[], 
	       int* LWORK, int* INFO);
  void dgeqpf_(int* M, int* N, double A[], int* LDA, int JPVT[],
	       double TAU[], double WORK[], int* INFO);
  void dsyev_(char* JOBZ, char* UPLO, int* N, double A[], int* LDA, double W[],
	      double WORK[], int* LWORK, int* INFO, long lengthA, 
	      long lengthB); 
}

CLOP_sub::CLOP_sub()
{
  A_sub = 0; Edof = 0; Phi = 0; jpvt = 0; 
}

CLOP_sub::~CLOP_sub()
{
  if (A_sub) delete A_sub;
  if (Edof) delete [] Edof;
  if (Phi) delete [] Phi;
  if (jpvt) delete [] jpvt;
}

void CLOP_sub::getmatrix_nnz(int subdofs_[], int ndof_sub, 
   const Epetra_CrsMatrix *A, int imap[], const Epetra_Comm* Comm_, int &nnz)
{
  int i, j, NumEntries, *Indices;
  double *Values;
  subdofs = subdofs_;
  ndof = ndof_sub;
  Comm = Comm_;
  MyPID = Comm->MyPID();
  for (i=0; i<ndof_sub; i++) imap[subdofs[i]] = i;
  nnz = 0;
  for (i=0; i<ndof_sub; i++) {
    A->ExtractMyRowView(subdofs[i], NumEntries, Values, Indices);
    for (j=0; j<NumEntries; j++) {
      if (imap[Indices[j]] > -1) nnz++;
    }
  }
 for (i=0; i<ndof_sub; i++) imap[subdofs[i]] = -1;
}

void CLOP_sub::factormatrix(const Epetra_CrsMatrix *A, int imap[], 
			    int rowbeg[], int colidx[], double K[])
{
  int i, j, NumEntries, *Indices, nnz;
  double *Values;
  for (i=0; i<ndof; i++) imap[subdofs[i]] = i;
  rowbeg[0] = 0;
  nnz = 0;
  for (i=0; i<ndof; i++) {
    A->ExtractMyRowView(subdofs[i], NumEntries, Values, Indices);
    for (j=0; j<NumEntries; j++) {
      if (imap[Indices[j]] > -1) nnz++;
    }
    rowbeg[i+1] = nnz;
  }
  nnz = 0;
  for (i=0; i<ndof; i++) {
    A->ExtractMyRowView(subdofs[i], NumEntries, Values, Indices);
    for (j=0; j<NumEntries; j++) {
      if (imap[Indices[j]] > -1) {
	colidx[nnz] = imap[Indices[j]];
	K[nnz] = Values[j];
	nnz++;
      }
    }
  }
  /*
  char filename[101];
  sprintf(filename,"%s%d","submat", MyPID);
  CRD_utils::spmat_datfile(ndof, rowbeg, colidx, K, filename);
  */
  A_sub = new sparse_lu();
  A_sub->factor(ndof, nnz, rowbeg, colidx, K);
  myzero(K, nnz);
  for (i=0; i<ndof; i++) imap[subdofs[i]] = -1;
}

void CLOP_sub::genpu(const Epetra_IntVector *LD, 
		     const Epetra_MultiVector *Coords, double rhs[], 
		     double sol[], double temp[], int atype_sub, int ndim_sub,
		     double WORK[], int LWORK, double Edof_sub[],
		     int & nneg)
{
  int i, j, k, nrhs, nnz_rhs, dof;
  double *coords;
  atype = atype_sub;
  ndim = ndim_sub;
  LD->ExtractView(&locdof);
  int ndof_proc = Coords->Stride();
  Coords->ExtractView(&coords, &ndof_proc);
  x = &coords[0];
  y = &coords[ndof_proc];
  z = &coords[2*ndof_proc];
  xcent = 0; ycent = 0; zcent = 0;
  for (i=0; i<ndof; i++) {
    dof = subdofs[i]; xcent += x[dof]; ycent += y[dof]; zcent += z[dof];
  }
  xcent /= ndof; ycent /= ndof; zcent /= ndof;
  if (atype == 1) {
    csdim_max = 1;
    nrhs = ndim + 1;
    nnz_rhs = ndof*nrhs;
    myzero(rhs, nnz_rhs);
    if (ndim == 2) {
      for (i=0; i<ndof; i++) {
	dof = subdofs[i];
	rhs[i+0*ndof] = 1;
	rhs[i+1*ndof] = x[dof] - xcent;
	rhs[i+2*ndof] = y[dof] - ycent;
      }
    }
    if (ndim == 3) {
      for (i=0; i<ndof; i++) {
	dof = subdofs[i];
	rhs[i+0*ndof] = 1;
	rhs[i+1*ndof] = x[dof] - xcent;
	rhs[i+2*ndof] = y[dof] - ycent;
	rhs[i+3*ndof] = z[dof] - zcent;
      }
    }
    //
    // check to make sure coordinate data is meaningful
    //
    int coord_flag = 0;
    for (i=0; i<ndof; i++) {
      dof = subdofs[i];
      if (x[dof] != 0) {
	coord_flag = 1;
        break;
      }
    }
    if (coord_flag == 0) nrhs = 1;
  }
  if (atype == 2) {
    if (ndim == 2) {
      csdim_max = 3;
      nrhs = 6;
      nnz_rhs = ndof*nrhs;
      myzero(rhs, nnz_rhs);
      for (i=0; i<ndof; i++) {
	dof = subdofs[i];
	if (abs(locdof[dof]) == 1) {
	  rhs[i+0*ndof] = 1;
	  rhs[i+1*ndof] = x[dof] - xcent;
	  rhs[i+2*ndof] = y[dof] - ycent;
	}
	if (abs(locdof[dof]) == 2) {
	  rhs[i+3*ndof] = 1;
	  rhs[i+4*ndof] = x[dof] - xcent;
	  rhs[i+5*ndof] = y[dof] - ycent;
	}
      }
    }
    if (ndim == 3) {
      csdim_max = 6;
      nrhs = 12;
      nnz_rhs = ndof*nrhs;
      if (ndof < nrhs) nnz_rhs = nrhs*nrhs;
      myzero(rhs, nnz_rhs);
      for (i=0; i<ndof; i++) {
	dof = subdofs[i];
	if (abs(locdof[dof]) == 1) {
	  rhs[i+0*ndof] = 1;
	  rhs[i+1*ndof] = x[dof] - xcent;
	  rhs[i+2*ndof] = y[dof] - ycent;
	  rhs[i+3*ndof] = z[dof] - zcent;
	}
	if (abs(locdof[dof]) == 2) {
	  rhs[i+4*ndof] = 1;
	  rhs[i+5*ndof] = x[dof] - xcent;
	  rhs[i+6*ndof] = y[dof] - ycent;
	  rhs[i+7*ndof] = z[dof] - zcent;
	}
	if (abs(locdof[dof]) == 3) {
	  rhs[i+8*ndof]  = 1;
	  rhs[i+9*ndof]  = x[dof] - xcent;
	  rhs[i+10*ndof] = y[dof] - ycent;
	  rhs[i+11*ndof] = z[dof] - zcent;
	}
      }
    }
  }
  if (atype == 3) {
    csdim_max = 3;
    nrhs = 9;
    nnz_rhs = ndof*nrhs;
    if (ndof < nrhs) nnz_rhs = nrhs*nrhs;
    myzero(rhs, nnz_rhs);
    for (i=0; i<ndof; i++) {
      dof = subdofs[i];
      if (abs(locdof[dof]) == 1) {
	rhs[i+0*ndof] = 1;
	rhs[i+1*ndof] = x[dof] - xcent;
	rhs[i+2*ndof] = y[dof] - ycent;
      }
      if (abs(locdof[dof]) == 2) {
	rhs[i+3*ndof] = 1;
	rhs[i+4*ndof] = x[dof] - xcent;
	rhs[i+5*ndof] = y[dof] - ycent;
      }
      if (abs(locdof[dof]) == 3) {
	rhs[i+6*ndof] = 1;
	rhs[i+7*ndof] = x[dof] - xcent;
	rhs[i+8*ndof] = y[dof] - ycent;
      }
    }
  }
  A_sub->sol(nrhs, rhs, sol, temp);
  //
  // form SOL' * K * SOL = SOL' * RHS
  //
  char TRANSA = 'T'; char TRANSB = 'N';
  double ALPHA = 1, BETA = 0;
  Epetra_BLAS EB;
  if (ndof == 0) {
    Edof = new double[0];
    int duma1(0), duma2;
    Comm->MaxAll(&duma1, &duma2, 1);
    return;
  }
  EB.GEMM(TRANSA, TRANSB, nrhs, nrhs, ndof, ALPHA,
	  sol, ndof, rhs, ndof, BETA, temp, nrhs);
  //
  // eigenvalues of SOL' * K * SOL are positive if K is positive, but
  // may be negative if K has one or more negative eigenvalues. check
  // the signs of eigenvalues of SOL' * K * SOL and make zero if
  // any are found to be negative
  //
  char JOBZ('V'), UPLO('U');
  double S[60];
  for (i=0; i<nrhs*nrhs; i++) rhs[i] = temp[i];
  dsyev_(&JOBZ, &UPLO, &nrhs, temp, &nrhs, S, WORK, &LWORK, &INFO, 1, 1);
  assert(INFO == 0);
  nneg = 0;
  double max_mag(0), tol_neg_eig(1e-8);
  for (i=0; i<nrhs; i++) if (fabs(S[i]) > max_mag) max_mag = fabs(S[i]);
  for (i=0; i<nrhs; i++) if (S[i] < -tol_neg_eig*max_mag) nneg++;
  if (nneg > 0) {
    myzero(rhs, nrhs*nrhs);
    for (i=0; i<nrhs; i++) {
      if (S[i] > tol_neg_eig*max_mag) {
	for (j=0; j<nrhs; j++) {
	  for (k=0; k<nrhs; k++) rhs[j+nrhs*k] += 
				   S[i]*temp[i*nrhs+j]*temp[i*nrhs+k];
	}
      }
    }
  }
  for (i=0; i<nrhs*nrhs; i++) temp[i] = rhs[i];
  //
  // calculate pseudo-inverse
  //
  double RCOND = 1e-8;
  int RANK;
  double rhs_sol[60], sol_sol[60];
  myzero(rhs, nrhs*nrhs);
  for (i=0; i<nrhs; i++) rhs[i*(nrhs+1)] = 1;
  dgelss_(&nrhs, &nrhs, &nrhs, temp, &nrhs, rhs, &nrhs, S, &RCOND, &RANK,
	  WORK, &LWORK, &INFO);
  assert(INFO == 0);
  //
  // calculate diagonal of flexibility matrix
  //
  Edof = new double[ndof];
  TRANSA = 'N';
  int INCX = 1, INCY = 1;
  for (i=0; i<ndof; i++) {
    for (j=0; j<nrhs; j++) sol_sol[j] = sol[i+j*ndof];
    EB.GEMV(TRANSA, nrhs, nrhs, ALPHA, rhs, nrhs,
	    sol_sol, BETA, rhs_sol);
    Edof[i] = EB.DOT(nrhs, sol_sol, rhs_sol);
    assert(Edof[i] >= 0);
  }
  for (i=0; i<ndof; i++) Edof_sub[subdofs[i]] += Edof[i];
  /*
  if (MyPID == 0) {
    cout << "xcent, ycent = " << xcent << " " << ycent << endl;
    cout << "singular values = " << endl;
    for (i=0; i<nrhs; i++) cout << S[i] << endl;
    //    cout << "x, y, Edof_orig" << endl;
    //    for (i=0; i<ndof; i++) {
    //      cout << x[dofa[i]] << " " << y[dofa[i]] << " " << Edof[i] << endl;
    //    }
  }
  */
}

void CLOP_sub::normalpu(double Edof_sub[], unsigned char nsubdof[])
{
  for (int i=0; i<ndof; i++) {
    if (Edof_sub[subdofs[i]] == 0)
      Edof[i] = double(1)/double(nsubdof[subdofs[i]]);
    else Edof[i] /= Edof_sub[subdofs[i]];
  }
}

void CLOP_sub::construct_coarse1(const Epetra_CrsMatrix *A, double rhs[], 
     double sol[], double temp[], int rowbeg[], int colidx[], double K[],
     int imap[], unsigned char nsubdof[], int & csdimP, int & ndof_rot_)
{
  int i, nnz_Phi, dof;
  nnz_Phi = ndof*csdim_max;
  Phi = new double[nnz_Phi];
  myzero(Phi, nnz_Phi);
  if (atype == 1) for (i=0; i<ndof; i++) Phi[i] = Edof[i];
  if (atype == 2) {
    if (ndim == 2) {
      for (i=0; i<ndof; i++) {
	dof = subdofs[i];
	if (abs(locdof[dof]) == 1) {
	  Phi[i+0*ndof] =  Edof[i];
	  Phi[i+2*ndof] = -Edof[i]*(y[dof] - ycent);
	}
	if (abs(locdof[dof]) == 2) {
	  Phi[i+1*ndof] =  Edof[i];
	  Phi[i+2*ndof] =  Edof[i]*(x[dof] - xcent);
	}
      }
    }
    if (ndim == 3) {
      for (i=0; i<ndof; i++) {
	dof = subdofs[i];
	if (abs(locdof[dof]) == 1) {
	  Phi[i+0*ndof] =  Edof[i];
	  Phi[i+4*ndof] =  Edof[i]*(z[dof] - zcent);
	  Phi[i+5*ndof] = -Edof[i]*(y[dof] - ycent);
	}
	if (abs(locdof[dof]) == 2) {
	  Phi[i+1*ndof] =  Edof[i];
	  Phi[i+5*ndof] =  Edof[i]*(x[dof] - xcent);
	  Phi[i+3*ndof] = -Edof[i]*(z[dof] - zcent);
	}
	if (abs(locdof[dof]) == 3) {
	  Phi[i+2*ndof] =  Edof[i];
	  Phi[i+3*ndof] =  Edof[i]*(y[dof] - ycent);
	  Phi[i+4*ndof] = -Edof[i]*(x[dof] - xcent);
	}
	if (abs(locdof[dof]) == 4) Phi[i+3*ndof] = Edof[i];
	if (abs(locdof[dof]) == 5) Phi[i+4*ndof] = Edof[i];
	if (abs(locdof[dof]) == 6) Phi[i+5*ndof] = Edof[i];
      }
    }
  }      
  if (atype == 3) {
    for (i=0; i<ndof; i++) {
      dof = subdofs[i];
      if (abs(locdof[dof]) == 1) {
	Phi[i+0*ndof] =  Edof[i];
	Phi[i+1*ndof] =  Edof[i]*(y[dof] - ycent);
	Phi[i+2*ndof] = -Edof[i]*(x[dof] - xcent);
      }
      if (abs(locdof[dof]) == 2) Phi[i+1*ndof] =  Edof[i];
      if (abs(locdof[dof]) == 3) Phi[i+2*ndof] =  Edof[i];
    }
  }
  //
  // use QR decomposition of Phi with column pivoting to determine
  // a linearly independent set of columns of Phi (Note: csdim is
  // reduced if needed)
  //
  csdim = csdim_max;
  jpvt = new int[csdim];
  double *tau  = new double[csdim];
  double *dwork = new double[3*csdim];
  memcpy(temp, Phi, nnz_Phi*sizeof(double));
  myzero(jpvt, csdim);
  if (ndof > 0) dgeqpf_(&ndof, &csdim, temp, &ndof, jpvt, tau, dwork, &INFO);
  else INFO = 0;
  assert(INFO == 0);
  if (csdim_max > ndof) csdim_max = ndof;
  csdim = 0;
  for (i=0; i<csdim_max; i++) {
    if (fabs(temp[i*(ndof+1)]) <= 1.0e-8*fabs(temp[0])) break;
    csdim++;
  }
  memcpy(temp, Phi, nnz_Phi*sizeof(double));
  for (i=0; i<csdim; i++) {
    int ii = ndof*i;
    int jj = ndof*(jpvt[i]-1);
    memcpy(&Phi[ii], &temp[jj], ndof*sizeof(double));
  }
  csdimP = csdim;
  //  cout << "csdim = " << csdim << endl;
  delete [] tau;
  delete [] dwork;
  /*
  if (MyPID == 1) {
    cout << "ndof = " << ndof << endl;
    cout << "Phi = " << endl;
    for (i=0; i<ndof; i++) {
      for (int j=0; j<csdim; j++) cout << Phi[i+j*ndof] << " ";
      cout << endl;
    }
  }
  */
  //
  // check for presence of rotational dofs for dofs shared by more than
  // one extended substructure
  //
  int ndof_rot(0);
  int *rot_dofs = new int[ndof];
  if (atype == 2) {
    for (i=0; i<ndof; i++) {
      dof = subdofs[i];
      if ((abs(locdof[dof]) >= 4) && (abs(locdof[dof]) <= 6) 
	                          && (nsubdof[dof] > 1)) {
	rot_dofs[ndof_rot] = i;
	imap[dof] = ndof_rot;
	ndof_rot++;
      }
    }
  }
  if (atype == 3) {
    for (i=0; i<ndof; i++) {
      dof = subdofs[i];
      if ((abs(locdof[dof]) >= 2) && (abs(locdof[dof]) <= 3) 
	                          && (nsubdof[dof] > 1)) {
	rot_dofs[ndof_rot] = i;
	imap[dof] = ndof_rot;
	ndof_rot++;
      }
    }
  }
  ndof_rot_ = ndof_rot;
  if (ndof_rot > 0) {
    //
    // form and factor stiffness matrix associated with rotational dofs
    //
    int NumEntries, *Indices, j, m, nnz, dof2, mm, mmm;
    double *Values;
    nnz = 0;
    rowbeg[0] = 0;
    for (i=0; i<ndof_rot; i++) {
      dof = subdofs[rot_dofs[i]];
      A->ExtractMyRowView(dof, NumEntries, Values, Indices);
      for (j=0; j<NumEntries; j++) {
	dof2 = Indices[j];
	if (imap[dof2] > -1) {
	  colidx[nnz] = imap[dof2];
	  K[nnz] = Values[j];
	  nnz++;
	}
      }
      rowbeg[i+1] = nnz;
    }
    sparse_lu *A_rot; A_rot = new sparse_lu();
    int colmax(-1), colmin(1);
    for (i=0; i<nnz; i++) {
      if (colidx[i] > colmax) colmax = colidx[i];
      if (colidx[i] < colmin) colmin = colidx[i];
    }
    assert (colmax < ndof_rot);
    assert (colmin == 0);
    A_rot->factor(ndof_rot, nnz, rowbeg, colidx, K);
    //
    // solve for rotational dofs via static condensation
    //
    for (i=0; i<ndof; i++) {
      if (imap[subdofs[i]] == -1) imap[subdofs[i]] = -10 - i;
    }
    myzero(rhs, ndof_rot*csdim);
    for (m=0; m<csdim; m++) {
      mm =  m*ndof;
      mmm = m*ndof_rot;
      for (i=0; i<ndof_rot; i++) {
	dof = subdofs[rot_dofs[i]];
	A->ExtractMyRowView(dof, NumEntries, Values, Indices);
	for (j=0; j<NumEntries; j++) {
	  dof2 = Indices[j];
	  if (imap[dof2] <= -10) {
	    dof2 = -imap[dof2] - 10;
	    rhs[i+mmm] -= Values[j]*Phi[dof2+mm];
	  }
	}
      }
    }
    A_rot->sol(csdim, rhs, sol, temp);
    delete A_rot;
    for (m=0; m<csdim; m++) {
      mm = m*ndof;
      mmm = m*ndof_rot;
      for (i=0; i<ndof_rot; i++) Phi[rot_dofs[i]+mm] = sol[i+mmm];
    }
    for (i=0; i<ndof; i++) imap[subdofs[i]] = -1;
    //
    // calculate residuals associated with Phi (should be zero for rotational
    // dofs)
    //
    /*
    for (i=0; i<ndof; i++) imap[subdofs[i]] = i;
    for (i=0; i<ndof; i++) {
      dof = subdofs[i];
      for (j=0; j<csdim; j++) {
	double dsum = 0;
	A->ExtractMyRowView(dof, NumEntries, Values, Indices);
	for (m=0; m<NumEntries; m++) {
	  dof2 = Indices[m];
	  if (imap[dof2] >= 0) dsum += Values[m]*Phi[j*ndof+imap[dof2]];
	}
	rhs[j*ndof+i] = dsum;
      }
    }
    for (i=0; i<ndof; i++) imap[subdofs[i]] = -1;
    if (MyPID == 1) {
      cout << "residuals for Phi" << endl;
      for (i=0; i<ndof; i++) {
	for (j=0; j<csdim; j++) cout << rhs[j*ndof+i] << " ";
	cout << endl;
      }
    }
    */
  }
  delete [] rot_dofs;
}

void CLOP_sub::get_cdof(int cs_local[], int & csdimP, double & xcent_,
			double & ycent_, double & zcent_)
{
  for (int i=0; i<csdim; i++) cs_local[i] = jpvt[i];
  csdimP = csdim;
  xcent_ = xcent; ycent_ = ycent; zcent_ = zcent;
}

void CLOP_sub::sum_scalar_multiply(double Edof_sub[], int rbm, double alpha)
{
  int i, j, ibeg, dof;
  for (i=0; i<csdim; i++) {
    if (jpvt[i] == rbm) {
      ibeg = i*ndof;
      for (j=0; j<ndof; j++) {
	dof = subdofs[j];
	Edof_sub[dof] += alpha*Phi[ibeg+j];
      }
    }
  }
}

void CLOP_sub::correct(double Edof_sub[], int rbm, unsigned char nsubdof[],
		       int ldof[], int lbound, int ubound)
{
  int i, j, ibeg, dof;
  for (i=0; i<csdim; i++) {
    if (jpvt[i] == rbm) {
      ibeg = i*ndof;
      for (j=0; j<ndof; j++) {
	dof = subdofs[j];
	if ((abs(ldof[dof]) > lbound) && (abs(ldof[dof]) < ubound)) {
	  if (abs(ldof[dof]) == rbm) {
	    Phi[ibeg+j] += (1 - Edof_sub[dof])/nsubdof[dof];
	  }
	  else {
	    Phi[ibeg+j] += (0 - Edof_sub[dof])/nsubdof[dof];
	  }
	}
      }
    }
  }
}

void CLOP_sub::statcond(unsigned char nsubdof[], unsigned char on_sub_bound[], 
			int imap[], int rowbeg[], int colidx[], double K[], 
			double rhs[], double sol[], double temp[], 
			const Epetra_CrsMatrix *A)
{
  int i, j, ndof_free(0), dof, dof2, NumEntries, *Indices, nnz(0);
  int m, mm, mmm;
  double *Values;
  for (i=0; i<ndof; i++) {
    dof = subdofs[i];
    if ((nsubdof[dof] > 1) && (on_sub_bound[dof] == 0)) ndof_free++;
  }
  int *free_dofs = new int[ndof_free];
  ndof_free = 0;
  for (i=0; i<ndof; i++) {
    dof = subdofs[i];
    if ((nsubdof[dof] > 1) && (on_sub_bound[dof] == 0)) {
      free_dofs[ndof_free] = i;
      imap[dof] = ndof_free;
      ndof_free++;
    }
  }
  /*
  if (MyPID == 1) {
    cout << "ndof = " << ndof << endl;
    cout << "Phi = " << endl;
    for (i=0; i<ndof; i++) {
      for (int j=0; j<csdim; j++) cout << Phi[i+j*ndof] << " ";
      cout << endl;
    }
  }
  */
  if ((ndof_free > 0) && (ndof_free != ndof)) {
    //
    // form and factor stiffness matrix associated with free dofs
    //
    rowbeg[0] = 0;
    for (i=0; i<ndof_free; i++) {
      dof = subdofs[free_dofs[i]];
      A->ExtractMyRowView(dof, NumEntries, Values, Indices);
      for (j=0; j<NumEntries; j++) {
	dof2 = Indices[j];
	if (imap[dof2] != -1) {
	  colidx[nnz] = imap[dof2];
	  K[nnz] = Values[j];
	  nnz++;
	}
      }
      rowbeg[i+1] = nnz;
    }
    sparse_lu *A_free; A_free = new sparse_lu();
    A_free->factor(ndof_free, nnz, rowbeg, colidx, K);
    //
    // solve for free dofs via static condensation
    //
    for (i=0; i<ndof; i++) {
      if (imap[subdofs[i]] == -1) imap[subdofs[i]] = -10 - i;
    }
    myzero(rhs, ndof_free*csdim);
    for (m=0; m<csdim; m++) {
      mm =  m*ndof;
      mmm = m*ndof_free;
      for (i=0; i<ndof_free; i++) {
	dof = subdofs[free_dofs[i]];
	A->ExtractMyRowView(dof, NumEntries, Values, Indices);
	for (j=0; j<NumEntries; j++) {
	  dof2 = Indices[j];
	  if (imap[dof2] <= -10) {
	    dof2 = -imap[dof2] - 10;
	    rhs[i+mmm] -= Values[j]*Phi[dof2+mm];
	  }
	}
      }
    }
    A_free->sol(csdim, rhs, sol, temp);
    delete A_free;
    for (m=0; m<csdim; m++) {
      mm = m*ndof;
      mmm = m*ndof_free;
      for (i=0; i<ndof_free; i++) Phi[free_dofs[i]+mm] = sol[i+mmm];
    }
    for (i=0; i<ndof; i++) imap[subdofs[i]] = -1;
    delete [] free_dofs;
  }
  //
  // warning: temporary code
  //
  if (xcent < 0.3) {
    //    myzero(Phi, csdim*ndof);
  }
  /*  
  if (MyPID == 1) {
    cout << "ndof = " << ndof << endl;
    cout << "Phi = " << endl;
    for (i=0; i<ndof; i++) {
      for (int j=0; j<csdim; j++) cout << Phi[i+j*ndof] << " ";
      cout << endl;
    }
  }
  */
}

void CLOP_sub::get_Phi_ptr(double* & Phi_ptr)
{
  Phi_ptr = Phi;
}

void CLOP_sub::subpre(double rr[], double zz[], double rhs_work[], 
		      double sol_work[], double tmp_work[])
{
  int i, nrhs(1);
  for (i=0; i<ndof; i++) rhs_work[i] = rr[subdofs[i]];
  A_sub->sol(nrhs, rhs_work, sol_work, tmp_work);
  for (i=0; i<ndof; i++) zz[subdofs[i]] += sol_work[i];
}
