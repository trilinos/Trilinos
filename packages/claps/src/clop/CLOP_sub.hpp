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
		     int imap[], const Epetra_Comm* Comm, int &nnz);
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
  int ndof, csdim, csdim_max, ndim, atype, INFO, MyPID;
  double xcent, ycent, zcent, *x, *y, *z, *Edof, *Phi;
  sparse_lu *A_sub;
  int *jpvt, *subdofs, *locdof;
};
  
