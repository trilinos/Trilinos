/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#ifndef __MLGGB_
#define __MLGGB_

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif



#define SHIFTS   0
#define MAX_ITRS 2
#define MODE     6

#include <stdio.h>
#include <stdlib.h>
#include "ml_common.h"
#include "ml_mat_formats.h"
#include "ml_lapack.h"
#include "ml_eigf2c.h"

#define DNAUPD_F77  F77_FUNC(dnaupd,DNAUPD)
#define PDNAUPD_F77  F77_FUNC(pdnaupd,PDNAUPD)


struct ML_Eigenvalue_Struct  {
  int     Max_Iter;                  /* User input from input file */
  int     Num_Eigenvalues;
  int     Arnoldi;
  double  Residual_Tol;
  int     Fattening;


  int     Nflag;          /* Flag to indicate the first Newton iteration */
  int     Pnconv;         /* Previous number of converged eigenvalues */
  double *Evec, *Eval;    /* eigenvectors and eigenvalues to be reused
			     with MGGB */

};

#ifndef ML_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif



void  ML_ARPACK_driver(char which[],
			 char bmat[], int iparam[], int mode,
			 int nev, int ncv, double tol,  ML *ml,
		       struct ML_CSR_MSRdata  *mydata, int Fattening,
		       struct ML_Eigenvalue_Struct *eigen_struct,
		       int Debug_Flag, int GGB_alp_flag);


  void ML_GGB2CSR (double *v, int nconv, int MatSize, int proc_id,
		   struct ML_CSR_MSRdata  *mydata, int Debug_Flag );


  void  ML_GGBalp (double *NewVec, int nconv, int nloc2, struct ML_Eigenvalue_Struct
		   *eigen_struct);

  extern double  ML_subspace (int nrows, double *inp1, int ncols1, double *inp2, int ncols2);



  extern void ML_ARPACK_GGB(
			    struct ML_Eigenvalue_Struct *eigen_struct,ML *ml,
			    struct ML_CSR_MSRdata *mydata, int Debug_Flag,
			    int GGB_alp_flag);

  extern int  ML_MGGB_angle(struct ML_Eigenvalue_Struct *eigen_struct,ML *ml,
			    struct ML_CSR_MSRdata *mydata);

  extern int  ML_Rayleigh (ML *ml, int nrows, double *q, int count);

  extern double *ML_complex_gdot(int leng, double *ureal, double *uimag, double *vreal, double *vimag,
				 ML_Comm *comm);

  extern double ML_normc(double *real, double *imag,  int leng );

  extern void ML_Eig_Destroy(void *data);

  extern int ML_OperatorGGB_Apply (double *densemat, int Nrows, int Ncols, double *din, double *dout, int Transpose);

void PREFIX DNAUPD_F77(int *, char *, int *, char *, int *, double *, double *,
		 int *, double *, int *, int *, int *, double *, double *,
	       int *, int *);

void PREFIX PDNAUPD_F77(int *, int *, char *, int *, char *, int *, double *, double *,
		int *, double *, int *, int *, int *, double *, double *,
		int *, int *);


#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif
