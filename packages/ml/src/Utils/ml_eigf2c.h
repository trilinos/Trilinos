#ifndef __MLEIGF2C__
#define __MLEIGF2C__

#include "ml_common.h"
#include "ml_lapack.h"
#include "ml_memory.h"
#include "ml_comm.h"

#define DMOUT_F77   F77_FUNC(dmout,DMOUT)
#define PDMOUT_F77  F77_FUNC(pdmout,PDMOUT)
#define DNEUPD_F77  F77_FUNC(dneupd,DNEUPD)

#ifndef ML_CPP
#ifdef __cplusplus
extern "C"
{
#endif
#endif



  extern int ml_pdmout__(USR_COMM *comm, int *lout, int *m, int *n, double *a, 
			   int *lda, int *idigit);
  
  extern int ml_pdneupc__(int *comm, 
			 int *ivec, char *howmny, int *celect, double *d__, 
			 double *v, int *ldv, double *workev,  char *bmat, int *n, 
			 char *which, int *nev, double *tol, double *resid, int *ncv, 
			 int *iparam, int *ipntr, double *workd, double *workl, 
			 int *lworkl, int *ierr, ftnlen howmny_len, ftnlen bmat_len, 
			 ftnlen which_len);
  /*
  extern int ml_pdneupc__(int *comm,
			  int *ivec, char *howmny, int *celect, double *d__, 
			  double *v, int *ldv, double *workev,  char *bmat, int *n, 
		       char *which, int *nev, double *tol, double *resid, int *ncv, 
			  int *iparam, int *ipntr, double *workd, double *workl, 
			  int *lworkl, int *ierr, ftnlen howmny_len, ftnlen bmat_len, 
			  ftnlen which_len);
  */
extern void PREFIX DMOUT_F77(int *, int *, int *, double *, int *lda, int *idigit,
		      char *, ftnlen);
  
extern void PREFIX  PDMOUT_F77(int *, 
		       int *, int *, int *, double *, int *lda, int *idigit,
		       char *, ftnlen);


extern void PREFIX DNEUPD_F77(int *, char *, int *, double *, double *, double *, 
		      int *, double *, double *, double *, char *bmat, 
		      int *n, char *which, int *nev, double *tol, double *, 
		      int *ncv, double *, int *, int *, int *, double *,
		      double *, int *, int *, ftnlen, ftnlen, ftnlen);
  
  
extern void PREFIX PDNEUPD_F77(int *, 
		       int *, char *, int *, double *, double *, double *, 
		       int *, double *, double *, double *, char *bmat, 
		       int *n, char *which, int *nev, double *tol, double *, 
		       int *ncv, double *, int *, int *, int *, double *,
		       double *, int *, int *, ftnlen, ftnlen, ftnlen);

#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif



#endif


