/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#ifndef __MLMLS__
#define __MLMLS__

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

#include "ml_common.h"

#define MLS_MAX_DEG  5 /* max. degree of MLS smoother        */

struct MLSthing {
 /*
  * set degree and
  * the precomputed coefficients used in MLS application
  */
  int     mlsDeg; /* degree for current level */
  double  mlsBoost;
  double  mlsOver;
  double  mlsOm[MLS_MAX_DEG];
  double  mlsOm2;
  double  mlsCf[MLS_MAX_DEG];
  double *pAux, *res, *y;  /* workarrays allocated in .... to be reused */
  double eig_ratio;
  double eig_boost;
  double beta_real, beta_img;

  ML_Sm_BGS_Data *block_scaling;/* these last arguments are used to */
  ML_Operator *unscaled_matrix;     /* implement block scaling instead  */
  ML_Operator *scaled_matrix;       /* of point scaling when doing      */
                                    /* Chebyshev smoothing. The basic   */
                                    /* idea is to turn off diagonal     */
                                    /* scaling by setting the diagonal  */
                                    /* to 1. Then to create a matrix    */
                                    /* wrapper that does Dinv*Amat*v    */
                                    /* when a matvec is requested (where*/
                                    /* D is some block diagonal */
};

#ifndef ML_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

int   ML_MLS_Smooth0( double b[], double vx[], double vy[], int deg,
		      double *om, double *cf, int nit, double over,
		      double wk);
int   ML_MLS_Smooth1( double b[], double vx[], double vy[], int deg,
		      double *om, double *cf, int nit, double over,
		      double wk);
int ML_MLS_SandwPres(void *sm, int inlen, double x[], int outlen, double y[]);
int ML_MLS_SandwPost(void *sm, int inlen, double x[], int outlen, double y[]);
int ML_MLS_SPrime_Apply(void *sm,int inlen,double x[],int outlen, double rhs[]);


#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif
