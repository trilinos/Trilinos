#ifndef __MLMLS__
#define __MLMLS__

#define MLS_MAX_DEG  5 /* max. degree of MLS smoother        */

struct MLSthing {
 /*
  * set degree and 
  * the precomputed coefficients used in MLS application
  */
  int     mlsReady;  
  int     mlsLevel;  
  int     mlsDeg; /* degree for current level */
  double  mlsBoost; 
  double  mlsOver; 
  double  mlsOm[MLS_MAX_DEG];
  double  mlsOm2;
  double  mlsCf[MLS_MAX_DEG];
  int     nIsolNod; 
  double *pIsolNod; 
  double *res0, *res, *y;  /* workarrays allocated in .... to be reused */
};


#ifdef __cplusplus
extern "C" {
#endif

int   ML_MLS_GetIsoln(int lv); 
int   ML_MLS_SetOm(int lv, int deg); 
int   ML_MLS_Smooth0(double b[], double vx[], double vy[], int deg, 
		     double *om, double *cf, int nit, double over, 
		     double wk);
int   ML_MLS_Smooth1(double b[], double vx[], double vy[], int deg, 
		     double *om, double *cf, int nit, double over, 
		     double wk);
int ML_MLS_SandwPres(void *sm, int inlen, double x[], int outlen, double y[]);
int ML_MLS_SandwPost(void *sm, int inlen, double x[], int outlen, double y[]);
int ML_MLS_SPrime_Apply(void *sm,int inlen,double x[],int outlen, double rhs[]);


#ifdef __cplusplus
}
#endif

#endif
