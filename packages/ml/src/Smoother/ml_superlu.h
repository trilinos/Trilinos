
/********************************************************************* */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/********************************************************************* */
/*          Utilities for Aztec/SuperLU users                          */
/********************************************************************* */


#ifndef __MLSUPERLU__
#define __MLSUPERLU__

#ifdef __cplusplus
extern "C" {
#endif

extern int ML_SuperLU_Solve(void *vsolver,int ilen,double *x,int olen,
			    double *rhs);

extern int ML_SuperLU_SolveLocal(void *vsolver, double *x, double *rhs);
extern int ML_CSolve_Clean_SuperLU( void *vsolver, ML_CSolveFunc *func);

#ifdef __cplusplus
}
#endif

#endif
