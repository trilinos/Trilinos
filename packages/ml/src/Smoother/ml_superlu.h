
/********************************************************************* */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/********************************************************************* */
/*          Utilities for Aztec/SuperLU users                          */
/********************************************************************* */


#ifndef __MLSUPERLU__
#define __MLSUPERLU__

typedef struct ML_Sm_Schwarz_Data_Struct ML_Sm_Schwarz_Data;

struct ML_Sm_Schwarz_Data_Struct 
{
   int           Nrows;
   int           **bmat_ia;
   int           **bmat_ja;
   double        **bmat_aa;
   int           **aux_bmat_ia;
   int           **aux_bmat_ja;
   double        **aux_bmat_aa;
   ML_CommInfoOP *getrow_comm;
   int           nblocks;
   int           *blk_info;
   int           *blk_size;
   int           **blk_indices;
   int           **perm_r;
   int           **perm_c;
#ifdef SUPERLU
   SuperMatrix   **slu_Amat;
   SuperMatrix   **slu_Lmat;
   SuperMatrix   **slu_Umat;
#endif
};

#ifdef __cplusplus
extern "C" {
#endif

extern int ML_SuperLU_Solve(void *vsolver,int ilen,double *x,int olen,
			    double *rhs);

extern int ML_SuperLU_SolveLocal(void *vsolver, double *x, double *rhs);
extern int ML_CSolve_Clean_SuperLU( void *vsolver, ML_CSolveFunc *func);
extern  int ML_Smoother_Create_Schwarz_Data(ML_Sm_Schwarz_Data **data);
extern  int ML_Smoother_VBlockSchwarzDecomposition(ML_Sm_Schwarz_Data *, 
                    ML_Operator *, ML_Comm *, int, int *,int*,double *,int *, 
                    int *,int);

#ifdef __cplusplus
}
#endif

#endif
