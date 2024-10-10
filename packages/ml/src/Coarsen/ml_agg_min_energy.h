/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#ifndef ML_AGG_MIN_ENERGY
#define ML_AGG_MIN_ENERGY

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

#include "ml_common.h"
#include "ml_include.h"

#ifndef ML_CPP
#ifdef __cplusplus
extern "C"
{
#endif
#endif

int ML_AGG_Gen_Prolongator_MinEnergy(ML *ml,int level, int clevel, void *data);
int ML_AGG_Gen_Restriction_MinEnergy(ML *ml,int level, int clevel, void *data);

void ML_multiply_all_vscale(ML_Operator* left, ML_Operator* right,
                               double* InnerProd, double* diagonal);
void ML_ImplicitAbs_Destroy(void *data);
int ML_ImplicitAbs_Getrow(ML_Operator *data, int N_requested_rows,
                       int requested_rows[], int allocated_space,
                       int columns[], double values[],
                       int row_lengths[]);
int ML_ImplicitAbs_Matvec(ML_Operator *Amat_in, int ilen, double p[],
                    int olen, double ap[]);
ML_Operator *ML_Operator_ImplicitAbs(ML_Operator *Amat, int OnDestroy_FreeChild);
void ML_Sort_Cols(struct ML_CSR_MSRdata * CSR_Data, int nRows);
double ML_MaxEntry(ML_Operator * A);
void ML_Enforce_Sparsity(ML_Operator * A, struct ML_CSR_MSRdata *Pattern);
void ML_print_mat(double * mat, int rows, int cols, char FileName[]);
void ML_Squeeze_Out_Zeros(ML_Operator *A);
void ML_Drop(ML_Operator * A, double droptol);
ML_Operator * ML_ODE_Strength_Matrix(ML_Operator * A, int num_steps, double t_final, double drop_tol);
void ML_Satisfy_Constraints(ML_Operator *Update, ML_Operator *Pattern, double *Bone, double *BtBinv, int * F, int numPDEs, int numDOFs, int numNodes, int NullDim);
int ML_AGG_Gen_Prolongator_MandelMinEnergy(ML *ml,int level, int clevel, void *data);




#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif
#endif
