#ifndef __MLXYT__
#define __MLXYT__


#ifdef __cplusplus
extern "C" {
#endif



extern void setup_henry(ML *my_ml, int grid0, int **imapper, int **separator,
        int **sep_size, int *Nseparators, int *Nlocal, int *Nghost,
        ML_Operator **matvec_data);

int CSR_submv(ML_Operator *Amat, double p[], double ap[]);
int CSR_submatvec(ML_Operator *Amat, double p[], double ap[], int mask);
int ML_submv(ML_Operator *Amat, double p[], double ap[]);
int ML_submatvec(ML_Operator *Amat, double p[], double ap[], int mask);
int ML_Gen_CoarseSolverXYT( ML *ml, int i);
void ML_XYTfree(void *temp);

#ifdef __cplusplus
}
#endif

#endif
