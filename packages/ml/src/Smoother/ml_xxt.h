#ifndef __MLXYT__
#define __MLXYT__
#include "ml_common.h"
extern void setup_henry(ML *my_ml, int grid0, int **imapper, int **separator,
        int **sep_size, int *Nseparators, int *Nlocal, int *Nghost,
        ML_Operator **matvec_data);

extern int CSR_submv(ML_Operator *Amat, double p[], double ap[], int mask);
int CSR_submatvec(ML_Operator *Amat, double p[], double ap[], int mask);
int ML_submv(ML_Operator *Amat, double p[], double ap[], int mask);
int ML_submatvec(ML_Operator *Amat, double p[], double ap[], int mask);

extern int ML_Comm_subGappendInt(ML_Comm *com_ptr, int *vals, int *cur_length, 
                    int total_length,int sub_mask);
#endif
