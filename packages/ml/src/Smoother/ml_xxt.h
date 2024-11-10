/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#ifndef __MLXYT__
#define __MLXYT__

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif
#include "ml_common.h"
extern void setup_henry_xxt(ML *my_ml, int grid0, int **imapper, int **separator,
        int **sep_size, int *Nseparators, int *Nlocal, int *Nghost,
        ML_Operator **matvec_data);

extern int CSR_submv(ML_Operator *Amat, double p[], double ap[], int mask);
int CSR_submatvec(ML_Operator *Amat, double p[], double ap[], int mask);
int ML_submv(ML_Operator *Amat, double p[], double ap[], int mask);
int ML_submatvec(ML_Operator *Amat, double p[], double ap[], int mask);

extern int ML_Comm_subGappendInt(ML_Comm *com_ptr, int *vals, int *cur_length,
                    int total_length,int sub_mask);
extern void ML_subexchange_bdry(double x[], ML_CommInfoOP *comm_info,
                      int start_location, int total_send, ML_Comm *comm,
                      int mask);
#endif
