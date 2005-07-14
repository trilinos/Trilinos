#include "ml_common.h"
#include "ml_include.h"
#include "ml_agg_min_energy.h"
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

using namespace std;

static double multiply(int row, struct ML_CSR_MSRdata* left, 
                       struct ML_CSR_MSRdata* right)
{
  double res = 0;

  int*    lrowptr = left->rowptr;
  int     llen    = lrowptr[row + 1] - lrowptr[row];
  int*    lbindx  = &(left->columns[lrowptr[row]]);
  double* lval    = &(left->values[lrowptr[row]]);

  int*    rrowptr = right->rowptr;
  int     rlen    = rrowptr[row + 1] - rrowptr[row];
  int*    rbindx  = &(right->columns[rrowptr[row]]);
  double* rval    = &(right->values[rrowptr[row]]);

  map<int, double> rmap;
  map<int, double>::iterator cur;

  for (int i = 0 ; i < llen ; ++i)
    rmap[rbindx[i]] = rval[i];

  for (int i = 0 ; i < llen ; ++i)
  {
    cur = rmap.find(lbindx[i]);
    if (cur != rmap.end())
      res += lval[i] * (cur->second);
  }

  return(res);
}

/* ************************************************************************* */
/* generate smooth prolongator by minimizing energy                          */
/* ------------------------------------------------------------------------- */

int ML_AGG_Gen_Prolongator_MinEnergy(ML *ml,int level, int clevel, void *data)
{
  int         Ncoarse, Nfine, gNfine, gNcoarse, jj;
  ML_Operator **prev_P_tentatives;
  ML_Aggregate * ag = (ML_Aggregate *) data;

#ifdef ML_TIMING
  double t0;
  t0 =  GetClock();
#endif

  ML_Operator* Amat = &(ml->Amat[level]);
  ML_Operator* Pmat = &(ml->Pmat[clevel]);

  if (Amat->num_PDEs < ag->num_PDE_eqns) Amat->num_PDEs = ag->num_PDE_eqns;
  Amat->num_PDEs = ag->num_PDE_eqns;
  prev_P_tentatives = ag->P_tentative;

  Nfine    = Amat->outvec_leng;
  gNfine   = ML_Comm_GsumInt(ml->comm, Nfine);
  ML_Aggregate_Set_CurrentLevel(ag, level);

  /* =============================================== */
  /* Creates the non-smoothed prolongator, P_0, then */
  /* Checks the dimension of the coarse problem.     */
  /* =============================================== */

  ML_Operator* P_0 = ML_Operator_Create(ml->comm);

  Ncoarse  = ML_Aggregate_Coarsen(ag,Amat,&P_0,ml->comm);
  gNcoarse = ML_Comm_GsumInt( ml->comm, Ncoarse);

  gNcoarse = gNcoarse / P_0->num_PDEs;
  gNfine = gNfine / Amat->num_PDEs;
  if (gNcoarse == 0 || ((1.0*gNfine) / (1.0*gNcoarse+0.1) < 1.05))
  {
    if (P_0 != NULL ) ML_Operator_Destroy(&P_0);
    return -1;
  }

  /* ============================================================ */
  /* Start construction Pmatrix to minimize the energy of         */
  /* each basis function. This requires two additional operators, */
  /* here called P_prime and P_second.                            */
  /* The procedure is as follows:                                 */
  /* 1) Extract the diagonal of A and store its inverse. This is  */
  /*    a point diagonal.                                         */
  /* 2) Create the operator DinvA with implicit scaling.          */
  /* 3) Create the operator P_prime = D A P_0                     */
  /* 4) Create the operator P_second = Z D A P_0                  */
  /* 5) Compute the damping parameters and create Omega           */
  /* 6) Create the operator OmegaP_second which scaled P_second   */
  /* 6) Compute Pmatrix = P_prime - OmegaP_second.                */
  /* ============================================================ */

  int row_length;
  int allocated = 128;
  int*    bindx = (int    *)  ML_allocate(allocated*sizeof(int   ));
  double* val   = (double *)  ML_allocate(allocated*sizeof(double));
  int n = Amat->getrow->Nrows;
  vector<double> Dinv(n);

  for (int i = 0 ; i < n ; i++) 
  {
    ML_get_matrix_row(Amat, 1, &i, &allocated, &bindx, &val,
                      &row_length, 0);
    for  (int j = 0; j < row_length; j++) {
      Dinv[i] = 0.0;
      if (i == bindx[j]) {
        if (val[j] == 0.0)
        {
          cerr << "*ML*ERR* zero diagonal element on row " << i << endl;
          exit(EXIT_FAILURE);
        }
        Dinv[i] = 1.0 / val[j];
        break;
      }
    }
  }

  ML_free(bindx);
  ML_free(val);

  ML_Operator* DinvA = 0;
  DinvA = ML_Operator_ImplicitlyVScale(Amat, &Dinv[0], 0);

  ML_Operator* P_prime = ML_Operator_Create(P_0->comm);
  ML_2matmult(DinvA, P_0, P_prime, ML_CSR_MATRIX);

  ML_Operator* DinvA_trans = ML_Operator_Create(P_0->comm);
  ML_Operator_Transpose_byrow(DinvA, DinvA_trans);

  ML_Operator* Z = ML_Operator_Create(P_0->comm);
  ML_Operator_Add(DinvA, DinvA_trans, Z, ML_CSR_MATRIX, +1.0);

  ML_Operator* P_second = ML_Operator_Create(P_0->comm);
  ML_2matmult(Z, P_prime, P_second, ML_CSR_MATRIX);

  ML_Operator* P_0_trans = ML_Operator_Create(P_0->comm);
  ML_Operator_Transpose_byrow(P_0, P_0_trans);

  ML_Operator* P_prime_trans = ML_Operator_Create(P_0->comm);
  ML_Operator_Transpose_byrow(P_prime, P_prime_trans);

  ML_Operator* P_second_trans = ML_Operator_Create(P_0->comm);
  ML_Operator_Transpose_byrow(P_second, P_second_trans);

  struct ML_CSR_MSRdata* P_0_data;
  P_0_data = (struct ML_CSR_MSRdata *)ML_Get_MyGetrowData(P_0_trans);

  struct ML_CSR_MSRdata* P_prime_data;
  P_prime_data = (struct ML_CSR_MSRdata *)ML_Get_MyGetrowData(P_prime_trans);

  struct ML_CSR_MSRdata* P_second_data;
  P_second_data = (struct ML_CSR_MSRdata *)ML_Get_MyGetrowData(P_second_trans);

  vector<double> ColOmega(P_0->outvec_leng);
  
  for (int i = 0 ; i < P_0->invec_leng ; ++i)
  {
    double num = multiply(i, P_0_data, P_second_data);
    double den = multiply(i, P_prime_data, P_second_data);
    ColOmega[i] = num / den;
  }

  /* scale P_0 */
  for (int row = 0 ; row < P_0->invec_leng ; ++row)
  {
    int*    rowptr = P_0_data->rowptr;
    int     len    = rowptr[row + 1] - rowptr[row];
    int*    bindx  = &(P_0_data->columns[rowptr[row]]);
    double* val    = &(P_0_data->values[rowptr[row]]);
    for (int i = 0 ; i < len ; ++i)
      val[i] *= ColOmega[bindx[i]];
  }

#if 0
  ML_Operator* OmegaP_second = NULL;
  OmegaP_second = ML_Operator_ImplicitlyVScale(P_second, &RowOmega[0], 0);

  ML_Operator_Add(P_prime, OmegaP_second, Pmat, ML_CSR_MATRIX, -1.0);
#endif

  struct ML_AGG_Matrix_Context widget;
  widget.near_bdry = NULL;
  widget.omega = 1.0;
  widget.drop_tol = ag->drop_tol_for_smoothing;
  widget.Amat   = Amat;
  widget.aggr_info = ag->aggr_info[level];
  ML_Operator* AGGsmoother = ML_Operator_Create(ml->comm);
  ML_Operator_Set_ApplyFuncData(AGGsmoother, widget.Amat->invec_leng,
                                widget.Amat->outvec_leng, &widget,
                                widget.Amat->matvec->Nrows, NULL, 0);
  ML_Operator_Set_Getrow(AGGsmoother, 
                         widget.Amat->getrow->Nrows, 
                         ML_AGG_JacobiSmoother_Getrows);
  ML_CommInfoOP_Clone(&(AGGsmoother->getrow->pre_comm),
                      widget.Amat->getrow->pre_comm);

  ML_2matmult(AGGsmoother, P_0, &(ml->Pmat[clevel]), ML_CSR_MATRIX);
     
  ML_Operator_Set_1Levels(&(ml->Pmat[clevel]), &(ml->SingleLevel[clevel]), 
                          &(ml->SingleLevel[level]));

  ML_Operator_Destroy(&P_0); // FIXME: if keep_P_tentative..
  ML_Operator_Destroy(&P_prime);
  // FIXME ML_Operator_Destroy(&P_second);
  ML_Operator_Destroy(&DinvA);
  ML_Operator_Destroy(&P_0_trans);
  ML_Operator_Destroy(&P_prime_trans);
  ML_Operator_Destroy(&P_second_trans);

  /* FIXME: scale back P_0 if P_tent is kept */
  /* FIXME: DELETE MEMORY!!!! */

#ifdef ML_TIMING
  ml->Pmat[clevel].build_time =  GetClock() - t0;
  ml->timing->total_build_time += ml->Pmat[clevel].build_time;
#endif

  return(0);
}
