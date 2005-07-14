#include "ml_common.h"
#include "ml_include.h"
#include "ml_agg_min_energy.h"
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include "float.h"

using namespace std;

// ============ //
// private data //
// ============ //

static int     Dinv_size = -1;
static double* Dinv      = 0;

// ====================================================================== 
inline static double multiply(int row, struct ML_CSR_MSRdata* left, 
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

  for (int i = 0 ; i < rlen ; ++i)
    rmap[rbindx[i]] = rval[i];

  for (int i = 0 ; i < llen ; ++i)
  {
    cur = rmap.find(lbindx[i]);
    if (cur != rmap.end())
      res += lval[i] * (cur->second);
  }

  return(res);
}

// ====================================================================== 
inline static double multiply_self(int row, struct ML_CSR_MSRdata* self)
{
  double res = 0;

  int*    rowptr = self->rowptr;
  int     len    = rowptr[row + 1] - rowptr[row];
  double* val    = &(self->values[rowptr[row]]);

  for (int i = 0 ; i < len ; ++i)
  {
    res += val[i] * val[i];
  }

  return(res);
}

// ====================================================================== 
// Note: Op must be an CSR matrix but it is easy to make it for
// a generic getrow()
// Note: off-processor connections are simply discarded
static void multiply_self_all(ML_Operator* Op, double* result)
{
  int n = Op->invec_leng;
  int n_rows = Op->getrow->Nrows;

  for (int i = 0 ; i < n ; ++i) result[i] = 0.0;

  struct ML_CSR_MSRdata* data = 0;
  data = (struct ML_CSR_MSRdata*)ML_Get_MyGetrowData(Op);
  int* rowptr = data->rowptr;
  int* columns = data->columns;
  double* values = data->values;

  for (int row = 0 ; row < n_rows ; ++row)
  {
    int     len    = rowptr[row + 1] - rowptr[row];
    int*    bindx  = &(columns[rowptr[row]]);
    double* val    = &(values[rowptr[row]]);

    for (int i = 0 ; i < len ; ++i)
      if (bindx[i] < n)
        result[bindx[i]] += val[i] * val[i];
  }
}

// ====================================================================== 
// Note: Op must be an CSR matrix but it is easy to make it for
// a generic getrow()
// Note: off-processor connections are simply discarded
inline static void multiply_all(ML_Operator* left, ML_Operator* right,
                                double* result)
{
  int n = left->invec_leng;
  int n_rows = left->getrow->Nrows;

  if (n != right->invec_leng || n_rows != right->getrow->Nrows)
  {
    cerr << "Error: non-comparible operators" << endl;
    cerr << "(file " << __FILE__ << ", line " << __LINE__ << ")" << endl;
    exit(EXIT_FAILURE);
  }

  for (int i = 0 ; i < n ; ++i) result[i] = 0.0;

  struct ML_CSR_MSRdata* left_data = 0;
  struct ML_CSR_MSRdata* right_data = 0;
  left_data = (struct ML_CSR_MSRdata*)ML_Get_MyGetrowData(left);
  right_data = (struct ML_CSR_MSRdata*)ML_Get_MyGetrowData(right);

  int*    lrowptr = left_data->rowptr;
  int*    rrowptr = right_data->rowptr;
  map<int, double>::iterator cur;

  for (int row = 0 ; row < n_rows ; ++row)
  {
    int     llen    = lrowptr[row + 1] - lrowptr[row];
    int*    lbindx  = &(left_data->columns[lrowptr[row]]);
    double* lval    = &(left_data->values[lrowptr[row]]);

    int     rlen    = rrowptr[row + 1] - rrowptr[row];
    int*    rbindx  = &(right_data->columns[rrowptr[row]]);
    double* rval    = &(right_data->values[rrowptr[row]]);

    map<int, double> lmap;

    for (int i = 0 ; i < llen ; ++i)
    {
      lmap[lbindx[i]] = lval[i];
    }

    for (int i = 0 ; i < rlen ; ++i)
    {
      int pos = rbindx[i];
      if (pos < n) 
      {
        cur = lmap.find(pos);
        if (cur != lmap.end())
          result[pos] += rval[i] * (cur->second);
      }
    }
  }
}

// ====================================================================== 
// generate smooth prolongator by minimizing energy                      
//
// \author Marzio Sala and Ray Tuminaro, 9214
//
// \date 14-Jul-05
// ====================================================================== 
int ML_AGG_Gen_Prolongator_MinEnergy(ML *ml,int level, int clevel, void *data)
{
  int         Ncoarse, Nfine, gNfine, gNcoarse;
  ML_Operator **prev_P_tentatives;
  ML_Aggregate * ag = (ML_Aggregate *) data;

#ifdef ML_TIMING
  double t0;
  t0 =  GetClock();
#endif

  ML_Operator* Amat = &(ml->Amat[level]); //already created and filled
  prev_P_tentatives = ag->P_tentative; // for keep_P_tentative

  if (Amat->num_PDEs < ag->num_PDE_eqns) Amat->num_PDEs = ag->num_PDE_eqns;
  Amat->num_PDEs = ag->num_PDE_eqns;

  Nfine    = Amat->outvec_leng;
  gNfine   = ML_Comm_GsumInt(ml->comm, Nfine);
  ML_Aggregate_Set_CurrentLevel(ag, level);

  /* =============================================== */
  /* Creates the non-smoothed prolongator, P_0, then */
  /* Checks the dimension of the coarse problem.     */
  /* Note that I always need P_0 to form R.          */
  /* =============================================== */

  ML_Operator* P_0 = ML_Operator_Create(ml->comm);

  if (prev_P_tentatives == 0) 
  {
    ag->P_tentative = ML_Operator_ArrayCreate(ag->max_levels);
    prev_P_tentatives = ag->P_tentative;
    for (int jj = 0 ; jj < ag->max_levels; jj++) prev_P_tentatives[jj] = 0;
  }

  if ((prev_P_tentatives != 0) && (prev_P_tentatives[clevel] != 0))
  {
    P_0 = prev_P_tentatives[clevel];
    Ncoarse = P_0->invec_leng;
  }
  else
  {
    Ncoarse  = ML_Aggregate_Coarsen(ag,Amat,&P_0,ml->comm);
    prev_P_tentatives[clevel] = P_0;
  }

  gNcoarse = ML_Comm_GsumInt( ml->comm, Ncoarse);
  gNcoarse = gNcoarse / P_0->num_PDEs;
  gNfine = gNfine / Amat->num_PDEs;
  if (gNcoarse == 0 || ((1.0*gNfine) / (1.0*gNcoarse+0.1) < 1.05))
  {
    ML_Operator_Destroy(&P_0);
    return -1;
  }

  // ============================== //
  // Start the construction of Pmat //
  // ============================== //
  
  int row_length;
  int allocated = 128;
  int*    bindx = (int    *)  ML_allocate(allocated*sizeof(int   ));
  double* val   = (double *)  ML_allocate(allocated*sizeof(double));
  int n = Amat->getrow->Nrows;

  if (Dinv != 0 || Dinv_size != -1) 
  {
    cerr << "Error: Static data Dinv is not null or Dinv_size is wrong!" << endl;
    cerr << "(file " << __FILE__ << ", line " << __LINE__ << ")" << endl;
    exit(EXIT_FAILURE);
  }

  // this vector is free'd in the generation of the restriction
  // so that I don't have to query for the diagonal.
  Dinv = (double*)ML_allocate(sizeof(double) * n);
  Dinv_size = n;

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

  ML_Operator* Scaled_A = 0;
  Scaled_A = ML_Operator_ImplicitlyVScale(Amat, Dinv, 0);
  ML_CommInfoOP_Clone(&(Scaled_A->getrow->pre_comm), Amat->getrow->pre_comm);

  ML_Operator* P_prime = ML_Operator_Create(P_0->comm);
  ML_2matmult(Scaled_A, P_0, P_prime, ML_CSR_MATRIX);

  ML_Operator* Z = 0;
  ML_Operator* Scaled_A_trans = 0;
  ML_Operator* P_second = 0;

  int n_0 = P_0->invec_leng;
  int n_0_tot = ML_gsum_int(n_0, P_0->comm);

  vector<double> num(n_0);
  vector<double> den(n_0);
  vector<double> ColOmega(n_0);

  switch (ag->minimizing_energy) {
  case 1:
    // Z_1 = I                          
    // This is simple, no need for P_second because == P_prime

    multiply_all(P_0, P_prime, &num[0]);
    multiply_self_all(P_prime, &den[0]);
    break;

  case 2:
    // Z_2 = A^T * A. Need to be smart here to avoid the construction of Z_2
    
    Scaled_A_trans = ML_Operator_Create(P_0->comm);
    ML_Operator_Transpose_byrow(Scaled_A, Scaled_A_trans);

    P_second = ML_Operator_Create(P_0->comm);
    ML_2matmult(Scaled_A, P_prime, P_second, ML_CSR_MATRIX);

    multiply_all(P_prime, P_second, &num[0]);
    multiply_all(P_second, P_second, &den[0]);
    break;

  case 3:
    // Z_3 = A^T + A
    // Need P_second, and to form A^T, then Z_3
    
    Scaled_A_trans = ML_Operator_Create(P_0->comm);
    ML_Operator_Transpose_byrow(Scaled_A, Scaled_A_trans);

    Z = ML_Operator_Create(P_0->comm);
    ML_Operator_Add(Scaled_A, Scaled_A_trans, Z, ML_CSR_MATRIX, 1.0);

    P_second = ML_Operator_Create(P_0->comm);
    ML_2matmult(Z, P_prime, P_second, ML_CSR_MATRIX);

    multiply_all(P_0, P_second, &num[0]);
    multiply_all(P_prime, P_second, &den[0]);
    break;

  default:
    // should never be here
    cerr << "Incorrect parameter (" << ag->minimizing_energy << ")" << endl;
    cerr << "(file " << __FILE__ << ", line " << __LINE__ << ")" << endl;
    exit(EXIT_FAILURE);
  }

  int zero_local = 0;
  double min_local = DBL_MAX;
  double max_local = DBL_MIN;

  for (int i = 0 ; i < n_0 ; ++i) 
  {
    ColOmega[i] = num[i] / den[i];
    double& val = ColOmega[i];
    if (val < 0.0) 
    {
      val = 0.0;
      ++zero_local;
    }
    if (val < min_local) min_local = val;
    if (val > max_local) max_local = val;
  }

  double min_all = ML_gmin_double(min_local, P_0->comm);
  double max_all = ML_gmax_double(max_local, P_0->comm);
  double zero_all = ML_gsum_int(zero_local, P_0->comm);

  if (ML_Get_PrintLevel() > 5 && P_0->comm->ML_mypid == 0)
  { 
    cout << endl;
    cout << "Prolongator Smoothing: Using energy minimization (scheme = " 
         << ag->minimizing_energy << ")" << endl;
    cout << "Damping parameter: min = " << min_all <<  ", max = " << max_all 
         << " (" << zero_all << " zeros out of " << n_0_tot << ")" << endl;
  }

  // convert the omega's from column-based to row-based

  vector<double> RowOmega(n);
  for (int i = 0 ; i < n ; i++) RowOmega[i] = DBL_MAX;

  int* aggr_info = ag->aggr_info[level];

  for (int row = 0 ; row < n ; row++) 
  {
    ML_get_matrix_row(Amat, 1, &row, &allocated, &bindx, &val,
                      &row_length, 0);

    for  (int j = 0; j < row_length; j++) 
    {
      if (bindx[j] < n) // this is white magic
      {
        int aggr = aggr_info[bindx[j]];
        double omega = ColOmega[aggr];
        if (omega < RowOmega[row]) RowOmega[row] = omega;
      }
    }
  }

  ML_Operator* Scaled_P_prime = 0;
  Scaled_P_prime = ML_Operator_ImplicitlyVScale(P_prime, &RowOmega[0], 0);
  ML_CommInfoOP_Clone(&(Scaled_P_prime->getrow->pre_comm), P_prime->getrow->pre_comm);

  ML_Operator_Add(P_0, Scaled_P_prime, &(ml->Pmat[clevel]), 
                  ML_CSR_MATRIX, -1.0);

  ML_Operator_Set_1Levels(&(ml->Pmat[clevel]), &(ml->SingleLevel[clevel]),
                          &(ml->SingleLevel[level]));

  ML_free(bindx);
  ML_free(val);

  if (Z)              ML_Operator_Destroy(&Z);
  if (P_prime)        ML_Operator_Destroy(&P_prime);
  if (P_second)       ML_Operator_Destroy(&P_second);
  if (Scaled_A)       ML_Operator_Destroy(&Scaled_A);
  if (Scaled_A_trans) ML_Operator_Destroy(&Scaled_A_trans);
  if (Scaled_P_prime) ML_Operator_Destroy(&Scaled_P_prime);

#ifdef ML_TIMING
  ml->Pmat[clevel].build_time =  GetClock() - t0;
  ml->timing->total_build_time += ml->Pmat[clevel].build_time;
#endif

  return(0);
}

// ====================================================================== 
// generate smoothed restriction by minimizing energy                      
// Note: in the symmetric case we don't really need this, a transpose
// of P suffices.
//
// \author Marzio Sala and Ray Tuminaro, 9214
//
// \date 14-Jul-05
// ====================================================================== 
int ML_AGG_Gen_Restriction_MinEnergy(ML *ml,int level, int clevel, void *data)
{
  ML_Operator **prev_P_tentatives;
  ML_Aggregate * ag = (ML_Aggregate *) data;

#ifdef ML_TIMING
  double t0;
  t0 =  GetClock();
#endif

  ML_Operator* Amat = &(ml->Amat[level]); //already created and filled
  prev_P_tentatives = ag->P_tentative; // for keep_P_tentative
  ML_Operator* P_0 = prev_P_tentatives[clevel]; // already created and filled

  /* ============================================================ */
  /* Start construction Pmatrix to minimize the energy of         */
  /* each basis function. This requires two additional operators, */
  /* here called P_prime and P_second.                            */
  /* ============================================================ */

  int row_length;
  int allocated = 128;
  int*    bindx = (int    *)  ML_allocate(allocated*sizeof(int   ));
  double* val   = (double *)  ML_allocate(allocated*sizeof(double));
  int n = Amat->getrow->Nrows;

  if (Dinv == 0 || Dinv_size != n) 
  {
    cerr << "Error: Static data Dinv is null or Dinv_size is wrong!" << endl;
    cerr << "(file " << __FILE__ << ", line " << __LINE__ << ")" << endl;
    exit(EXIT_FAILURE);
  }

  ML_Operator* Scaled_A = 0;
  // FIXME: scale by rows or scale by columns??
  Scaled_A = ML_Operator_ImplicitlyVScale(Amat, Dinv, 0);
  ML_CommInfoOP_Clone(&(Scaled_A->getrow->pre_comm), Amat->getrow->pre_comm);

  ML_Operator* Scaled_A_trans = ML_Operator_Create(P_0->comm);
  ML_Operator_Transpose_byrow(Scaled_A, Scaled_A_trans);

  ML_Operator* P_prime = ML_Operator_Create(P_0->comm);
  ML_2matmult(Scaled_A_trans, P_0, P_prime, ML_CSR_MATRIX);

  ML_Operator* Z = 0;
  ML_Operator* P_second = 0;

  int n_0 = P_0->invec_leng;
  int n_0_tot = ML_gsum_int(n_0, P_0->comm);
  vector<double> num(n_0);
  vector<double> den(n_0);
  vector<double> ColOmega(n_0);

  switch (ag->minimizing_energy) {
  case 1:
    // Z_1 = I                          
    // This is simple, no need for P_second because == P_prime
    multiply_all(P_0, P_prime, &num[0]);
    multiply_self_all(P_prime, &den[0]);
    break;

  case 2:
    // Z_2 = A^T * A. Need to be smart here to avoid the construction of Z_2
    
    P_second = ML_Operator_Create(P_0->comm);
    ML_2matmult(Scaled_A_trans, P_prime, P_second, ML_CSR_MATRIX);

    multiply_all(P_prime, P_second, &num[0]);
    multiply_all(P_second, P_second, &den[0]);
    break;

  case 3:
    // Z_3 = A^T + A
    // Need P_second, and to form A^T, then Z_3
    
    Z = ML_Operator_Create(P_0->comm);
    ML_Operator_Add(Scaled_A, Scaled_A_trans, Z, ML_CSR_MATRIX, 1.0);

    P_second = ML_Operator_Create(P_0->comm);
    ML_2matmult(Z, P_prime, P_second, ML_CSR_MATRIX);

    multiply_all(P_0, P_second, &num[0]);
    multiply_all(P_prime, P_second, &den[0]);
    break;

  default:
    // should never be here
    cerr << "Incorrect parameter (" << ag->minimizing_energy << ")" << endl;
    cerr << "(file " << __FILE__ << ", line " << __LINE__ << ")" << endl;
    exit(EXIT_FAILURE);
  }

  int zero_local = 0;

  for (int i = 0 ; i < n_0 ; ++i) 
  {
    ColOmega[i] = num[i] / den[i];
    if (ColOmega[i] < 0.) 
    {
      ColOmega[i] = 0;
      ++zero_local;
    }
  }

  double min_local = DBL_MAX;
  double max_local = DBL_MIN;

  for (int i = 0 ; i < n_0 ; ++i)
  {
    if (ColOmega[i] < min_local) min_local = ColOmega[i];
    if (ColOmega[i] > max_local) max_local = ColOmega[i];
  }

  double min_all = ML_gmin_double(min_local, P_0->comm);
  double max_all = ML_gmax_double(max_local, P_0->comm);
  double zero_all = ML_gsum_int(zero_local, P_0->comm);

  if (ML_Get_PrintLevel() > 5 && P_0->comm->ML_mypid == 0)
  { 
    cout << endl;
    cout << "Restriction Smoothing: Using energy minimization (scheme = " 
         << ag->minimizing_energy << ")" << endl;
    cout << "Damping parameter: min = " << min_all <<  ", max = " << max_all 
         << " (" << zero_all << " zeros out of " << n_0_tot << ")" << endl;
    cout << endl;
  }


  // convert the omega's from column-based to row-based

  vector<double> RowOmega(n);
  for (int i = 0 ; i < n ; i++) RowOmega[i] = DBL_MAX;

  int* aggr_info = ag->aggr_info[level];

  for (int row = 0 ; row < n ; row++) 
  {
    ML_get_matrix_row(Amat, 1, &row, &allocated, &bindx, &val,
                      &row_length, 0);

    for  (int j = 0; j < row_length; j++) 
    {
      int pos = bindx[j];
      if (pos < n) // another big of magic
      {
        int aggr = aggr_info[pos];
        double omega = ColOmega[aggr];
        if (omega < RowOmega[row]) RowOmega[row] = omega;
      }
    }
  }

  ML_Operator* Scaled_P_prime = 0;
  Scaled_P_prime = ML_Operator_ImplicitlyVScale(P_prime, &RowOmega[0], 0);
  ML_CommInfoOP_Clone(&(Scaled_P_prime->getrow->pre_comm), 
                      P_prime->getrow->pre_comm);

  ML_Operator* temp = ML_Operator_Create(P_0->comm);
  ML_Operator_Add(P_0, Scaled_P_prime, temp, ML_CSR_MATRIX, -1.0);

  ML_Operator_Transpose_byrow(temp, &(ml->Rmat[level]));

  ML_Operator_Set_1Levels(&(ml->Rmat[level]), &(ml->SingleLevel[level]), 
                          &(ml->SingleLevel[clevel]));

  ML_free(bindx);
  ML_free(val);

  if (Scaled_A)       ML_Operator_Destroy(&Scaled_A);
  if (Scaled_A_trans) ML_Operator_Destroy(&Scaled_A_trans);
  if (P_prime)        ML_Operator_Destroy(&P_prime);
  if (Scaled_P_prime) ML_Operator_Destroy(&Scaled_P_prime);
  if (Z)              ML_Operator_Destroy(&Z);
  if (P_second)       ML_Operator_Destroy(&P_second);
  if (temp)           ML_Operator_Destroy(&temp);

  ML_free(Dinv); Dinv = 0;
  Dinv_size = -1;

#ifdef ML_TIMING
  ml->Rmat[level].build_time =  GetClock() - t0;
  ml->timing->total_build_time += ml->Rmat[level].build_time;
#endif

  return(0);
}
