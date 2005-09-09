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

void peek(char *str, ML_Operator *mat, int row ) 
{
  int allocated = 0, *bindx = NULL, row_length, i;
  double *val = NULL;

  ML_get_matrix_row(mat, 1, &row, &allocated, &bindx, &val, &row_length, 0);

  for (i = 0; i < row_length; i++) 
    printf("%s(%d,%d) = %20.13e;\n",str,row+1,bindx[i]+1,val[i]);

  ML_free(bindx); ML_free(val);
}

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
  double t0,t1;

#ifdef ML_TIMING
  double t0;
  t0 =  GetClock();
#endif
  t0 =  GetClock();


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

// rst: comment out 2 lines above and put
// rst: Scaled_A = Amat;

  ML_Operator* DinvAmat;
  DinvAmat = ML_Operator_ImplicitlyBlockDinvScale(Amat);

  ML_Operator *Pprime_gz = NULL, *temp1 = NULL, *temp2 = NULL;
  int *root_pts = NULL, Nroots;
  int count, i, j;
  struct ML_CSR_MSRdata *csr_data;
  int *row_ptr, *columns;
  double *values;
  int *newrowptr, *newcolumns, max_nz_per_row, nnz_in_row;
  double *newvalues;


#define newstuff
  int compress = 1;


  ML_Operator* P_prime = ML_Operator_Create(P_0->comm);
  t1 =  GetClock(); printf("start time = %e\n",t1-t0);  t0 = t1;
  ML_2matmult(Amat, P_0, P_prime, ML_CSR_MATRIX);
  t1 =  GetClock(); printf("A*P0 time = %e\n",t1-t0);  t0 = t1;
  ML_CSR_DropSmall(P_prime, 1e-4, 1e-4, 1e-4);
  t1 =  GetClock(); printf("drop time = %e\n",t1-t0);  t0 = t1;
  ML_AGG_DinvP(P_prime, (MLSthing *) DinvAmat->data, Amat->num_PDEs,0,1);
  t1 =  GetClock(); printf("D^-1 A*P0 = %e\n",t1-t0);  t0 = t1;

  Pprime_gz = P_prime;
  if (compress == 1) {
    temp1     =  ML_Operator_Create(Amat->comm);
    temp2     =  ML_Operator_Create(Amat->comm);
    Pprime_gz = ML_Operator_Create(Amat->comm);

    ML_Operator_Transpose(P_0,temp1);
    t1 =  GetClock(); printf("P_0 transpose = %e\n",t1-t0);  t0 = t1;
    ML_2matmult(temp1, P_prime, temp2, ML_CSR_MATRIX);
    t1 =  GetClock(); printf("Acoarse = %e\n",t1-t0);  t0 = t1;
    Nroots = ML_Operator_MisRootPts( temp2,  Amat->num_PDEs, &root_pts);
    t1 =  GetClock(); printf("Misroot time = %e\n",t1-t0);  t0 = t1;

    csr_data = (struct ML_CSR_MSRdata *) P_prime->data;
    row_ptr = csr_data->rowptr;
    columns = csr_data->columns;
    values  = csr_data->values;
    count = 0;
    for (i = 0; i < n ; i++) {
      for (j = row_ptr[i]; j < row_ptr[i+1]; j++) {
	if ( root_pts[columns[j]] != -1) count++;
      }
    }
    newrowptr  = (int    *)  ML_allocate(sizeof(int)*(n+1));
    newcolumns = (int    *)  ML_allocate(sizeof(int)*(count+1));
    newvalues  = (double *)  ML_allocate(sizeof(double)*(count+1));
    count      = 0;
    newrowptr[0]   = 0;
    max_nz_per_row = 0;
    for (i = 0; i < n ; i++) {
      nnz_in_row = 0;
      for (j = row_ptr[i]; j < row_ptr[i+1]; j++) {
	if ( root_pts[columns[j]] != -1) {
	  nnz_in_row++;
	  newvalues[count   ] = values[j];
	  newcolumns[count++] = root_pts[columns[j]];
	}
      }
      if (nnz_in_row > max_nz_per_row) max_nz_per_row = nnz_in_row;
      newrowptr[i+1] = count;
    }
    csr_data = (struct ML_CSR_MSRdata *) ML_allocate(sizeof(struct ML_CSR_MSRdata));
    if (csr_data == NULL) pr_error("no space for csr_data\n");
    csr_data->columns = newcolumns;
    csr_data->values  = newvalues;
    csr_data->rowptr  = newrowptr;

    ML_Operator_Set_ApplyFuncData(Pprime_gz, Nroots, n, csr_data, n, NULL, 0);
    ML_Operator_Set_Getrow(Pprime_gz, n, CSR_getrow);
    ML_Operator_Set_ApplyFunc (Pprime_gz, CSR_matvec);
    //ML_CommInfoOP_Clone( &(Amat->getrow->pre_comm),  
    Pprime_gz->data_destroy = ML_CSR_MSRdata_Destroy;
    Pprime_gz->max_nz_per_row = max_nz_per_row;
    Pprime_gz->N_nonzeros     = count;
    t1 =  GetClock(); printf("Compress time = %e\n",t1-t0);  t0 = t1;
    //    ML_CSR_DropSmall(Pprime_gz, 1.e-9);

  }

  ML_Operator* Z = 0;
  ML_Operator* P_second = 0;

  int n_0     = P_0->invec_leng;
  int n_0_tot = ML_gsum_int(n_0, P_0->comm);

  vector<double> num(n_0);
  vector<double> tmp(n_0);
  vector<double> den(n_0);
  vector<double> ColOmega(n_0);
#ifdef newstuff
  vector<double> tnum(Nroots);
  vector<double> tden(Nroots);
#endif
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
    ML_2matmult(Amat, Pprime_gz, P_second, ML_CSR_MATRIX);
//  ML_CSR_DropSmall(P_prime, 1e-4);
    t1 =  GetClock(); printf("second 2mat = %e\n",t1-t0);  t0 = t1;
    ML_AGG_DinvP(P_second, (MLSthing *) DinvAmat->data, Amat->num_PDEs,0,1);
    t1 =  GetClock(); printf("second scale = %e\n",t1-t0);  t0 = t1;
    if (compress == 1) {
      multiply_all(Pprime_gz, P_second, &tnum[0]);
      t1 =  GetClock(); printf("first mult_all = %e\n",t1-t0);  t0 = t1;
      multiply_all(P_second, P_second, &tden[0]);
      t1 =  GetClock(); printf("second mult_all = %e\n",t1-t0);  t0 = t1;
    }
    else {
      multiply_all(Pprime_gz, P_second, &num[0]);
      t1 =  GetClock(); printf("first mult_all = %e\n",t1-t0);  t0 = t1;
      multiply_all(P_second, P_second, &den[0]);
      t1 =  GetClock(); printf("second mult_all = %e\n",t1-t0);  t0 = t1;
    }
    break;

  case 3:
    //             diag( P0' ( A'D' + DA) D A P0)
    //   omega =   -----------------------------
    //             diag( P0'A'D' ( A'D' + DA) D A P0)
    //
    //             diag( Prime'Prime + P0'Psecond)
    //         =   -----------------------------
    //                2*diag( Psecond'Prime)
    //
    //    where Prime = D A P0  and Psecond = D A Pprime
    //
    P_second = ML_Operator_Create(P_0->comm);
    
  // if (Amat->num_PDEs == 1) {
  //    Do this Marzio's way as it avoids a couple of transposes.
  //    
  //    ML_2matmult(Scaled_A, P_prime, P_second, ML_CSR_MATRIX);
  // }
  // else {
        ML_2matmult(Amat, P_prime, P_second, ML_CSR_MATRIX);
t1 =  GetClock(); printf("next matmat = %e\n",t1-t0);  t0 = t1;
 ML_AGG_DinvP(P_second, (MLSthing *) DinvAmat->data, Amat->num_PDEs,0,1);
t1 =  GetClock(); printf("next scale = %e\n",t1-t0);  t0 = t1;
  // }   


    multiply_all(P_0, P_second,  &num[0]);
t1 =  GetClock(); printf("first mult_all = %e\n",t1-t0);  t0 = t1;
    multiply_all(P_prime, P_prime, &tmp[0]);
t1 =  GetClock(); printf("second mult_all = %e\n",t1-t0);  t0 = t1;
#ifdef newstuff
    for (int i = 0 ; i < n_0 ; ++i) tnum[i] += tmp[i];
    multiply_all(P_prime, P_second, &tden[0]);
t1 =  GetClock(); printf("third mult_all = %e\n",t1-t0);  t0 = t1;
    for (int i = 0 ; i < n_0 ; ++i) tden[i] *= 2.;
#else
    for (int i = 0 ; i < n_0 ; ++i) num[i] += tmp[i];
    multiply_all(P_prime, P_second, &den[0]);
t1 =  GetClock(); printf("third mult_all = %e\n",t1-t0);  t0 = t1;
    for (int i = 0 ; i < n_0 ; ++i) den[i] *= 2.;
#endif


    //
    //  The old way of calculating this was to first create a
    //  matrix Z such that Z = A^T + A. I have left the old way
    //  commented out as it avoids an addtional 'multiply_all'. 
    //  However, it does this at the expense of transposing the 
    //  entire matrix. I would guess that the transpose costs 
    //  more than any savings. -rst
    // 
    //     Scaled_A_trans = ML_Operator_Create(P_0->comm);
    //     ML_Operator_Transpose_byrow(Scaled_A, Scaled_A_trans);
    //     Z = ML_Operator_Create(P_0->comm);
    //     ML_Operator_Add(Scaled_A, Scaled_A_trans, Z, ML_CSR_MATRIX, 1.0);
    //     P_second = ML_Operator_Create(P_0->comm);
    //     ML_2matmult(Z, P_prime, P_second, ML_CSR_MATRIX);
    //     multiply_all(P_0, P_second, &den[0]);
    //     multiply_all(P_prime, P_second, &den[0]);
    break;

  default:
    // should never be here
    cerr << "Incorrect parameter (" << ag->minimizing_energy << ")" << endl;
    cerr << "(file " << __FILE__ << ", line " << __LINE__ << ")" << endl;
    exit(EXIT_FAILURE);
  }

  /* convert from root points to all the coarse points */

  double temp_omega, temp_den, temp_num;
  int flag;
  if (compress == 1) {
t1 =  GetClock(); printf("start expand = %e\n",t1-t0);  t0 = t1;

    for (int row = 0 ; row < temp2->invec_leng ; row++) {
      flag = -1;
      temp_omega = 1e+20; temp_num = -1.e10; temp_den = -1.;
      ML_get_matrix_row(temp2, 1, &row, &allocated, &bindx, &val,
			&row_length, 0); 
      for  (int j = 0; j < row_length; j++)     {
	if ( (flag != 1) && ( root_pts[bindx[j]] != -1 )) {
	  if ( bindx[j] == row ) {
	    temp_num = tnum[root_pts[bindx[j]]];
	    temp_den = tden[root_pts[bindx[j]]];
	    flag = 1;
	  }
	  else {
	    if ( tnum[root_pts[bindx[j]]]/tden[root_pts[bindx[j]]] <
		 temp_num/temp_den) {
	      temp_num = tnum[root_pts[bindx[j]]];
	      temp_den = tden[root_pts[bindx[j]]];
	    }
	  }
	}
      }
      if (temp_den == -1.) { printf("problems\n"); exit(-1); }
      num[row] = temp_num; den[row] = temp_den;
    }
t1 =  GetClock(); printf("end expand = %e\n",t1-t0);  t0 = t1;
  }
  int zero_local   = 0;
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

  double min_all  = ML_gmin_double(min_local, P_0->comm);
  double max_all  = ML_gmax_double(max_local, P_0->comm);
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

  //  vector<double> RowOmega(n);
  double *RowOmega;
  RowOmega = (double *) ML_allocate(sizeof(double)*(n+1));
  ag->old_RowOmegas = RowOmega;
  for (int i = 0 ; i < n ; i++) RowOmega[i] = DBL_MAX;

  for (int row = 0 ; row < n ; row++) 
  {
    ML_get_matrix_row(P_prime, 1, &row, &allocated, &bindx, &val,
                      &row_length, 0); // should really be P0 + Pprime

    for  (int j = 0; j < row_length; j++) 
    {
      if (bindx[j] < n) // this is white magic
      {
	double omega = ColOmega[bindx[j]];
        if (omega < RowOmega[row]) RowOmega[row] = omega;
      }
    }
  }

  ML_Operator* Scaled_P_prime = 0;

  Scaled_P_prime = ML_Operator_ImplicitlyVScale(P_prime, &RowOmega[0], 0);
  ML_CommInfoOP_Clone(&(Scaled_P_prime->getrow->pre_comm), P_prime->getrow->pre_comm);

  ML_Operator_Add(P_0, Scaled_P_prime, &(ml->Pmat[clevel]), 
                  ML_CSR_MATRIX, -1.0);
  //  ML_CSR_DropSmall(&(ml->Pmat[clevel]), 1.e-2);

  ML_Operator_Set_1Levels(&(ml->Pmat[clevel]), &(ml->SingleLevel[clevel]),
                          &(ml->SingleLevel[level]));

  if (bindx != NULL) ML_free(bindx);
  if (val   != NULL) ML_free(val);

  if (DinvAmat != NULL) ML_Operator_Destroy(&DinvAmat);
  if (Pprime_gz != P_prime) ML_Operator_Destroy(&Pprime_gz);
  if (temp2)          ML_Operator_Destroy(&temp2);
  if (temp1)          ML_Operator_Destroy(&temp1);
  if (Z)              ML_Operator_Destroy(&Z);
  if (P_prime)        ML_Operator_Destroy(&P_prime);
  if (P_second)       ML_Operator_Destroy(&P_second);
  if (Scaled_P_prime) ML_Operator_Destroy(&Scaled_P_prime);
  if (root_pts != NULL) ML_free(root_pts);

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

  double *RowOmega = ag->old_RowOmegas;

  /* ============================================================ */
  /* Start construction Pmatrix to minimize the energy of         */
  /* each basis function. This requires two additional operators, */
  /* here called P_prime and P_second.                            */
  /* ============================================================ */

  ML_Operator* P0_trans = ML_Operator_Create(P_0->comm);
  ML_Operator* P0TA     = ML_Operator_Create(P_0->comm);

  ML_Operator_Transpose_byrow(P_0, P0_trans);
  ML_2matmult(P0_trans, Amat, P0TA, ML_CSR_MATRIX);
  ML_CSR_DropSmall(P0TA, 1.e-5, 1e-5, 1e-5);
  ML_Operator* ttt = ML_Operator_ImplicitlyBlockDinvScale(Amat);

  ML_AGG_DinvP(P0TA, (MLSthing *) ttt->data, Amat->num_PDEs,1,0);
  ML_Operator* Scaled_P0TA = 0;
  Scaled_P0TA = ML_Operator_ImplicitlyVCScale(P0TA, &(RowOmega[0]), 0);

#ifdef out
  ML_CommInfoOP_Clone(&(Scaled_P_prime->getrow->pre_comm), 
                      P_prime->getrow->pre_comm);
#endif

  ML_Operator_Add(P0_trans,Scaled_P0TA,&(ml->Rmat[level]),ML_CSR_MATRIX,-1.0);

  //  ML_CSR_DropSmall(&(ml->Rmat[level]), 1.e-2);


  ML_Operator_Set_1Levels(&(ml->Rmat[level]), &(ml->SingleLevel[level]), 
                          &(ml->SingleLevel[clevel]));



  if ( ag->old_RowOmegas != NULL) ML_free(ag->old_RowOmegas);
  ag->old_RowOmegas = NULL;

  Dinv_size = -1;

#ifdef ML_TIMING
  ml->Rmat[level].build_time =  GetClock() - t0;
  ml->timing->total_build_time += ml->Rmat[level].build_time;
#endif


  if (ttt != NULL) ML_Operator_Destroy(&ttt);
  if (P0_trans != NULL) ML_Operator_Destroy(&P0_trans);  
  if (P0TA != NULL) ML_Operator_Destroy(&P0TA);
  if (Scaled_P0TA != NULL) ML_Operator_Destroy(&Scaled_P0TA);

  return(0);
}
