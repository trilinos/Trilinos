#ifndef ML_OPERATOR_UTILS_H
#define ML_OPERATOR_UTILS_H
#include "ml_include.h"
#include <iostream>
#include "ml_operator.h"
#include "ml_epetra.h"
#include "ml_amesos.h"
#include "ml_epetra_utils.h"
#include "ml_amesos_wrap.h"
#ifdef HAVE_ML_ANASAZI
#include "ml_anasazi.h"
#endif
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "MLAPI_Space.h"
#include "MLAPI_MultiVector.h"
#include "MLAPI_Operator.h"

/*!
\file MLAPI_Operator_Utils

\brief Functions to create and compose MLAPI::Operator's.

\author Marzio Sala, SNL 9214.

\date Last updated on Feb-05.
*/

using namespace std;

namespace Teuchos {
  class ParameterList;
}

namespace MLAPI {

//! Performs a triple matrix-matrix product, res = R * A *P.
Operator RAP(const Operator& R, const Operator& A, 
             const Operator& P)
{
  ML_Operator* Rmat = R.GetML_Operator();
  ML_Operator* Amat = A.GetML_Operator();
  ML_Operator* Pmat = P.GetML_Operator();
  ML_Operator* result = 0;

  result = ML_Operator_Create (Rmat->comm);

  ML_rap(Rmat, Amat, Pmat, result, MatrixType);

  Operator op(P.GetDomainSpace(),P.GetDomainSpace(), result);
  return(op);
}

#include "ml_aggregate.h"
#include "ml_agg_METIS.h"

// DELETE THIS CRAP
#if 0
//! Builds the tentative prolongator using aggregation.
Operator BuildP(const Operator& A, Teuchos::ParameterList& List)
{
  ML_Aggregate* agg_object;
  ML_Aggregate_Create(&agg_object);
  ML_Aggregate_Set_MaxLevels(agg_object,2);
  ML_Aggregate_Set_StartLevel(agg_object,0);
  ML_Aggregate_Set_Threshold(agg_object,0.0);
  agg_object->curr_threshold = 0.0;
  
  ML_Operator* ML_Ptent = 0;
  ML_Ptent = ML_Operator_Create(GetML_Comm());
  string CoarsenType = List.get("aggregation: type","MIS");

  // FIXME: as in MLP
  int NullSpaceDim = List.get("nullspace dimension", 1);
  double* NullSpace = List.get("nullspace", (double*)0);

  int size = A.DomainSpace().NumMyElements();
  if (NullSpace)
  {
    ML_memory_alloc((void **)&(agg_object->nullspace_vect), 
                    sizeof(double) * size, "ns");
    for (int i = 0 ; i < size * NullSpaceDim ; ++i)
      agg_object->nullspace_vect[i] = NullSpace[i];
  }
  agg_object->nullspace_dim = NullSpaceDim;

  int NewSize;
  
  if (CoarsenType == "Uncoupled") {
    NewSize = ML_Aggregate_CoarsenUncoupled(agg_object, A.GetData(),
                                            &ML_Ptent, GetML_Comm());
  }
  else if (CoarsenType == "MIS") {
    NewSize = ML_Aggregate_CoarsenMIS(agg_object, A.GetData(),
                                      &ML_Ptent, GetML_Comm());
  }
  else if (CoarsenType == "METIS") {
    int NodesPerAggr = List.get("aggregation: nodes per aggregate",1);
    ML ml_object;
    ml_object.ML_num_levels = 1;
    ML_Aggregate_Set_NodesPerAggr(&ml_object,agg_object,0,NodesPerAggr);
    NewSize = ML_Aggregate_CoarsenMETIS(agg_object, A.GetData(),
                                      &ML_Ptent, GetML_Comm());
  }
  else {
    cerr << "ERROR: In BuildP()" << endl;
    cerr << "ERROR: (file " << __FILE__ << ", line " << __LINE__ << ")" << endl;
    cerr << "ERROR: Requested aggregation schem (" << CoarsenType
         << ") not recognized." << endl;
    throw(-1);
  }

  int NumMyElements = NewSize;
  Space CoarseSpace(-1,NumMyElements);
  Operator Ptent(CoarseSpace,A.RangeSpace(),ML_Ptent,true);

  if (NullSpace == 0)
  {
    NullSpace = new double[NewSize];
    List.set("nullspace",NullSpace);
  }

  // store next-level nullspace vector in list
  List.set("nullspace dimension", NullSpaceDim);
  for (int i = 0 ; i < NewSize ; ++i)
    NullSpace[i] = agg_object->nullspace_vect[i];

  ML_Aggregate_Destroy(&agg_object);

  return(Ptent);

}
#endif

//! Returns a newly created transpose of \c A.
Operator Transpose(const Operator& A) 
{
  ML_Operator* ML_transp;
  ML_transp = ML_Operator_Create(GetML_Comm());
  ML_Operator_Transpose_byrow(A.GetML_Operator(),ML_transp);

  Operator transp(A.GetRangeSpace(),A.GetDomainSpace(), ML_transp,true);
  return(transp);
}

//! Returns the identity matrix.
Operator Identity(const Space& DomainSpace, const Space& RangeSpace)
{
  ML_Operator* ML_eye = ML_Operator_Create(GetML_Comm());
  int size = DomainSpace.GetNumMyElements();
  ML_Operator_Set_ApplyFuncData(ML_eye, size, size,
            NULL, size, eye_matvec, 0);
  ML_Operator_Set_Getrow(ML_eye, size, eye_getrows);
  Operator* eye = new Operator(DomainSpace,DomainSpace,ML_eye,true);
  return(*eye);
}

int diag_matvec(ML_Operator *Amat_in, int ilen, double p[], 
                int olen, double ap[])
{
  MultiVector* D = (MultiVector*)Amat_in->data;
  
  for (int i = 0; i < olen; i++) ap[i] = (*D)(i) * p[i];

  return(1);
}

int diag_getrows(ML_Operator *data, int N_requested_rows, int requested_rows[],
                 int allocated_space, int columns[], double values[],
                 int row_lengths[])
{
  MultiVector* D = (MultiVector*)data->data;

  if (allocated_space < N_requested_rows) {
    ML_avoid_unused_param(data);
    return(0);
  }

  for (int i = 0; i < N_requested_rows; i++) {
    row_lengths[i] = 1;
    columns[i]     = requested_rows[i];
    values[i]      = (*D)(i);
  }
  return(1);
}

//! Returns a vector containing the diagonal elements of \c A.
MultiVector Diagonal(const Operator& A)
{
  // FIXME
  if (A.GetDomainSpace() != A.GetRangeSpace()) {
    cerr << "ERROR: In Diagonal()" << endl;
    cerr << "ERROR: (file " << __FILE__ << ", line " << __LINE__ << ")" << endl;
    cerr << "ERROR: Currently only square matrices are supported." << endl;
    throw(-1);
  }

  MultiVector D(A.GetDomainSpace());
  D = 0.0;
  
  ML_Operator* matrix = A.GetML_Operator();

  if (matrix->getrow == NULL) {
    cerr << "ERROR: In Diagonal()" << endl;
    cerr << "ERROR: (file " << __FILE__ << ", line " << __LINE__ << ")" << endl;
    cerr << "ERROR: getrow() not set!" << endl;
    throw(-1);
  }

  int row_length;
  int allocated = 128;
  int*    bindx = (int    *)  ML_allocate(allocated*sizeof(int   ));
  double* val   = (double *)  ML_allocate(allocated*sizeof(double));

  for (int i = 0 ; i < matrix->getrow->Nrows; i++) {
    ML_get_matrix_row(matrix, 1, &i, &allocated, &bindx, &val,
                      &row_length, 0);
    for  (int j = 0; j < row_length; j++) {
      if (bindx[j] == i)
        D(i) = val[j];
    }
  }

  ML_free(val);
  ML_free(bindx);
  return (D);

}

//! Returns a newly created operator, containing D on the diagonal.
Operator Diagonal(const MultiVector& D)
{
  ML_Operator* MLDiag = ML_Operator_Create(GetML_Comm());
  int size = D.GetMyLength();
  // FIXME: this is a memory leak!
  MultiVector* D2 = new MultiVector(D);
  ML_Operator_Set_ApplyFuncData(MLDiag, size, size,
            (void*)D2, size, diag_matvec, 0);
  ML_Operator_Set_Getrow(MLDiag, size, diag_getrows);
  Operator Diag(D.GetVectorSpace(),D.GetVectorSpace(),MLDiag,true);
  return(Diag);
}

//! Returns an operator defined as (I - Damping A).
Operator JacobiIterationOperator(const Operator& Amat, double Damping,
                                 struct ML_AGG_Matrix_Context* widget)
{

  widget->Amat = Amat.GetML_Operator();
  widget->omega  = Damping;
  ML_Operator* AGGsmoother = ML_Operator_Create(GetML_Comm());
  ML_Operator_Set_ApplyFuncData(AGGsmoother, widget->Amat->invec_leng,
                                widget->Amat->outvec_leng, widget,
                                widget->Amat->matvec->Nrows, NULL, 0);
  ML_Operator_Set_Getrow(AGGsmoother, 
                         widget->Amat->getrow->Nrows, 
                         ML_AGG_JacobiSmoother_Getrows);
  ML_CommInfoOP_Clone(&(AGGsmoother->getrow->pre_comm),
                      widget->Amat->getrow->pre_comm);

  Operator tmp(Amat.GetDomainSpace(), Amat.GetRangeSpace(), AGGsmoother, false);

  return(tmp);
}

//! Creates C = scalarA * A + scalarB * B.
int ML_Operator_Add2(ML_Operator *A, ML_Operator *B, ML_Operator *C,
		    int matrix_type, double scalarA, double scalarB)
{
  int A_allocated = 0, *A_bindx = NULL, B_allocated = 0, *B_bindx = NULL;
  double *A_val = NULL, *B_val = NULL, *hashed_vals;
  int i, A_length, B_length, *hashed_inds;
  int max_nz_per_row = 0, j;
  int hash_val, index_length;
  int *columns, *rowptr, nz_ptr, hash_used, global_col;
  double *values;
  struct ML_CSR_MSRdata *temp;
  int *A_gids, *B_gids;
  int max_per_proc;
#ifdef ML_WITH_EPETRA
  int count;
#endif

  if (A->getrow == NULL) 
    pr_error("ML_Operator_Add: A does not have a getrow function.\n");

  if (B->getrow == NULL) 
    pr_error("ML_Operator_Add: B does not have a getrow function.\n");

  if (A->getrow->Nrows != B->getrow->Nrows) {
    printf("ML_Operator_Add: Can not add, two matrices do not have the same");
    printf(" number of rows %d vs %d",A->getrow->Nrows,B->getrow->Nrows);
    exit(1);
  }

  if (A->invec_leng != B->invec_leng) {
    printf("ML_Operator_Add: Can not add, two matrices do not have the same");
    printf(" number of columns %d vs %d",A->getrow->Nrows,B->getrow->Nrows);
    exit(1);
  }

  /* let's just count some things */
  index_length = A->invec_leng + 1;
  if (A->getrow->pre_comm != NULL) {
    ML_CommInfoOP_Compute_TotalRcvLength(A->getrow->pre_comm);
    index_length += A->getrow->pre_comm->total_rcv_length;
  }
  if (B->getrow->pre_comm != NULL) {
    ML_CommInfoOP_Compute_TotalRcvLength(B->getrow->pre_comm);
    index_length += B->getrow->pre_comm->total_rcv_length;
  }

  ML_create_unique_col_id(A->invec_leng, &A_gids, A->getrow->pre_comm,
			  &max_per_proc,A->comm);
  ML_create_unique_col_id(B->invec_leng, &B_gids, B->getrow->pre_comm,
			  &max_per_proc,B->comm);


  hashed_inds = (int *) ML_allocate(sizeof(int)*index_length);
  hashed_vals = (double *) ML_allocate(sizeof(double)*index_length);

  for (i = 0; i < index_length; i++) hashed_inds[i] = -1;
  for (i = 0; i < index_length; i++) hashed_vals[i] = 0.;

  nz_ptr = 0;
  for (i = 0 ; i < A->getrow->Nrows; i++) {
    hash_used = 0;
      ML_get_matrix_row(A, 1, &i, &A_allocated, &A_bindx, &A_val,
                        &A_length, 0);
      for (j = 0; j < A_length; j++) {
	global_col = A_gids[A_bindx[j]];
	hash_val = ML_hash_it(global_col, hashed_inds, index_length,&hash_used);
        hashed_inds[hash_val] = global_col;
        hashed_vals[hash_val] += scalarA * A_val[j];
	A_bindx[j] = hash_val;
      }

      ML_get_matrix_row(B, 1, &i, &B_allocated, &B_bindx, &B_val,
                        &B_length, 0);
      for (j = 0; j < B_length; j++) {
	global_col = B_gids[B_bindx[j]];
	hash_val = ML_hash_it(global_col, hashed_inds, index_length,&hash_used);
        hashed_inds[hash_val] = global_col;
        hashed_vals[hash_val] += scalarB*B_val[j];
        B_bindx[j] = hash_val;
      }

      for (j = 0; j < A_length; j++) {
        nz_ptr++;
	hashed_inds[A_bindx[j]] = -1;
	hashed_vals[A_bindx[j]] = 0.;
      }
      for (j = 0; j < B_length; j++) {
        if (hashed_inds[B_bindx[j]] != -1) {
	  nz_ptr++;
	  hashed_inds[B_bindx[j]] = -1;
	  hashed_vals[B_bindx[j]] = 0.;
	}
      }
  }
  nz_ptr++;

  rowptr = (int    *) ML_allocate(sizeof(int)*(A->outvec_leng+1));
  if (matrix_type == ML_CSR_MATRIX) {
    columns= (int    *) ML_allocate(sizeof(int)*nz_ptr);
    values = (double *) ML_allocate(sizeof(double)*nz_ptr);
  }
#ifdef ML_WITH_EPETRA
  else if (matrix_type == ML_EpetraCRS_MATRIX) {
    columns= (int    *) ML_allocate(sizeof(int)*(index_length+1));
    values = (double *) ML_allocate(sizeof(double)*(index_length+1));
  }
#endif
  else {
    pr_error("ML_Operator_Add: Unknown matrix type\n");
  }


  columns= (int    *) ML_allocate(sizeof(int)*nz_ptr);
  values = (double *) ML_allocate(sizeof(double)*nz_ptr);
  if (values == NULL) pr_error("ML_Operator_Add: out of space\n");


  nz_ptr = 0;
  rowptr[0] = 0;
  for (i = 0 ; i < A->getrow->Nrows; i++) {
    hash_used = 0;
      ML_get_matrix_row(A, 1, &i, &A_allocated, &A_bindx, &A_val,
                        &A_length, 0);
      for (j = 0; j < A_length; j++) {
	global_col = A_gids[A_bindx[j]];
	hash_val = ML_hash_it(global_col, hashed_inds, index_length,&hash_used);
        hashed_inds[hash_val] = global_col;
        hashed_vals[hash_val] += scalarA * A_val[j];
	A_bindx[j] = hash_val;
      }

      ML_get_matrix_row(B, 1, &i, &B_allocated, &B_bindx, &B_val,
                        &B_length, 0);
      for (j = 0; j < B_length; j++) {
	global_col = B_gids[B_bindx[j]];
	hash_val = ML_hash_it(global_col, hashed_inds, index_length,&hash_used);
        hashed_inds[hash_val] = global_col;
        hashed_vals[hash_val] += scalarB*B_val[j];
        B_bindx[j] = hash_val;
      }
#ifdef ML_WITH_EPETRA
      if (matrix_type == ML_EpetraCRS_MATRIX) {
	for (j = 0; j < A_length; j++) {
	  columns[j] = hashed_inds[A_bindx[j]];
	  values[j]  = hashed_vals[A_bindx[j]];
	  nz_ptr++;
	  hashed_inds[A_bindx[j]] = -1;
	  hashed_vals[A_bindx[j]] = 0.;
	}
	count = A_length;
	for (j = 0; j < B_length; j++) {
	  if (hashed_inds[B_bindx[j]] != -1) {
	    columns[count] = hashed_inds[B_bindx[j]];
	    values[count++]  = hashed_vals[B_bindx[j]];
	    nz_ptr++;
	    hashed_inds[B_bindx[j]] = -1;
	    hashed_vals[B_bindx[j]] = 0.;
	  }
	}
	ML_Epetra_CRSinsert(C,i,columns,values,count);
      }
      else {
#endif
	for (j = 0; j < A_length; j++) {
	  columns[nz_ptr] = hashed_inds[A_bindx[j]];
	  values[nz_ptr]  = hashed_vals[A_bindx[j]];
	  nz_ptr++;
	  hashed_inds[A_bindx[j]] = -1;
	  hashed_vals[A_bindx[j]] = 0.;
	}
	for (j = 0; j < B_length; j++) {
	  if (hashed_inds[B_bindx[j]] != -1) {
	    columns[nz_ptr] = hashed_inds[B_bindx[j]];
	    values[nz_ptr]  = hashed_vals[B_bindx[j]];
	    nz_ptr++;
	    hashed_inds[B_bindx[j]] = -1;
	    hashed_vals[B_bindx[j]] = 0.;
	  }
	}
#ifdef ML_WITH_EPETRA
      }
#endif
      rowptr[i+1] = nz_ptr;
      if (rowptr[i+1] - rowptr[i] > max_nz_per_row)
	max_nz_per_row = rowptr[i+1] - rowptr[1];
  }
  if (matrix_type == ML_CSR_MATRIX) {
    temp = (struct ML_CSR_MSRdata *) ML_allocate(sizeof(struct ML_CSR_MSRdata));
    if (temp == NULL) pr_error("ML_Operator_Add: no space for temp\n");
    temp->columns = columns;
    temp->values  = values;
    temp->rowptr   = rowptr;

    ML_Operator_Set_ApplyFuncData(C, B->invec_leng, A->outvec_leng, 
				  temp,A->outvec_leng, NULL,0);
    ML_Operator_Set_Getrow(C, A->outvec_leng, CSR_getrow);
    ML_Operator_Set_ApplyFunc (C, CSR_matvec);
    ML_globalcsr2localcsr(C, max_per_proc);
    C->data_destroy = ML_CSR_MSRdata_Destroy;

    C->max_nz_per_row = max_nz_per_row;
    C->N_nonzeros     = nz_ptr;
  }
#ifdef ML_WITH_EPETRA
  else {
    ML_free(rowptr); 
    ML_free(columns);
    ML_free(values);
  }
#endif

  ML_free(A_gids);
  ML_free(B_gids);
  ML_free(hashed_vals);
  ML_free(hashed_inds);
  ML_free(A_val);
  ML_free(A_bindx);
  ML_free(B_val);
  ML_free(B_bindx);

  return 1;

}

void AnalyzeCheap(const Operator& A) 
{
  ML_Operator_Analyze(A.GetML_Operator(), const_cast<char*>(A.GetLabel().c_str()));
}

void PrintSparsity(const Operator& A, char* title,
                   char* FileName, int NumPDEEquations)
{
  ML_Operator_PrintSparsity(A.GetML_Operator(), title,
                            FileName, ML_TRUE, NumPDEEquations = 1);
}

} // namespace MLAPI
#endif // ML_OPERATOR_UTILS_H
