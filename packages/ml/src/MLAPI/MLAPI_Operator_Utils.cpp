#include "ml_common.h"
#if defined(HAVE_ML_MLAPI)

#include "ml_include.h"
#include <iostream>
#include "ml_operator.h"
#include "ml_epetra.h"
#include "ml_amesos.h"
#include "ml_epetra_utils.h"
#include "ml_amesos_wrap.h"
#ifdef HAVE_ML_ANASAxI
#include "ml_anasazi.h"
#endif
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "MLAPI_Error.h"
#include "MLAPI_Space.h"
#include "MLAPI_MultiVector.h"
#include "MLAPI_Operator.h"
#include "Ifpack_Utils.h"
#ifdef MB_MODIF_QR
#include "ml_qr_fix.h"
#endif

using namespace std;

namespace Teuchos {
  class ParameterList;
}

namespace MLAPI {

// ====================================================================== 

Operator GetRAP(const Operator& R, const Operator& A, 
                const Operator& P)
{
  ML_Operator* Rmat = R.GetML_Operator();
  ML_Operator* Amat = A.GetML_Operator();
  ML_Operator* Pmat = P.GetML_Operator();
  ML_Operator* result = 0;

  result = ML_Operator_Create (Rmat->comm);

/* The fixing of coarse matrix only works if it is in MSR format */
   int myMatrixType = ML_MSR_MATRIX;
   ML_rap(Rmat, Amat, Pmat, result, myMatrixType);
   result->num_PDEs = Pmat->num_PDEs;
#ifdef  MB_MODIF_QR
   ML_fixCoarseMtx(result, myMatrixType);
#endif/*MB_MODIF_QR*/

  Operator op(P.GetDomainSpace(),P.GetDomainSpace(), result);
  return(op);
}

#include "ml_aggregate.h"
#include "ml_agg_METIS.h"

// ====================================================================== 
Operator GetTranspose(const Operator& A, const bool byrow = true) 
{
  ML_Operator* ML_transp;
  ML_transp = ML_Operator_Create(GetML_Comm());
  if (byrow)
    ML_Operator_Transpose_byrow(A.GetML_Operator(),ML_transp);
  else
    ML_Operator_Transpose(A.GetML_Operator(),ML_transp);

  Operator transp(A.GetRangeSpace(),A.GetDomainSpace(), ML_transp,true);
  return(transp);
}

// ====================================================================== 
Operator GetIdentity(const Space& DomainSpace, const Space& RangeSpace)
{
  ML_Operator* ML_eye = ML_Operator_Create(GetML_Comm());
  int size = DomainSpace.GetNumMyElements();
  ML_Operator_Set_ApplyFuncData(ML_eye, size, size,
            NULL, size, eye_matvec, 0);
  ML_Operator_Set_Getrow(ML_eye, size, eye_getrows);
  Operator* eye = new Operator(DomainSpace,DomainSpace,ML_eye,true);
  return(*eye);
}

// ====================================================================== 
MultiVector GetDiagonal(const Operator& A)
{
  // FIXME
  if (A.GetDomainSpace() != A.GetRangeSpace())
    ML_THROW("Currently only square matrices are supported", -1);

  MultiVector D(A.GetDomainSpace());
  D = 0.0;
  
  ML_Operator* matrix = A.GetML_Operator();

  if (matrix->getrow == NULL)
    ML_THROW("getrow() not set!", -1);

  int row_length;
  int allocated = 128;
  int*    bindx = (int    *)  ML_allocate(allocated*sizeof(int   ));
  double* val   = (double *)  ML_allocate(allocated*sizeof(double));

  for (int i = 0 ; i < matrix->getrow->Nrows; i++) {
    ML_get_matrix_row(matrix, 1, &i, &allocated, &bindx, &val,
                      &row_length, 0);
    for  (int j = 0; j < row_length; j++) {
      D(i) = 0.0;
      if (i == bindx[j]) {
        D(i) = val[j];
        break;
      }
    }
  }

  ML_free(val);
  ML_free(bindx);
  return (D);

}

// ====================================================================== 
MultiVector GetDiagonal(const Operator& A, const int offset)
{
  // FIXME
  if (A.GetDomainSpace() != A.GetRangeSpace())
    ML_THROW("Currently only square matrices are supported", -1);

  MultiVector D(A.GetDomainSpace());
  D = 0.0;
  
  ML_Operator* matrix = A.GetML_Operator();

  if (matrix->getrow == NULL)
    ML_THROW("getrow() not set!", -1);

  int row_length;
  int allocated = 128;
  int*    bindx = (int    *)  ML_allocate(allocated*sizeof(int   ));
  double* val   = (double *)  ML_allocate(allocated*sizeof(double));

  for (int i = 0 ; i < matrix->getrow->Nrows; i++) {
    int GlobalRow = A.GetGRID(i);
    ML_get_matrix_row(matrix, 1, &i, &allocated, &bindx, &val,
                      &row_length, 0);
    for  (int j = 0; j < row_length; j++) {
      D(i) = 0.0;
      if (A.GetGCID(bindx[j]) == GlobalRow + offset) {
        D(i) = val[j];
        break;
      }
    }
  }

  ML_free(val);
  ML_free(bindx);
  return (D);

}

// ====================================================================== 
// DIAGONAL OPERATOR
// - takes a MultiVector (w/ only one vector) in input
// - allocates memory to store the diagonal as a double vector
// - sets the pointers
// - the destructor deletes the double vector
// ====================================================================== 

static int diag_matvec(ML_Operator *Amat_in, int ilen, double p[], 
                int olen, double ap[])
{
  double* D = (double*)Amat_in->data;
  
  for (int i = 0; i < olen; i++) ap[i] = D[i] * p[i];

  return(1);
}

static int diag_getrows(ML_Operator *data, int N_requested_rows, int requested_rows[],
                 int allocated_space, int columns[], double values[],
                 int row_lengths[])
{
  double* D = (double*)data->data;

  if (allocated_space < N_requested_rows) {
    ML_avoid_unused_param(data);
    return(0);
  }

  for (int i = 0; i < N_requested_rows; i++) {
    row_lengths[i] = 1;
    columns[i]     = requested_rows[i];
    values[i]      = D[requested_rows[i]];
  }
  return(1);
}

static void diag_destroy(void* data)
{
  double* D = (double*)data;
  delete[] D;
}

// ====================================================================== 
Operator GetDiagonal(const MultiVector& D)
{
  if (D.GetNumVectors() != 1)
    ML_THROW("D.GetNumVectors() != 1", -1);

  int size = D.GetMyLength();
  if (size == 0)
    ML_THROW("empty diagonal vector in input", -1);

  double* diag = new double[size];
  for (int i = 0 ; i < size ; ++i) 
    diag[i] = D(i);

  // creates the ML operator and store the diag pointer,
  // as well as the function pointers
  ML_Operator* MLDiag = ML_Operator_Create(GetML_Comm());

  MLDiag->invec_leng = size;
  MLDiag->outvec_leng = size;
  MLDiag->data = (void*)diag;
  MLDiag->matvec->func_ptr = diag_matvec;

  MLDiag->matvec->ML_id = ML_NONEMPTY;
  MLDiag->matvec->Nrows = size;
  MLDiag->from_an_ml_operator = 0;
  MLDiag->data_destroy = diag_destroy;

  MLDiag->getrow->func_ptr = diag_getrows;

  MLDiag->getrow->ML_id = ML_NONEMPTY;
  MLDiag->getrow->Nrows = size;

  // creates the MLAPI wrapper
  Operator Diag(D.GetVectorSpace(),D.GetVectorSpace(),MLDiag,true);
  return(Diag);
}

// ====================================================================== 
static void widget_destroy(void* data)
{
  struct ML_AGG_Matrix_Context* widget = (struct ML_AGG_Matrix_Context*)data;
  delete widget;
}

// ====================================================================== 
Operator GetJacobiIterationOperator(const Operator& Amat, double Damping)
{

  struct ML_AGG_Matrix_Context* widget = new struct ML_AGG_Matrix_Context;
  widget->near_bdry = 0;
  widget->aggr_info = 0;
  widget->drop_tol  = 0.0;

  widget->Amat = Amat.GetML_Operator();
  widget->omega = Damping;

  ML_Operator* tmp_ML = ML_Operator_Create(GetML_Comm());
  ML_Operator_Set_ApplyFuncData(tmp_ML, widget->Amat->invec_leng,
                                widget->Amat->outvec_leng, widget,
                                widget->Amat->matvec->Nrows, NULL, 0);

  tmp_ML->data_destroy = widget_destroy;

  ML_Operator_Set_Getrow(tmp_ML, widget->Amat->getrow->Nrows, 
                         ML_AGG_JacobiSmoother_Getrows);

  // Creates a new copy of pre_comm, so that the old pre_comm
  // can be destroyed without worry
  ML_CommInfoOP_Clone(&(tmp_ML->getrow->pre_comm),
                      widget->Amat->getrow->pre_comm);

  Operator tmp(Amat.GetDomainSpace(), Amat.GetRangeSpace(), tmp_ML, true,
               Amat.GetRCPOperatorBox());

  return(tmp);
}

// ====================================================================== 
static int Ptent1D_matvec(ML_Operator *Amat_in, int ilen, double p[], 
                int olen, double ap[])
{
  double* D = (double*)Amat_in->data;
  
  for (int i = 0; i < olen; i++) ap[i] = D[i] * p[i / 3];

  return(1);
}

// ====================================================================== 
static int Ptent1D_getrows(ML_Operator *data, int N_requested_rows, int requested_rows[],
                 int allocated_space, int columns[], double values[],
                 int row_lengths[])
{
  double* D = (double*)data->data;

  if (allocated_space < N_requested_rows) {
    ML_avoid_unused_param(data);
    return(0);
  }

  for (int i = 0; i < N_requested_rows; i++) {
    row_lengths[i] = 1;
    columns[i]     = requested_rows[i] / 3;
    values[i]      = D[i];
  }
  return(1);
}

static void Ptent1D_destroy(void* data)
{
  double* D = (double*)data;
  delete[] D;
}

// ====================================================================== 
Operator GetPtent1D(const MultiVector& D, const int offset = 0)
{
  if (D.GetNumVectors() != 1)
    ML_THROW("D.GetNumVectors() != 1", -1);

  int size = D.GetMyLength();
  if (size == 0)
    ML_THROW("empty diagonal vector in input", -1);

  double* diag = new double[size];
  for (int i = 0 ; i < size ; ++i)
    diag[i] = D(i);

  // creates the ML operator and store the diag pointer,
  // as well as the function pointers
  ML_Operator* MLDiag = ML_Operator_Create(GetML_Comm());

  int invec_leng = size / 3 + size % 3;
  int outvec_leng = size;

  MLDiag->invec_leng = invec_leng;
  MLDiag->outvec_leng = outvec_leng;
  MLDiag->data = (void*)diag;
  MLDiag->data_destroy = Ptent1D_destroy;
  MLDiag->matvec->func_ptr = Ptent1D_matvec;

  MLDiag->matvec->ML_id = ML_NONEMPTY;
  MLDiag->matvec->Nrows = outvec_leng;
  MLDiag->from_an_ml_operator = 0;

  MLDiag->getrow->func_ptr = Ptent1D_getrows;

  MLDiag->getrow->ML_id = ML_NONEMPTY;
  MLDiag->getrow->Nrows = outvec_leng;

  // creates the domain space
  vector<int> MyGlobalElements(invec_leng);
  for (int i = 0 ; i < invec_leng ; ++i) 
    MyGlobalElements[i] = D.GetVectorSpace()(i * 3) / 3;
  Space DomainSpace(invec_leng, -1, &MyGlobalElements[0]);
  Space RangeSpace = D.GetVectorSpace();

  // creates the MLAPI wrapper
  Operator Diag(DomainSpace,RangeSpace,MLDiag,true);
  return(Diag);
}

// ====================================================================== 
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
	ML_hash_it(global_col, hashed_inds, index_length,&hash_used,&hash_val);
        hashed_inds[hash_val] = global_col;
        hashed_vals[hash_val] += scalarA * A_val[j];
	A_bindx[j] = hash_val;
      }

      ML_get_matrix_row(B, 1, &i, &B_allocated, &B_bindx, &B_val,
                        &B_length, 0);
      for (j = 0; j < B_length; j++) {
	global_col = B_gids[B_bindx[j]];
	ML_hash_it(global_col, hashed_inds, index_length,&hash_used, &hash_val);
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

  columns = 0;
  values = 0;

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

  nz_ptr = 0;
  rowptr[0] = 0;
  for (i = 0 ; i < A->getrow->Nrows; i++) {
    hash_used = 0;
      ML_get_matrix_row(A, 1, &i, &A_allocated, &A_bindx, &A_val,
                        &A_length, 0);
      for (j = 0; j < A_length; j++) {
	global_col = A_gids[A_bindx[j]];
	ML_hash_it(global_col, hashed_inds, index_length,&hash_used, &hash_val);
        hashed_inds[hash_val] = global_col;
        hashed_vals[hash_val] += scalarA * A_val[j];
	A_bindx[j] = hash_val;
      }

      ML_get_matrix_row(B, 1, &i, &B_allocated, &B_bindx, &B_val,
                        &B_length, 0);
      for (j = 0; j < B_length; j++) {
	global_col = B_gids[B_bindx[j]];
	ML_hash_it(global_col, hashed_inds, index_length,&hash_used, &hash_val);
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

// ====================================================================== 
void AnalyzeCheap(const Operator& A) 
{
  Ifpack_Analyze(*(A.GetRowMatrix()));
}

// ====================================================================== 
void PrintSparsity(const Operator& A, int NumPDEEquations)
{
  string FileName = A.GetLabel() + ".ps";
  Ifpack_PrintSparsity(*(A.GetRowMatrix()), FileName.c_str(),
                       NumPDEEquations);
}

// ====================================================================== 
Operator GetScaledOperator(const Operator& A, const double alpha) 
{
  ML_Operator* ScaledA = 0;
  ScaledA = ML_Operator_ExplicitlyScale(A.GetML_Operator(),
                                        (double)alpha);

  if (ScaledA == 0)
    ML_THROW("ML_Operator_ExplicitlyScale returned 0", -1);

  Operator res;
  res.Reshape(A.GetDomainSpace(), A.GetRangeSpace(), ScaledA,
              true, A.GetRCPOperatorBox());
    
  return(res);
}

// ====================================================================== 
Operator Duplicate(const Operator& A)
{
  return(GetScaledOperator(A, 1.0));
}

} // namespace MLAPI

#endif // HAVE_ML_MLAPI
