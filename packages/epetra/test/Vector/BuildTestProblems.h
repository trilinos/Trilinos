#include "Epetra_Map.h"
#include "Epetra_Vector.h"

int  BuildMatrixTests (Epetra_Vector & C,
			     const char transa, const char transb, 
			     const double alpha, 
			     Epetra_Vector& A, 
			     Epetra_Vector& B,
			     const double beta,
			     Epetra_Vector& C_GEMM );

  
int  BuildVectorTests (Epetra_Vector & C,
				const double alpha, 
				Epetra_Vector& A, 
				Epetra_Vector& sqrtA,
				Epetra_Vector& B,
				Epetra_Vector& C_alphaA,
				Epetra_Vector& C_alphaAplusB,
				Epetra_Vector& C_plusB,
				double* const dotvec_AB,
				double* const norm1_A,
				double* const norm2_sqrtA,
				double* const norminf_A,
				double* const normw_A,
				Epetra_Vector& Weights,
				double* const minval_A,
				double* const maxval_A,
				double* const meanval_A );  
  
