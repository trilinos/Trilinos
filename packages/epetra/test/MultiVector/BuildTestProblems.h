#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"

int  BuildMatrixTests (Epetra_MultiVector & C,
			     const char transa, const char transb, 
			     const double alpha, 
			     Epetra_MultiVector& A, 
			     Epetra_MultiVector& B,
			     const double beta,
			     Epetra_MultiVector& C_GEMM );

  
int  BuildMultiVectorTests (Epetra_MultiVector & C,
				const double alpha, 
				Epetra_MultiVector& A, 
				Epetra_MultiVector& sqrtA,
				Epetra_MultiVector& B,
				Epetra_MultiVector& C_alphaA,
				Epetra_MultiVector& C_alphaAplusB,
				Epetra_MultiVector& C_plusB,
				double* const dotvec_AB,
				double* const norm1_A,
				double* const norm2_sqrtA,
				double* const norminf_A,
				double* const normw_A,
				Epetra_MultiVector& Weights,
				double* const minval_A,
				double* const maxval_A,
				double* const meanval_A );  
  
