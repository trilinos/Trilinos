#include "Petra_Petra.h"
#include "Petra_Map.h"
#include "Petra_RDP_MultiVector.h"

int  BuildRDP_MatrixTests (Petra_RDP_MultiVector & C,
			     const char transa, const char transb, 
			     const double alpha, 
			     Petra_RDP_MultiVector& A, 
			     Petra_RDP_MultiVector& B,
			     const double beta,
			     Petra_RDP_MultiVector& C_GEMM );

  
int  BuildRDP_MultiVectorTests (Petra_RDP_MultiVector & C,
				const double alpha, 
				Petra_RDP_MultiVector& A, 
				Petra_RDP_MultiVector& sqrtA,
				Petra_RDP_MultiVector& B,
				Petra_RDP_MultiVector& C_alphaA,
				Petra_RDP_MultiVector& C_alphaAplusB,
				Petra_RDP_MultiVector& C_plusB,
				double* const dotvec_AB,
				double* const norm1_A,
				double* const norm2_sqrtA,
				double* const norminf_A,
				double* const normw_A,
				Petra_RDP_MultiVector& Weights,
				double* const minval_A,
				double* const maxval_A,
				double* const meanval_A );  
  
