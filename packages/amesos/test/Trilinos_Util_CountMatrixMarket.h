#include <vector>
#include "Epetra_Comm.h"
#if 1
void Trilinos_Util_CountMatrixMarket( const char *data_file, 
				      vector<int> &non_zeros,
				      int &N_rows, int &nnz, 
				      const Epetra_Comm  &comm) ;

#else
void Trilinos_Util_CountMatrixMarket( const char *data_file,   vector<int> &non_zeros,
				      int &N_rows, int &nnz 
 );
#endif
