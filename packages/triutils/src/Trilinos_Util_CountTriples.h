#include <vector>
#include "Epetra_Comm.h"

void Trilinos_Util_CountTriples( const char *data_file, 
				 bool symmetric, 
				 vector<int> &non_zeros,
				 int &N_rows, int &nnz, 
				 const Epetra_Comm  &comm) ;
