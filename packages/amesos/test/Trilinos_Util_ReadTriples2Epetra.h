#include <vector>
#include "Epetra_Comm.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

int Trilinos_Util_ReadTriples2Epetra( char *data_file,
				 bool symmetric, 
				 const Epetra_Comm &comm, 
				 Epetra_Map *& map, 
				 Epetra_CrsMatrix *& A, 
				 Epetra_Vector *& x, 
				 Epetra_Vector *& b,
				 Epetra_Vector *&xexact ) ;
