#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "../epetra_test_err.h"
int checkmap(Epetra_Map & Map, int NumGlobalElements, int NumMyElements,
	     int * MyGlobalElements, int IndexBase, Epetra_Comm & Comm,
          bool DistributedGlobal);

