#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "../epetra_test_err.h"
int checkmap(Epetra_BlockMap & Map, int NumGlobalElements, int NumMyElements,
	     int * MyGlobalElements, int ElementSize, int * ElementSizeList, 
	     int NumGlobalEquations, int NumMyEquations,
          int IndexBase, Epetra_Comm & Comm,
          bool DistributedGlobal);

