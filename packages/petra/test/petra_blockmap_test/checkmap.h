#include "Petra_Comm.h"
#include "Petra_Map.h"
int checkmap(Petra_BlockMap & Map, int NumGlobalElements, int NumMyElements,
	     int * MyGlobalElements, int ElementSize, int * ElementSizeList, 
	     int NumGlobalEquations, int NumMyEquations,
          int IndexBase, Petra_Comm & Comm,
          bool DistributedGlobal);

