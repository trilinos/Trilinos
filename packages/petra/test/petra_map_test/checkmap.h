#include "Petra_Comm.h"
#include "Petra_Map.h"
int checkmap(Petra_Map & Map, int NumGlobalElements, int NumMyElements,
	     int * MyGlobalElements, int IndexBase, Petra_Comm & Comm,
          bool DistributedGlobal);

