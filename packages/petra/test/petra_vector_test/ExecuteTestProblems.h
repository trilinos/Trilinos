#include "Petra_Petra.h"
#include "Petra_BlockMap.h"
#include "Petra_Map.h"
#include "Petra_LocalMap.h"
#include "Petra_RDP_Vector.h"
#include "BuildTestProblems.h"
int RDP_MatrixTests(const Petra_BlockMap & map, const Petra_LocalMap & LocalMap,
		    bool verbose);

int RDP_VectorTests(const Petra_BlockMap & Map, bool verbose);

int BadResidual(double * Residual);
