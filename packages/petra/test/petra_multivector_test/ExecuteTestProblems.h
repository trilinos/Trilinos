#include "Petra_Petra.h"
#include "Petra_BlockMap.h"
#include "Petra_Map.h"
#include "Petra_LocalMap.h"
#include "Petra_RDP_MultiVector.h"
#include "BuildTestProblems.h"
int RDP_MatrixTests(const Petra_BlockMap & map, const Petra_LocalMap & LocalMap, int NumVectors,
		    bool verbose);

int RDP_MultiVectorTests(const Petra_BlockMap & Map, int NumVectors, bool verbose);

int BadResidual(double * Residual, int NumVectors);
