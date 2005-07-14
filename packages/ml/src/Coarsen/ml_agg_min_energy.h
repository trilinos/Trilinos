#ifndef ML_AGG_MIN_ENERGY
#define ML_AGG_MIN_ENERGY

#include "ml_common.h"
#include "ml_include.h"

#ifndef ML_CPP
#ifdef __cplusplus
extern "C" 
{
#endif
#endif

int ML_AGG_Gen_Prolongator_MinEnergy(ML *ml,int level, int clevel, void *data);
int ML_AGG_Gen_Restriction_MinEnergy(ML *ml,int level, int clevel, void *data);

#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif
#endif
