/********************************************************************* */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */   
/* ******************************************************************** */

/********************************************************************* */
/*          Decomposition with METIS                                  */
/********************************************************************* */

#ifndef __MLAGGPARMETIS__
#define __MLAGGPARMETIS__

extern int ML_Aggregate_Set_ReqLocalCoarseSize( ML *ml, ML_Aggregate *ag, 
					 int level,
					 int desired_aggre_per_proc );
					 
#ifndef ML_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

#include "ml_comm.h"

  extern ML_Operator * ML_BuildQt( int StartingNumElements,
				   int ReorderedNumElements,
				   int reordered_decomposition[],
				   USR_COMM mpi_communicator,
				   ML_Comm *ml_communicator );
  
#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif /* #ifndef __MLAGGPARMETIS__ */
