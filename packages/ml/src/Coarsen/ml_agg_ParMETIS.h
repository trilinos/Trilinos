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
  
  extern int ML_Aggregate_Set_ReqLocalCoarseSize( ML *ml, ML_Aggregate *ag, 
						  int level,
						  int desired_aggre_per_proc );
  extern int ML_Aggregate_CoarsenParMETIS( ML_Aggregate *ml_ag,
					   ML_Operator *Amatrix, 
					   ML_Operator **Pmatrix,
					   ML_Comm *comm);
  

  
#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif /* #ifndef __MLAGGPARMETIS__ */
