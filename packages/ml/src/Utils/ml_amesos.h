/* ************************************************************************* */
/* See the file COPYRIGHT for a complete copyright notice, contact person,   */
/* and disclaimer.                                                           */
/* ************************************************************************* */

#ifndef __ML_AMESOS_H__
#define __ML_AMESOS_H__

#include <stdlib.h>
#include <assert.h>
#include "ml_struct.h"

#define ML_AMESOS_KLU            0
#define ML_AMESOS_UMFPACK        1
#define ML_AMESOS_SUPERLUDIST    2
#define ML_AMESOS_MUMPS          3
#define ML_AMESOS_SCALAPACK      4


#ifndef ML_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

  /** Applies the Amesos direct solver to the given vector. */
  extern int ML_Smoother_Amesos(ML_Smoother *sm,int inlen,double x[],int outlen,
				double rhs[]);

  /** Clean memory associated to Amesos_Handle. */
  void ML_Smoother_Clean_Amesos(void *Amesos_Handle);
  
  /** Generates the coarse solver using Amesos. */
  /*! Generates the coarse solver using one of the Amesos supported
    direct solvers.

    \param \inout ml: ML_Structure;

    \param \in nl: level for which we have to define the coarse solver;

    \param \in AmesosSolver: integer variable, that can be ML_AMESOS_KLU,
    ML_AMESOS_UMFPACK, ML_AMESOS_SUPERLUDIST, ML_AMESOS_MUMPS,
    ML_AMESOS_SCALAPACK;

    \param \in MaxProcs: integer defining the maximum number of
    processes to use in the coarse solution (only for some of the
    supported Amesos solvers).
  */
  
int ML_Gen_Smoother_Amesos(ML *ml, int nl, int AmesosSolver, int MaxProcs);


#ifndef ML_CPP
#ifdef __cplusplus
   }
#endif
#endif

#endif /* #ifndef __ML_AMESOS_H__ */
