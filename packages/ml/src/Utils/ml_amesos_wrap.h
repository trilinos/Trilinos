

#ifndef _MLAMESOSWRAP_
#define _MLAMESOSWRAP_

#include "ml_include.h"

#ifndef ML_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

  /** Generates the direct solver. */
  /*! This function performs several operations:
    - wrap the ML_Operator for the given level to an Epetra_CrsMatrix
    - creates the Amesos object;
    - compute the symbolic and numeric factorzation.

    \param \inout ml: ML_Structure

    \param \in curr_level: specifies the level for which we have to
    define the direct solver;

    \param \in choice: integer variable, that can be ML_AMESOS_KLU,
    ML_AMESOS_UMFPACK, ML_AMESOS_SUPERLUDIST, ML_AMESOS_MUMPS,
    ML_AMESOS_SCALAPACK;

    \param \in MaxProcs: integer defining the maximum number of
    processes to use in the coarse solution (only for some of the
    supported Amesos solvers);

    \param \out: it will contain a pointer to the Amesos object (casted
    to void *).
  */
  int ML_Amesos_Gen(ML *ml, int curr_level, int choice,
		    int MaxProcs, void **Amesos_Handle);

  /** Solves using Amesos, and the factorization computed by ML_Amesos_Gen. */
  int ML_Amesos_Solve( void *Amesos_Handle, double x[], double rhs[] ) ;

  /** Destroy the Amesos object. Prints out some timing. */
  void ML_Amesos_Destroy(void *Amesos_Handle);

#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif
