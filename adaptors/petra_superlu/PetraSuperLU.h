#ifndef _PETRASUPERLU_H_
#define _PETRASUPERLU_H_

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#ifdef EPETRA_MPI
#include "mpi.h"
#endif
#ifndef __cplusplus
#define __cplusplus
#endif
#include "Epetra_Comm.h"
#include "Epetra_Time.h"
#include "Epetra_Map.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Vector.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_LinearProblem.h"
#include "dsp_defs.h"
#include "util.h"


/*! \file 
\brief Aztec2Petra:  A function that converts an Aztec linear problem to a Petra linear problem.

    Aztec2Petra takes the Aztec proc_config, Amat, az_x and az_b objects and converts them into
    corresponding Petra equivalents: comm, map (layout information), A, x, and b.  This function
    is used by AZOO_iterate, but can be used independently by someone making a transistion from 
    Aztec to Trilinos/AztecOO.
*/
/*! \fn int Epetra_SuperLU( Epetra_LinearProblem * Problem)

\brief Converts from an Aztec linear problem to a Petra linear problem.

\param proc_config (In)
       Aztec array containing information about the parallel machine.
\param Amat (In)
       An Aztec AZ_MATRIX structure.  Must be an MSR or VBR matrix at this time.
\param az_x (In)
       The Aztec initial guess/solution vector.  Must be of adequate length on each processor
       for any ghost values (unnecessary in uniprocessor mode).
\param az_b (In)
       The Aztec right hand side vector.  .
\param comm (Out)
       A pointer to a Epetra_Comm object.  Must be deleted by the caller of this function.
\param map (Out)
       A pointer to a Epetra_BlockMap object.  Must be deleted by the caller of this function.  
       Note:  This object may actually be a Epetra_Map object, but Epetra_BlockMap is a base 
       clase for Epetra_Map.
\param A (Out)
       A pointer to a Epetra_RowMatrix object containing a \bf deep copy of the matrix in Amat.  
       Must be deleted by the caller of this function.  Note:  This pointer will actually point to a 
       Epetra_CrsMatrix or a Epetra_VbrMatrix.  We cast the pointer to a Epetra_RowMatrix
       because it is the abstract base class used by AztecOO.
\param x (Out)
       A pointer to a Epetra_Vector object containing a \bf shallow copy (view) of az_x.  
       Must be deleted by the caller of this function.
\param b (Out)
       A pointer to a Epetra_Vector object containing a \bf shallow copy (view) of az_b.  
       Must be deleted by the caller of this function.
*/
int PetraSuperLU( Epetra_LinearProblem * Problem);

#endif /* _PETRASUPERLU_H_ */
