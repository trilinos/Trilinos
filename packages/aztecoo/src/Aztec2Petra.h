
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

#ifndef _AZTEC2PETRA_H_
#define _AZTEC2PETRA_H_

#ifndef __cplusplus
#define __cplusplus
#endif

#include "az_aztec.h"

#ifdef AZTEC_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_CrsMatrix.h"


/*! \file 
\brief Aztec2Petra:  A function that converts an Aztec linear problem to a Petra linear problem.

    Aztec2Petra takes the Aztec proc_config, Amat, az_x and az_b objects and converts them into
    corresponding Petra equivalents: comm, map (layout information), A, x, and b.  This function
    is used by AZOO_iterate, but can be used independently by someone making a transistion from 
    Aztec to Trilinos/AztecOO.
*/
/*! \fn int Aztec2Petra(int * proc_config,
          AZ_MATRIX * Amat, double * az_x,double * az_b,
          Epetra_Comm * &comm,
          Epetra_BlockMap * & map,
          Epetra_RowMatrix * &A,
          Epetra_Vector * & x,
          Epetra_Vector * & b,
	  int ** global_indices)

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
       A pointer to a Epetra_RowMatrix object containing a \bf deep copy of the matrix in Amat, if
       the user matrix is an Msr matrix.  It is a \bf shallow copy of the matrix if the user matrix
       is in Vbr format.
       Must be deleted by the caller of this function.  Note:  This pointer will actually point to a 
       Epetra_CrsMatrix or a Epetra_VbrMatrix.  We cast the pointer to a Epetra_RowMatrix
       because it is the abstract base class used by AztecOO.
\param x (Out)
       A pointer to a Epetra_Vector object containing a \bf shallow copy (view) of az_x.  
       Must be deleted by the caller of this function.
\param b (Out)
       A pointer to a Epetra_Vector object containing a \bf shallow copy (view) of az_b.  
       Must be deleted by the caller of this function.
\param global_indices (Out)
       A pointer to an internally created integer array.  If the user matrix is in Vbr format,
       this array contains a copy of the column index values in global mode.  By using this
       array, we can avoid a deep copy of the user matrix in this case. SPECIAL NOTE:  This array
       must be delete using the special Aztec function as follows:
       if (global_indices!=0) AZ_free((void *) global_indices);

*/
int Aztec2Petra(int * proc_config,
          AZ_MATRIX * Amat, double * az_x,double * az_b,
          Epetra_Comm * &comm,
          Epetra_BlockMap * & map,
          Epetra_RowMatrix * &A,
          Epetra_Vector * & x,
          Epetra_Vector * & b,
	  int ** global_indices);

#endif /* _AZTEC2PETRA_H_ */
