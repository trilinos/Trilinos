
/*
//@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER
*/

#ifndef _AZTEC2PETRA_H_
#define _AZTEC2PETRA_H_

#if defined(AztecOO_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The AztecOO package is deprecated"
#endif
#endif

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
