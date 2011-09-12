//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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

#ifndef EPETRAEXT_READ_EPETRA_LINEAR_SYSTEM_H
#define EPETRAEXT_READ_EPETRA_LINEAR_SYSTEM_H

#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Teuchos_RefCountPtr.hpp"

namespace EpetraExt {

/** \brief Read in an Epetra linear system from a file.
 *
 * \param  fileName   [in] Name of the file to read in the linear system (see file formats below).
 * \param  comm       [in] The communicator
 * \param  map        [out] The createe map.  <tt>map==NULL</tt> is allowed on input in
 *                    which case this will not be returned.
 * \param A           [out] The created matrix.  <tt>A==NULL</tt> is allowed on input in
 *                    which case this will not be returned.
 * \param x           [out] The created LHS vector.  <tt>x==NULL</tt> is allowed on input in
 *                    which case this will not be returned.
 * \param b           [out] The created RHS vector.  <tt>b==NULL</tt> is allowed on input in
 *                    which case this will not be returned.
 * \param xExact      [out] The created exact LHS vector (if known).  <tt>xExact==NULL</tt> is allowed on input in
 *                    which case this will not be returned.
 *
 * This function reads from a number file formats (*.triU, *.triS, *.mtx,
 * *.hb)
 *
 * ToDo: Finish documentation!
 *
 * ToDo: Put this in EpetraExt after the release is finished.
 */
void readEpetraLinearSystem(
  const std::string                               &fileName
  ,const Epetra_Comm                              &comm
  ,Teuchos::RefCountPtr<Epetra_CrsMatrix>         *A        = NULL
  ,Teuchos::RefCountPtr<Epetra_Map>               *map      = NULL
  ,Teuchos::RefCountPtr<Epetra_Vector>            *x        = NULL
  ,Teuchos::RefCountPtr<Epetra_Vector>            *b        = NULL
  ,Teuchos::RefCountPtr<Epetra_Vector>            *xExact   = NULL
  );

} // namespace EpetraExt

#endif // EPETRAEXT_READ_EPETRA_LINEAR_SYSTEM_H
