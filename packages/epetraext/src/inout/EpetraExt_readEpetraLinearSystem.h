//@HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
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
