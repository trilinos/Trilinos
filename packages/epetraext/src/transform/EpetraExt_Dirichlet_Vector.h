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
                                                                                                    
#ifndef EpetraExt_VECTOR_DIRICHLET_H
#define EpetraExt_VECTOR_DIRICHLET_H

#include <EpetraExt_Transform.h>

#include <Epetra_IntVector.h>
#include <Epetra_Vector.h>

namespace EpetraExt {

///
/** Given an input Epetra_Vector, apply given dirichlet conditions
 */
class Vector_Dirichlet : public InPlaceTransform<Epetra_Vector>
{
 public:

  ///
  /** Destructor
   */
  ~Vector_Dirichlet() {}

  ///
  /** Constructor
   @param Locations Integer Vector containing 1's for Dirichlet rows and 0's if not
   @param Values Vector containing values of the Dirichlet BC's
   */
  Vector_Dirichlet( const Epetra_IntVector & Locations,
                    const Epetra_Vector & Values )
  : locations_( Locations ),
    values_( Values )
  {}

  ///
  /** Applies Dirichlet BC's
   */
  bool fwd();

  ///
  /** NoOp
   */
  bool rvs();

 private:

  const Epetra_IntVector locations_; 
  const Epetra_Vector values_;
};

} //namespace EpetraExt

#endif //EpetraExt_VECTOR_DIRICHLET_H

