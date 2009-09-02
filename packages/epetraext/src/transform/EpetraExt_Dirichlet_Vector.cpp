#if 0 // Dead source?
// @HEADER
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
// @HEADER

#include <EpetraExt_Dirichlet_Vector.h>

namespace EpetraExt {

bool
Vector_Dirichlet::
fwd()
{
  Epetra_Vector & Vec = *origObj_;
  int MyLength = Vec.MyLength();

  for( int i = 0 ; i < MyLength; ++i )
    if( locations_[i] ) Vec[i] = values_[i];

  return true;
}

bool
Vector_Dirichlet::
rvs()
{
  return true;
}

} //namespace EpetraExt

#endif // if 0
