/*@HEADER
// ***********************************************************************
// 
//       Tifpack: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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
*/

#ifndef _TIFPACK_CONFIGDEFS_HPP_
#define _TIFPACK_CONFIGDEFS_HPP_

#include <Tifpack_config.h>
#include <Teuchos_ScalarTraits.hpp>

//The sgn function isn't well defined for complex.
//Is it correct to operate on the real part of x as is done below?
template<class Scalar>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
TIFPACK_SGN(const Scalar& x)
{
  static const typename Teuchos::ScalarTraits<Scalar>::magnitudeType one = Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::one());
  return Teuchos::ScalarTraits<Scalar>::real(x) > 0.0 ? -one : one;
}

#include <Tpetra_ConfigDefs.hpp>
namespace Tifpack {
  typedef Tpetra::global_size_t global_size_t;
}

#endif /*_TIFPACK_CONFIGDEFS_HPP_*/
