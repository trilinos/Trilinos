// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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

// ////////////////////////////////////////////////////////////////////////
// Teuchos_PrimitiveTypeTraits.hpp

#ifndef TEUCHOS_PRIMITIVE_TYPE_TRAITS_H
#define TEUCHOS_PRIMITIVE_TYPE_TRAITS_H

#include "Teuchos_TestForException.hpp"

/** \file Teuchos_PrimitiveTypeTraits.hpp
 */

namespace Teuchos {

///
/** Declaration of a traits class for an decomposing object into an
 * array of primitive objects.
 *
 * The idea behind this traits class it that it allows an object of
 * semi-complex structure to be externalized into an array of
 * primitive data types.
 *
 * This default traits class works just fine for types that are
 * already primitive.
 */
template <class T> class PrimitiveTypeTraits {
public:
  ///
  typedef T  primitiveType;
  ///
  static int numPrimitiveObjs() { return 1; }
  ///
  static void extractPrimitiveObjs(
    const T                &obj
    ,const int             numPrimitiveObjs
    ,primitiveType         primitiveObjs[]
    )
    {
#ifdef _DEBUG
      TEST_FOR_EXCEPTION( numPrimitiveObjs!=1 || primitiveObjs==NULL, std::invalid_argument, "Error!" );
#endif
      primitiveObjs[0] = obj;
    }
  ///
  static void loadPrimitiveObjs(
    const int              numPrimitiveObjs
    ,const primitiveType   primitiveObjs[]
    ,T                     *obj
    )
    {
#ifdef _DEBUG
      TEST_FOR_EXCEPTION( numPrimitiveObjs!=1 || primitiveObjs==NULL, std::invalid_argument, "Error!" );
#endif
      *obj = primitiveObjs[0];
    }
};

#ifdef HAVE_COMPLEX

///
/** Partial specialization of <tt>PrimitiveTypeTraits</tt> for <tt>std::complex</tt>.
 */
template <class T> class PrimitiveTypeTraits< std::complex<T> > {
public:
  ///
  typedef T  primitiveType;
  ///
  static int numPrimitiveObjs() { return 2; }
  ///
  static void extractPrimitiveObjs(
    const std::complex<T>  &obj
    ,const int             numPrimitiveObjs
    ,primitiveType         primitiveObjs[]
    )
    {
#ifdef _DEBUG
      TEST_FOR_EXCEPTION( numPrimitiveObjs!=2 || primitiveObjs==NULL, std::invalid_argument, "Error!" );
#endif
      primitiveObjs[0] = obj.real();
      primitiveObjs[1] = obj.imag();
    }
  ///
  static void loadPrimitiveObjs(
    const int              numPrimitiveObjs
    ,const primitiveType   primitiveObjs[]
    ,std::complex<T>       *obj
    )
    {
#ifdef _DEBUG
      TEST_FOR_EXCEPTION( numPrimitiveObjs!=2 || primitiveObjs==NULL, std::invalid_argument, "Error!" );
#endif
      *obj = std::complex<T>( primitiveObjs[0], primitiveObjs[1] );
    }
};

#endif // HAVE_COMPLEX

} // namespace Teuchos

#endif // TEUCHOS_PRIMITIVE_TYPE_TRAITS_H
