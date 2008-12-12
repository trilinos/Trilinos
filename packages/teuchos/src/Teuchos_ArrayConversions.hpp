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

#ifndef TEUCHOS_ARRAY_CONVERSIONS_H
#define TEUCHOS_ARRAY_CONVERSIONS_H

/*! \file Teuchos_ArrayConversions.hpp
  \brief Templated conversions between Array<RCP<T> > and ArrayView<Ptr<T> >
*/

#include "Teuchos_Ptr.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"


namespace Teuchos {


/** \brief Utility function to convert an Array[View,RCP]<[const] PTR<T_in> >
 * object to an Array<Ptr<T_out> > object.
 *
 * Note that the input array type can be any valid array object that supports
 * size(), operator[](const int), and who's elemnets support the member
 * function ptr() that returns an object that is implicitly convertible into a
 * Ptr<T_out> object.  Comparable input array classes include Array,
 * ArrayView, ArrayRCP, std::vector, and any other type that support the basic
 * operations.
 *
 * NOTE: This function will convert from base type to derived types as well.
 *
 * NOTE: This function will perform dynamic memory allocation in the creation
 * of the returned Array<Ptr<T_out> > object.  However, this should not be a
 * significant performance disadvantage because typically any class type that
 * you would manipulate through an RCP or Ptr object would tend to be larger
 * an more complex so that such overhead would be small and insignificant.
 */
template<class T_out, class ArrayPtrT_in>
Array<Ptr<T_out> > arrayPtrConv(const ArrayPtrT_in &a_in)
{
  using Teuchos::as;
  Array<Ptr<T_out> > a_out(a_in.size());
  for (Teuchos_Ordinal i = 0; i < as<Teuchos_Ordinal>(a_in.size()); ++i) {
    a_out[i] = a_in[i].ptr();
  }
  return a_out;
}


/** \brief Utility function to convert any Array[View,RCP]<[const] RCP<T_in> >
 * object to an Array<RCP<T_out> > object.
 *
 * Similar to the function arrayPtrConv() except this function returns
 * Array<RCP<T_out> > objects instead of Array<Ptr<T_out> > objects.
 */
template<class T_out, class ArrayPtrT_in>
Array<RCP<T_out> > arrayRcpConv(const ArrayPtrT_in &a_in)
{
  using Teuchos::as;
  Array<RCP<T_out> > a_out(a_in.size());
  for (Teuchos_Ordinal i = 0; i < as<Teuchos_Ordinal>(a_in.size()); ++i) {
    a_out[i] = a_in[i];
  }
  return a_out;
}


} // namespace Teuchos

#endif // TEUCHOS_ARRAY_CONVERSIONS_H

