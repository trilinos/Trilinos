// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// @HEADER

#ifndef TEUCHOS_ARRAY_CONVERSIONS_H
#define TEUCHOS_ARRAY_CONVERSIONS_H

/*! \file Teuchos_ArrayConversions.hpp
  \brief Templated conversions between Array<RCP<T> > and ArrayView<Ptr<T> >
*/

#include "Teuchos_Ptr.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Assert.hpp"


namespace Teuchos {


/** \brief Utility function to convert from an an input
 * Array[View,RCP]<[const] PTR<T_in> > object to an output
 * ArrayView<Ptr<T_out> > object.
 */
template<class ArrayPtrT_in, class T_out>
void arrayViewPtrConv( const ArrayPtrT_in &a_in,
  const ArrayView<Ptr<T_out> > &a_out )
{
  using Teuchos::as;
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(as<Teuchos_Ordinal>(a_in.size()),
    as<Teuchos_Ordinal>(a_out.size()));
#endif
  for (Teuchos_Ordinal i = 0; i < as<Teuchos_Ordinal>(a_in.size()); ++i) {
    a_out[i] = a_in[i].ptr();
  }
}


/** \brief Utility function to convert from an input Array[View,RCP]<[const]
 * RCP<T_in> > object to an output ArrayView<RCP<T_out> > object.
 */
template<class ArrayPtrT_in, class T_out>
void arrayViewRcpConv( const ArrayPtrT_in &a_in,
  const ArrayView<RCP<T_out> > &a_out )
{
  using Teuchos::as;
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(as<Teuchos_Ordinal>(a_in.size()),
    as<Teuchos_Ordinal>(a_out.size()));
#endif
  for (Teuchos_Ordinal i = 0; i < as<Teuchos_Ordinal>(a_in.size()); ++i) {
    a_out[i] = a_in[i];
  }
}


/** \brief Utility function to convert an Array[View,RCP]<[const] PTR<T_in> >
 * object to an Array<Ptr<T_out> > object.
 *
 * The input array type can be Array, ArrayView, ArrayRCP, or any other array
 * object that supports size() and operator[](const int).
 *
 * The reason that the while array input type is a template argument instead
 * of just passing in ArrayView templated on the pointer type is that implicit
 * conversions of the array input type will not take place unless you specify
 * the type of the input.  That is painful and verbose and we can avoid this
 * by just templating the entire array type as is done here.
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
  arrayViewPtrConv(a_in, a_out());
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
  arrayViewRcpConv(a_in, a_out());
  return a_out;
}


/** \brief Utility function that does a reinterpret case to convert an
 * ArrayView<const Ptr<T> > object to an ArrayView<const Ptr<const T> >
 * object.
 *
 * Making this conversion requires an reinterpret case since Ptr<T> is not
 * exactly the same thing as a raw pointer but this conversion should be 100
 * percent portable and safe.
 *
 *
 */
template<class T>
ArrayView<const Ptr<const T> >
arrayConstPtrConstCast(const ArrayView<const Ptr<T> > &a_in)
{
  return av_reinterpret_cast<const Ptr<const T> >(a_in);
}


/** \brief Utility function that does a reinterpret case to convert an
 * ArrayView<const RCP<T> > object to an ArrayView<const RCP<const T> >
 * object.
 *
 * Making this conversion requires an reinterpret case since RCP<T> is not
 * exactly the same thing as a raw pointer but this conversion should be 100
 * percent portable and safe.
 *
 *
 */
template<class T>
ArrayView<const RCP<const T> >
arrayConstRcpConstCast(const ArrayView<const RCP<T> > &a_in)
{
  return av_reinterpret_cast<const RCP<const T> >(a_in);
}


} // namespace Teuchos


#endif // TEUCHOS_ARRAY_CONVERSIONS_H

