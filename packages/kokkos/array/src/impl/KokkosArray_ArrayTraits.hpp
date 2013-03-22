/*
//@HEADER
// ************************************************************************
//
//   KokkosArray: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOSARRAYTRAITS_HPP
#define KOKKOSARRAYTRAITS_HPP

namespace KokkosArray {
namespace Impl {

/* TR1 conformal compile-time array traits utilities.
 * Prefer to use TR1 when portably available.
 */
//----------------------------------------------------------------------------

template< bool , class = void >
struct enable_if ;

template< class T >
struct enable_if< true , T > { typedef T type ; };

template< class , class T = void >
struct enable_if_type { typedef T type ; };

//----------------------------------------------------------------------------

template <class T, T v>
struct integral_constant
{
  static const T value = v;
  typedef T value_type;
  typedef integral_constant<T,v> type;
};

typedef integral_constant<bool,false> false_type ;
typedef integral_constant<bool,true>  true_type ;

template< bool B >
struct bool_ : public integral_constant<bool,B> {};

template< unsigned I >
struct unsigned_ : public integral_constant<unsigned,I> {};

template< int I >
struct int_ : public integral_constant<int,I> {};

//----------------------------------------------------------------------------
/** \brief  If is the same type */
template< class X , class Y > struct is_same : public false_type {};
template< class X >           struct is_same<X,X> : public true_type {};

//----------------------------------------------------------------------------

template <size_t N>
struct is_power_of_two
{
  enum type { value = (N > 0) && !(N & (N-1)) };
};

template < size_t N , bool OK = is_power_of_two<N>::value >
struct power_of_two ;

template < size_t N >
struct power_of_two<N,true>
{
  enum type { value = 1+ power_of_two<(N>>1),true>::value };
};

template <>
struct power_of_two<2,true>
{
  enum type { value = 1 };
};

template <>
struct power_of_two<1,true>
{
  enum type { value = 0 };
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOSARRAYTRAITS_HPP */

