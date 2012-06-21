/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_ARRAYTRAITS_HPP
#define KOKKOS_ARRAYTRAITS_HPP

namespace KokkosArray {
namespace Impl {

//----------------------------------------------------------------------------
/* TR1 conformal compile-time array traits utilities.
 * Prefer to use TR1 when portably available.
 */

/** \brief  Remove 'const' from type */
template< class T > struct remove_const ;

/** \brief  If is the same type */
template< class X , class Y > struct is_same ;

//----------------------------------------------------------------------------

template< class T > struct remove_const< const T > { typedef T type ; };

template< class T > struct remove_const { typedef T type ; };

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

template< class X , class Y >
struct is_same : public false_type {};

template< class X >
struct is_same<X,X> : public true_type {};

//----------------------------------------------------------------------------
/** \brief  Whether an empty extent is different from a zero extent
 *          is compiler-dependent.
 */
typedef bool_< ! is_same<int[],int[0]>::value >::type
  extent_empty_not_same_extent_zero ;

/** \brief  Rank of an array: rank<T>::value */
template< class T , class = true_type > struct rank ;

template< class T > struct rank<T,true_type> : public unsigned_<0> {};

template< class T , unsigned N >
struct rank< T[N] , true_type >
  : public unsigned_< rank<T,true_type>::value + 1 > {};

template< class T >
struct rank< T[0] , true_type >
  : public unsigned_< rank<T,true_type>::value + 1 > {};

template< class T >
struct rank< T[] , extent_empty_not_same_extent_zero >
  : public unsigned_<1> {};

//----------------------------------------------------------------------------

/** \brief  Remove the first (left most) dimension of the type. */
template< class T , class = true_type > struct remove_extent ;

template< class T >
struct remove_extent<T,true_type> { typedef T type ; };

template< class T >
struct remove_extent<T[0],true_type> { typedef T type ; };

template< class T , unsigned N >
struct remove_extent<T[N],true_type> { typedef T type ; };

template< class T >
struct remove_extent<T[], extent_empty_not_same_extent_zero >
  { typedef T type ; };

//----------------------------------------------------------------------------

template< class T , class = true_type > struct remove_all_extents ;

/** \brief  Remove all dimensions of the array */
template< class T >
struct remove_all_extents<T,true_type> { typedef T type ; };

template< class T , unsigned N >
struct remove_all_extents<T[N],true_type>
  { typedef typename remove_all_extents<T>::type type ; };

template< class T >
struct remove_all_extents<T[0],true_type>
  { typedef typename remove_all_extents<T>::type type ; };

template< class T >
struct remove_all_extents<T[], extent_empty_not_same_extent_zero >
  { typedef typename remove_all_extents<T>::type type ; };

//----------------------------------------------------------------------------

/** \brief  Extent (a.k.a. dimension) of an array: extent<T,I>::value.
 *          If rank I is not a dimension of the array then value == 0.
 *          If rank I is the leading unspecified dimension [] then value == 0.
 */
template< class T , unsigned I = 0 , class = true_type > struct extent ;

template< class T , unsigned I >
struct extent<T,I,true_type> : public unsigned_< 0> {};

/* To prevent roll-over */
template< class T >
struct extent<T, ~0u, true_type> : public unsigned_< 0> {};

template< class T , unsigned N >
struct extent<T[N],0u,true_type> : public unsigned_< N> {};

template< class T , unsigned N , unsigned I >
struct extent<T[N],I,true_type>
  : public unsigned_< extent<T,I-1>::value > {};

template< class T >
struct extent<T[0],0u,true_type> : public unsigned_< 0> {};

template< class T , unsigned I >
struct extent<T[0],I,true_type>
  : public unsigned_< extent<T,I-1>::value > {};

template< class T >
struct extent<T[], 0u, extent_empty_not_same_extent_zero >
  : public unsigned_< 0> {};

template< class T , unsigned I >
struct extent<T[], I, extent_empty_not_same_extent_zero >
  : public unsigned_< extent<T,I-1>::value > {};

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_ARRAYTRAITS_HPP */

