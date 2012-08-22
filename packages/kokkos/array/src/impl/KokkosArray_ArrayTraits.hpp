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

#ifndef KOKKOSARRAYTRAITS_HPP
#define KOKKOSARRAYTRAITS_HPP

namespace KokkosArray {
namespace Impl {

/* TR1 conformal compile-time array traits utilities.
 * Prefer to use TR1 when portably available.
 */
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
/** \brief  Rank of an array: rank<T>::value.
 *
 *  Equating, or not, the types 'T[]' and 'T[0]' is compiler dependent.
 *  However, the type 'T[]' is specified by the C++ standard to be used
 *  only for the first array specification.
 */

template< class T , unsigned J = 0 >
struct rank : public unsigned_< J > {};

template< class T >
struct rank< T[] , 0 > : public unsigned_< rank<T,1>::value > {};

template< class T , unsigned J >
struct rank< T[0] , J > : public unsigned_< rank<T,J+1>::value > {};

template< class T , unsigned N , unsigned J >
struct rank< T[N] , J > : public unsigned_< rank<T,J+1>::value > {};

//----------------------------------------------------------------------------
/** \brief  Add 'const' to the value type of an array.
 *
 *  The precedence of the 'const' and array type specifications
 *  is compiler dependent.  As such the template recurses to
 *  the value type for the array before adding the const qualifier.
 *
 *  An empty extent '[]' is replaced with a zero extent '[0]'
 *  as for Kokkos Arrays this is the general case.
 */

template< class T , unsigned J = 0 >
struct add_const { typedef const T type ; };

template< class T , unsigned J >
struct add_const< const T , J > { typedef const T type ; };


template< class T >
struct add_const< T[] , 0 >
{ typedef typename add_const<T,1>::type type [0] ; };

template< class T >
struct add_const< const T[] , 0 >
{ typedef typename add_const<T,1>::type type [0] ; };


template< class T , unsigned J >
struct add_const< T[0] , J >
{ typedef typename add_const<T,J+1>::type type [0] ; };

template< class T , unsigned J >
struct add_const< const T[0] , J >
{ typedef typename add_const<T,J+1>::type type [0] ; };


template< class T , unsigned N , unsigned J >
struct add_const< T[N] , J >
{ typedef typename add_const<T,J+1>::type type [N] ; };

template< class T , unsigned N , unsigned J >
struct add_const< const T[N] , J >
{ typedef typename add_const<T,J+1>::type type [N] ; };

//----------------------------------------------------------------------------
/** \brief  Remove 'const' from the value type of an array.
 *
 *  An empty extent '[]' is replaced with a zero extent '[0]'
 *  as for Kokkos Arrays this is the general case.
 */

template< class T , unsigned J = 0 >
struct remove_const { typedef T type ; };

template< class T , unsigned J >
struct remove_const< const T , J > { typedef T type ; };


template< class T >
struct remove_const< T[] , 0 >
{ typedef typename remove_const<T,1>::type type [0] ; };

template< class T >
struct remove_const< const T[] , 0 >
{ typedef typename remove_const<T,1>::type type [0] ; };


template< class T , unsigned J >
struct remove_const< T[0] , J >
{ typedef typename remove_const<T,J+1>::type type [0] ; };

template< class T , unsigned J >
struct remove_const< const T[0] , J >
{ typedef typename remove_const<T,J+1>::type type [0] ; };


template< class T , unsigned N , unsigned J >
struct remove_const< T[N] , J >
{ typedef typename remove_const<T,J+1>::type type [N] ; };

template< class T , unsigned N , unsigned J >
struct remove_const< const T[N] , J >
{ typedef typename remove_const<T,J+1>::type type [N] ; };

//----------------------------------------------------------------------------
/** \brief  Remove all dimensions of the array */

template< class T , unsigned J = 0 >
struct remove_all_extents { typedef T type ; };

template< class T >
struct remove_all_extents< T[] , 0 >
{ typedef typename remove_all_extents<T,1>::type type ; };

template< class T , unsigned J >
struct remove_all_extents< T[0] , J >
{ typedef typename remove_all_extents<T,J+1>::type type ; };

template< class T , unsigned N , unsigned J >
struct remove_all_extents< T[N] , J >
{ typedef typename remove_all_extents<T,J+1>::type type ; };

//----------------------------------------------------------------------------

/** \brief  Extent (a.k.a. dimension) of an array: extent<T,I>::value.
 *          If rank I is not a dimension of the array then value == 0.
 *          If rank I is the leading unspecified dimension [] then value == 0.
 */

template< class T , unsigned I , unsigned J = 0 >
struct extent : public unsigned_< 0 > {};


template< class T >
struct extent< T[] , 0 , 0 > : public unsigned_< 0 > {};

template< class T , unsigned J >
struct extent< T[0] , 0 , J > : public unsigned_< 0 > {};

template< class T , unsigned N , unsigned J >
struct extent< T[N] , 0 , J > : public unsigned_< N > {};


template< class T , unsigned I >
struct extent< T[] , I , 0 >
  : public unsigned_< extent<T,I-1,1>::value > {};

template< class T , unsigned I , unsigned J >
struct extent< T[0] , I , J >
  : public unsigned_< extent<T,I-1,J+1>::value > {};

template< class T , unsigned N , unsigned I , unsigned J >
struct extent< T[N] , I , J >
  : public unsigned_< extent<T,I-1,J+1>::value > {};

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOSARRAYTRAITS_HPP */

