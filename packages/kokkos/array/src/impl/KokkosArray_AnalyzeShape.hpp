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

#ifndef KOKKOSARRAY_ANALYZESHAPE_HPP
#define KOKKOSARRAY_ANALYZESHAPE_HPP

#include <impl/KokkosArray_Shape.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

//----------------------------------------------------------------------------

template< class ShapeType , unsigned N ,
          unsigned R = ShapeType::rank_dynamic >
struct ShapeInsert ;

template< class ShapeType , unsigned N >
struct ShapeInsert< ShapeType , N , 0 >
{
  typedef Shape< typename ShapeType::array_layout ,
                 ShapeType::scalar_size ,
                 ShapeType::rank + 1 ,
                 N ,
                 ShapeType::StaticN0 ,
                 ShapeType::StaticN1 ,
                 ShapeType::StaticN2 ,
                 ShapeType::StaticN3 ,
                 ShapeType::StaticN4 ,
                 ShapeType::StaticN5 ,
                 ShapeType::StaticN6 > type ;
};

template< class ShapeType , unsigned N >
struct ShapeInsert< ShapeType , N , 1 >
{
  typedef Shape< typename ShapeType::array_layout ,
                 ShapeType::scalar_size ,
                 ShapeType::rank + 1 ,
                 ShapeType::StaticN0 ,
                 N ,
                 ShapeType::StaticN1 ,
                 ShapeType::StaticN2 ,
                 ShapeType::StaticN3 ,
                 ShapeType::StaticN4 ,
                 ShapeType::StaticN5 ,
                 ShapeType::StaticN6 > type ;
};

template< class ShapeType , unsigned N >
struct ShapeInsert< ShapeType , N , 2 >
{
  typedef Shape< typename ShapeType::array_layout ,
                 ShapeType::scalar_size ,
                 ShapeType::rank + 1 ,
                 ShapeType::StaticN0 ,
                 ShapeType::StaticN1 ,
                 N ,
                 ShapeType::StaticN2 ,
                 ShapeType::StaticN3 ,
                 ShapeType::StaticN4 ,
                 ShapeType::StaticN5 ,
                 ShapeType::StaticN6 > type ;
};

template< class ShapeType , unsigned N >
struct ShapeInsert< ShapeType , N , 3 >
{
  typedef Shape< typename ShapeType::array_layout ,
                 ShapeType::scalar_size ,
                 ShapeType::rank + 1 ,
                 ShapeType::StaticN0 ,
                 ShapeType::StaticN1 ,
                 ShapeType::StaticN2 ,
                 N ,
                 ShapeType::StaticN3 ,
                 ShapeType::StaticN4 ,
                 ShapeType::StaticN5 ,
                 ShapeType::StaticN6 > type ;
};

template< class ShapeType , unsigned N >
struct ShapeInsert< ShapeType , N , 4 >
{
  typedef Shape< typename ShapeType::array_layout ,
                 ShapeType::scalar_size ,
                 ShapeType::rank + 1 ,
                 ShapeType::StaticN0 ,
                 ShapeType::StaticN1 ,
                 ShapeType::StaticN2 ,
                 ShapeType::StaticN3 ,
                 N ,
                 ShapeType::StaticN4 ,
                 ShapeType::StaticN5 ,
                 ShapeType::StaticN6 > type ;
};

template< class ShapeType , unsigned N >
struct ShapeInsert< ShapeType , N , 5 >
{
  typedef Shape< typename ShapeType::array_layout ,
                 ShapeType::scalar_size ,
                 ShapeType::rank + 1 ,
                 ShapeType::StaticN0 ,
                 ShapeType::StaticN1 ,
                 ShapeType::StaticN2 ,
                 ShapeType::StaticN3 ,
                 ShapeType::StaticN4 ,
                 N ,
                 ShapeType::StaticN5 ,
                 ShapeType::StaticN6 > type ;
};

template< class ShapeType , unsigned N >
struct ShapeInsert< ShapeType , N , 6 >
{
  typedef Shape< typename ShapeType::array_layout ,
                 ShapeType::scalar_size ,
                 ShapeType::rank + 1 ,
                 ShapeType::StaticN0 ,
                 ShapeType::StaticN1 ,
                 ShapeType::StaticN2 ,
                 ShapeType::StaticN3 ,
                 ShapeType::StaticN4 ,
                 ShapeType::StaticN5 ,
                 N ,
                 ShapeType::StaticN6 > type ;
};

template< class ShapeType , unsigned N >
struct ShapeInsert< ShapeType , N , 7 >
{
  typedef Shape< typename ShapeType::array_layout ,
                 ShapeType::scalar_size ,
                 ShapeType::rank + 1 ,
                 ShapeType::StaticN0 ,
                 ShapeType::StaticN1 ,
                 ShapeType::StaticN2 ,
                 ShapeType::StaticN3 ,
                 ShapeType::StaticN4 ,
                 ShapeType::StaticN5 ,
                 ShapeType::StaticN6 ,
                 N > type ;
};

//----------------------------------------------------------------------------

/** \brief  Analyze the array shape defined by a KokkosArray::View data type.
 *
 *  It is presumed that the data type can be mapped down to a multidimensional
 *  array of an intrinsic scalar numerical type (double, float, int, ... ).
 *  The 'value_type' of an array may be an embedded aggregate type such
 *  as a fixed length array 'Array<T,N>'.  In this case the 'scalar_type'
 *  is 'T' and the 'value_type' is 'Array<T,N>' to enable data layout 
 *  according to shape and scalar_type AND data access by value_type.
 *
 *  The embedded aggregate type must have an AnalyzeShape specialization
 *  to map it down to a shape and intrinsic scalar numerical type.
 */

template< class T , class Layout = void >
struct AnalyzeShape {
  typedef       T  scalar_type ;
  typedef       T  array_type ;
  typedef       T  value_type ;
  typedef       T  type ;
  typedef const T  const_scalar_type ;
  typedef const T  const_array_type ;
  typedef const T  const_value_type ;
  typedef const T  const_type ;
  typedef       T  non_const_scalar_type ;
  typedef       T  non_const_array_type ;
  typedef       T  non_const_value_type ;
  typedef       T  non_const_type ;

  typedef Shape< Layout, sizeof(T), 0 >  shape ;

  static const unsigned rank = shape::rank ;
};

template< class T , class Layout >
struct AnalyzeShape< const T , Layout > {
private:
  typedef AnalyzeShape<T,Layout> nested ;
  typedef typename nested::shape nested_shape ;
public:

  typedef typename nested::const_scalar_type scalar_type ;
  typedef typename nested::const_array_type  array_type ;
  typedef typename nested::const_value_type  value_type ;
  typedef typename nested::const_type        type ;

  typedef typename nested::const_scalar_type const_scalar_type ;
  typedef typename nested::const_array_type  const_array_type ;
  typedef typename nested::const_value_type  const_value_type ;
  typedef typename nested::const_type        const_type ;

  typedef typename nested::non_const_scalar_type non_const_scalar_type ;
  typedef typename nested::non_const_array_type  non_const_array_type ;
  typedef typename nested::non_const_value_type  non_const_value_type ;
  typedef typename nested::non_const_type        non_const_type ;

  typedef nested_shape shape ;

  static const unsigned rank = shape::rank ;
};

template< class T , class Layout >
struct AnalyzeShape< T * , Layout >
{
private:
  typedef AnalyzeShape<T, Layout > nested ;
  typedef typename nested::shape nested_shape ;
public:

  typedef typename nested::scalar_type  scalar_type ;
  typedef typename nested::array_type * array_type ;
  typedef typename nested::value_type   value_type ;
  typedef typename nested::type       * type ;

  typedef typename nested::const_scalar_type  const_scalar_type ;
  typedef typename nested::const_array_type * const_array_type ;
  typedef typename nested::const_value_type   const_value_type ;
  typedef typename nested::const_type       * const_type ;

  typedef typename nested::non_const_scalar_type  non_const_scalar_type ;
  typedef typename nested::non_const_array_type * non_const_array_type ;
  typedef typename nested::non_const_value_type   non_const_value_type ;
  typedef typename nested::non_const_type       * non_const_type ;

  typedef typename ShapeInsert< nested_shape , 0 >::type shape ;

  static const unsigned rank = shape::rank ;
};

template< class T , class Layout >
struct AnalyzeShape< T[] , Layout >
{
private:
  typedef AnalyzeShape<T,Layout> nested ;
  typedef typename nested::shape nested_shape ;
public:

  typedef typename nested::scalar_type scalar_type ;
  typedef typename nested::array_type  array_type [] ;
  typedef typename nested::value_type  value_type ;
  typedef typename nested::type        type [] ;

  typedef typename nested::const_scalar_type const_scalar_type ;
  typedef typename nested::const_array_type  const_array_type [] ;
  typedef typename nested::const_value_type  const_value_type ;
  typedef typename nested::const_type        const_type [] ;

  typedef typename nested::non_const_scalar_type non_const_scalar_type ;
  typedef typename nested::non_const_array_type  non_const_array_type [] ;
  typedef typename nested::non_const_value_type  non_const_value_type ;
  typedef typename nested::non_const_type        non_const_type [] ;

  typedef typename ShapeInsert< nested_shape , 0 >::type shape ;

  static const unsigned rank = shape::rank ;
};

template< class T , class Layout >
struct AnalyzeShape< const T[] , Layout >
{
private:
  typedef AnalyzeShape< const T,Layout> nested ;
  typedef typename nested::shape nested_shape ;
public:

  typedef typename nested::scalar_type scalar_type ;
  typedef typename nested::array_type  array_type [] ;
  typedef typename nested::value_type  value_type ;
  typedef typename nested::type        type [] ;

  typedef typename nested::const_scalar_type const_scalar_type ;
  typedef typename nested::const_array_type  const_array_type [] ;
  typedef typename nested::const_value_type  const_value_type ;
  typedef typename nested::const_type        const_type [] ;

  typedef typename nested::non_const_scalar_type non_const_scalar_type ;
  typedef typename nested::non_const_array_type  non_const_array_type [] ;
  typedef typename nested::non_const_value_type  non_const_value_type ;
  typedef typename nested::non_const_type        non_const_type [] ;

  typedef typename ShapeInsert< nested_shape , 0 >::type shape ;

  static const unsigned rank = shape::rank ;
};

template< class T , unsigned N , class Layout >
struct AnalyzeShape< T[N] , Layout >
{
private:
  typedef AnalyzeShape<T, Layout > nested ;
  typedef typename nested::shape nested_shape ;
public:

  typedef typename nested::scalar_type scalar_type ;
  typedef typename nested::array_type  array_type [N] ;
  typedef typename nested::value_type  value_type ;
  typedef typename nested::type        type [N] ;

  typedef typename nested::const_scalar_type const_scalar_type ;
  typedef typename nested::const_array_type  const_array_type [N] ;
  typedef typename nested::const_value_type  const_value_type ;
  typedef typename nested::const_type        const_type [N] ;

  typedef typename nested::non_const_scalar_type non_const_scalar_type ;
  typedef typename nested::non_const_array_type  non_const_array_type [N] ;
  typedef typename nested::non_const_value_type  non_const_value_type ;
  typedef typename nested::non_const_type        non_const_type [N] ;

  typedef typename ShapeInsert< nested_shape , N >::type shape ;

  static const unsigned rank = shape::rank ;
};

template< class T , unsigned N , class Layout >
struct AnalyzeShape< const T[N] , Layout >
{
private:
  typedef AnalyzeShape< const T, Layout > nested ;
  typedef typename nested::shape nested_shape ;
public:

  typedef typename nested::scalar_type scalar_type ;
  typedef typename nested::array_type  array_type [N] ;
  typedef typename nested::value_type  value_type ;
  typedef typename nested::type        type [N] ;

  typedef typename nested::const_scalar_type const_scalar_type ;
  typedef typename nested::const_array_type  const_array_type [N] ;
  typedef typename nested::const_value_type  const_value_type ;
  typedef typename nested::const_type        const_type [N] ;

  typedef typename nested::non_const_scalar_type non_const_scalar_type ;
  typedef typename nested::non_const_array_type  non_const_array_type [N] ;
  typedef typename nested::non_const_value_type  non_const_value_type ;
  typedef typename nested::non_const_type        non_const_type [N] ;

  typedef typename ShapeInsert< nested_shape , N >::type shape ;

  static const unsigned rank = shape::rank ;
};

} // namespace Impl
} // namespace KokkosArray

#endif /* #ifndef KOKKOSARRAY_ANALYZESHAPE_HPP */

