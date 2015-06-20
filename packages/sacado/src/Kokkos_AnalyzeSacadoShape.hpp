// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef KOKKOS_ANALYZE_SACADO_SHAPE_HPP
#define KOKKOS_ANALYZE_SACADO_SHAPE_HPP

#include <impl/Kokkos_Shape.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------

/** \brief  Analyze the array shape defined by a Kokkos::View data type.
 *
 *  It is presumed that the data type can be mapped down to a multidimensional
 *  array of an intrinsic scalar numerical type (double, float, int, ... ).
 *  The 'value_type' of an array may be an embedded aggregate type such
 *  as a fixed length array 'Array<T,N>'.  In this case the 'array_instrinsic_type'
 *  represents the underlying array of intrinsic numeric type.
 *
 *  The embedded aggregate type must have an AnalyzeSacadoShape specialization
 *  to map it down to a shape and intrinsic scalar numerical type.
 *
 *  This is a slight variation of the original AnalyzeShape in that takes the
 *  the Layout as a template parameter.  This allows the sacado specializations
 *  to put the sacado dimension in different places depending on the layout.
 */

template< class T, class Layout >
struct AnalyzeSacadoShape : public Shape< sizeof(T) , 0 >
{
  typedef void specialize ;

  typedef Shape< sizeof(T), 0 >  shape ;

  typedef       T  array_intrinsic_type ;
  typedef       T  flat_array_type ;
  typedef       T  value_type ;
  typedef       T  type ;
  typedef const T  const_array_intrinsic_type ;
  typedef const T  const_flat_array_type ;
  typedef const T  const_value_type ;
  typedef const T  const_type ;
  typedef       T  non_const_array_intrinsic_type ;
  typedef       T  non_const_flat_array_type ;
  typedef       T  non_const_value_type ;
  typedef       T  non_const_type ;
};

template< class Layout>
struct AnalyzeSacadoShape<void, Layout> : public Shape< 0 , 0 >
{
  typedef void specialize ;

  typedef Shape< 0 , 0 >  shape ;

  typedef       void  array_intrinsic_type ;
  typedef       void  flat_array_type ;
  typedef       void  value_type ;
  typedef       void  type ;
  typedef const void  const_array_intrinsic_type ;
  typedef const void  const_flat_array_type ;
  typedef const void  const_value_type ;
  typedef const void  const_type ;
  typedef       void  non_const_array_intrinsic_type ;
  typedef       void  non_const_flat_array_type ;
  typedef       void  non_const_value_type ;
  typedef       void  non_const_type ;
};

template< class T, class Layout >
struct AnalyzeSacadoShape< const T, Layout > :
    public AnalyzeSacadoShape<T, Layout>::shape
{
private:
  typedef AnalyzeSacadoShape<T, Layout> nested ;
public:

  typedef typename nested::specialize specialize ;

  typedef typename nested::shape shape ;

  typedef typename nested::const_array_intrinsic_type  array_intrinsic_type ;
  typedef typename nested::const_flat_array_type       flat_array_type ;
  typedef typename nested::const_value_type            value_type ;
  typedef typename nested::const_type                  type ;

  typedef typename nested::const_array_intrinsic_type   const_array_intrinsic_type ;
  typedef typename nested::const_flat_array_type        const_flat_array_type ;
  typedef typename nested::const_value_type             const_value_type ;
  typedef typename nested::const_type                   const_type ;

  typedef typename nested::non_const_array_intrinsic_type  non_const_array_intrinsic_type ;
  typedef typename nested::non_const_flat_array_type       non_const_flat_array_type ;
  typedef typename nested::non_const_value_type            non_const_value_type ;
  typedef typename nested::non_const_type                  non_const_type ;
};

template< class T, class Layout >
struct AnalyzeSacadoShape< T *, Layout >
  : public ShapeInsert< typename AnalyzeSacadoShape<T,Layout>::shape , 0 >::type
{
private:
  typedef AnalyzeSacadoShape<T, Layout> nested ;
public:

  typedef typename nested::specialize specialize ;

  typedef typename ShapeInsert< typename nested::shape , 0 >::type shape ;

  typedef typename nested::array_intrinsic_type * array_intrinsic_type ;
  typedef typename nested::flat_array_type      * flat_array_type ;
  typedef typename nested::value_type             value_type ;
  typedef typename nested::type                 * type ;

  typedef typename nested::const_array_intrinsic_type * const_array_intrinsic_type ;
  typedef typename nested::const_flat_array_type      * const_flat_array_type ;
  typedef typename nested::const_value_type             const_value_type ;
  typedef typename nested::const_type                 * const_type ;

  typedef typename nested::non_const_array_intrinsic_type * non_const_array_intrinsic_type ;
  typedef typename nested::non_const_flat_array_type      * non_const_flat_array_type ;
  typedef typename nested::non_const_value_type             non_const_value_type ;
  typedef typename nested::non_const_type                 * non_const_type ;
};

template< class T, class Layout >
struct AnalyzeSacadoShape< T[], Layout >
  : public ShapeInsert< typename AnalyzeSacadoShape<T, Layout>::shape , 0 >::type
{
private:
  typedef AnalyzeSacadoShape<T, Layout> nested ;
public:

  typedef typename nested::specialize specialize ;

  typedef typename ShapeInsert< typename nested::shape , 0 >::type shape ;

  typedef typename nested::array_intrinsic_type  array_intrinsic_type [] ;
  typedef typename nested::flat_array_type       flat_array_type [] ;
  typedef typename nested::value_type            value_type ;
  typedef typename nested::type                  type [] ;

  typedef typename nested::const_array_intrinsic_type  const_array_intrinsic_type [] ;
  typedef typename nested::const_flat_array_type       const_flat_array_type [] ;
  typedef typename nested::const_value_type            const_value_type ;
  typedef typename nested::const_type                  const_type [] ;

  typedef typename nested::non_const_array_intrinsic_type non_const_array_intrinsic_type [] ;
  typedef typename nested::non_const_flat_array_type      non_const_flat_array_type [] ;
  typedef typename nested::non_const_value_type           non_const_value_type ;
  typedef typename nested::non_const_type                 non_const_type [] ;
};

template< class T, class Layout >
struct AnalyzeSacadoShape< const T[], Layout >
  : public ShapeInsert< typename AnalyzeSacadoShape< const T, Layout >::shape , 0 >::type
{
private:
  typedef AnalyzeSacadoShape< const T, Layout > nested ;
public:

  typedef typename nested::specialize specialize ;

  typedef typename ShapeInsert< typename nested::shape , 0 >::type shape ;

  typedef typename nested::array_intrinsic_type  array_intrinsic_type [] ;
  typedef typename nested::flat_array_type       flat_array_type [] ;
  typedef typename nested::value_type            value_type ;
  typedef typename nested::type                  type [] ;

  typedef typename nested::const_array_intrinsic_type  const_array_intrinsic_type [] ;
  typedef typename nested::const_flat_array_type       const_flat_array_type [] ;
  typedef typename nested::const_value_type            const_value_type ;
  typedef typename nested::const_type                  const_type [] ;

  typedef typename nested::non_const_array_intrinsic_type non_const_array_intrinsic_type [] ;
  typedef typename nested::non_const_flat_array_type      non_const_flat_array_type [] ;
  typedef typename nested::non_const_value_type           non_const_value_type ;
  typedef typename nested::non_const_type                 non_const_type [] ;
};

template< class T, class Layout , unsigned N >
struct AnalyzeSacadoShape< T[N], Layout >
  : public ShapeInsert< typename AnalyzeSacadoShape<T, Layout>::shape , N >::type
{
private:
  typedef AnalyzeSacadoShape<T, Layout> nested ;
public:

  typedef typename nested::specialize specialize ;

  typedef typename ShapeInsert< typename nested::shape , N >::type shape ;

  typedef typename nested::array_intrinsic_type  array_intrinsic_type [N] ;
  typedef typename nested::flat_array_type       flat_array_type [N] ;
  typedef typename nested::value_type            value_type ;
  typedef typename nested::type                  type [N] ;

  typedef typename nested::const_array_intrinsic_type  const_array_intrinsic_type [N] ;
  typedef typename nested::const_flat_array_type       const_flat_array_type [N] ;
  typedef typename nested::const_value_type            const_value_type ;
  typedef typename nested::const_type                  const_type [N] ;

  typedef typename nested::non_const_array_intrinsic_type non_const_array_intrinsic_type [N] ;
  typedef typename nested::non_const_flat_array_type      non_const_flat_array_type [N] ;
  typedef typename nested::non_const_value_type           non_const_value_type ;
  typedef typename nested::non_const_type                 non_const_type [N] ;
};

template< class T, class Layout , unsigned N >
struct AnalyzeSacadoShape< const T[N], Layout >
  : public ShapeInsert< typename AnalyzeSacadoShape< const T, Layout >::shape , N >::type
{
private:
  typedef AnalyzeSacadoShape< const T, Layout > nested ;
public:

  typedef typename nested::specialize specialize ;

  typedef typename ShapeInsert< typename nested::shape , N >::type shape ;

  typedef typename nested::array_intrinsic_type  array_intrinsic_type [N] ;
  typedef typename nested::flat_array_type       flat_array_type [N] ;
  typedef typename nested::value_type            value_type ;
  typedef typename nested::type                  type [N] ;

  typedef typename nested::const_array_intrinsic_type  const_array_intrinsic_type [N] ;
  typedef typename nested::const_flat_array_type       const_flat_array_type [N] ;
  typedef typename nested::const_value_type            const_value_type ;
  typedef typename nested::const_type                  const_type [N] ;

  typedef typename nested::non_const_array_intrinsic_type non_const_array_intrinsic_type [N] ;
  typedef typename nested::non_const_flat_array_type      non_const_flat_array_type [N] ;
  typedef typename nested::non_const_value_type           non_const_value_type ;
  typedef typename nested::non_const_type                 non_const_type [N] ;
};

} // namespace Impl
} // namespace Kokkos

#endif /* #ifndef KOKKOS_ANALYZE_SACADO_SHAPE_HPP */
