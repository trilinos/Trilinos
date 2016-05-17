// @HEADER
// ***********************************************************************
//
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
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
