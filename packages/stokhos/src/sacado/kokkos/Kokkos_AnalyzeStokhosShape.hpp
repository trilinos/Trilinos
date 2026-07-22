// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#if false
#ifndef KOKKOS_ANALYZE_STOKHOS_SHAPE_HPP
#define KOKKOS_ANALYZE_STOKHOS_SHAPE_HPP

//#include <impl/Kokkos_Shape.hpp>

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
 *  The embedded aggregate type must have an AnalyzeStokhosShape specialization
 *  to map it down to a shape and intrinsic scalar numerical type.
 *
 *  This is a slight variation of the original AnalyzeShape in that takes the
 *  the Layout as a template parameter.  This allows the sacado specializations
 *  to put the sacado dimension in different places depending on the layout.
 */

template< class T, class Layout >
struct AnalyzeStokhosShape : public Shape< sizeof(T) , 0 >
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
struct AnalyzeStokhosShape<void, Layout> : public Shape< 0 , 0 >
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
struct AnalyzeStokhosShape< const T, Layout > :
    public AnalyzeStokhosShape<T, Layout>::shape
{
private:
  typedef AnalyzeStokhosShape<T, Layout> nested ;
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
struct AnalyzeStokhosShape< T *, Layout >
  : public ShapeInsert< typename AnalyzeStokhosShape<T,Layout>::shape , 0 >::type
{
private:
  typedef AnalyzeStokhosShape<T, Layout> nested ;
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
struct AnalyzeStokhosShape< T[], Layout >
  : public ShapeInsert< typename AnalyzeStokhosShape<T, Layout>::shape , 0 >::type
{
private:
  typedef AnalyzeStokhosShape<T, Layout> nested ;
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
struct AnalyzeStokhosShape< const T[], Layout >
  : public ShapeInsert< typename AnalyzeStokhosShape< const T, Layout >::shape , 0 >::type
{
private:
  typedef AnalyzeStokhosShape< const T, Layout > nested ;
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
struct AnalyzeStokhosShape< T[N], Layout >
  : public ShapeInsert< typename AnalyzeStokhosShape<T, Layout>::shape , N >::type
{
private:
  typedef AnalyzeStokhosShape<T, Layout> nested ;
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
struct AnalyzeStokhosShape< const T[N], Layout >
  : public ShapeInsert< typename AnalyzeStokhosShape< const T, Layout >::shape , N >::type
{
private:
  typedef AnalyzeStokhosShape< const T, Layout > nested ;
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
#endif
#endif /* #ifndef KOKKOS_ANALYZE_STOKHOS_SHAPE_HPP */
