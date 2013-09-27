/*
//@HEADER
// ************************************************************************
//
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
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

#ifndef KOKKOS_VIEWSUPPORT_HPP
#define KOKKOS_VIEWSUPPORT_HPP

#include <impl/Kokkos_Shape.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/** \brief  Evaluate if LHS = RHS view assignment is allowed. */
template< class ViewLHS , class ViewRHS >
struct ViewAssignable
{
  // Same memory space.
  // Same value type.
  // Compatible 'const' qualifier
  // Cannot assign managed = unmannaged
  enum { assignable_value =
    ( is_same< typename ViewLHS::value_type ,
               typename ViewRHS::value_type >::value
      ||
      is_same< typename ViewLHS::value_type ,
               typename ViewRHS::const_value_type >::value )
    &&
    is_same< typename ViewLHS::memory_space ,
             typename ViewRHS::memory_space >::value
    &&
    ( ! ( ViewLHS::is_managed && ! ViewRHS::is_managed ) )
  };

  enum { assignable_shape =
    // Compatible shape and matching layout:
    ( ShapeCompatible< typename ViewLHS::shape_type ,
                       typename ViewRHS::shape_type >::value
      &&
      is_same< typename ViewLHS::array_layout ,
               typename ViewRHS::array_layout >::value )
    ||
    // Matching layout, same rank, and LHS dynamic rank
    ( is_same< typename ViewLHS::array_layout ,
               typename ViewRHS::array_layout >::value
      &&
      int(ViewLHS::rank) == int(ViewRHS::rank)
      &&
      int(ViewLHS::rank) == int(ViewLHS::rank_dynamic) )
    ||
    // Both rank-0, any shape and layout
    ( int(ViewLHS::rank) == 0 && int(ViewRHS::rank) == 0 )
    ||
    // Both rank-1 and LHS is dynamic rank-1, any shape and layout
    ( int(ViewLHS::rank) == 1 && int(ViewRHS::rank) == 1 &&
      int(ViewLHS::rank_dynamic) == 1 )
    };

  enum { value = assignable_value && assignable_shape };
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class ShapeType , class LayoutType , class Enable = void >
class LayoutStride ;

/* Arrays with rank <= 1 have no stride */
template< class ShapeType , class LayoutType >
class LayoutStride< ShapeType , LayoutType ,
                    typename enable_if< ShapeType::rank <= 1 >::type >
{
public:

  enum { dynamic = false };
  enum { value = 0 };

  KOKKOS_INLINE_FUNCTION static
  void assign( LayoutStride & , const unsigned ) {}

  KOKKOS_INLINE_FUNCTION static
  void assign_no_padding( LayoutStride & , const ShapeType & ) {}

  KOKKOS_INLINE_FUNCTION static
  void assign_with_padding( LayoutStride & , const ShapeType & ) {}
};

/* Array with LayoutLeft and 0 == rank_dynamic have static stride that are is not padded. */
template< class ShapeType >
class LayoutStride< ShapeType , LayoutLeft ,
                    typename enable_if<(
                      ( 1 <  ShapeType::rank ) &&
                      ( 0 == ShapeType::rank_dynamic )
                    )>::type >
{
public:

  enum { dynamic = false };
  enum { value   = ShapeType::N0 };

  KOKKOS_INLINE_FUNCTION static
  void assign( LayoutStride & , const unsigned ) {}

  KOKKOS_INLINE_FUNCTION static
  void assign_no_padding( LayoutStride & , const ShapeType & ) {}

  KOKKOS_INLINE_FUNCTION static
  void assign_with_padding( LayoutStride & , const ShapeType & ) {}
};

/* Array with LayoutRight and 1 >= rank_dynamic have static stride that is not padded */
template< class ShapeType >
class LayoutStride< ShapeType , LayoutRight ,
                    typename enable_if<(
                      ( 1 <  ShapeType::rank ) &&
                      ( 1 >= ShapeType::rank_dynamic )
                    )>::type >
{
public:

  enum { dynamic = false };
  enum { value   = ShapeType::N1 * ShapeType::N2 * ShapeType::N3 *
                   ShapeType::N4 * ShapeType::N5 * ShapeType::N6 * ShapeType::N7 };

  KOKKOS_INLINE_FUNCTION static
  void assign( LayoutStride & , const unsigned ) {}

  KOKKOS_INLINE_FUNCTION static
  void assign_no_padding( LayoutStride & , const ShapeType & ) {}

  KOKKOS_INLINE_FUNCTION static
  void assign_with_padding( LayoutStride & , const ShapeType & ) {}
};


/* Otherwise array has runtime stride that is padded. */
template< class ShapeType , class LayoutType , class Enable >
class LayoutStride
{
public:

  enum { dynamic = true };

  unsigned value ;

  KOKKOS_INLINE_FUNCTION static
  void assign( LayoutStride & stride , const unsigned n ) { stride.value = n ; }

  KOKKOS_INLINE_FUNCTION static
  void assign_no_padding( LayoutStride & vs , const ShapeType & sh )
    {
      enum { left = is_same< LayoutType , LayoutLeft >::value };

      // Left  layout arrays are aligned on the first dimension.
      // Right layout arrays are aligned on blocks of the 2-8th dimensions.
      vs.value = ShapeType::rank <= 1 ? 0 : (
                 left ? sh.N0
                      : sh.N1 * sh.N2 * sh.N3 * sh.N4 * sh.N5 * sh.N6 * sh.N7 );
    }

  KOKKOS_INLINE_FUNCTION static
  void assign_with_padding( LayoutStride & vs , const ShapeType & sh )
    {
      enum { div   = MEMORY_ALIGNMENT / ShapeType::scalar_size };
      enum { mod   = MEMORY_ALIGNMENT % ShapeType::scalar_size };
      enum { align = 0 == mod ? div : 0 };

      assign_no_padding( vs , sh );

      if ( align && MEMORY_ALIGNMENT_THRESHOLD * align < vs.value ) {

        const unsigned count_mod = vs.value % ( div ? div : 1 );

        if ( count_mod ) { vs.value += align - count_mod ; }
      }
    }
};

template< class ShapeType , class LayoutType >
KOKKOS_INLINE_FUNCTION
size_t capacity( const ShapeType & shape ,
                 const LayoutStride< ShapeType , LayoutType > & stride )
{
  enum { left = is_same< LayoutType , LayoutLeft >::value };

  return ShapeType::rank <= 1 ? size_t(shape.N0) : (
         left ? size_t( stride.value * shape.N1 * shape.N2 * shape.N3 * shape.N4 * shape.N5 * shape.N6 * shape.N7 )
              : size_t( stride.value * shape.N0 ));
}

template< typename iType , class ShapeType , class LayoutType >
KOKKOS_INLINE_FUNCTION
void stride( iType * const s , const ShapeType & shape ,
                               const LayoutStride< ShapeType , LayoutType > & stride )
{
  enum { rank = ShapeType::rank };
  enum { left = is_same< LayoutType , LayoutLeft >::value };

  if ( 0 < rank ) {
    if ( 1 == rank ) {
      s[0] = 1 ;
    }
    else if ( left ) {
      s[0] = 1 ;
      s[1] = stride.value ;
      if ( 2 < rank ) { s[2] = s[1] * shape.N1 ; }
      if ( 3 < rank ) { s[3] = s[2] * shape.N2 ; }
      if ( 4 < rank ) { s[4] = s[3] * shape.N3 ; }
      if ( 5 < rank ) { s[5] = s[4] * shape.N4 ; }
      if ( 6 < rank ) { s[6] = s[5] * shape.N5 ; }
      if ( 7 < rank ) { s[7] = s[6] * shape.N6 ; }
    }
    else {
      s[rank-1] = 1 ;
      if ( 7 < rank ) { s[6] = s[7] * shape.N7 ; }
      if ( 6 < rank ) { s[5] = s[6] * shape.N6 ; }
      if ( 5 < rank ) { s[4] = s[5] * shape.N5 ; }
      if ( 4 < rank ) { s[3] = s[4] * shape.N4 ; }
      if ( 3 < rank ) { s[2] = s[3] * shape.N3 ; }
      if ( 2 < rank ) { s[1] = s[2] * shape.N2 ; }
      s[0] = stride.value ;
    }
  }
}

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/** \brief  View tracking increment/decrement only happens when
 *          view memory is managed and executing in the host space.
 */
template< class ViewTraits , class Enable = void >
struct ViewTracking {
  KOKKOS_INLINE_FUNCTION static void increment( const void * ) {}
  KOKKOS_INLINE_FUNCTION static void decrement( const void * ) {}
};

template< class ViewTraits >
struct ViewTracking< ViewTraits ,
                     typename enable_if<(
                       ViewTraits::is_managed &&
                       Impl::is_same< HostSpace , ExecutionSpace >::value
                     )>::type >
{
  typedef typename ViewTraits::memory_space memory_space ;

  KOKKOS_INLINE_FUNCTION static void increment( const void * ptr )
    { memory_space::increment( ptr ); }

  KOKKOS_INLINE_FUNCTION static void decrement( const void * ptr )
    { memory_space::decrement( ptr ); }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class DstMemorySpace , class SrcMemorySpace >
struct DeepCopy ;

template< class OutputView , unsigned Rank = OutputView::Rank >
struct ViewInit
{
  typedef typename OutputView::device_type device_type ;
  typedef typename OutputView::scalar_type scalar_type ;
  typedef typename device_type::size_type  size_type ;

  const OutputView output ;

  explicit ViewInit( const OutputView & arg_out ) : output( arg_out )
    { parallel_for( output.dimension_0() , *this ); }

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type i0 ) const
  {
    const scalar_type default_value = scalar_type();

    for ( size_type i1 = 0 ; i1 < output.dimension_1() ; ++i1 ) {
    for ( size_type i2 = 0 ; i2 < output.dimension_2() ; ++i2 ) {
    for ( size_type i3 = 0 ; i3 < output.dimension_3() ; ++i3 ) {
    for ( size_type i4 = 0 ; i4 < output.dimension_4() ; ++i4 ) {
    for ( size_type i5 = 0 ; i5 < output.dimension_5() ; ++i5 ) {
    for ( size_type i6 = 0 ; i6 < output.dimension_6() ; ++i6 ) {
    for ( size_type i7 = 0 ; i7 < output.dimension_7() ; ++i7 ) {
      new (&output.at(i0,i1,i2,i3,i4,i5,i6,i7)) scalar_type(default_value) ;
    }}}}}}}
  }
};

template< class OutputView >
struct ViewInit< OutputView , 1 >
{
  typedef typename OutputView::device_type device_type ;
  typedef typename OutputView::value_type  value_type ;
  typedef typename device_type::size_type  size_type ;

  const OutputView output ;

  explicit ViewInit( const OutputView & arg_out ) : output( arg_out )
    { parallel_for( output.dimension_0() , *this ); }

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type i0 ) const
  {
    value_type default_value = value_type();
    new (&output(i0)) value_type(default_value) ;
  }
};

template< class OutputView >
struct ViewInit< OutputView , 0 >
{
  typedef typename OutputView::device_type device_type ;
  typedef typename OutputView::value_type  value_type ;
  typedef typename device_type::size_type  size_type ;

  const OutputView output ;

  explicit ViewInit( const OutputView & arg_out ) : output( arg_out )
    { parallel_for( 1 , *this ); }

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type /*i0*/ ) const
  {
    value_type default_value = value_type();
    new (&(*output)) value_type(default_value) ;
  }
};

template< class Device >
struct ViewInitialize
{
  template< class ViewType >
  inline explicit ViewInitialize( const ViewType & view )
    { ViewInit<ViewType> init( view ); }
};

template< class OutputView , class InputView  , unsigned Rank = OutputView::Rank >
struct ViewRemap
{
  typedef typename OutputView::device_type device_type ;
  typedef typename device_type::size_type  size_type ;

  const OutputView output ;
  const InputView  input ;
  const size_type n0 ;
  const size_type n1 ;
  const size_type n2 ;
  const size_type n3 ;
  const size_type n4 ;
  const size_type n5 ;
  const size_type n6 ;
  const size_type n7 ;

  ViewRemap( const OutputView & arg_out , const InputView & arg_in )
    : output( arg_out ), input( arg_in )
    , n0( std::min( (size_t)arg_out.dimension_0() , (size_t)arg_in.dimension_0() ) )
    , n1( std::min( (size_t)arg_out.dimension_1() , (size_t)arg_in.dimension_1() ) )
    , n2( std::min( (size_t)arg_out.dimension_2() , (size_t)arg_in.dimension_2() ) )
    , n3( std::min( (size_t)arg_out.dimension_3() , (size_t)arg_in.dimension_3() ) )
    , n4( std::min( (size_t)arg_out.dimension_4() , (size_t)arg_in.dimension_4() ) )
    , n5( std::min( (size_t)arg_out.dimension_5() , (size_t)arg_in.dimension_5() ) )
    , n6( std::min( (size_t)arg_out.dimension_6() , (size_t)arg_in.dimension_6() ) )
    , n7( std::min( (size_t)arg_out.dimension_7() , (size_t)arg_in.dimension_7() ) )
    {
      parallel_for( n0 , *this );
    }

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type i0 ) const
  {
    for ( size_type i1 = 0 ; i1 < n1 ; ++i1 ) {
    for ( size_type i2 = 0 ; i2 < n2 ; ++i2 ) {
    for ( size_type i3 = 0 ; i3 < n3 ; ++i3 ) {
    for ( size_type i4 = 0 ; i4 < n4 ; ++i4 ) {
    for ( size_type i5 = 0 ; i5 < n5 ; ++i5 ) {
    for ( size_type i6 = 0 ; i6 < n6 ; ++i6 ) {
    for ( size_type i7 = 0 ; i7 < n7 ; ++i7 ) {
      output.at(i0,i1,i2,i3,i4,i5,i6,i7) = input.at(i0,i1,i2,i3,i4,i5,i6,i7);
    }}}}}}}
  }
};

template< class OutputView , class InputView  >
struct ViewRemap< OutputView ,  InputView , 0 >
{
  typedef typename OutputView::value_type   value_type ;
  typedef typename OutputView::memory_space dst_space ;
  typedef typename InputView ::memory_space src_space ;

  ViewRemap( const OutputView & arg_out , const InputView & arg_in )
  {
    DeepCopy< dst_space , src_space >( arg_out.ptr_on_device() ,
                                       arg_in.ptr_on_device() ,
                                       sizeof(value_type) );
  }
};

template< class OutputView , unsigned Rank = OutputView::Rank >
struct ViewFill
{
  typedef typename OutputView::device_type       device_type ;
  typedef typename OutputView::const_value_type  const_value_type ;
  typedef typename device_type::size_type        size_type ;

  const OutputView output ;
  const_value_type input ;

  ViewFill( const OutputView & arg_out , const_value_type & arg_in )
    : output( arg_out ), input( arg_in )
    {
      parallel_for( output.dimension_0() , *this );
    }

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type i0 ) const
  {
    for ( size_type i1 = 0 ; i1 < output.dimension_1() ; ++i1 ) {
    for ( size_type i2 = 0 ; i2 < output.dimension_2() ; ++i2 ) {
    for ( size_type i3 = 0 ; i3 < output.dimension_3() ; ++i3 ) {
    for ( size_type i4 = 0 ; i4 < output.dimension_4() ; ++i4 ) {
    for ( size_type i5 = 0 ; i5 < output.dimension_5() ; ++i5 ) {
    for ( size_type i6 = 0 ; i6 < output.dimension_6() ; ++i6 ) {
    for ( size_type i7 = 0 ; i7 < output.dimension_7() ; ++i7 ) {
      output.at(i0,i1,i2,i3,i4,i5,i6,i7) = input ;
    }}}}}}}
  }
};

template< class OutputView >
struct ViewFill< OutputView , 0 >
{
  typedef typename OutputView::device_type       device_type ;
  typedef typename OutputView::const_value_type  const_value_type ;
  typedef typename OutputView::memory_space      dst_space ;

  ViewFill( const OutputView & arg_out , const_value_type & arg_in )
  {
    DeepCopy< dst_space , dst_space >( arg_out.ptr_on_device() , & arg_in ,
                                       sizeof(const_value_type) );
  }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_VIEWSUPPORT_HPP */


