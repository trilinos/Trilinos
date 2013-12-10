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

#ifndef KOKKOS_EXAMPLE_LINALG_BLAS_HPP
#define KOKKOS_EXAMPLE_LINALG_BLAS_HPP

#include <cmath>
#include <utility>
#include <Kokkos_View.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Example {

template< class ViewTypeX , class ViewTypeY , class TagType = void >
struct Dot { enum { exists = false }; };

template< class ViewTypeX , class ViewTypeY >
inline
typename Kokkos::Example::Dot< ViewTypeX , ViewTypeY >::value_type
dot( const ViewTypeX & x ,
     const ViewTypeY & y )
          
{
  Kokkos::Example::Dot< ViewTypeX , ViewTypeY > functor( x , y );

  return functor.result ;
}

template< class ViewTypeX , class ViewTypeY , class TagType >
typename Kokkos::Example::Dot< ViewTypeX , ViewTypeY , TagType >::value_type
dot( const ViewTypeX & x ,
     const ViewTypeY & y ,
     TagType  tag )
          
{
  Kokkos::Example::Dot< ViewTypeX , ViewTypeY , TagType > functor( x , y );

  return functor.result ;
}

} // namespace Example
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Example {

// w(i) = alpha * x(i) + beta * y(i)
template< typename ScalarTypeA ,
          class    ViewTypeX ,
          typename ScalarTypeB ,
          class    ViewTypeY ,
          class    ViewTypeW ,
          class    TagType = void >
struct WAXPBY ;

template< typename ScalarTypeA , class ViewTypeX ,
          typename ScalarTypeB , class ViewTypeY ,
                                 class ViewTypeW >
void waxpby( const ScalarTypeA & alpha ,
             const ViewTypeX   & x ,
             const ScalarTypeB & beta ,
             const ViewTypeY   & y ,
             const ViewTypeW   & w ,
             typename WAXPBY< ScalarTypeA , ViewTypeX ,
                              ScalarTypeB , ViewTypeY ,
                                            ViewTypeW >::tag_type * = 0 )
{
  WAXPBY< ScalarTypeA , ViewTypeX , ScalarTypeB , ViewTypeY , ViewTypeW >( alpha , x , beta , y , w );
}

template< typename ScalarTypeA , class ViewTypeX ,
          typename ScalarTypeB , class ViewTypeY ,
                                 class ViewTypeW ,
          class TagType >
void waxpby( const ScalarTypeA & alpha ,
             const ViewTypeX   & x ,
             const ScalarTypeB & beta ,
             const ViewTypeY   & y ,
             const ViewTypeW   & w ,
             const TagType     & tag )
{
  WAXPBY< ScalarTypeA , ViewTypeX , ScalarTypeB , ViewTypeY , ViewTypeW , TagType >( alpha , x , beta , y , w );
}

template< typename ScalarTypeA , class ViewTypeX ,
                                 class ViewTypeY >
void axpy( const ScalarTypeA & alpha ,
           const ViewTypeX   & x ,
           const ViewTypeY   & y ,
           typename WAXPBY< ScalarTypeA , ViewTypeX , void , ViewTypeY , void >::tag_type * = 0 )
{
  WAXPBY< ScalarTypeA , ViewTypeX , void , ViewTypeY , void >( alpha , x , y );
}

template< class ViewTypeX , typename ScalarTypeB , class ViewTypeY >
void xpby( const ViewTypeX   & x ,
           const ScalarTypeB & beta ,
           const ViewTypeY   & y ,
           typename WAXPBY< void , ViewTypeX , ScalarTypeB , ViewTypeY , void >::tag_type * = 0 )
{
  WAXPBY< void , ViewTypeX , ScalarTypeB , ViewTypeY , void >( x , beta , y );
}

} // namespace Example
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Example {

template< class ViewTypeX , class ViewTypeY >
struct Dot<
  ViewTypeX ,
  ViewTypeY ,
  typename Kokkos::Impl::enable_if<(
    // Must be a view type:
    Kokkos::is_view< ViewTypeX >::value &&
    Kokkos::is_view< ViewTypeY >::value &&
    // Must be a view to scalar type:
    Kokkos::Impl::is_same< typename ViewTypeX::value_type ,
                           typename ViewTypeX::scalar_type >::value &&
    Kokkos::Impl::is_same< typename ViewTypeY::value_type ,
                           typename ViewTypeY::scalar_type >::value &&
    // Must have same device:
    Kokkos::Impl::is_same< typename ViewTypeX::device_type ,
                           typename ViewTypeY::device_type >::value &&
    // Must be rank one:
    ViewTypeX::rank == 1 &&
    ViewTypeY::rank == 1
 )>::type >
{
private:

  const View< typename ViewTypeX::const_data_type ,
              typename ViewTypeX::array_layout ,
              typename ViewTypeX::device_type ,
              MemoryUnmanaged >  x ;

  const View< typename ViewTypeY::const_data_type ,
              typename ViewTypeY::array_layout ,
              typename ViewTypeY::device_type ,
              MemoryUnmanaged >  y ;

public:

  typedef void tag_type ;
  typedef typename ViewTypeX::device_type           device_type ;
  typedef typename ViewTypeX::non_const_value_type  value_type ;

  value_type result ;

  inline
  Dot( const ViewTypeX & arg_x ,
       const ViewTypeY & arg_y )
    : x( arg_x ), y( arg_y ), result(0)
  {
    parallel_reduce( std::min( x.dimension_0() , y.dimension_0() ) , *this , result );
  }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  void operator()( const iType & i , value_type & update ) const
    { update += x(i) * y(i); }

  KOKKOS_INLINE_FUNCTION
  static void join( volatile       value_type & update ,
                    volatile const value_type & source )
    { update += source ; }

  KOKKOS_INLINE_FUNCTION
  static void init( value_type & update ) { update = 0 ; }
}; // Dot

#if 0

template< class ViewTypeX , class ViewTypeY , class TagType , class Tail >
struct Dot<
  ViewTypeX ,
  ViewTypeY ,
  TagList< TagType , Tail > >
{
  enum { exists = false };

  // Advance the tag list to an entry that exists:
  // Dot<?,?,void>::tag_list = TagList< void , void >
  typedef typename
    Kokkos::Impl::if_c<
      Dot< ViewTypeX , ViewTypeY , TagType >::exists ,
      TagList< TagType , Tail > ,
      Dot< ViewTypeX , ViewTypeY , Tail >::tag_list >::type 
    tag_list ;

  inline
  Dot( const ViewTypeX & arg_x ,
       const ViewTypeY & arg_y ,
       const TagList< TagType , Tail > & arg_tag )
  {
    typedef typename tag_list::head head ;
    typedef typename tag_list::tail tail ;
    const tag_list & tag = arg_tag.advance( head );
    if ( tag.enabled() ) {
      Dot< ViewTypeX , ViewTypeY , head > s( arg_x , arg_y );
    }
    else {
      Dot< ViewTypeX , ViewTypeY , tail > s( arg_x , arg_y , tag.next() );
    }
  }
}; // Dot

template< class ViewTypeX , class ViewTypeY ,
          class Tag0 , class Tag1 , class Tag2 , class Tag3 >
struct Dot<
  ViewTypeX ,
  ViewTypeY ,
  TypeVector< Tag0 , Tag1 , Tag2 , Tag3 > >
{
  inline
  Dot( const ViewTypeX & arg_x ,
       const ViewTypeY & arg_y ,
       const 
  {
  }
}; // Dot

#endif

} // namespace Example
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Example {

template< typename ScalarTypeA ,
          class    ViewTypeX ,
          typename ScalarTypeB ,
          class    ViewTypeY ,
          class    ViewTypeW >
struct WAXPBY< ScalarTypeA ,
               ViewTypeX ,
               ScalarTypeB ,
               ViewTypeY ,
               ViewTypeW ,
  typename Kokkos::Impl::enable_if<(
    Kokkos::is_view< ViewTypeX >::value &&
    Kokkos::is_view< ViewTypeY >::value &&
    Kokkos::is_view< ViewTypeW >::value &&
    // Same device:
    Kokkos::Impl::is_same< typename ViewTypeX::device_type ,
                           typename ViewTypeW::device_type >::value &&
    Kokkos::Impl::is_same< typename ViewTypeY::device_type ,
                           typename ViewTypeW::device_type >::value &&
    // Non-const output
    ! Kokkos::Impl::is_const< typename ViewTypeW::value_type >::value &&
    // Scalar type:
    Kokkos::Impl::is_same< typename ViewTypeX::value_type ,
                           typename ViewTypeX::scalar_type >::value &&
    Kokkos::Impl::is_same< typename ViewTypeY::value_type ,
                           typename ViewTypeY::scalar_type >::value &&
    Kokkos::Impl::is_same< typename ViewTypeW::value_type ,
                           typename ViewTypeW::scalar_type >::value &&
    // Rank one:
    ViewTypeX::rank == 1 &&
    ViewTypeY::rank == 1 &&
    ViewTypeW::rank == 1
  )>::type >
{
private:

  const Kokkos::View< typename ViewTypeW::data_type ,
                      typename ViewTypeW::array_layout ,
                      typename ViewTypeW::device_type ,
                      MemoryUnmanaged > w ;

  const Kokkos::View< typename ViewTypeX::data_type ,
                      typename ViewTypeX::array_layout ,
                      typename ViewTypeX::device_type ,
                      MemoryUnmanaged > x ;

  const Kokkos::View< typename ViewTypeY::data_type ,
                      typename ViewTypeY::array_layout ,
                      typename ViewTypeY::device_type ,
                      MemoryUnmanaged > y ;

  const ScalarTypeA  alpha ;
  const ScalarTypeB  beta ;

public:

  typedef void tag_type ;
  typedef typename ViewTypeW::device_type  device_type ;

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  void operator()( const iType & i ) const
  {
    w(i) = alpha * x(i) + beta * y(i);
  }

  inline
  WAXPBY( const ScalarTypeA & arg_alpha ,
          const ViewTypeX   & arg_x ,
          const ScalarTypeB & arg_beta ,
          const ViewTypeY   & arg_y ,
          const ViewTypeW   & arg_w )
    : w( arg_w ), x( arg_x ), y( arg_y )
    , alpha( arg_alpha ), beta( arg_beta )
  {
    parallel_for( std::min( w.dimension_0() ,
                  std::min( y.dimension_0() ,
                            x.dimension_0() ) ) , *this );
  }
}; // WAXPBY

// y(i) += alpha * x(i)
template< typename ScalarTypeA ,
          class    ViewTypeX ,
          class    ViewTypeY >
struct WAXPBY< ScalarTypeA ,
               ViewTypeX ,
               void ,
               ViewTypeY ,
               void ,
  typename Kokkos::Impl::enable_if<(
    Kokkos::is_view< ViewTypeX >::value &&
    Kokkos::is_view< ViewTypeY >::value &&
    // Same device:
    Kokkos::Impl::is_same< typename ViewTypeX::device_type ,
                           typename ViewTypeY::device_type >::value &&
    // Non-const output
    ! Kokkos::Impl::is_const< typename ViewTypeY::value_type >::value &&
    // Scalar type:
    Kokkos::Impl::is_same< typename ViewTypeX::value_type ,
                           typename ViewTypeX::scalar_type >::value &&
    Kokkos::Impl::is_same< typename ViewTypeY::value_type ,
                           typename ViewTypeY::scalar_type >::value &&
    // Rank one:
    ViewTypeX::rank == 1 &&
    ViewTypeY::rank == 1
  )>::type >
{
private:

  const Kokkos::View< typename ViewTypeX::data_type ,
                      typename ViewTypeX::array_layout ,
                      typename ViewTypeX::device_type ,
                      MemoryUnmanaged > x ;

  const Kokkos::View< typename ViewTypeY::data_type ,
                      typename ViewTypeY::array_layout ,
                      typename ViewTypeY::device_type ,
                      MemoryUnmanaged > y ;

  const ScalarTypeA  alpha ;

public:

  typedef void tag_type ;
  typedef typename ViewTypeY::device_type  device_type ;

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  void operator()( const iType & i ) const
  {
    y(i) += alpha * x(i);
  }

  inline
  WAXPBY( const ScalarTypeA & arg_alpha ,
          const ViewTypeX   & arg_x ,
          const ViewTypeY   & arg_y )
    : x( arg_x ), y( arg_y ), alpha( arg_alpha )
  {
    parallel_for( std::min( x.dimension_0() , y.dimension_0() ) , *this );
  }
}; // WAXPBY

// y(i) = x(i) + beta * y(i);
template< class    ViewTypeX ,
          typename ScalarTypeB ,
          class    ViewTypeY >
struct WAXPBY< void ,
               ViewTypeX ,
               ScalarTypeB ,
               ViewTypeY ,
               void ,
  typename Kokkos::Impl::enable_if<(
    Kokkos::is_view< ViewTypeX >::value &&
    Kokkos::is_view< ViewTypeY >::value &&
    // Same device:
    Kokkos::Impl::is_same< typename ViewTypeX::device_type ,
                           typename ViewTypeY::device_type >::value &&
    // Non-const output
    ! Kokkos::Impl::is_const< typename ViewTypeY::value_type >::value &&
    // Scalar type:
    Kokkos::Impl::is_same< typename ViewTypeX::value_type ,
                           typename ViewTypeX::scalar_type >::value &&
    Kokkos::Impl::is_same< typename ViewTypeY::value_type ,
                           typename ViewTypeY::scalar_type >::value &&
    // Rank one:
    ViewTypeX::rank == 1 &&
    ViewTypeY::rank == 1
  )>::type >
{
private:

  const Kokkos::View< typename ViewTypeX::data_type ,
                      typename ViewTypeX::array_layout ,
                      typename ViewTypeX::device_type ,
                      MemoryUnmanaged > x ;

  const Kokkos::View< typename ViewTypeY::data_type ,
                      typename ViewTypeY::array_layout ,
                      typename ViewTypeY::device_type ,
                      MemoryUnmanaged > y ;

  const ScalarTypeB  beta ;

public:

  typedef void tag_type ;
  typedef typename ViewTypeY::device_type  device_type ;

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  void operator()( const iType & i ) const
    { y(i) = x(i) + beta * y(i); }

  inline
  WAXPBY( const ViewTypeX   & arg_x ,
          const ScalarTypeB & arg_beta ,
          const ViewTypeY   & arg_y )
    : x( arg_x ), y( arg_y ), beta( arg_beta )
  {
    parallel_for( std::min( x.dimension_0() , y.dimension_0() ) , *this );
  }
}; // WAXPBY

} // namespace Example
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_EXAMPLE_LINALG_BLAS_HPP */


