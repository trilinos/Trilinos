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

#ifndef KOKKOS_HOST_VIEW_HPP
#define KOKKOS_HOST_VIEW_HPP

#include <algorithm>

#include <Kokkos_Host.hpp>
#include <Kokkos_View.hpp>

#include <Host/Kokkos_Host_Parallel.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class OutputView , class InputView  , unsigned Rank = OutputView::Rank >
struct HostViewRemap
{
  typedef Host device_type ;

  const OutputView output ;
  const InputView  input ;
  const Host::size_type n0 ;
  const Host::size_type n1 ;
  const Host::size_type n2 ;
  const Host::size_type n3 ;
  const Host::size_type n4 ;
  const Host::size_type n5 ;
  const Host::size_type n6 ;
  const Host::size_type n7 ;

  HostViewRemap( const OutputView & arg_out , const InputView & arg_in )
    : output( arg_out ), input( arg_in )
    , n0( std::min( arg_out.dimension_0() , arg_in.dimension_0() ) )
    , n1( std::min( arg_out.dimension_1() , arg_in.dimension_1() ) )
    , n2( std::min( arg_out.dimension_2() , arg_in.dimension_2() ) )
    , n3( std::min( arg_out.dimension_3() , arg_in.dimension_3() ) )
    , n4( std::min( arg_out.dimension_4() , arg_in.dimension_4() ) )
    , n5( std::min( arg_out.dimension_5() , arg_in.dimension_5() ) )
    , n6( std::min( arg_out.dimension_6() , arg_in.dimension_6() ) )
    , n7( std::min( arg_out.dimension_7() , arg_in.dimension_7() ) )
    {
      parallel_for( n0 , *this );
    }

  inline
  void operator()( const Host::size_type i0 ) const
  {
    for ( Host::size_type i1 = 0 ; i1 < n1 ; ++i1 ) {
    for ( Host::size_type i2 = 0 ; i2 < n2 ; ++i2 ) {
    for ( Host::size_type i3 = 0 ; i3 < n3 ; ++i3 ) {
    for ( Host::size_type i4 = 0 ; i4 < n4 ; ++i4 ) {
    for ( Host::size_type i5 = 0 ; i5 < n5 ; ++i5 ) {
    for ( Host::size_type i6 = 0 ; i6 < n6 ; ++i6 ) {
    for ( Host::size_type i7 = 0 ; i7 < n7 ; ++i7 ) {
      output(i0,i1,i2,i3,i4,i5,i6,i7) = input(i0,i1,i2,i3,i4,i5,i6,i7);
    }}}}}}}
  }
};

template< class OutputView , class InputView >
struct HostViewRemap< OutputView , InputView , 1 >
{
  typedef Host device_type ;

  const OutputView output ;
  const InputView  input ;

  inline
  void operator()( const Host::size_type i0 ) const
    { output(i0) = input(i0); }

  HostViewRemap( const OutputView & arg_out , const InputView & arg_in )
    : output( arg_out ), input( arg_in )
    { parallel_for( std::min( arg_out.dimension_0() , arg_in.dimension_0() ) , *this ); }
};

template< class OutputView , class InputView >
struct HostViewRemap< OutputView , InputView , 0 >
{
  HostViewRemap( const OutputView & arg_out , const InputView & arg_in )
    { *arg_out = *arg_in ; }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

/** \brief Deep copy equal dimension arrays in the host space which
 *         have different layouts or specializations.
 */

template< class DT , class DL , class DM , class DS ,
          class ST , class SL , class SM , class SS >
inline
void deep_copy( const View< DT, DL, Host, DM, DS> & dst ,
                const View< ST, SL, Host, SM, SS> & src ,
                const typename Impl::enable_if<(
                  // Destination is not constant:
                  Impl::is_same< typename ViewTraits<DT,DL,Host,DM>::value_type ,
                                 typename ViewTraits<DT,DL,Host,DM>::non_const_value_type >::value
                  &&
                  // Same rank
                  ( unsigned( ViewTraits<DT,DL,Host,DM>::rank ) ==
                    unsigned( ViewTraits<ST,SL,Host,SM>::rank ) )
                  &&
                  // Different layout or different specialization:
                  ( ( ! Impl::is_same< typename DL::array_layout ,
                                       typename SL::array_layout >::value )
                    ||
                    ( ! Impl::is_same< DS , SS >::value )
                  )
                )>::type * = 0 )
{
  typedef View< DT, DL, Host, DM, DS> dst_type ;
  typedef View< ST, SL, Host, SM, SS> src_type ;

  assert_shapes_equal_dimension( dst.shape() , src.shape() );

  Impl::HostViewRemap< dst_type , src_type >( dst , src );
}

/** \brief  Deep copy of scalar value */

template< typename ValueType , class LayoutSrc , class MemoryTraits >
inline
void deep_copy( ValueType & dst ,
                const View< ValueType , LayoutSrc , Host , MemoryTraits , Impl::LayoutScalar > & src )
{ dst = src ; }

template< typename ValueType , class LayoutDst , class MemoryTraits >
inline
void deep_copy( const View< ValueType , LayoutDst , Host , MemoryTraits , Impl::LayoutScalar > & dst ,
                const ValueType & src )
{ dst = src ; }

} // namespace Kokkos

#endif /* #ifndef KOKKOS_HOST_VIEW_HPP */


