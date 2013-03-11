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

#ifndef KOKKOSARRAY_HOST_VIEW_HPP
#define KOKKOSARRAY_HOST_VIEW_HPP

#include <KokkosArray_Host.hpp>
#include <KokkosArray_View.hpp>

#include <Host/KokkosArray_Host_Parallel.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template< class DataType , class LayoutType , class ManagedType , class Specialize >
struct ViewInitialize< View< DataType , LayoutType , Host , ManagedType , Specialize > >
{
  typedef View< DataType , LayoutType , Host , ManagedType , Specialize > view_type ;
  typedef typename view_type::scalar_type scalar_type ;

  static void apply( const view_type & view )
  {
    const size_t count = ViewAssignment< Specialize >::allocation_count( view );
    Impl::HostParallelFill< scalar_type >( view.ptr_on_device() , 0 , count );
  }
};

template< class DataType , class LayoutType , class ManagedType >
struct ViewInitialize< View< DataType , LayoutType , Host , ManagedType , Impl::LayoutScalar > >
{
  typedef View< DataType , LayoutType , Host , ManagedType , Impl::LayoutScalar > view_type ;
  typedef typename view_type::scalar_type scalar_type ;

  static void apply( const view_type & view )
  {
    Impl::HostParallelFill< scalar_type >( view.ptr_on_device() , 0 , 1 );
  }
};

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< class OutputView , class InputView  , unsigned Rank >
struct HostViewRemap ;

template< class OutputView , class InputView >
struct HostViewRemap< OutputView , InputView , 8 >
{
  const OutputView output ;
  const InputView  input ;

  HostViewRemap( const OutputView & arg_out , const InputView & arg_in )
    : output( arg_out ), input( arg_in )
    { HostParallelLaunch< HostViewRemap >( *this ); }

  void operator()( HostThread & this_thread ) const
  {
    std::pair<Host::size_type,Host::size_type> range =
      this_thread.work_range( output.dimension_0() );

    for ( Host::size_type i0 = range.first ; i0 < range.second ; ++i0 ) {
    for ( Host::size_type i1 = 0 ; i1 < output.dimension_1() ; ++i1 ) {
    for ( Host::size_type i2 = 0 ; i2 < output.dimension_2() ; ++i2 ) {
    for ( Host::size_type i3 = 0 ; i3 < output.dimension_3() ; ++i3 ) {
    for ( Host::size_type i4 = 0 ; i4 < output.dimension_4() ; ++i4 ) {
    for ( Host::size_type i5 = 0 ; i5 < output.dimension_5() ; ++i5 ) {
    for ( Host::size_type i6 = 0 ; i6 < output.dimension_6() ; ++i6 ) {
    for ( Host::size_type i7 = 0 ; i7 < output.dimension_7() ; ++i7 ) {
      output(i0,i1,i2,i3,i4,i5,i6,i7) = input(i0,i1,i2,i3,i4,i5,i6,i7);
    }}}}}}}}
  }
};

template< class OutputView , class InputView >
struct HostViewRemap< OutputView , InputView , 7 >
{
  const OutputView output ;
  const InputView  input ;

  HostViewRemap( const OutputView & arg_out , const InputView & arg_in )
    : output( arg_out ), input( arg_in )
    { HostParallelLaunch< HostViewRemap >( *this ); }

  void operator()( HostThread & this_thread ) const
  {
    std::pair<Host::size_type,Host::size_type> range =
      this_thread.work_range( output.dimension_0() );

    for ( Host::size_type i0 = range.first ; i0 < range.second ; ++i0 ) {
    for ( Host::size_type i1 = 0 ; i1 < output.dimension_1() ; ++i1 ) {
    for ( Host::size_type i2 = 0 ; i2 < output.dimension_2() ; ++i2 ) {
    for ( Host::size_type i3 = 0 ; i3 < output.dimension_3() ; ++i3 ) {
    for ( Host::size_type i4 = 0 ; i4 < output.dimension_4() ; ++i4 ) {
    for ( Host::size_type i5 = 0 ; i5 < output.dimension_5() ; ++i5 ) {
    for ( Host::size_type i6 = 0 ; i6 < output.dimension_6() ; ++i6 ) {
      output(i0,i1,i2,i3,i4,i5,i6) = input(i0,i1,i2,i3,i4,i5,i6);
    }}}}}}}
  }
};

template< class OutputView , class InputView >
struct HostViewRemap< OutputView , InputView , 6 >
{
  const OutputView output ;
  const InputView  input ;

  HostViewRemap( const OutputView & arg_out , const InputView & arg_in )
    : output( arg_out ), input( arg_in )
    { HostParallelLaunch< HostViewRemap >( *this ); }

  void operator()( HostThread & this_thread ) const
  {
    std::pair<Host::size_type,Host::size_type> range =
      this_thread.work_range( output.dimension_0() );

    for ( Host::size_type i0 = range.first ; i0 < range.second ; ++i0 ) {
    for ( Host::size_type i1 = 0 ; i1 < output.dimension_1() ; ++i1 ) {
    for ( Host::size_type i2 = 0 ; i2 < output.dimension_2() ; ++i2 ) {
    for ( Host::size_type i3 = 0 ; i3 < output.dimension_3() ; ++i3 ) {
    for ( Host::size_type i4 = 0 ; i4 < output.dimension_4() ; ++i4 ) {
    for ( Host::size_type i5 = 0 ; i5 < output.dimension_5() ; ++i5 ) {
      output(i0,i1,i2,i3,i4,i5) = input(i0,i1,i2,i3,i4,i5);
    }}}}}}
  }
};

template< class OutputView , class InputView >
struct HostViewRemap< OutputView , InputView , 5 >
{
  const OutputView output ;
  const InputView  input ;

  HostViewRemap( const OutputView & arg_out , const InputView & arg_in )
    : output( arg_out ), input( arg_in )
    { HostParallelLaunch< HostViewRemap >( *this ); }

  void operator()( HostThread & this_thread ) const
  {
    std::pair<Host::size_type,Host::size_type> range =
      this_thread.work_range( output.dimension_0() );

    for ( Host::size_type i0 = range.first ; i0 < range.second ; ++i0 ) {
    for ( Host::size_type i1 = 0 ; i1 < output.dimension_1() ; ++i1 ) {
    for ( Host::size_type i2 = 0 ; i2 < output.dimension_2() ; ++i2 ) {
    for ( Host::size_type i3 = 0 ; i3 < output.dimension_3() ; ++i3 ) {
    for ( Host::size_type i4 = 0 ; i4 < output.dimension_4() ; ++i4 ) {
      output(i0,i1,i2,i3,i4) = input(i0,i1,i2,i3,i4);
    }}}}}
  }
};

template< class OutputView , class InputView >
struct HostViewRemap< OutputView , InputView , 4 >
{
  const OutputView output ;
  const InputView  input ;

  HostViewRemap( const OutputView & arg_out , const InputView & arg_in )
    : output( arg_out ), input( arg_in )
    { HostParallelLaunch< HostViewRemap >( *this ); }

  void operator()( HostThread & this_thread ) const
  {
    std::pair<Host::size_type,Host::size_type> range =
      this_thread.work_range( output.dimension_0() );

    for ( Host::size_type i0 = range.first ; i0 < range.second ; ++i0 ) {
    for ( Host::size_type i1 = 0 ; i1 < output.dimension_1() ; ++i1 ) {
    for ( Host::size_type i2 = 0 ; i2 < output.dimension_2() ; ++i2 ) {
    for ( Host::size_type i3 = 0 ; i3 < output.dimension_3() ; ++i3 ) {
      output(i0,i1,i2,i3) = input(i0,i1,i2,i3);
    }}}}
  }
};

template< class OutputView , class InputView >
struct HostViewRemap< OutputView , InputView , 3 >
{
  const OutputView output ;
  const InputView  input ;

  HostViewRemap( const OutputView & arg_out , const InputView & arg_in )
    : output( arg_out ), input( arg_in )
    { HostParallelLaunch< HostViewRemap >( *this ); }

  void operator()( HostThread & this_thread ) const
  {
    std::pair<Host::size_type,Host::size_type> range =
      this_thread.work_range( output.dimension_0() );

    for ( Host::size_type i0 = range.first ; i0 < range.second ; ++i0 ) {
    for ( Host::size_type i1 = 0 ; i1 < output.dimension_1() ; ++i1 ) {
    for ( Host::size_type i2 = 0 ; i2 < output.dimension_2() ; ++i2 ) {
      output(i0,i1,i2) = input(i0,i1,i2);
    }}}
  }
};

template< class OutputView , class InputView >
struct HostViewRemap< OutputView , InputView , 2 >
{
  const OutputView output ;
  const InputView  input ;

  HostViewRemap( const OutputView & arg_out , const InputView & arg_in )
    : output( arg_out ), input( arg_in )
    { HostParallelLaunch< HostViewRemap >( *this ); }

  void operator()( HostThread & this_thread ) const
  {
    std::pair<Host::size_type,Host::size_type> range =
      this_thread.work_range( output.dimension_0() );

    for ( Host::size_type i0 = range.first ; i0 < range.second ; ++i0 ) {
    for ( Host::size_type i1 = 0 ; i1 < output.dimension_1() ; ++i1 ) {
      output(i0,i1) = input(i0,i1);
    }}
  }
};

template< class OutputView , class InputView >
struct HostViewRemap< OutputView , InputView , 1 >
{
  const OutputView output ;
  const InputView  input ;

  HostViewRemap( const OutputView & arg_out , const InputView & arg_in )
    : output( arg_out ), input( arg_in )
    { HostParallelLaunch< HostViewRemap >( *this ); }

  void operator()( HostThread & this_thread ) const
  {
    std::pair<Host::size_type,Host::size_type> range =
      this_thread.work_range( output.dimension_0() );

    for ( Host::size_type i0 = range.first ; i0 < range.second ; ++i0 ) {
      output(i0) = input(i0);
    }
  }
};

template< class OutputView , class InputView >
struct HostViewRemap< OutputView , InputView , 0 >
{
  HostViewRemap( const OutputView & arg_out , const InputView & arg_in )
    { *arg_out = *arg_in ; }
};

} // namespace Impl

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
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
                  ( ! ViewTraits<DT,DL,Host,DM>::is_const )
                  &&
                  // Different layout or different specialization:
                  ( ( ! Impl::is_same< typename DL::array_layout ,
                                       typename SL::array_layout >::value )
                    ||
                    ( ! Impl::is_same< DS , SS >::value )
                  )
                  &&
                  // Same rank
                  ( unsigned( ViewTraits<DT,DL,Host,DM>::rank ) ==
                    unsigned( ViewTraits<ST,SL,Host,SM>::rank ) )
                )>::type * = 0 )
{
  typedef View< DT, DL, Host, DM, DS> dst_type ;
  typedef View< ST, SL, Host, SM, SS> src_type ;

  assert_shapes_equal_dimension( dst.shape() , src.shape() );

  Impl::HostViewRemap< dst_type , src_type , dst_type::rank >( dst , src );
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

} // namespace KokkosArray

#endif /* #ifndef KOKKOSARRAY_HOST_VIEW_HPP */


