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

#ifndef KOKKOS_HOST_VIEW_HPP
#define KOKKOS_HOST_VIEW_HPP

#include <KokkosArray_View.hpp>
#include <Host/KokkosArray_Host_Parallel.hpp>

#include <KokkosArray_Host_macros.hpp>
#include <impl/KokkosArray_ViewOperLeft_macros.hpp>
#include <impl/KokkosArray_ViewOperRight_macros.hpp>
#include <impl/KokkosArray_View_macros.hpp>
#include <KokkosArray_Clear_macros.hpp>

#if 1

namespace KokkosArray {

template< class DataType , class LayoutType >
void View< DataType , LayoutType , Host >::create(
  const std::string & label ,
  const View< DataType , LayoutType , Host >::shape_type shape )
{
  const size_t count = Impl::allocation_count( shape );

  oper_type::m_shape = shape ;
  oper_type::m_ptr_on_device = (value_type *)
    memory_space::allocate( label ,
                            typeid(value_type) ,
                            sizeof(value_type) ,
                            count );

  Impl::HostParallelFill<value_type>( oper_type::m_ptr_on_device , 0 , count );
}

}

#endif


namespace KokkosArray {
namespace Impl {

//----------------------------------------------------------------------------

/** \brief  Deep copy a single value */
template< typename DataType , class LayoutType >
struct Factory< View< DataType , LayoutType , Host > , DataType >
{
  typedef View< DataType , LayoutType , Host >  output_type ;
  typedef DataType                              input_type ;
  typedef typename output_type::shape_type      shape_type ;
  typedef typename output_type::value_type      value_type ;

  typedef typename
    assert_shape_is_rank_zero< shape_type >::type ok_rank ;

  typedef typename
    StaticAssertAssignable< value_type , DataType >::type ok_assign ;

  static inline
  void deep_copy( const output_type & output , const input_type & input )
    { *output = input ; }
};

/** \brief  Deep copy a single value */
template< typename DataType , class LayoutType >
struct Factory< DataType , View< DataType , LayoutType , Host > >
{
  typedef DataType                              output_type ;
  typedef View< DataType , LayoutType , Host >  input_type ;
  typedef typename input_type::shape_type       shape_type ;
  typedef typename output_type::value_type      value_type ;

  typedef typename
    assert_shape_is_rank_zero< shape_type >::type ok_rank ;

  typedef typename
    StaticAssertAssignable< DataType , value_type >::type ok_assign ;

  static inline
  void deep_copy( output_type & output , const input_type & input )
    { output = *input ; }
};

//----------------------------------------------------------------------------
/** \brief  Identical arrays */
template< class DataType , class LayoutType >
struct Factory< View< DataType , LayoutType , Host > ,
                View< DataType , LayoutType , Host > >
{
public:
  typedef View< DataType , LayoutType , Host > output_type ;
  typedef View< DataType , LayoutType , Host > input_type ;

  static inline
  void deep_copy( const output_type & output ,
                  const input_type  & input )
  {
    typedef typename output_type::value_type value_type ;

    if ( output != input ) {

      assert_shapes_are_equal( output.m_shape , input.m_shape );

      const size_t count = allocation_count( output.m_shape );

      HostParallelCopy<value_type,value_type>( output.ptr_on_device() ,
                                               input. ptr_on_device() ,
                                               count );
    }
  }

  static inline
  void deep_copy( const output_type & output ,
                  const input_type  & input ,
                  const size_t count )
  {
    typedef typename output_type::value_type value_type ;

    // Only for rank-1 arrays, or arrays where higher ranks are one

    assert_shape_effective_rank1_at_leastN( output.m_shape , count );
    assert_shape_effective_rank1_at_leastN( input.m_shape , count );

    HostParallelCopy<value_type,value_type>( output.ptr_on_device() ,
                                             input. ptr_on_device() ,
                                             count );
  }

  // Called by create_mirror
  static inline
  output_type create( const input_type & input )
  {
    return output_type("mirror",input.m_shape);
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< class OutputView , unsigned OutputRank ,
          class InputView  , unsigned InputRank >
struct HostViewRemap ;

template< class OutputView , class InputView >
struct HostViewRemap< OutputView , 8 , InputView , 8 >
  : public HostThreadWorker<void>
{
  const OutputView output ;
  const InputView  input ;

  HostViewRemap( const OutputView & arg_out , const InputView & arg_in )
    : output( arg_out ), input( arg_in )
    { HostThreadWorker<void>::execute( *this ); }

  void execute_on_thread( HostThread & this_thread ) const
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

    this_thread.barrier();
  }
};

template< class OutputView , class InputView >
struct HostViewRemap< OutputView , 7 , InputView , 7 >
  : public HostThreadWorker<void>
{
  const OutputView output ;
  const InputView  input ;

  HostViewRemap( const OutputView & arg_out , const InputView & arg_in )
    : output( arg_out ), input( arg_in )
    { HostThreadWorker<void>::execute( *this ); }

  void execute_on_thread( HostThread & this_thread ) const
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

    this_thread.barrier();
  }
};

template< class OutputView , class InputView >
struct HostViewRemap< OutputView , 6 , InputView , 6 >
  : public HostThreadWorker<void>
{
  const OutputView output ;
  const InputView  input ;

  HostViewRemap( const OutputView & arg_out , const InputView & arg_in )
    : output( arg_out ), input( arg_in )
    { HostThreadWorker<void>::execute( *this ); }

  void execute_on_thread( HostThread & this_thread ) const
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

    this_thread.barrier();
  }
};

template< class OutputView , class InputView >
struct HostViewRemap< OutputView , 5 , InputView , 5 >
  : public HostThreadWorker<void>
{
  const OutputView output ;
  const InputView  input ;

  HostViewRemap( const OutputView & arg_out , const InputView & arg_in )
    : output( arg_out ), input( arg_in )
    { HostThreadWorker<void>::execute( *this ); }

  void execute_on_thread( HostThread & this_thread ) const
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

    this_thread.barrier();
  }
};

template< class OutputView , class InputView >
struct HostViewRemap< OutputView , 4 , InputView , 4 >
  : public HostThreadWorker<void>
{
  const OutputView output ;
  const InputView  input ;

  HostViewRemap( const OutputView & arg_out , const InputView & arg_in )
    : output( arg_out ), input( arg_in )
    { HostThreadWorker<void>::execute( *this ); }

  void execute_on_thread( HostThread & this_thread ) const
  {
    std::pair<Host::size_type,Host::size_type> range =
      this_thread.work_range( output.dimension_0() );

    for ( Host::size_type i0 = range.first ; i0 < range.second ; ++i0 ) {
    for ( Host::size_type i1 = 0 ; i1 < output.dimension_1() ; ++i1 ) {
    for ( Host::size_type i2 = 0 ; i2 < output.dimension_2() ; ++i2 ) {
    for ( Host::size_type i3 = 0 ; i3 < output.dimension_3() ; ++i3 ) {
      output(i0,i1,i2,i3) = input(i0,i1,i2,i3);
    }}}}

    this_thread.barrier();
  }
};

template< class OutputView , class InputView >
struct HostViewRemap< OutputView , 3 , InputView , 3 >
  : public HostThreadWorker<void>
{
  const OutputView output ;
  const InputView  input ;

  HostViewRemap( const OutputView & arg_out , const InputView & arg_in )
    : output( arg_out ), input( arg_in )
    { HostThreadWorker<void>::execute( *this ); }

  void execute_on_thread( HostThread & this_thread ) const
  {
    std::pair<Host::size_type,Host::size_type> range =
      this_thread.work_range( output.dimension_0() );

    for ( Host::size_type i0 = range.first ; i0 < range.second ; ++i0 ) {
    for ( Host::size_type i1 = 0 ; i1 < output.dimension_1() ; ++i1 ) {
    for ( Host::size_type i2 = 0 ; i2 < output.dimension_2() ; ++i2 ) {
      output(i0,i1,i2) = input(i0,i1,i2);
    }}}

    this_thread.barrier();
  }
};

template< class OutputView , class InputView >
struct HostViewRemap< OutputView , 2 , InputView , 2 >
  : public HostThreadWorker<void>
{
  const OutputView output ;
  const InputView  input ;

  HostViewRemap( const OutputView & arg_out , const InputView & arg_in )
    : output( arg_out ), input( arg_in )
    { HostThreadWorker<void>::execute( *this ); }

  void execute_on_thread( HostThread & this_thread ) const
  {
    std::pair<Host::size_type,Host::size_type> range =
      this_thread.work_range( output.dimension_0() );

    for ( Host::size_type i0 = range.first ; i0 < range.second ; ++i0 ) {
    for ( Host::size_type i1 = 0 ; i1 < output.dimension_1() ; ++i1 ) {
      output(i0,i1) = input(i0,i1);
    }}

    this_thread.barrier();
  }
};

template< class OutputView , class InputView >
struct HostViewRemap< OutputView , 1 , InputView , 1 >
  : public HostThreadWorker<void>
{
  const OutputView output ;
  const InputView  input ;

  HostViewRemap( const OutputView & arg_out , const InputView & arg_in )
    : output( arg_out ), input( arg_in )
    { HostThreadWorker<void>::execute( *this ); }

  void execute_on_thread( HostThread & this_thread ) const
  {
    std::pair<Host::size_type,Host::size_type> range =
      this_thread.work_range( output.dimension_0() );

    for ( Host::size_type i0 = range.first ; i0 < range.second ; ++i0 ) {
      output(i0) = input(i0);
    }

    this_thread.barrier();
  }
};

template< class OutputView , class InputView >
struct HostViewRemap< OutputView , 0 , InputView , 0 >
  : public HostThreadWorker<void>
{
  const OutputView output ;
  const InputView  input ;

  HostViewRemap( const OutputView & arg_out , const InputView & arg_in )
    : output( arg_out ), input( arg_in )
    { *arg_out = *arg_in ; }
};

//----------------------------------------------------------------------------

template< class DataTypeOutput , class LayoutOutput ,
          class DataTypeInput ,  class LayoutInput >
struct Factory< View< DataTypeOutput , LayoutOutput , Host > ,
                View< DataTypeInput ,  LayoutInput ,  Host > >
{
  typedef View< DataTypeOutput , LayoutOutput , Host > output_type ;
  typedef View< DataTypeInput ,  LayoutInput ,  Host > input_type ;

  static inline
  void deep_copy( const output_type & output ,
                  const input_type  & input )
  {
    typedef typename output_type::value_type output_value_type ;
    typedef typename input_type ::value_type input_value_type ;

    typedef typename output_type::shape_type output_shape ;
    typedef typename input_type ::shape_type input_shape ;

    assert_shapes_equal_dimension( output.m_shape , input.m_shape );

    if ( output.m_shape == input.m_shape ) {
      HostParallelCopy<output_value_type,input_value_type>(
        output.ptr_on_device() ,
        input. ptr_on_device() ,
        allocation_count( output.m_shape ) );
    }
    else {
      HostViewRemap< output_type , output_shape::rank ,
                     input_type  , input_shape::rank >( output , input );
    }
  }
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray

#endif /* #ifndef KOKKOS_HOST_VIEW_HPP */


