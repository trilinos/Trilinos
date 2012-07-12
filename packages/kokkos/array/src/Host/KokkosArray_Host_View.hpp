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
#include <impl/KokkosArray_Shape_macros.hpp>
#include <impl/KokkosArray_View_macros.hpp>
#include <KokkosArray_Clear_macros.hpp>


namespace KokkosArray {
namespace Impl {

//----------------------------------------------------------------------------
/** \brief  Primary creation / allocation factory based upon shape */

template< typename DataType , class LayoutType >
struct Factory< View< DataType , LayoutType , Host > , void >
{
  typedef View< DataType , LayoutType , Host >  output_type ;
  typedef typename output_type::shape_type      shape_type ;
  typedef typename output_type::memory_space    memory_space ;
  typedef typename output_type::value_type      value_type ;

  inline static
  output_type create( const std::string & label , const shape_type & shape )
  {
    const size_t count = allocation_count( shape );

    output_type output ;

    output.m_shape = shape ;
    output.m_ptr_on_device = (value_type *)
      memory_space::allocate( label ,
                              typeid(value_type) ,
                              sizeof(value_type) ,
                              count );

    HostParallelFill<value_type>( output.m_ptr_on_device , 0 , count );

    return output ;
  }
};

//----------------------------------------------------------------------------

/** \brief  Deep copy a single value */
template< typename DataType , class LayoutType >
struct Factory< View< DataType , LayoutType , Host > , DataType >
{
  typedef View< DataType , LayoutType , Host >  output_type ;
  typedef DataType                              input_type ;
  typedef typename output_type::shape_type      shape_type ;
  typedef typename output_type::value_type      value_type ;

  typedef typename assert_shape_is_rank< shape_type , 0 >::type ok_rank ;

  typedef typename StaticAssertAssignable< value_type , DataType >::type ok_assign ;

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

  typedef typename assert_shape_is_rank< shape_type , 0 >::type ok_rank ;

  typedef typename StaticAssertAssignable< DataType , value_type >::type ok_assign ;

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
    typedef Host::memory_space            memory_space ;
    typedef typename output_type::value_type  value_type ;

    output_type output ;
    output.m_shape = input.m_shape ;
    output.m_ptr_on_device = (value_type *)
      memory_space::allocate( std::string("mirror") ,
                              typeid(value_type) ,
                              sizeof(value_type) ,
                              allocation_count( output.m_shape ) );

    return output ;
  }
};

/** \brief  Deep different but compatible arrays */
template< class OutDataType , class OutLayoutType ,
          class InDataType ,  class InLayoutType >
struct Factory< View< OutDataType , OutLayoutType , Host > ,
                View< InDataType ,  InLayoutType ,  Host > >
{
public:
  typedef View< OutDataType , OutLayoutType , Host > output_type ;
  typedef View< InDataType ,  InLayoutType ,  Host > input_type ;

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
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#if 0

template< typename DataType , class LayoutOutput , class LayoutInput >
struct Factory< View< DataType , LayoutOutput , Host > ,
                View< DataType , LayoutInput , Host > >
{
  typedef View< DataType , LayoutOutput , Host >   output_type ;
  typedef View< DataType , LayoutInput , Host > input_type ;

  inline static
  void deep_copy( const output_type & output , const input_type & input )
  {
    typedef typename output_type::value_type value_type ;

    HostIndexMapDeepCopy< value_type,
                          typename output_type::index_map ,
                          typename input_type ::index_map >
      ::deep_copy( output.m_data , output.m_index_map ,
                   input .m_data , input .m_index_map );
  }
};

#endif

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray

#endif /* #ifndef KOKKOS_HOST_VIEW_HPP */

