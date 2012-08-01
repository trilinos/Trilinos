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

#ifndef KOKKOS_CUDA_VIEW_HPP
#define KOKKOS_CUDA_VIEW_HPP

#include <KokkosArray_Host.hpp>
#include <KokkosArray_View.hpp>

#include <KokkosArray_Cuda_macros.hpp>
#include <impl/KokkosArray_ViewOperLeft_macros.hpp>
#include <impl/KokkosArray_ViewOperRight_macros.hpp>
#include <impl/KokkosArray_View_macros.hpp>
#include <KokkosArray_Clear_macros.hpp>

//----------------------------------------------------------------------------

namespace KokkosArray {

template< class DataType , class LayoutType >
void View< DataType , LayoutType , Cuda >::create(
  const std::string & label ,
  const View< DataType , LayoutType , Cuda >::shape_type shape )
{
  const size_t count = Impl::allocation_count( shape );

  oper_type::m_shape = shape ;
  oper_type::m_ptr_on_device = (value_type *)
    memory_space::allocate( label , 
                            typeid(value_type) , 
                            sizeof(value_type) , 
                            count );
}

}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template< typename DataType , class LayoutType >
struct Factory< View< DataType , LayoutType , Cuda > , DataType >
{
  typedef View< DataType , LayoutType , Cuda >  output_type ;
  typedef DataType                              input_type ;
  typedef typename output_type::value_type      value_type ;

  static inline
  void deep_copy( const output_type & output , const input_type & input )
  {
    typedef typename output_type::shape_type shape_type ;

    typedef typename
      assert_shape_is_rank_zero< shape_type >::type ok_rank ;

    typedef typename
      StaticAssertAssignable< value_type , DataType >::type ok_assign ;

    typedef Cuda::memory_space  memory_space ;

    memory_space::copy_to_device_from_device( output.ptr_on_device() ,
                                              & input ,
                                              sizeof(value_type) );
  }
};

template< typename DataType , class LayoutType >
struct Factory< DataType , View< DataType , LayoutType , Cuda > >
{
  typedef DataType                              output_type ;
  typedef View< DataType , LayoutType , Cuda >  input_type ;
  typedef typename output_type::value_type      value_type ;

  static inline
  void deep_copy( output_type & output , const input_type & input )
  {
    typedef typename output_type::shape_type shape_type ;

    typedef typename
      assert_shape_is_rank_zero< shape_type >::type ok ;

    typedef typename
      StaticAssertAssignable< DataType , value_type >::type ok_assign ;

    typedef Cuda::memory_space  memory_space ;

    memory_space::copy_to_device_from_device( & output ,
                                              input.ptr_on_device() ,
                                              sizeof(value_type) );
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< class DataType , class LayoutType >
struct Factory< View< DataType , LayoutType , Cuda > ,
                View< DataType , LayoutType , Cuda > >
{
public:
  typedef View< DataType , LayoutType , Cuda > output_type ;
  typedef View< DataType , LayoutType , Cuda > input_type ;

  static inline
  void deep_copy( const output_type & output ,
                  const input_type  & input )
  {
    typedef Cuda::memory_space                memory_space ;
    typedef typename output_type::value_type  value_type ;

    if ( output != input ) {

      assert_shapes_are_equal( output.m_shape , input.m_shape );


      const size_t size = output.m_shape.value_size *
                          allocation_count( output.m_shape );

      memory_space::copy_to_device_from_device( output.ptr_on_device() ,
                                                input. ptr_on_device() ,
                                                size );
    }
  }

  static inline
  void deep_copy( const output_type & output ,
                  const input_type  & input ,
                  const size_t count )
  {
    typedef Cuda::memory_space            memory_space ;
    typedef typename output_type::value_type value_type ;

    // Only for rank-1 arrays, or arrays where higher ranks are one

    assert_shape_effective_rank1_at_leastN( output.m_shape , count );
    assert_shape_effective_rank1_at_leastN( input.m_shape , count );

    const size_t size = output.m_shape.value_size * count ;

    memory_space::copy_to_device_from_device( output.ptr_on_device() ,
                                              input. ptr_on_device() ,
                                              size );
  }
};

template< class OutDataType , class OutLayoutType ,
          class InDataType , class InLayoutType >
struct Factory< View< OutDataType , OutLayoutType , Cuda > ,
                View< InDataType , InLayoutType , Cuda > >
{
public:
  typedef View< OutDataType , OutLayoutType , Cuda > output_type ;
  typedef View< InDataType , InLayoutType , Cuda > input_type ;

  static inline
  void deep_copy( const output_type & output ,
                  const input_type  & input )
  {
    typedef Cuda::memory_space            memory_space ;
    typedef typename output_type::value_type  value_type ;

    if ( output != input ) {

      assert_shapes_are_equal( output.m_shape , input.m_shape );

      const size_t size = output.m_shape.value_size *
                          allocation_count( output.m_shape );

      memory_space::copy_to_device_from_device( output.ptr_on_device() ,
                                                input. ptr_on_device() ,
                                                size );
    }
  }

  static inline
  void deep_copy( const output_type & output ,
                  const input_type  & input ,
                  const size_t count )
  {
    typedef Cuda::memory_space            memory_space ;
    typedef typename output_type::value_type value_type ;

    // Only for rank-1 arrays, or arrays where higher ranks are one

    assert_shape_effective_rank1_at_leastN( output.m_shape , count );
    assert_shape_effective_rank1_at_leastN( input.m_shape , count );

    const size_t size = output.m_shape.value_size * count ;

    memory_space::copy_to_device_from_device( output.ptr_on_device() ,
                                              input. ptr_on_device() ,
                                              size );
  }
};

template< class OutDataType , class OutLayoutType ,
          class InDataType , class InLayoutType >
struct Factory< View< OutDataType , OutLayoutType , Host > ,
                View< InDataType , InLayoutType , Cuda > >
{
public:
  typedef View< OutDataType , OutLayoutType , Host > output_type ;
  typedef View< InDataType , InLayoutType , Cuda > input_type ;

  static inline
  void deep_copy( const output_type & output ,
                  const input_type  & input )
  {
    typedef Cuda::memory_space            memory_space ;
    typedef typename output_type::value_type  value_type ;

    assert_shapes_are_equal( output.m_shape , input.m_shape );

    const size_t size = output.m_shape.value_size *
                        allocation_count( output.m_shape );

    memory_space::copy_to_host_from_device( output.ptr_on_device() ,
                                            input. ptr_on_device() ,
                                            size );
  }

  static inline
  void deep_copy( const output_type & output ,
                  const input_type  & input ,
                  const size_t count )
  {
    typedef Cuda::memory_space  memory_space ;
    typedef typename output_type::value_type value_type ;

    // Only for rank-1 arrays, or arrays where higher ranks are one

    assert_shape_effective_rank1_at_leastN( output.m_shape , count );
    assert_shape_effective_rank1_at_leastN( input.m_shape , count );

    const size_t size = output.m_shape.value_size * count ;

    memory_space::copy_to_host_from_device( output.ptr_on_device() ,
                                            input. ptr_on_device() ,
                                            size );
  }

  // Called by create_mirror
  static inline
  output_type create( const input_type & input )
  {
    return output_type("mirror",input.m_shape);
  }
};

/* Copy data to Cuda from HostMirror */
template< class OutDataType , class OutLayoutType ,
          class InDataType , class InLayoutType >
struct Factory< View< OutDataType , OutLayoutType , Cuda > ,
                View< InDataType , InLayoutType , Host > >
{
public:
  typedef View< OutDataType , OutLayoutType , Cuda > output_type ;
  typedef View< InDataType , InLayoutType , Host > input_type ;

  static inline
  void deep_copy( const output_type & output ,
                  const input_type  & input )
  {
    typedef Cuda::memory_space            memory_space ;
    typedef typename output_type::value_type  value_type ;

    assert_shapes_are_equal( output.m_shape , input.m_shape );

    const size_t size = output.m_shape.value_size *
                        allocation_count( output.m_shape );

    memory_space::copy_to_device_from_host( output.ptr_on_device() ,
                                            input. ptr_on_device() ,
                                            size );
  }

  static inline
  void deep_copy( const output_type & output ,
                  const input_type  & input ,
                  const size_t count )
  {
    typedef Cuda::memory_space  memory_space ;
    typedef typename output_type::value_type value_type ;

    // Only for rank-1 arrays, or arrays where higher ranks are one

    assert_shape_effective_rank1_at_leastN( output.m_shape , count );
    assert_shape_effective_rank1_at_leastN( input.m_shape , count );

    const size_t size = output.m_shape.value_size * count ;

    memory_space::copy_to_device_from_host( output.ptr_on_device() ,
                                            input. ptr_on_device() ,
                                            size );
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray

#endif /* #ifndef KOKKOS_CUDA_VIEW_HPP */

