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

#include <KokkosArray_Cuda_macros.hpp>
#include <impl/KokkosArray_Shape_macros.hpp>
#include <impl/KokkosArray_View_macros.hpp>
#include <KokkosArray_Clear_macros.hpp>


namespace KokkosArray {
namespace Impl {

//----------------------------------------------------------------------------

template< typename DataType , class LayoutSpec >
struct Factory< View< DataType , LayoutSpec , Cuda > , unsigned_<1> >
{
  typedef View< DataType , LayoutSpec , Cuda >  output_type ;
  typedef typename output_type::shape_type      shape_type ;

  typedef typename StaticAssertSame<
                      unsigned_< shape_type::rank_dynamic > ,
                      unsigned_< 1 > >::type ok_rank ;

  static output_type create( const std::string & label , size_t n0 )
  {
    typedef Cuda::memory_space_new  memory_space ;
    typedef typename output_type::value_type value_type ;
    typedef typename output_type::shape_type shape_type ;

    output_type output ;

    output.m_shape = Factory< shape_type , memory_space >::create(n0);
    output.m_ptr_on_device = (value_type *)
      memory_space::allocate( label ,
                              typeid(value_type) ,
                              sizeof(value_type) ,
                              allocation_count( output.m_shape ) );
    return output ;
  }
};

template< typename DataType , class LayoutSpec >
struct Factory< View< DataType , LayoutSpec , Cuda > , unsigned_<2> >
{
  typedef View< DataType , LayoutSpec , Cuda >  output_type ;
  typedef typename output_type::shape_type      shape_type ;

  typedef typename StaticAssertSame<
                      unsigned_< shape_type::rank_dynamic > ,
                      unsigned_< 2 > >::type ok_rank ;

  static output_type create( const std::string & label ,
                             size_t n0 , size_t n1 )
  {
    typedef Cuda::memory_space_new  memory_space ;
    typedef typename output_type::value_type value_type ;
    typedef typename output_type::shape_type shape_type ;

    output_type output ;

    output.m_shape = Factory< shape_type , memory_space >::create(n0,n1);
    output.m_ptr_on_device = (value_type *)
      memory_space::allocate( label ,
                              typeid(value_type) ,
                              sizeof(value_type) ,
                              allocation_count( output.m_shape ) );
    return output ;
  }
};

template< typename DataType , class LayoutSpec >
struct Factory< View< DataType , LayoutSpec , Cuda > , unsigned_<3> >
{
  typedef View< DataType , LayoutSpec , Cuda >  output_type ;
  typedef typename output_type::shape_type      shape_type ;

  typedef typename StaticAssertSame<
                      unsigned_< shape_type::rank_dynamic > ,
                      unsigned_< 3 > >::type ok_rank ;

  static output_type create( const std::string & label ,
                             size_t n0 , size_t n1 , size_t n2 )
  {
    typedef Cuda::memory_space_new  memory_space ;
    typedef typename output_type::value_type value_type ;
    typedef typename output_type::shape_type shape_type ;

    output_type output ;

    output.m_shape = Factory< shape_type , memory_space >::create(n0,n1,n2);
    output.m_ptr_on_device = (value_type *)
      memory_space::allocate( label ,
                              typeid(value_type) ,
                              sizeof(value_type) ,
                              allocation_count( output.m_shape ) );
    return output ;
  }
};

template< typename DataType , class LayoutSpec >
struct Factory< View< DataType , LayoutSpec , Cuda > , unsigned_<4> >
{
  typedef View< DataType , LayoutSpec , Cuda >  output_type ;
  typedef typename output_type::shape_type      shape_type ;

  typedef typename StaticAssertSame<
                      unsigned_< shape_type::rank_dynamic > ,
                      unsigned_< 4 > >::type ok_rank ;

  static output_type create( const std::string & label ,
                             size_t n0 , size_t n1 , size_t n2 , size_t n3 )
  {
    typedef Cuda::memory_space_new  memory_space ;
    typedef typename output_type::value_type value_type ;
    typedef typename output_type::shape_type shape_type ;

    output_type output ;

    output.m_shape = Factory< shape_type , memory_space >::create(n0,n1,n2,n3);
    output.m_ptr_on_device = (value_type *)
      memory_space::allocate( label ,
                              typeid(value_type) ,
                              sizeof(value_type) ,
                              allocation_count( output.m_shape ) );
    return output ;
  }
};

template< typename DataType , class LayoutSpec >
struct Factory< View< DataType , LayoutSpec , Cuda > , unsigned_<5> >
{
  typedef View< DataType , LayoutSpec , Cuda >  output_type ;
  typedef typename output_type::shape_type      shape_type ;

  typedef typename StaticAssertSame<
                      unsigned_< shape_type::rank_dynamic > ,
                      unsigned_< 5 > >::type ok_rank ;

  static output_type create( const std::string & label ,
                             size_t n0 , size_t n1 , size_t n2 , size_t n3 ,
                             size_t n4 )
  {
    typedef Cuda::memory_space_new  memory_space ;
    typedef typename output_type::value_type value_type ;
    typedef typename output_type::shape_type shape_type ;

    output_type output ;

    output.m_shape = Factory< shape_type , memory_space >
                       ::create(n0,n1,n2,n3,n4);
    output.m_ptr_on_device = (value_type *)
      memory_space::allocate( label ,
                              typeid(value_type) ,
                              sizeof(value_type) ,
                              allocation_count( output.m_shape ) );
    return output ;
  }
};

template< typename DataType , class LayoutSpec >
struct Factory< View< DataType , LayoutSpec , Cuda > , unsigned_<6> >
{
  typedef View< DataType , LayoutSpec , Cuda >  output_type ;
  typedef typename output_type::shape_type      shape_type ;

  typedef typename StaticAssertSame<
                      unsigned_< shape_type::rank_dynamic > ,
                      unsigned_< 6 > >::type ok_rank ;

  static output_type create( const std::string & label ,
                             size_t n0 , size_t n1 , size_t n2 , size_t n3 ,
                             size_t n4 , size_t n5 )
  {
    typedef Cuda::memory_space_new  memory_space ;
    typedef typename output_type::value_type value_type ;
    typedef typename output_type::shape_type shape_type ;

    output_type output ;

    output.m_shape = Factory< shape_type , memory_space >
                       ::create(n0,n1,n2,n3,n4,n5);
    output.m_ptr_on_device = (value_type *)
      memory_space::allocate( label ,
                              typeid(value_type) ,
                              sizeof(value_type) ,
                              allocation_count( output.m_shape ) );
    return output ;
  }
};

template< typename DataType , class LayoutSpec >
struct Factory< View< DataType , LayoutSpec , Cuda > , unsigned_<7> >
{
  typedef View< DataType , LayoutSpec , Cuda >  output_type ;
  typedef typename output_type::shape_type      shape_type ;

  typedef typename StaticAssertSame<
                      unsigned_< shape_type::rank_dynamic > ,
                      unsigned_< 7 > >::type ok_rank ;

  static output_type create( const std::string & label ,
                             size_t n0 , size_t n1 , size_t n2 , size_t n3 ,
                             size_t n4 , size_t n5 , size_t n6 )
  {
    typedef Cuda::memory_space_new  memory_space ;
    typedef typename output_type::value_type value_type ;
    typedef typename output_type::shape_type shape_type ;

    output_type output ;

    output.m_shape = Factory< shape_type , memory_space >
                       ::create(n0,n1,n2,n3,n4,n5,n6);
    output.m_ptr_on_device = (value_type *)
      memory_space::allocate( label ,
                              typeid(value_type) ,
                              sizeof(value_type) ,
                              allocation_count( output.m_shape ) );
    return output ;
  }
};

template< typename DataType , class LayoutSpec >
struct Factory< View< DataType , LayoutSpec , Cuda > , unsigned_<8> >
{
  typedef View< DataType , LayoutSpec , Cuda >  output_type ;
  typedef typename output_type::shape_type      shape_type ;

  typedef typename StaticAssertSame<
                      unsigned_< shape_type::rank_dynamic > ,
                      unsigned_< 8 > >::type ok_rank ;

  static output_type create( const std::string & label ,
                             size_t n0 , size_t n1 , size_t n2 , size_t n3 ,
                             size_t n4 , size_t n5 , size_t n6 , size_t n7 )
  {
    typedef Cuda::memory_space_new  memory_space ;
    typedef typename output_type::value_type value_type ;
    typedef typename output_type::shape_type shape_type ;

    output_type output ;

    output.m_shape = Factory< shape_type , memory_space >
                       ::create(n0,n1,n2,n3,n4,n5,n6,n7);
    output.m_ptr_on_device = (value_type *)
      memory_space::allocate( label ,
                              typeid(value_type) ,
                              sizeof(value_type) ,
                              allocation_count( output.m_shape ) );
    return output ;
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< class DataType , class LayoutSpec >
struct Factory< View< DataType , LayoutSpec , Cuda > ,
                View< DataType , LayoutSpec , Cuda > >
{
public:
  typedef View< DataType , LayoutSpec , Cuda > output_type ;
  typedef View< DataType , LayoutSpec , Cuda > input_type ;

  static inline
  void deep_copy( const output_type & output ,
                  const input_type  & input )
  {
    typedef Cuda::memory_space_new            memory_space ;
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

  // Called by create_mirror
  static inline
  output_type create( const input_type & input )
  {
    typedef Cuda::memory_space_new            memory_space ;
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

template< class OutDataType , class OutLayoutSpec ,
          class InDataType , class InLayoutSpec >
struct Factory< View< OutDataType , OutLayoutSpec , Cuda > ,
                View< InDataType , InLayoutSpec , Cuda > >
{
public:
  typedef View< OutDataType , OutLayoutSpec , Cuda > output_type ;
  typedef View< InDataType , InLayoutSpec , Cuda > input_type ;

  static inline
  void deep_copy( const output_type & output ,
                  const input_type  & input )
  {
    typedef Cuda::memory_space_new            memory_space ;
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
};

template< class OutDataType , class OutLayoutSpec ,
          class InDataType , class InLayoutSpec >
struct Factory< View< OutDataType , OutLayoutSpec , Host > ,
                View< InDataType , InLayoutSpec , Cuda > >
{
public:
  typedef View< OutDataType , OutLayoutSpec , Host > output_type ;
  typedef View< InDataType , InLayoutSpec , Cuda > input_type ;

  static inline
  void deep_copy( const output_type & output ,
                  const input_type  & input )
  {
    typedef Cuda::memory_space_new            memory_space ;
    typedef typename output_type::value_type  value_type ;

    assert_shapes_are_equal( output.m_shape , input.m_shape );

    const size_t size = output.m_shape.value_size *
                        allocation_count( output.m_shape );

    memory_space::copy_to_host_from_device( output.ptr_on_device() ,
                                            input. ptr_on_device() ,
                                            size );
  }

  // Called by create_mirror
  static inline
  output_type create( const input_type & input )
  {
    typedef Host::memory_space_new            memory_space ;
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

/* Copy data to Cuda from HostMirror */
template< class OutDataType , class OutLayoutSpec ,
          class InDataType , class InLayoutSpec >
struct Factory< View< OutDataType , OutLayoutSpec , Cuda > ,
                View< InDataType , InLayoutSpec , Host > >
{
public:
  typedef View< OutDataType , OutLayoutSpec , Cuda > output_type ;
  typedef View< InDataType , InLayoutSpec , Host > input_type ;

  static inline
  void deep_copy( const output_type & output ,
                  const input_type  & input )
  {
    typedef Cuda::memory_space_new            memory_space ;
    typedef typename output_type::value_type  value_type ;

    assert_shapes_are_equal( output.m_shape , input.m_shape );

    const size_t size = output.m_shape.value_size *
                        allocation_count( output.m_shape );

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

