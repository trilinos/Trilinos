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

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#ifndef KOKKOS_HOST_MULTIVECTOR_HPP
#define KOKKOS_HOST_MULTIVECTOR_HPP

#include <string>

#include <KokkosArray_Host_macros.hpp>
#include <impl/KokkosArray_MultiVector_macros.hpp>
#include <KokkosArray_Clear_macros.hpp>

namespace KokkosArray {
namespace Impl {

template< typename ValueType >
struct Factory< MultiVector< ValueType , Host > , void >
{
  //------------------------------------

  template< class Device >
  struct Initialize : public HostThreadWorker<void> {
  private:

    typedef MultiVector< ValueType , Device > type ;

    type dst ;

    Initialize( const type & arg_dst ) : dst( arg_dst ) {}

    void execute_on_thread( Impl::HostThread & this_thread ) const
    {
      typedef Host::size_type size_type ;

      const std::pair<size_type,size_type> range =
        this_thread.work_range( dst.length() );

      for ( size_type i = 0 ; i < dst.count() ; ++i ) {
        ValueType *       x     = & dst( range.first , i );
        ValueType * const x_end = & dst( range.second , i );

        while ( x_end != x ) { *x++ = 0 ; }
      }

      this_thread.barrier();
    }

  public:

    static void run( const type & arg )
    {
      Initialize driver( arg );

      HostThreadWorker<void>::execute( driver );
    }
  };

  //------------------------------------

  typedef MultiVector< ValueType , Host > output_type ;

  static inline
  output_type 
  create( const std::string & label , size_t length , size_t count )
  {
    typedef Host::memory_space_new  memory_space ;
    typedef typename output_type::value_type value_type ;
    typedef typename output_type::view_type  view_type ;
    typedef typename view_type::shape_type   shape_type ;

    output_type output ;

    // Want different 'first touch' and the View factory will give

    output.m_memory.m_shape = Factory< shape_type , memory_space >
                                ::create(length,count);

    output.m_memory.m_ptr_on_device = (value_type *)
      memory_space::allocate( label ,
                              typeid(value_type) ,
                              sizeof(value_type) ,
                              allocation_count( output.m_memory.m_shape ) );

    Initialize< Host >::run( output );

    return output ;
  }
};


template< typename ValueType , class Device >
struct Factory< MultiVector< ValueType , HostMapped< Device > > , 
                MultiVector< ValueType , Device > >
{
  typedef MultiVector< ValueType , HostMapped< Device > > output_type ;
  typedef MultiVector< ValueType , Device >               input_type ;

  static inline
  output_type create( const input_type & input )
  {
    typedef Host::memory_space_new  memory_space ;
    typedef typename output_type::value_type value_type ;
    typedef typename output_type::shape_type shape_type ;

    output_type output ;

    // Want different 'first touch' and the View factory will give

    output.m_memory.m_shape = input.m_memory.m_shape ;

    output.m_memory.m_ptr_on_device = (value_type *)
      memory_space::allocate( std::string("mirror") ,
                              typeid(value_type) ,
                              sizeof(value_type) ,
                              allocation_count( output.m_shape ) );

    typedef typename
      Factory< MultiVector< ValueType , Host > , void >
        ::template Initialize< HostMapped< Device > > Initialize ;

    Initialize::run( output );

    return output ;
  }
};

//----------------------------------------------------------------------------

template< typename ValueType >
struct Factory< MultiVector< ValueType , Host > ,
                MultiVector< ValueType , Host > >
  : public HostThreadWorker<void>
{
  typedef Host::size_type  size_type ;

  typedef MultiVector< ValueType , Host > output_type ;
  typedef MultiVector< ValueType , Host > input_type ;

private:

  output_type  output ;
  input_type   input ;

  Factory( const output_type & arg_output , const input_type & arg_input ,
           const size_t arg_length , const size_t arg_count )
    : output( arg_output )
    , input( arg_input )
    {}

  void execute_on_thread( HostThread & this_thread ) const
  {
    const size_type count = output.m_memory.dimension_1();

    const std::pair<size_type,size_type> range =
      this_thread.work_range( output.m_memory.dimension_0() );

    for ( size_type i = 0 ; i < count ; ++i ) {
      const output_type x( output , i );
      const input_type  y( input , i );

      ValueType * const x_end = x.ptr_on_device() + range.second ;
      ValueType *       x_ptr = x.ptr_on_device() + range.first ;
      const ValueType * y_ptr = y.ptr_on_device() + range.first ;
      while ( x_end != x_ptr ) *x_ptr++ = *y_ptr++ ;
    }

    this_thread.barrier();
  }

public:

  static inline
  void deep_copy( const output_type & output , const input_type & input )
  {
    typedef MemoryManager< Host::memory_space > HostMemorySpace ;
    Factory driver( output , input , output.length() , output.count() );
    HostThreadWorker<void>::execute( driver );
  }

  static inline
  void deep_copy( const output_type & output , const input_type & input ,
                  const size_t length )
  {
    typedef MemoryManager< Host::memory_space > HostMemorySpace ;
    Factory driver( output , input , length , 1 );
    HostThreadWorker<void>::execute( driver );
  }

  static inline
  output_type create( const input_type & input )
  {
    output_type output ;
    output.m_memory = Factory< typename output_type::view_type ,
                               typename input_type::view_type >
                        ::create( input.m_memory );
    return output ;
  }
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray

#endif /* #ifndef KOKKOS_HOST_MULTIVECTOR_HPP */

