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

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------

template< typename ValueType >
struct Factory< MDArray< ValueType , Host > , void >
  : public HostThreadWorker<void>
{
public:
  typedef Host::size_type              size_type ;
  typedef MDArray< ValueType , Host >  output_type ;

private:

  output_type dst ;

  Factory( const output_type & arg_dst ) : dst( arg_dst ) {}

  /** \brief  First touch the allocated memory with the proper thread.
   *
   *  My memory is dense, contiguous, and slowest stride on parallel dim.
   */
  void execute_on_thread( Impl::HostThread & this_thread ) const
  {
    const std::pair<size_type,size_type> range =
      this_thread.work_range( dst.m_map.allocation_size() );

    ValueType * const x_end = dst.ptr_on_device() + range.second ;
    ValueType *       x     = dst.ptr_on_device() + range.first ;
    while ( x_end != x ) { *x++ = 0 ; }

    this_thread.barrier();
  }

public:

  static output_type create( const std::string & label ,
                             size_t nP , size_t n1 , size_t n2 , size_t n3 ,
                             size_t n4 , size_t n5 , size_t n6 , size_t n7 )
  {
    typedef MemoryManager< Host > memory_manager ;

    output_type array ;

    array.m_map.template assign< ValueType >(nP,n1,n2,n3,n4,n5,n6,n7);
    array.m_memory.allocate( array.m_map.allocation_size() , label );

    memory_manager::disable_memory_view_tracking();

    Factory driver( array );

    memory_manager::enable_memory_view_tracking();
    
    HostThreadWorker<void>::execute( driver );

    return array ;
  }
};

//----------------------------------------------------------------------------

template< typename ValueType >
struct Factory< MDArray< ValueType , Host > ,
                MDArray< ValueType , Host > > : HostThreadWorker<void>
{
public:
  typedef Host::size_type             size_type ;
  typedef MDArray< ValueType , Host > output_type ;
  typedef MDArray< ValueType , Host > input_type ;

private:

  output_type  output ;
  input_type  input ;

  Factory( const output_type & arg_output , const input_type & arg_input )
    : output( arg_output ) , input( arg_input ) {} 

  void execute_on_thread( Impl::HostThread & this_thread ) const
  {
    const std::pair<size_type,size_type> range =
      this_thread.work_range( output.m_map.allocation_size() );

    ValueType * const x_end = output.ptr_on_device() + range.second ;
    ValueType *       x     = output.ptr_on_device() + range.first ;
    const ValueType * y     = input.ptr_on_device() + range.first ;
    while ( x_end != x ) { *x++ = *y++ ; }

    this_thread.barrier();
  }

public:

  static inline
  void deep_copy( const output_type & arg_output , const input_type & arg_input )
  {
    typedef MemoryManager< Host::memory_space > memory_manager ;
    memory_manager::disable_memory_view_tracking();
    Factory driver( arg_output , arg_input );
    memory_manager::enable_memory_view_tracking();
    HostThreadWorker<void>::execute( driver );
  }

  static inline
  output_type create( const input_type & input )
  {
    return Factory< output_type , void >::create(
      std::string(),
      ( 0 < input.rank() ? input.dimension(0) : 0 ),
      ( 1 < input.rank() ? input.dimension(1) : 0 ),
      ( 2 < input.rank() ? input.dimension(2) : 0 ),
      ( 3 < input.rank() ? input.dimension(3) : 0 ),
      ( 4 < input.rank() ? input.dimension(4) : 0 ),
      ( 5 < input.rank() ? input.dimension(5) : 0 ),
      ( 6 < input.rank() ? input.dimension(6) : 0 ),
      ( 7 < input.rank() ? input.dimension(7) : 0 ) );
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace Impl
} // namespace Kokkos

