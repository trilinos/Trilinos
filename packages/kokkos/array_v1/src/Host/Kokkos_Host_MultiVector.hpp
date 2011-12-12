/** \HEADER
 *************************************************************************
 *
 *                            Kokkos
 *                 Copyright 2010 Sandia Corporation
 *
 *  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 *  the U.S. Government retains certain rights in this software.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are
 *  met:
 *
 *  1. Redistributions of source code must retain the above copyright
 *  notice, this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the Corporation nor the names of the
 *  contributors may be used to endorse or promote products derived from
 *  this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 *  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 *  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 *  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *************************************************************************
 */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< typename ValueType >
class Initialize< MultiVector<ValueType,Host> >
  : public HostThreadWorker<void> {
public:
  typedef Host::size_type  size_type ;

  typedef MultiVector< ValueType , Host > dst_type ;

  void execute_on_thread( HostThread & this_thread ) const
  {
    const std::pair<size_type,size_type> range =
      this_thread.work_range( dst.length() );

    for ( size_type i = 0 ; i < dst.count() ; ++i ) {
      ValueType * const x_end = & dst(0,i) + range.second ;
      ValueType *       x_ptr = & dst(0,i) + range.first ;
      while ( x_end != x_ptr ) *x_ptr++ = 0 ;
    }

    this_thread.barrier();
  }

  Initialize( const dst_type & arg_dst ) : dst( arg_dst ) {}

  static void run( const dst_type & arg_dst )
  {
    typedef MemoryManager< Host::memory_space > memory_manager ;
    memory_manager::disable_memory_view_tracking();
    Initialize driver( arg_dst );
    memory_manager::enable_memory_view_tracking();
    HostThreadWorker<void>::execute( driver );
  }

private:
  dst_type dst ;
};

//----------------------------------------------------------------------------

template< typename ValueType >
class DeepCopy< MultiVector< ValueType , Host > ,
                MultiVector< ValueType , Host > >
  : public HostThreadWorker<void> {
public:
  typedef Host::size_type  size_type ;

  typedef MultiVector< ValueType , Host > dst_type ;
  typedef MultiVector< ValueType , Host > src_type ;

  void execute_on_thread( HostThread & this_thread ) const
  {
    const std::pair<size_type,size_type> range =
      this_thread.work_range( dst.length() );

    for ( size_type i = 0 ; i < dst.count() ; ++i ) {
      const dst_type x( dst , i );
      const src_type y( src , i );
      ValueType * const x_end = x.ptr_on_device() + range.second ;
      ValueType *       x_ptr = x.ptr_on_device() + range.first ;
      const ValueType * y_ptr = y.ptr_on_device() + range.first ;
      while ( x_end != x_ptr ) *x_ptr++ = *y_ptr++ ;
    }

    this_thread.barrier();
  }

  static void run( const dst_type & arg_dst , const src_type & arg_src )
  {
    typedef MemoryManager< Host::memory_space > memory_manager ;
    memory_manager::disable_memory_view_tracking();
    DeepCopy driver( arg_dst , arg_src );
    memory_manager::enable_memory_view_tracking();
    HostThreadWorker<void>::execute( driver );
  }

private:

  dst_type dst ;
  src_type src ;

  DeepCopy( const dst_type & arg_dst , const src_type & arg_src )
    : dst( arg_dst ), src( arg_src ) {}
};

} // namespace Impl
} // namespace Kokkos

