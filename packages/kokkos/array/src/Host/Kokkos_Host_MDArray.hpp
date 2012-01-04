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

/** \brief  This is a device type hidden from the user to handle
 *          host memory space with a different mdarray map.
 *
 *  Use case: 
 *  class MDArray<T,Cuda> {
 *     typedef MDArray<T,HostMapped< Cuda::mdarray_map > > HostMirror ;
 *  };
 */
template< class MDArrayMapType >
class HostMapped {
public:
  typedef Host::size_type     size_type ;
  typedef Host::memory_space  memory_space ;
  typedef MDArrayMapType      mdarray_map ;
};

template< typename ValueType >
class MDArrayHostMirror< ValueType , Host::mdarray_map > {
public:
  typedef MDArray< ValueType , Host > type ;
};

template< typename ValueType , class MDArrayMapType >
class MDArrayHostMirror {
public:
  typedef MDArray< ValueType , HostMapped< MDArrayMapType > > type ;
};

//----------------------------------------------------------------------------

template< typename ValueType >
class Initialize< MDArray< ValueType , Host > >
  : public HostThreadWorker<void>
{
public:
  typedef Host::size_type              size_type ;
  typedef MDArray< ValueType , Host >  dst_type ;

  dst_type   dst ;
  size_type  chunk ;

  Initialize( const dst_type & arg_dst )
    : dst( arg_dst )
    , chunk( 1 )
    {
      for ( size_type r = 1 ; r < dst.rank() ; ++r ) {
        chunk *= dst.dimension(r);
      }
    }

  /** \brief  First touch the allocated memory with the proper thread.
   *
   *  My memory is dense, contiguous, and slowest stride on parallel dim.
   */
  void execute_on_thread( Impl::HostThread & this_thread ) const
  {
    const std::pair<size_type,size_type> range =
      this_thread.work_range( dst.dimension(0) );

    ValueType * const x_end = dst.ptr_on_device() + chunk * range.second ;
    ValueType *       x     = dst.ptr_on_device() + chunk * range.first ;
    while ( x_end != x ) { *x++ = 0 ; }

    this_thread.barrier();
  }

  static void run( const dst_type & arg_dst )
  {
    typedef MemoryManager< Host::memory_space > memory_manager ;
    memory_manager::disable_memory_view_tracking();
    Initialize driver( arg_dst );
    memory_manager::enable_memory_view_tracking();
    HostThreadWorker<void>::execute( driver );
  }
};

//----------------------------------------------------------------------------

template< typename DstType , unsigned N >
class MDArrayInitializeHostPermute : public HostThreadWorker<void>
{
private:
  MDArrayInitializeHostPermute();
  MDArrayInitializeHostPermute( const MDArrayInitializeHostPermute & );
  MDArrayInitializeHostPermute & operator = ( const MDArrayInitializeHostPermute & );
public:
  typedef Host::size_type size_type ;

  DstType m_dst ;
  Rank<N> rank ;

  explicit
  MDArrayInitializeHostPermute( const DstType & arg_dst )
  : m_dst( arg_dst ), rank() {}

  static void run( const DstType & arg_dst )
  {
    typedef MemoryManager< Host::memory_space > memory_manager ;
    memory_manager::disable_memory_view_tracking();
    MDArrayInitializeHostPermute driver( arg_dst );
    memory_manager::enable_memory_view_tracking();
    HostThreadWorker<void>::execute( driver );
  }

  void execute_on_thread( Impl::HostThread & this_thread ) const
  {
    const typename DstType::value_type zero = 0 ;

    const std::pair<size_type,size_type> range =
      this_thread.work_range( m_dst.dimension(0) );

    for ( size_type iP = range.first ; iP < range.second ; ++iP ) {
      m_dst.fill( rank , iP , zero );
    }

    this_thread.barrier();
  }
};

template< typename ValueType , class MDArrayMap >
class Initialize< MDArray< ValueType , HostMapped< MDArrayMap > > > {
public:
  typedef MDArray< ValueType , HostMapped< MDArrayMap > > dst_type ;

  static void run( const dst_type & dst )
  {
    switch( dst.rank() ) {
    case 1 : MDArrayInitializeHostPermute< dst_type,1>::run(dst); break ;
    case 2 : MDArrayInitializeHostPermute< dst_type,2>::run(dst); break ;
    case 3 : MDArrayInitializeHostPermute< dst_type,3>::run(dst); break ;
    case 4 : MDArrayInitializeHostPermute< dst_type,4>::run(dst); break ;
    case 5 : MDArrayInitializeHostPermute< dst_type,5>::run(dst); break ;
    case 6 : MDArrayInitializeHostPermute< dst_type,6>::run(dst); break ;
    case 7 : MDArrayInitializeHostPermute< dst_type,7>::run(dst); break ;
    case 8 : MDArrayInitializeHostPermute< dst_type,8>::run(dst); break ;
    }
  }
};

//----------------------------------------------------------------------------

template< typename ValueType >
class DeepCopy< MDArray< ValueType , Host > ,
                MDArray< ValueType , Host > > : HostThreadWorker<void>
{
public:
  typedef Host::size_type             size_type ;
  typedef MDArray< ValueType , Host > dst_type ;
  typedef MDArray< ValueType , Host > src_type ;

private:

  dst_type  dst ;
  src_type  src ;
  size_type chunk ;

  DeepCopy( const dst_type & arg_dst , const src_type & arg_src )
    : dst( arg_dst )
    , src( arg_src )
    , chunk( 1 )
    {
      for ( size_type r = 1 ; r < dst.rank() ; ++r ) {
        chunk *= dst.dimension(r);
      }
    }

public:

  void execute_on_thread( Impl::HostThread & this_thread ) const
  {
    const std::pair<size_type,size_type> range =
      this_thread.work_range( dst.dimension(0) );

    ValueType * const x_end = dst.ptr_on_device() + chunk * range.second ;
    ValueType *       x     = dst.ptr_on_device() + chunk * range.first ;
    const ValueType * y     = src.ptr_on_device() + chunk * range.first ;
    while ( x_end != x ) { *x++ = *y++ ; }

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
};

//----------------------------------------------------------------------------

template< class DstType , class SrcType, unsigned N >
class MDArrayDeepCopyHostPermute : HostThreadWorker<void> {
public:
  const DstType dst ;
  const SrcType src ;
  const Rank<N> rank ;

  MDArrayDeepCopyHostPermute( const DstType & arg_dst ,
                              const SrcType & arg_src )
    : dst( arg_dst ) , src( arg_src ), rank()
    { HostThreadWorker<void>::execute( *this ); }

  void execute_on_thread( Impl::HostThread & this_thread ) const
  {
    typedef Host::size_type size_type ;

    const std::pair<size_type,size_type> range =
      this_thread.work_range( dst.dimension(0) );

    for ( size_type i0 = range.first ; i0 < range.second ; ++i0 ) {
      dst.assign( rank , i0 , src );
    }

    this_thread.barrier();
  }
};

template< typename ValueType , class MDArrayMap >
class DeepCopy< MDArray< ValueType , Host > ,
                MDArray< ValueType , HostMapped< MDArrayMap > > >
  : HostThreadWorker<void>
{
public:
  typedef MDArray< ValueType , Host >                     dst_type ;
  typedef MDArray< ValueType , HostMapped< MDArrayMap > > src_type ;

  static void run( const dst_type & dst , const src_type & src )
  {
    typedef MemoryManager< Host::memory_space > memory_manager ;
    memory_manager::disable_memory_view_tracking();
    switch( dst.rank() ) {
    case 1 : MDArrayDeepCopyHostPermute<dst_type,src_type,1>(dst,src); break ;
    case 2 : MDArrayDeepCopyHostPermute<dst_type,src_type,2>(dst,src); break ;
    case 3 : MDArrayDeepCopyHostPermute<dst_type,src_type,3>(dst,src); break ;
    case 4 : MDArrayDeepCopyHostPermute<dst_type,src_type,4>(dst,src); break ;
    case 5 : MDArrayDeepCopyHostPermute<dst_type,src_type,5>(dst,src); break ;
    case 6 : MDArrayDeepCopyHostPermute<dst_type,src_type,6>(dst,src); break ;
    case 7 : MDArrayDeepCopyHostPermute<dst_type,src_type,7>(dst,src); break ;
    case 8 : MDArrayDeepCopyHostPermute<dst_type,src_type,8>(dst,src); break ;
    }
    memory_manager::enable_memory_view_tracking();
  }
};

template< typename ValueType , class MDArrayMap >
class DeepCopy< MDArray< ValueType , HostMapped< MDArrayMap > > ,
                MDArray< ValueType , Host > >
  : HostThreadWorker<void>
{
public:
  typedef MDArray< ValueType , HostMapped< MDArrayMap > > dst_type ;
  typedef MDArray< ValueType , Host >                     src_type ;

  static void run( const dst_type & dst , const src_type & src )
  {
    typedef MemoryManager< Host::memory_space > memory_manager ;
    memory_manager::disable_memory_view_tracking();
    switch( dst.rank() ) {
    case 1 : MDArrayDeepCopyHostPermute<dst_type,src_type,1>(dst,src); break ;
    case 2 : MDArrayDeepCopyHostPermute<dst_type,src_type,2>(dst,src); break ;
    case 3 : MDArrayDeepCopyHostPermute<dst_type,src_type,3>(dst,src); break ;
    case 4 : MDArrayDeepCopyHostPermute<dst_type,src_type,4>(dst,src); break ;
    case 5 : MDArrayDeepCopyHostPermute<dst_type,src_type,5>(dst,src); break ;
    case 6 : MDArrayDeepCopyHostPermute<dst_type,src_type,6>(dst,src); break ;
    case 7 : MDArrayDeepCopyHostPermute<dst_type,src_type,7>(dst,src); break ;
    case 8 : MDArrayDeepCopyHostPermute<dst_type,src_type,8>(dst,src); break ;
    }
    memory_manager::enable_memory_view_tracking();
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace Impl
} // namespace Kokkos

