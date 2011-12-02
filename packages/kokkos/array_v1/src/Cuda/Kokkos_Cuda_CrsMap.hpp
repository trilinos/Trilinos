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

template< typename SizeType >
class CreateCrsMap< CrsMap< Cuda , SizeType > , 
                    CrsMap< Cuda , SizeType > > {
public:
  typedef CrsMap< Cuda , SizeType > type ;

  typedef MemoryManager< Cuda > memory_manager ;

  template< typename IteratorType >
  static
  type create( const std::string & label ,
               const IteratorType second_count_begin ,
               const IteratorType second_count_end ) 
  {
    type m ;
    
    SizeType first_count = 0 ;
    SizeType total_count = 0 ;
    
    for ( IteratorType i = second_count_begin ; i != second_count_end ; ++i ) {
      ++first_count ;
      total_count += *i ;
    }

    m.m_first_count = first_count ;
    m.m_total_count = total_count ;
    m.m_memory.allocate( first_count + 1 , label );

    SizeType * const offset = new SizeType[ first_count + 1 ];

    first_count = 0 ;
    total_count = 0 ;
    offset[0] = 0 ; 
    
    for ( IteratorType i = second_count_begin ; i != second_count_end ; ++i ) {
      offset[ ++first_count ] = ( total_count += *i );
    }

    memory_manager::copy_to_device_from_host(
      m.m_memory.ptr_on_device(), offset ,
      ( first_count + 1 ) * sizeof(SizeType) );

    delete[] offset ;

    return m ;
  }

  static
  type create( const type & rhs ) { return rhs ; }
};


template< typename DstSizeType , typename SrcSizeType >
class CreateCrsMap< CrsMap< Cuda , DstSizeType > , 
                    CrsMap< Cuda , SrcSizeType > > {
public:
  typedef CrsMap< Cuda , DstSizeType > type ;

  static
  type create( const CrsMap<Cuda,SrcSizeType> & rhs )
  {
    type m ;

    m.m_first_count = (DstSizeType) rhs.m_first_count ;
    m.m_total_count = (DstSizeType) rhs.m_total_count ;
    m.m_memory.allocate( m.m_first_count + 1 , std::string() );

    DstSizeType * const dst = m  .m_memory.ptr_on_device();
    SrcSizeType * const src = rhs.m_memory.ptr_on_device();

    CudaParallelCopy<DstSizeType,SrcSizeType>(
      dst , src , m.m_first_count + 1 );

    return m ;
  }
};

//----------------------------------------------------------------------------

template<>
class CreateCrsMap< Cuda , Cuda > {
public:
  typedef CrsMap< Cuda , Cuda::size_type > type ;

  template< typename IteratorType >
  static
  type create( const std::string & label ,
               const IteratorType second_count_begin ,
               const IteratorType second_count_end ) 
  {
    return CreateCrsMap<type,type>::create( label , second_count_begin ,
                                                    second_count_end );
  }
};

//----------------------------------------------------------------------------

template< typename SizeType >
class CreateCrsMap< CrsMap< Host , SizeType > , 
                    CrsMap< Cuda , SizeType > > {
public:
  typedef CrsMap< Host , SizeType > type ;

  static
  type create( const CrsMap<Cuda,SizeType> & rhs )
  {
    type m ;

    m.m_first_count = rhs.m_first_count ;
    m.m_total_count = rhs.m_total_count ;
    m.m_memory.allocate( m.m_first_count + 1 , std::string() );

    SizeType * const dst = m  .m_memory.ptr_on_device();
    SizeType * const src = rhs.m_memory.ptr_on_device();

    MemoryManager< Cuda >::copy_to_host_from_device(
      dst , src , ( m.m_first_count + 1 ) * sizeof(SizeType) );

    return m ;
  }
};

template< typename DstSizeType , typename SrcSizeType >
class CreateCrsMap< CrsMap< Host , DstSizeType > , 
                    CrsMap< Cuda , SrcSizeType > > {
public:
  typedef CrsMap< Host , DstSizeType > type ;

  static
  type create( const CrsMap<Cuda,SrcSizeType> & rhs )
  {
    type m ;

    m.m_first_count = (DstSizeType) rhs.m_first_count ;
    m.m_total_count = (DstSizeType) rhs.m_total_count ;
    m.m_memory.allocate( m.m_first_count + 1 , std::string() );

    DstSizeType * const dst = m  .m_memory.ptr_on_device();
    SrcSizeType * const src = rhs.m_memory.ptr_on_device();

    Impl::MemoryView< DstSizeType , Cuda > tmp ;

    tmp.allocate( m.m_first_count + 1 , std::string() );

    CudaParallelCopy<DstSizeType,SrcSizeType>(
      tmp.ptr_on_device() , src , ( m.m_first_count + 1 ) );

    MemoryManager< Cuda >::copy_to_host_from_device(
      dst , tmp.ptr_on_device() ,
      ( m.m_first_count + 1 ) * sizeof(DstSizeType) );

    return m ;
  }
};

template< typename SizeType >
class CreateCrsMap< Host , CrsMap< Cuda , SizeType > > {
public:
  typedef CrsMap< Host , Host::size_type > type ;

  static type create( const CrsMap<Cuda,SizeType> & rhs )
  { return CreateCrsMap<type,CrsMap<Cuda,SizeType> >::create( rhs ); }
};

//----------------------------------------------------------------------------

template< typename SizeType >
class CreateCrsMap< CrsMap< Cuda , SizeType > , 
                    CrsMap< Host , SizeType > > {
public:
  typedef CrsMap< Cuda , SizeType > type ;

  static
  type create( const CrsMap<Host,SizeType> & rhs )
  {
    type m ;

    m.m_first_count = rhs.m_first_count ;
    m.m_total_count = rhs.m_total_count ;
    m.m_memory.allocate( m.m_first_count + 1 , std::string() );

    SizeType * const dst = m  .m_memory.ptr_on_device();
    SizeType * const src = rhs.m_memory.ptr_on_device();

    MemoryManager< Cuda >::copy_to_device_from_host(
      dst , src , ( m.m_first_count + 1 ) * sizeof(SizeType) );

    return m ;
  }
};

template< typename DstSizeType , typename SrcSizeType >
class CreateCrsMap< CrsMap< Cuda , DstSizeType > , 
                    CrsMap< Host , SrcSizeType > > {
public:
  typedef CrsMap< Cuda , DstSizeType > type ;

  static
  type create( const CrsMap<Host,SrcSizeType> & rhs )
  {
    type m ;

    m.m_first_count = (DstSizeType) rhs.m_first_count ;
    m.m_total_count = (DstSizeType) rhs.m_total_count ;
    m.m_memory.allocate( m.m_first_count + 1 , std::string() );

    DstSizeType * const dst = m  .m_memory.ptr_on_device();
    SrcSizeType * const src = rhs.m_memory.ptr_on_device();

    Impl::MemoryView< SrcSizeType , Cuda > tmp ;

    tmp.allocate( m.m_first_count + 1 , std::string() );

    MemoryManager< Cuda >::copy_to_device_from_host(
      tmp.ptr_on_device() , src ,
      ( m.m_first_count + 1 ) * sizeof(SrcSizeType) );

    CudaParallelCopy<DstSizeType,SrcSizeType>(
      dst , tmp.ptr_on_device() , m.m_first_count + 1 );

    return m ;
  }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

