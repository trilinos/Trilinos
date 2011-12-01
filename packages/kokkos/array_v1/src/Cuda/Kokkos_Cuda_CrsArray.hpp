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

template< typename ValueType , typename SizeType >
class CreateCrsArray< ValueType , Cuda , SizeType > {
public:

  typedef MemoryManager< Cuda > memory_manager ;

  template< typename IteratorType >
  static
  CrsArray<ValueType,Cuda,SizeType>
    create( const std::string & label ,
            const IteratorType row_count_begin ,
            const IteratorType row_count_end ) 
  {
    CrsArray<ValueType,Cuda,SizeType> array ;
    
    SizeType row_count = 0 ;
    SizeType value_count = 0 ;
    
    for ( IteratorType i = row_count_begin ; i != row_count_end ; ++i ) {
      ++row_count ;
      value_count += *i ;
    }
    
    array.m_row_count = row_count ;
    array.m_value_count = value_count ;
    array.m_offset.allocate( row_count + 1 , label );
    array.m_values.allocate( value_count , label );

    SizeType * const offset = new SizeType[ row_count + 1 ];

    row_count = 0 ;
    value_count = 0 ;
    offset[0] = 0 ; 
    
    for ( IteratorType i = row_count_begin ; i != row_count_end ; ++i ) {
      offset[ ++row_count ] = ( value_count += *i );
    }

    memory_manager::copy_to_device_from_host(
      array.m_offset.ptr_on_device(), offset ,
      ( row_count + 1 ) * sizeof(SizeType) );

    delete[] offset ;

    return array ;
  }

  static
  CrsArray<ValueType,Host,SizeType>
    create_mirror( const CrsArray<ValueType,Cuda,SizeType> & rhs )
  {
    CrsArray<ValueType,Host,SizeType> array ;

    const size_t row_count   = rhs.m_row_count ;
    const size_t value_count = rhs.m_value_count ;

    array.m_row_count   = row_count ;
    array.m_value_count = value_count ;
    array.m_offset.allocate( row_count + 1 , std::string() );
    array.m_values.allocate( value_count , std::string() );

    SizeType * const dst = array.m_offset.ptr_on_device();
    SizeType * const src = rhs  .m_offset.ptr_on_device();

    memory_manager::copy_to_host_from_device(
      array.m_offset.ptr_on_device() ,
      rhs  .m_offset.ptr_on_device() ,
      ( row_count + 1 ) * sizeof(SizeType) );

    return array ;
  }
};

//----------------------------------------------------------------------------

template< typename ValueType >
class CreateMirror< CrsArray< ValueType , Cuda > , false /* copy */ >
{
public:
  typedef  CrsArray< ValueType , Cuda >  View ;
  typedef  typename View::HostView       HostView ;

  static
  HostView create( const View & v )
  { return CreateCrsArray<ValueType,Cuda,Cuda::size_type>::create_mirror( v ); }
};

//----------------------------------------------------------------------------

template< typename ValueType >
class DeepCopy< CrsArray< ValueType , Cuda > ,
                CrsArray< ValueType , Cuda > > {
public:
  static void run( const CrsArray< ValueType , Cuda > & dst ,
                   const CrsArray< ValueType , Cuda > & src )
  {
    const Cuda::size_type offset_size =
      ( src.m_row_count + 1 ) * sizeof(Cuda::size_type);

    const Cuda::size_type values_size =
      src.m_value_count * sizeof(ValueType);

    MemoryManager< Cuda >::
      copy_to_device_from_device( dst.m_offset.ptr_on_device(),
                                  src.m_offset.ptr_on_device(),
                                  offset_size );

    MemoryManager< Cuda >::
      copy_to_device_from_device( dst.m_values.ptr_on_device(),
                                  src.m_values.ptr_on_device(),
                                  values_size );
  }
};

template< typename ValueType >
class DeepCopy< CrsArray< ValueType , Cuda > ,
                typename CrsArray< ValueType , Cuda >::HostView > {
public:
  typedef CrsArray< ValueType , Cuda >                    dst_type ;
  typedef typename CrsArray< ValueType , Cuda >::HostView src_type ;

  static void run( const dst_type & dst , const src_type & src )
  {
    const Cuda::size_type offset_size =
      ( src.m_row_count + 1 ) * sizeof(Cuda::size_type);

    const Cuda::size_type values_size =
      src.m_value_count * sizeof(ValueType);

    MemoryManager< Cuda >::
      copy_to_device_from_host( dst.m_offset.ptr_on_device(),
                                src.m_offset.ptr_on_device(),
                                offset_size );

    MemoryManager< Cuda >::
      copy_to_device_from_host( dst.m_values.ptr_on_device(),
                                src.m_values.ptr_on_device(),
                                values_size );
  }
};

template< typename ValueType >
class DeepCopy< typename CrsArray< ValueType , Cuda >::HostView ,
                CrsArray< ValueType , Cuda > > {
public:
  typedef typename CrsArray< ValueType , Cuda >::HostView dst_type ;
  typedef CrsArray< ValueType , Cuda >                    src_type ;

  static void run( const dst_type & dst , const src_type & src )
  {
    const Cuda::size_type offset_size =
      ( src.m_row_count + 1 ) * sizeof(Cuda::size_type);

    const Cuda::size_type values_size =
      src.m_value_count * sizeof(ValueType);

    MemoryManager< Cuda >::
      copy_to_host_from_device( dst.m_offset.ptr_on_device(),
                                src.m_offset.ptr_on_device(),
                                offset_size );

    MemoryManager< Cuda >::
      copy_to_host_from_device( dst.m_values.ptr_on_device(),
                                src.m_values.ptr_on_device(),
                                values_size );
  }
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace Kokkos

