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

#ifndef KOKKOS_HOST_CRSARRAY_HPP
#define KOKKOS_HOST_CRSARRAY_HPP

namespace Kokkos {
namespace Impl {

template< typename ValueType >
class CreateCrsArray< ValueType , Host > {
public:

  template< typename IteratorType >
  static
  CrsArray<ValueType,Host>
    create( const std::string & label ,
            const IteratorType row_count_begin ,
            const IteratorType row_count_end )
  {
    typedef Host::size_type size_type ;

    CrsArray<ValueType,Host> array ;

    size_type row_count = 0 ;
    size_type value_count = 0 ;

    for ( IteratorType i = row_count_begin ; i != row_count_end ; ++i ) {
      ++row_count ;
      value_count += *i ;
    }

    array.m_row_count = row_count ;
    array.m_offset.allocate( row_count + 1 , label );
    array.m_values.allocate( value_count , label );

    size_type * const offset = array.m_offset.ptr_on_device();

    row_count = 0 ;
    value_count = 0 ;
    offset[0] = 0 ;

    for ( IteratorType i = row_count_begin ; i != row_count_end ; ++i ) {
      offset[ ++row_count ] = ( value_count += *i );
    }

    return array ;
  }

  static
  CrsArray<ValueType,Host> create( const CrsArray<ValueType,Host> & rhs )
  {
    typedef Host::size_type size_type ;

    CrsArray<ValueType,Host> array ;

    const size_t row_count = rhs.row_dimension();
    const size_t size      = rhs.size();

    array.m_row_count = row_count ;
    array.m_offset.allocate( row_count + 1 , std::string() );
    array.m_values.allocate( size , std::string() );

    size_type * const dst = array.m_offset.ptr_on_device();
    size_type * const src = rhs  .m_offset.ptr_on_device();

    for ( size_type i = 0 ; i <= row_count ; ++i ) { dst[i] = src[i] ; }

    return array ;
  }

};

//----------------------------------------------------------------------------

template< typename ValueType >
class CreateMirror< CrsArray< ValueType , Host > , false /* copy */ >
{
public:
  typedef  CrsArray< ValueType , Host >  View ;
  typedef  typename View::HostView       HostView ;

  static
  HostView create( const View & v )
  { return CreateCrsArray<ValueType,Host>::create( v ); }
};

//----------------------------------------------------------------------------

template< typename ValueType >
class DeepCopy< CrsArray< ValueType , Host > ,
                CrsArray< ValueType , Host > >
  : public HostThreadWorker<void>
{
public:
  typedef Host::size_type  size_type ;

  typedef CrsArray< ValueType , Host > dst_type ;
  typedef CrsArray< ValueType , Host > src_type ;

  void execute_on_thread( HostThread & this_thread ) const
  {
    {
      const std::pair<size_type,size_type> range =
        this_thread.work_range( src.row_dimension() + 1 );

      size_type * const x_end = dst.m_offset.ptr_on_device() + range.second ;
      size_type *       x_ptr = dst.m_offset.ptr_on_device() + range.first ;
      const size_type * y_ptr = src.m_offset.ptr_on_device() + range.first ;
      while ( x_end != x_ptr ) *x_ptr++ = *y_ptr++ ;
    }

    {
      const std::pair<size_type,size_type> range =
        this_thread.work_range( src.size() );

      ValueType * const x_end = dst.m_values.ptr_on_device() + range.second ;
      ValueType *       x_ptr = dst.m_values.ptr_on_device() + range.first ;
      const ValueType * y_ptr = src.m_values.ptr_on_device() + range.first ;
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

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_HOST_CRSARRAY_HPP */

