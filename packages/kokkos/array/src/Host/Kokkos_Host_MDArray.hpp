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
{
  typedef Host::size_type              size_type ;
  typedef MDArray< ValueType , Host >  output_type ;

  static output_type create( const std::string & label ,
                             size_t nP , size_t n1 , size_t n2 , size_t n3 ,
                             size_t n4 , size_t n5 , size_t n6 , size_t n7 )
  {
    typedef MemoryManager< Host > memory_manager ;

    output_type array ;

    array.m_map.template assign< ValueType >(nP,n1,n2,n3,n4,n5,n6,n7);
    array.m_memory.allocate( array.m_map.allocation_size() , label );

    HostParallelFill<ValueType>( array.m_memory.ptr_on_device() , 0 ,
                                 array.m_map.allocation_size() );

    return array ;
  }
};

//----------------------------------------------------------------------------

template< typename ValueType , class Device >
struct Factory< MDArray< ValueType , HostMapped< Device > > , void >
{
  typedef Host::size_type                              size_type ;
  typedef MDArray< ValueType , HostMapped< Device > >  output_type ;

  static output_type create( const std::string & label ,
                             size_t nP , size_t n1 , size_t n2 , size_t n3 ,
                             size_t n4 , size_t n5 , size_t n6 , size_t n7 )
  {
    typedef MemoryManager< Host > memory_manager ;

    output_type array ;

    array.m_map.template assign< ValueType >(nP,n1,n2,n3,n4,n5,n6,n7);
    array.m_memory.allocate( array.m_map.allocation_size() , label );

    HostParallelFill<ValueType>( array.m_memory.ptr_on_device() , 0 ,
                                 array.m_map.allocation_size() );

    return array ;
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< typename ValueType >
struct Factory< MDArray< ValueType , Host > , MDArray< ValueType , Host > >
{
public:
  typedef Host::size_type             size_type ;
  typedef MDArray< ValueType , Host > output_type ;
  typedef MDArray< ValueType , Host > input_type ;

  static inline
  void deep_copy( const output_type & output ,
                  const input_type  & input )
  {
    HostParallelCopy<ValueType,ValueType>( output.ptr_on_device() ,
                                           input. ptr_on_device() ,
                                           output.m_map.allocation_size() );
  }

  // Called by create_mirror
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

template< typename ValueType , class DeviceOutput , class DeviceInput >
struct HostMDArrayMappedDeepCopy : public HostThreadWorker<void>
{
  typedef Host::size_type                     size_type ;
  typedef MDArray< ValueType , DeviceOutput > output_type ;
  typedef MDArray< ValueType , DeviceInput >  input_type ;

private:

  output_type dst ;
  input_type  src ;
  const size_type N1 ;
  const size_type N2 ;
  const size_type N3 ;
  const size_type N4 ;
  const size_type N5 ;
  const size_type N6 ;
  const size_type N7 ;

  void copy8( const std::pair<size_type,size_type> range ) const
  {
    for ( size_type i0 = range.first ; i0 < range.second ; ++i0 ) {
    for ( size_type i1 = 0 ; i1 < N1 ; ++i1 ) {
    for ( size_type i2 = 0 ; i2 < N2 ; ++i2 ) {
    for ( size_type i3 = 0 ; i3 < N3 ; ++i3 ) {
    for ( size_type i4 = 0 ; i4 < N4 ; ++i4 ) {
    for ( size_type i5 = 0 ; i5 < N5 ; ++i5 ) {
    for ( size_type i6 = 0 ; i6 < N6 ; ++i6 ) {
    for ( size_type i7 = 0 ; i7 < N7 ; ++i7 ) {
      dst(i0,i1,i2,i3,i4,i5,i6,i7) = src(i0,i1,i2,i3,i4,i5,i6,i7);
    }}}}}}}}
  }

  void copy7( const std::pair<size_type,size_type> range ) const
  {
    for ( size_type i0 = range.first ; i0 < range.second ; ++i0 ) {
    for ( size_type i1 = 0 ; i1 < N1 ; ++i1 ) {
    for ( size_type i2 = 0 ; i2 < N2 ; ++i2 ) {
    for ( size_type i3 = 0 ; i3 < N3 ; ++i3 ) {
    for ( size_type i4 = 0 ; i4 < N4 ; ++i4 ) {
    for ( size_type i5 = 0 ; i5 < N5 ; ++i5 ) {
    for ( size_type i6 = 0 ; i6 < N6 ; ++i6 ) {
      dst(i0,i1,i2,i3,i4,i5,i6) = src(i0,i1,i2,i3,i4,i5,i6);
    }}}}}}}
  }

  void copy6( const std::pair<size_type,size_type> range ) const
  {
    for ( size_type i0 = range.first ; i0 < range.second ; ++i0 ) {
    for ( size_type i1 = 0 ; i1 < N1 ; ++i1 ) {
    for ( size_type i2 = 0 ; i2 < N2 ; ++i2 ) {
    for ( size_type i3 = 0 ; i3 < N3 ; ++i3 ) {
    for ( size_type i4 = 0 ; i4 < N4 ; ++i4 ) {
    for ( size_type i5 = 0 ; i5 < N5 ; ++i5 ) {
      dst(i0,i1,i2,i3,i4,i5) = src(i0,i1,i2,i3,i4,i5);
    }}}}}}
  }

  void copy5( const std::pair<size_type,size_type> range ) const
  {
    for ( size_type i0 = range.first ; i0 < range.second ; ++i0 ) {
    for ( size_type i1 = 0 ; i1 < N1 ; ++i1 ) {
    for ( size_type i2 = 0 ; i2 < N2 ; ++i2 ) {
    for ( size_type i3 = 0 ; i3 < N3 ; ++i3 ) {
    for ( size_type i4 = 0 ; i4 < N4 ; ++i4 ) {
      dst(i0,i1,i2,i3,i4) = src(i0,i1,i2,i3,i4);
    }}}}}
  }

  void copy4( const std::pair<size_type,size_type> range ) const
  {
    for ( size_type i0 = range.first ; i0 < range.second ; ++i0 ) {
    for ( size_type i1 = 0 ; i1 < N1 ; ++i1 ) {
    for ( size_type i2 = 0 ; i2 < N2 ; ++i2 ) {
    for ( size_type i3 = 0 ; i3 < N3 ; ++i3 ) {
      dst(i0,i1,i2,i3) = src(i0,i1,i2,i3);
    }}}}
  }

  void copy3( const std::pair<size_type,size_type> range ) const
  {
    for ( size_type i0 = range.first ; i0 < range.second ; ++i0 ) {
    for ( size_type i1 = 0 ; i1 < N1 ; ++i1 ) {
    for ( size_type i2 = 0 ; i2 < N2 ; ++i2 ) {
      dst(i0,i1,i2) = src(i0,i1,i2);
    }}}
  }

  void copy2( const std::pair<size_type,size_type> range ) const
  {
    for ( size_type i0 = range.first ; i0 < range.second ; ++i0 ) {
    for ( size_type i1 = 0 ; i1 < N1 ; ++i1 ) {
      dst(i0,i1) = src(i0,i1);
    }}
  }

  void copy1( const std::pair<size_type,size_type> range ) const
  {
    for ( size_type i0 = range.first ; i0 < range.second ; ++i0 ) {
      dst(i0) = src(i0);
    }
  }

  void execute_on_thread( HostThread & this_thread ) const
  {
    const std::pair<size_type,size_type> range =
      this_thread.work_range( dst.dimension(0) );

    switch( dst.rank() ) {
    case 1: copy1( range ); break ;
    case 2: copy2( range ); break ;
    case 3: copy3( range ); break ;
    case 4: copy4( range ); break ;
    case 5: copy5( range ); break ;
    case 6: copy6( range ); break ;
    case 7: copy7( range ); break ;
    case 8: copy8( range ); break ;
    }
  }

  HostMDArrayMappedDeepCopy( const output_type & output ,
                             const input_type  & input )
    : dst( output ), src( input )
    , N1( 1 < output.rank() ? output.dimension(1) : 0 )
    , N2( 2 < output.rank() ? output.dimension(2) : 0 )
    , N3( 3 < output.rank() ? output.dimension(3) : 0 )
    , N4( 4 < output.rank() ? output.dimension(4) : 0 )
    , N5( 5 < output.rank() ? output.dimension(5) : 0 )
    , N6( 6 < output.rank() ? output.dimension(6) : 0 )
    , N7( 7 < output.rank() ? output.dimension(7) : 0 )
    {}

public:

  static inline
  void deep_copy( const output_type & output ,
                  const input_type  & input )
  {
    typedef MemoryManager< Host > memory_manager ;

    memory_manager::disable_memory_view_tracking();

    const HostMDArrayMappedDeepCopy driver( output , input );

    memory_manager::enable_memory_view_tracking();

    HostThreadWorker<void>::execute( driver );
  }
};


template< typename ValueType , class Device >
struct Factory< MDArray< ValueType , Host > ,
                MDArray< ValueType , HostMapped< Device > > >
{
  typedef MDArray< ValueType , Host >                 output_type ;
  typedef MDArray< ValueType , HostMapped< Device > > input_type ;

  inline static
  void deep_copy( const output_type & output , const input_type & input )
  {
    HostMDArrayMappedDeepCopy<ValueType, Host , HostMapped< Device > >
      ::deep_copy( output , input );
  }
};

template< typename ValueType , class Device >
struct Factory< MDArray< ValueType , HostMapped< Device > > ,
                MDArray< ValueType , Host > >
{
  typedef MDArray< ValueType , HostMapped< Device > > output_type ;
  typedef MDArray< ValueType , Host >                 input_type ;

  inline static
  void deep_copy( const output_type & output , const input_type & input )
  {
    HostMDArrayMappedDeepCopy<ValueType,HostMapped< Device > , Host >
      ::deep_copy( output , input );
  }
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace Kokkos

