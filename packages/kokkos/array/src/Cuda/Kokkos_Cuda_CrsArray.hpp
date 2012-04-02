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

#ifndef KOKKOS_CUDA_CRSARRAY_HPP
#define KOKKOS_CUDA_CRSARRAY_HPP

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------

/** \brief  The hostview is identically mapped */

template< typename SizeType >
struct Factory< CrsArray< void , Cuda , SizeType > ,
                CrsArray< void , HostMapped< Cuda > , SizeType > >
{
  typedef CrsArray< void, Cuda >                        output_type ;
  typedef CrsArray< void, HostMapped<Cuda>, SizeType >  input_type ;

  static inline
  void deep_copy( const output_type & output , const input_type & input )
  {
    typedef typename output_type::value_type value_type ;

    const size_t size_row = sizeof(SizeType)*(output.m_row_count + 1 );

    MemoryManager< Cuda >::
      copy_to_device_from_host( output.m_row_map.ptr_on_device(),
                                input.m_row_map.ptr_on_device(),
                                size_row );
  }
};

//----------------------------------------------------------------------------

template< class ArrayType , typename SizeType >
struct Factory< CrsArray< ArrayType , Cuda , SizeType > ,
                CrsArray< ArrayType , HostMapped< Cuda > , SizeType > >
{
  typedef CrsArray< ArrayType, Cuda >                        output_type ;
  typedef CrsArray< ArrayType, HostMapped<Cuda>, SizeType >  input_type ;

  static inline
  void deep_copy( const output_type & output , const input_type & input )
  {
    typedef typename output_type::value_type value_type ;

    const size_t size_row = sizeof(SizeType)*(output.m_row_count + 1 );
    const size_t size_data = sizeof(value_type)*output.m_map.allocation_size();

    MemoryManager< Cuda >::
      copy_to_device_from_host( output.m_row_map.ptr_on_device(),
                                input.m_row_map.ptr_on_device(),
                                size_row );
    MemoryManager< Cuda >::
      copy_to_device_from_host( output.m_data.ptr_on_device(),
                                input.m_data.ptr_on_device(),
                                size_data );
  }
};

//----------------------------------------------------------------------------

template< typename SizeType >
struct Factory< CrsArray< void , HostMapped< Cuda > , SizeType > ,
                CrsArray< void , Cuda , SizeType > >
{
  typedef CrsArray< void, HostMapped<Cuda>, SizeType > output_type ;
  typedef CrsArray< void, Cuda >                       input_type ;

  static void deep_copy( const output_type & output , const input_type & input )
  {
    const size_t size_row = sizeof(SizeType)*(output.m_row_count + 1 );

    MemoryManager< Cuda >::
      copy_to_host_from_device( output.m_row_map.ptr_on_device(),
                                input.m_row_map.ptr_on_device(),
                                size_row );
  }

  static inline
  output_type create( const input_type & input )
  {
    output_type output ;

    output.m_row_count = input.m_row_count ;
    output.m_row_map.allocate( output.m_row_count + 1 , std::string() );

    deep_copy( output , input );

    return output ;
  }
};

template< class ArrayType , typename SizeType >
struct Factory< CrsArray< ArrayType , HostMapped< Cuda > , SizeType > ,
                CrsArray< ArrayType , Cuda , SizeType > >
{
  typedef CrsArray< ArrayType, HostMapped<Cuda>, SizeType > output_type ;
  typedef CrsArray< ArrayType, Cuda >                       input_type ;
  typedef typename output_type::value_type value_type ;

  static void deep_copy( const output_type & output , const input_type & input )
  {
    const size_t size_row = sizeof(SizeType)*(output.m_row_count + 1 );
    const size_t size_data = sizeof(value_type)*output.m_index_map.allocation_size();

    MemoryManager< Cuda >::
      copy_to_host_from_device( output.m_row_map.ptr_on_device(),
                                input.m_row_map.ptr_on_device(),
                                size_row );

    MemoryManager< Cuda >::
      copy_to_host_from_device( output.m_data.ptr_on_device(),
                                input.m_data.ptr_on_device(),
                                size_data );
  }

  static inline
  output_type create( const input_type & input )
  {
    output_type output ;

    output.m_row_count = input.m_row_count ;
    output.m_index_map.template assign< value_type >( input.entry_dimension(0) );
    output.m_row_map.allocate( output.m_row_count + 1 , std::string() );
    output.m_data.allocate( output.m_index_map.allocation_size(), std::string());

    deep_copy( output , input );

    return output ;
  }
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace Kokkos

#endif /* #ifndef KOKKOS_CUDA_CRSARRAY_HPP */

