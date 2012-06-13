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

#ifndef KOKKOS_CUDA_MULTIVECTOR_HPP
#define KOKKOS_CUDA_MULTIVECTOR_HPP

#include <string>

#include <KokkosArray_Cuda_macros.hpp>
#include <impl/KokkosArray_MultiVector_macros.hpp>
#include <KokkosArray_Clear_macros.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template< typename ValueType >
struct Factory< MultiVector< ValueType , Cuda > ,
                MultiVector< ValueType , Cuda > > {
public:
  static void deep_copy( const MultiVector< ValueType , Cuda > & dst ,
                         const MultiVector< ValueType , Cuda > & src )
  {
    const Cuda::size_type alloc_size =
      dst.m_count == 1 ? 
        dst.m_length * sizeof(ValueType) : 
        dst.m_stride * dst.m_count * sizeof(ValueType );

    MemoryManager< Cuda >::
      copy_to_device_from_device( dst.m_ptr_on_device ,
                                  src.m_ptr_on_device ,
                                  alloc_size );
  }

  static void deep_copy( const MultiVector< ValueType , Cuda > & dst ,
                         const MultiVector< ValueType , Cuda > & src ,
                         const size_t length )
  {
    MemoryManager< Cuda >::
      copy_to_device_from_device( dst.m_ptr_on_device ,
                                  src.m_ptr_on_device ,
                                  length * sizeof(ValueType) );
  }
};


template< typename ValueType >
struct Factory< MultiVector< ValueType , Cuda > ,
                MultiVector< ValueType , Host > >
{
  typedef MultiVector< ValueType , Cuda >  dst_type ;
  typedef MultiVector< ValueType , Host >  src_type ;

  static void deep_copy( const dst_type & dst , const src_type & src )
  {
    const Cuda::size_type size = 
      dst.m_count == 1 ? 
        dst.m_length * sizeof(ValueType) : 
        dst.m_stride * dst.m_count * sizeof(ValueType);

    // Require src.m_stride == dst.m_stride
    // Require src.m_count  == dst.m_count 

    MemoryManager< Cuda >::
      copy_to_device_from_host( dst.m_ptr_on_device ,
                                src.m_ptr_on_device ,
                                size );
  }

  static void deep_copy( const dst_type & dst , const src_type & src ,
                         const size_t length )
  {
    MemoryManager< Cuda >::
      copy_to_device_from_host( dst.m_ptr_on_device ,
                                src.m_ptr_on_device ,
                                length * sizeof(ValueType) );
  }
};

template< typename ValueType >
struct Factory< MultiVector< ValueType , Host > ,
                MultiVector< ValueType , Cuda > >
{
  typedef MultiVector< ValueType , Host >  output_type ;
  typedef MultiVector< ValueType , Cuda >  input_type ;

  static void deep_copy( const output_type & output , const input_type & input )
  {
    const Cuda::size_type size =
      output.m_stride * output.m_count * sizeof(ValueType);

    // Require src.m_stride == dst.m_stride
    // Require src.m_count  == dst.m_count 

    MemoryManager< Cuda >::
      copy_to_host_from_device( output.m_ptr_on_device ,
                                input.m_ptr_on_device ,
                                size );
  }

  static void deep_copy( const output_type & output , const input_type & input ,
                         const size_t length )
  {
    MemoryManager< Cuda >::
      copy_to_host_from_device( output.m_ptr_on_device ,
                                input.m_ptr_on_device ,
                                length * sizeof(ValueType) );
  }

  // Mirror creation:
  static inline
  output_type create( const input_type & input )
  {
    return Factory< output_type , void >
             ::create( std::string(),
                       input.m_length , input.m_count , input.m_stride );
  }
};

} // namespace Impl
} // namespace KokkosArray

#endif /* #ifndef KOKKOS_CUDA_MULTIVECTOR_HPP */

