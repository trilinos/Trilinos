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

#ifndef KOKKOS_CUDA_PREFIXSUM_HPP
#define KOKKOS_CUDA_PREFIXSUM_HPP

#include <string>

#include <Cuda/KokkosArray_Cuda_IndexMap.hpp>

#include <KokkosArray_Cuda_macros.hpp>
#include <impl/KokkosArray_PrefixSum_macros.hpp>
#include <KokkosArray_Clear_macros.hpp>

// For the host mirror:

#include <KokkosArray_Host_macros.hpp>
#undef KOKKOS_MACRO_DEVICE
#define KOKKOS_MACRO_DEVICE HostMapped< Cuda >
#include <impl/KokkosArray_PrefixSum_macros.hpp>
#include <KokkosArray_Clear_macros.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

//----------------------------------------------------------------------------

/** \brief  The hostview is identically mapped */

template< typename IntType >
struct Factory< PrefixSum< IntType , Cuda > ,
                PrefixSum< IntType , HostMapped< Cuda > > >
{
  typedef PrefixSum< IntType, Cuda >              output_type ;
  typedef PrefixSum< IntType, HostMapped<Cuda> >  input_type ;

  static inline
  void deep_copy( output_type & output , const input_type & input )
  {
    const size_t size_data = sizeof(IntType)*(output.m_length + 1 );

    MemoryManager< Cuda >::
      copy_to_device_from_host( output.m_data.ptr_on_device(),
                                input.m_data.ptr_on_device(),
                                size_data );
    output.m_sum = input.m_sum ;
  }
};

//----------------------------------------------------------------------------

template< typename IntType >
struct Factory< PrefixSum< IntType , HostMapped< Cuda > > ,
                PrefixSum< IntType , Cuda > >
{
  typedef PrefixSum< IntType, HostMapped<Cuda> > output_type ;
  typedef PrefixSum< IntType, Cuda >             input_type ;

  static void deep_copy( output_type & output , const input_type & input )
  {
    const size_t size_data = sizeof(IntType)*(output.m_length + 1 );

    MemoryManager< Cuda >::
      copy_to_host_from_device( output.m_data.ptr_on_device(),
                                input.m_data.ptr_on_device(),
                                size_data );
    output.m_sum = input.m_sum ;
  }

  static inline
  output_type create( const input_type & input )
  {
    output_type output ;

    output.m_length = input.m_length ;
    output.m_sum    = input.m_sum ;
    output.m_data.allocate( output.m_length + 1 , std::string() );

    deep_copy( output , input );

    return output ;
  }
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray

#endif /* #ifndef KOKKOS_CUDA_PREFIXSUM_HPP */

