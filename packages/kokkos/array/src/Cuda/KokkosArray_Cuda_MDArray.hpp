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

#ifndef KOKKOS_CUDA_MDARRAY_HPP
#define KOKKOS_CUDA_MDARRAY_HPP

#include <string>

#include <Cuda/KokkosArray_Cuda_IndexMap.hpp>

#include <KokkosArray_Cuda_macros.hpp>
#include <impl/KokkosArray_MDArray_macros.hpp>
#include <KokkosArray_Clear_macros.hpp>

// For the host-mapped memory view:

#include <KokkosArray_Host_macros.hpp>
#undef KOKKOS_MACRO_DEVICE
#define KOKKOS_MACRO_DEVICE HostMapped< Cuda >
#include <impl/KokkosArray_MDArray_macros.hpp>
#include <KokkosArray_Clear_macros.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template< typename ValueType >
struct Factory< MDArray< ValueType , Cuda > , void >
{
  typedef MDArray< ValueType , Cuda > output_type ;

  static output_type create( const std::string & label ,
                             size_t nP , size_t n1 , size_t n2 , size_t n3 ,
                             size_t n4 , size_t n5 , size_t n6 , size_t n7 )
  {
    output_type array ;

    array.m_map.template assign< ValueType >(nP,n1,n2,n3,n4,n5,n6,n7);
    array.m_memory.allocate( array.m_map.allocation_size() , label );

    // Cuda 'allocate' initializes to zero 

    return array ;
  }
};

template< typename ValueType >
struct Factory< MDArray< ValueType , Cuda > ,
                MDArray< ValueType , Cuda > >
{
  typedef MDArray< ValueType , Cuda > output_type ;

  static inline
  void deep_copy( const output_type & dst ,
                  const output_type & src )
  {
    const size_t size = dst.m_map.allocation_size() * sizeof(ValueType);

    MemoryManager< Cuda >::
      copy_to_device_from_device( dst.m_memory.ptr_on_device(),
                                  src.m_memory.ptr_on_device(),
                                  size );
  }

  static inline
  output_type create( const output_type & input )
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

/** \brief  The hostview is identically mapped */
template< typename ValueType >
struct Factory< MDArray< ValueType , Cuda > ,
                MDArray< ValueType , HostMapped< Cuda > > >
{
  typedef MDArray< ValueType , Cuda >                output_type ;
  typedef MDArray< ValueType , HostMapped< Cuda > >  input_type ;

  static inline
  void deep_copy( const output_type & output , const input_type & input )
  {
    const size_t size = output.m_map.allocation_size() * sizeof(ValueType);

    MemoryManager< Cuda >::
      copy_to_device_from_host( output.m_memory.ptr_on_device(),
                                input.m_memory.ptr_on_device(),
                                size );
  }
};

template< typename ValueType >
struct Factory< MDArray< ValueType , HostMapped< Cuda > > ,
                MDArray< ValueType , Cuda > >
{
  typedef MDArray< ValueType , HostMapped< Cuda > > output_type ;
  typedef MDArray< ValueType , Cuda >               input_type ;

  static void deep_copy( const output_type & output , const input_type & input )
  {
    const size_t size = input.m_map.allocation_size() * sizeof(ValueType);

    MemoryManager< Cuda >::
      copy_to_host_from_device( output.m_memory.ptr_on_device(),
                                input.m_memory.ptr_on_device(),
                                size );
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

} // namespace Impl
} // namespace KokkosArray

#endif /* #ifndef KOKKOS_CUDA_MDARRAY_HPP */

