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

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/** \brief No-op because memory is initialized at allocation */
template< typename ValueType >
class Initialize< MDArray< ValueType , Cuda > > {
public:
  static void run( const MDArray< ValueType , Cuda > & ) {}
};

template< typename ValueType >
class DeepCopy< MDArray< ValueType , Cuda > ,
                MDArray< ValueType , Cuda > > {
public:
  static void run( const MDArray< ValueType , Cuda > & dst ,
                   const MDArray< ValueType , Cuda > & src )
  {
    const size_t size = dst.m_map.allocation_size() * sizeof(ValueType);

    MemoryManager< Cuda >::
      copy_to_device_from_device( dst.m_memory.ptr_on_device(),
                                  src.m_memory.ptr_on_device(),
                                  size );
  }
};

/** \brief  The hostview is identically mapped */
template< typename ValueType >
class DeepCopy< MDArray< ValueType , Cuda > ,
                typename MDArray< ValueType , Cuda >::HostMirror > {
public:
  typedef MDArray< ValueType , Cuda >                    dst_type ;
  typedef typename MDArray< ValueType , Cuda >::HostMirror src_type ;

  static void run( const dst_type & dst , const src_type & src )
  {
    const size_t size = dst.m_map.allocation_size() * sizeof(ValueType);

    MemoryManager< Cuda >::
      copy_to_device_from_host( dst.m_memory.ptr_on_device(),
                                src.m_memory.ptr_on_device(),
                                 size );
  }
};

template< typename ValueType >
class DeepCopy< typename MDArray< ValueType , Cuda >::HostMirror ,
                MDArray< ValueType , Cuda > > {
public:
  typedef typename MDArray< ValueType , Cuda >::HostMirror dst_type ;
  typedef MDArray< ValueType , Cuda >                    src_type ;

  static void run( const dst_type & dst , const src_type & src )
  {
    const size_t size = src.m_map.allocation_size() * sizeof(ValueType);

    MemoryManager< Cuda >::
      copy_to_host_from_device( dst.m_memory.ptr_on_device(),
                                src.m_memory.ptr_on_device(),
                                size );
  }
};

} // namespace Impl
} // namespace Kokkos

