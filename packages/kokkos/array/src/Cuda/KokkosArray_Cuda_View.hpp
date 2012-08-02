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

#ifndef KOKKOS_CUDA_VIEW_HPP
#define KOKKOS_CUDA_VIEW_HPP

#include <KokkosArray_Host.hpp>
#include <KokkosArray_View.hpp>

#include <KokkosArray_Cuda_macros.hpp>
#include <impl/KokkosArray_ViewOperLeft_macros.hpp>
#include <impl/KokkosArray_ViewOperRight_macros.hpp>
#include <impl/KokkosArray_View_macros.hpp>
#include <KokkosArray_Clear_macros.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {

template< class DataType , class LayoutType >
void View< DataType , LayoutType , Cuda >::create(
  const std::string & label ,
  const typename View< DataType , LayoutType , Cuda >::shape_type shape )
{
  const size_t count = Impl::allocation_count( shape );

  oper_type::m_shape = shape ;
  oper_type::m_ptr_on_device = (value_type *)
    memory_space::allocate( label , 
                            typeid(value_type) , 
                            sizeof(value_type) , 
                            count );
}

}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template< class DeviceDst , class DeviceSrc >
struct CudaDeepCopy ;

template<>
struct CudaDeepCopy<Cuda,Cuda>
{
  template< typename ValueType >
  CudaDeepCopy( ValueType * dst ,
                ValueType * src ,
                size_t count )
  {
    Cuda::memory_space
        ::copy_to_device_from_device( dst , src , sizeof(ValueType) * count );
  }
};

template<>
struct CudaDeepCopy<Cuda,Host>
{
  template< typename ValueType >
  CudaDeepCopy( ValueType * dst ,
                ValueType * src ,
                size_t count )
  {
    Cuda::memory_space
        ::copy_to_device_from_host( dst , src , sizeof(ValueType) * count );
  }
};

template<>
struct CudaDeepCopy<Host,Cuda>
{
  template< typename ValueType >
  CudaDeepCopy( ValueType * dst ,
                ValueType * src ,
                size_t count )
  {
    Cuda::memory_space
        ::copy_to_host_from_device( dst , src , sizeof(ValueType) * count );
  }
};

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_CUDA_VIEW_HPP */

