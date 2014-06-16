/*
//@HEADER
// ************************************************************************
//
//                             Kokkos
//         Manycore Performance-Portable Multidimensional Arrays
//
//              Copyright (2012) Sandia Corporation
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
// Questions?  Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_TAGS_HPP
#define KOKKOS_TAGS_HPP

#include <impl/Kokkos_Traits.hpp>

namespace Kokkos {
#if   defined ( KOKKOS_HAVE_CUDA )
class Cuda ;
#endif
#if   defined ( KOKKOS_HAVE_OPENMP )
class OpenMP ;
#endif
#if   defined ( KOKKOS_HAVE_PTHREAD )
class Threads ;
#endif
class Serial ;
} // namespace Kokkos


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/// Define DefaultDeviceType
/// The DefaultDeviceType can either be set externally via
/// KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_**** or is chosen according to the active
/// Device types with following priority: Cuda,OpenMP,Threads,Serial

namespace Kokkos {
namespace Impl {
  #if   defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_CUDA ) && \
       !defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_OPENMP ) && \
       !defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_THREADS ) && \
       !defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_SERIAL )
    typedef Cuda DefaultDeviceType;
  #elif defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_OPENMP ) && \
       !defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_CUDA ) && \
       !defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_THREADS ) && \
       !defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_SERIAL )
    typedef OpenMP DefaultDeviceType;
  #elif defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_THREADS ) && \
       !defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_OPENMP ) && \
       !defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_CUDA ) && \
       !defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_SERIAL )
    typedef Threads DefaultDeviceType;
  #elif defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_SERIAL ) && \
       !defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_OPENMP ) && \
       !defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_THREADS ) && \
       !defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_CUDA )
    typedef Serial DefaultDeviceType;
  #else
    #if   defined ( KOKKOS_HAVE_CUDA )
      typedef Kokkos::Cuda DefaultDeviceType;
      #define KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_CUDA
    #elif defined ( KOKKOS_HAVE_OPENMP )
      typedef Kokkos::OpenMP DefaultDeviceType;
      #define KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_OPENMP
    #elif defined ( KOKKOS_HAVE_PTHREAD )
      typedef Kokkos::Threads DefaultDeviceType;
      #define KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_THREADS
    #else
      typedef Kokkos::Serial DefaultDeviceType;
      #define KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_SERIAL
    #endif
  #endif
}
}

namespace Kokkos {
namespace Impl {

struct LayoutTag {};
struct DeviceTag {};
struct MemoryTraitsTag {};
struct ViewTag {};
}

template< class C , class Enable = void >
struct is_layout : public Impl::false_type {};

template<class C>
struct is_layout<C,typename Impl::enable_if< ! Impl::is_same<typename C::kokkos_tag,int>::value>::type > {
  enum {value=bool(Impl::is_same<Impl::LayoutTag,typename C::kokkos_tag>::value)};
};

template< class C , class Enable = void >
struct is_device : public Impl::false_type {};

template<class C>
struct is_device<C,typename Impl::enable_if< ! Impl::is_same<typename C::kokkos_tag,int>::value>::type > {
  enum {value=bool(Impl::is_same<Impl::DeviceTag,typename C::kokkos_tag>::value)};
};

template< class C , class Enable = void >
struct is_memorytraits : public Impl::false_type {};

template<class C>
struct is_memorytraits<C,typename Impl::enable_if< ! Impl::is_same<typename C::kokkos_tag,int>::value>::type > {
  enum {value=bool(Impl::is_same<Impl::MemoryTraitsTag,typename C::kokkos_tag>::value)};
};

template< class C , class Enable = void >
struct is_view : public Impl::false_type {};

template<class C>
struct is_view<C,typename Impl::enable_if< ! Impl::is_same<typename C::kokkos_tag,int>::value>::type > {
  enum {value=bool(Impl::is_same<Impl::ViewTag,typename C::kokkos_tag>::value)};
};

}

#endif
