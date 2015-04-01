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

#include <Kokkos_Blas1_impl.hpp>
#include <climits>

namespace KokkosBlas {
namespace Impl {

#ifdef KOKKOS_HAVE_SERIAL
#define IMPL_EXEC_BLAS Kokkos::Serial
#define IMPL_MEM_BLAS Kokkos::HostSpace
double
Dot<const double*,Kokkos::LayoutLeft,Kokkos::Device<IMPL_EXEC_BLAS,IMPL_MEM_BLAS>,
           Kokkos::MemoryTraits<Kokkos::Unmanaged>,Kokkos::Impl::ViewDefault,
           const double*,Kokkos::LayoutLeft,Kokkos::Device<IMPL_EXEC_BLAS,IMPL_MEM_BLAS>,
           Kokkos::MemoryTraits<Kokkos::Unmanaged>,Kokkos::Impl::ViewDefault>::
  dot (const XVector& x, const XVector& y) {
    double result;

    // With Intel 15 in a serial test with 100000 elements for 1000 trials
    // using int instead of size_t is 2x faster
    if(x.dimension_0() < typename XVector::size_type(INT_MAX)) {
      Kokkos::RangePolicy<typename XVector::execution_space, int> policy(0,x.dimension_0());
      Kokkos::parallel_reduce(policy,
        DotFunctor<XVector,XVector,int>(x,y),
        result);
      return result;
    } else {
      Kokkos::RangePolicy<typename XVector::execution_space> policy(0,x.dimension_0());
      Kokkos::parallel_reduce(policy,
        DotFunctor<XVector,XVector>(x,y),
        result);
      return result;
    }
  }
#undef IMPL_EXEC_BLAS
#undef IMPL_MEM_BLAS
#endif

#ifdef KOKKOS_HAVE_OPENMP
#define IMPL_EXEC_BLAS Kokkos::OpenMP
#define IMPL_MEM_BLAS Kokkos::HostSpace
double
Dot<const double*,Kokkos::LayoutLeft,Kokkos::Device<IMPL_EXEC_BLAS,IMPL_MEM_BLAS>,
           Kokkos::MemoryTraits<Kokkos::Unmanaged>,Kokkos::Impl::ViewDefault,
           const double*,Kokkos::LayoutLeft,Kokkos::Device<IMPL_EXEC_BLAS,IMPL_MEM_BLAS>,
           Kokkos::MemoryTraits<Kokkos::Unmanaged>,Kokkos::Impl::ViewDefault>::
  dot (const XVector& x, const XVector& y) {
    double result;

    // With Intel 15 in a serial test with 100000 elements for 1000 trials
    // using int instead of size_t is 2x faster
    if(x.dimension_0() < typename XVector::size_type(INT_MAX)) {
      Kokkos::RangePolicy<typename XVector::execution_space, int> policy(0,x.dimension_0());
      Kokkos::parallel_reduce(policy,
        DotFunctor<XVector,XVector,int>(x,y),
        result);
      return result;
    } else {
      Kokkos::RangePolicy<typename XVector::execution_space> policy(0,x.dimension_0());
      Kokkos::parallel_reduce(policy,
        DotFunctor<XVector,XVector>(x,y),
        result);
      return result;
    }
  }
#undef IMPL_EXEC_BLAS
#undef IMPL_MEM_BLAS
#endif

#ifdef KOKKOS_HAVE_PTHREAD
#define IMPL_EXEC_BLAS Kokkos::Threads
#define IMPL_MEM_BLAS Kokkos::HostSpace
double
Dot<const double*,Kokkos::LayoutLeft,Kokkos::Device<IMPL_EXEC_BLAS,IMPL_MEM_BLAS>,
           Kokkos::MemoryTraits<Kokkos::Unmanaged>,Kokkos::Impl::ViewDefault,
           const double*,Kokkos::LayoutLeft,Kokkos::Device<IMPL_EXEC_BLAS,IMPL_MEM_BLAS>,
           Kokkos::MemoryTraits<Kokkos::Unmanaged>,Kokkos::Impl::ViewDefault>::
  dot (const XVector& x, const XVector& y) {
    double result;

    // With Intel 15 in a serial test with 100000 elements for 1000 trials
    // using int instead of size_t is 2x faster
    if(x.dimension_0() < typename XVector::size_type(INT_MAX)) {
      Kokkos::RangePolicy<typename XVector::execution_space, int> policy(0,x.dimension_0());
      Kokkos::parallel_reduce(policy,
        DotFunctor<XVector,XVector,int>(x,y),
        result);
      return result;
    } else {
      Kokkos::RangePolicy<typename XVector::execution_space> policy(0,x.dimension_0());
      Kokkos::parallel_reduce(policy,
        DotFunctor<XVector,XVector>(x,y),
        result);
      return result;
    }
  }
#undef IMPL_EXEC_BLAS
#undef IMPL_MEM_BLAS
#endif

#ifdef KOKKOS_HAVE_CUDA
#define IMPL_EXEC_BLAS Kokkos::Cuda
#define IMPL_MEM_BLAS Kokkos::CudaSpace
double
Dot<const double*,Kokkos::LayoutLeft,Kokkos::Device<IMPL_EXEC_BLAS,IMPL_MEM_BLAS>,
           Kokkos::MemoryTraits<Kokkos::Unmanaged>,Kokkos::Impl::ViewDefault,
           const double*,Kokkos::LayoutLeft,Kokkos::Device<IMPL_EXEC_BLAS,IMPL_MEM_BLAS>,
           Kokkos::MemoryTraits<Kokkos::Unmanaged>,Kokkos::Impl::ViewDefault>::
  dot (const XVector& x, const XVector& y) {
    double result;

    // With Intel 15 in a serial test with 100000 elements for 1000 trials
    // using int instead of size_t is 2x faster
    if(x.dimension_0() < typename XVector::size_type(INT_MAX)) {
      Kokkos::RangePolicy<typename XVector::execution_space, int> policy(0,x.dimension_0());
      Kokkos::parallel_reduce(policy,
        DotFunctor<XVector,XVector,int>(x,y),
        result);
      return result;
    } else {
      Kokkos::RangePolicy<typename XVector::execution_space> policy(0,x.dimension_0());
      Kokkos::parallel_reduce(policy,
        DotFunctor<XVector,XVector>(x,y),
        result);
      return result;
    }
  }
#undef IMPL_EXEC_BLAS
#undef IMPL_MEM_BLAS

#define IMPL_EXEC_BLAS Kokkos::Cuda
#define IMPL_MEM_BLAS Kokkos::CudaUVMSpace
double
Dot<const double*,Kokkos::LayoutLeft,Kokkos::Device<IMPL_EXEC_BLAS,IMPL_MEM_BLAS>,
           Kokkos::MemoryTraits<Kokkos::Unmanaged>,Kokkos::Impl::ViewDefault,
           const double*,Kokkos::LayoutLeft,Kokkos::Device<IMPL_EXEC_BLAS,IMPL_MEM_BLAS>,
           Kokkos::MemoryTraits<Kokkos::Unmanaged>,Kokkos::Impl::ViewDefault>::
  dot (const XVector& x, const XVector& y) {
    double result;

    // With Intel 15 in a serial test with 100000 elements for 1000 trials
    // using int instead of size_t is 2x faster
    if(x.dimension_0() < typename XVector::size_type(INT_MAX)) {
      Kokkos::RangePolicy<typename XVector::execution_space, int> policy(0,x.dimension_0());
      Kokkos::parallel_reduce(policy,
        DotFunctor<XVector,XVector,int>(x,y),
        result);
      return result;
    } else {
      Kokkos::RangePolicy<typename XVector::execution_space> policy(0,x.dimension_0());
      Kokkos::parallel_reduce(policy,
        DotFunctor<XVector,XVector>(x,y),
        result);
      return result;
    }
  }
#undef IMPL_EXEC_BLAS
#undef IMPL_MEM_BLAS
#endif

}
}
