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
#ifndef KOKKOS_BLAS1_IMPL_HPP_
#define KOKKOS_BLAS1_IMPL_HPP_

#include <TpetraKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_InnerProductSpaceTraits.hpp>

namespace KokkosBlas {
namespace Impl {

template<class XVector, class YVector, typename SizeType = typename XVector::size_type>
struct DotFunctor
{
  typedef typename XVector::execution_space        execution_space;
  typedef SizeType                                 size_type;
  typedef typename XVector::non_const_value_type   xvalue_type;
  typedef Kokkos::Details::InnerProductSpaceTraits<xvalue_type>  IPT;
  typedef typename IPT::dot_type                   value_type;
  XVector  m_x ;
  YVector  m_y ;

  //--------------------------------------------------------------------------
  DotFunctor(const XVector& x,const YVector& y):
    m_x(x),m_y(y)
  { }

  KOKKOS_FORCEINLINE_FUNCTION
  void operator()( const size_type &i, value_type &sum ) const
  {
    sum += IPT::dot( m_x(i), m_y(i) );  // m_x(i) * m_y(i)
  }

  KOKKOS_INLINE_FUNCTION
  void init (volatile value_type& update) const
  {
    update = Kokkos::Details::ArithTraits<value_type>::zero ();
  }

  KOKKOS_INLINE_FUNCTION
  void join( volatile value_type &update ,
             const volatile value_type &source ) const
  {
    update += source ;
  }
};

template<class XT, class XL, class XD, class XM, class XS,
         class YT, class YL, class YD, class YM, class YS>
struct Dot {
  static typename Kokkos::Details::InnerProductSpaceTraits<typename Kokkos::View<XT,XL,XD,XM,XS>::non_const_value_type>::dot_type
  dot (const Kokkos::View<XT,XL,XD,XM,XS>& x, const Kokkos::View<YT,YL,YD,YM,YS>& y) {
    typedef Kokkos::View<XT,XL,XD,XM,XS> XVector;
    typedef Kokkos::View<YT,YL,YD,YM,YS> YVector;

    typename Kokkos::Details::InnerProductSpaceTraits<typename XVector::non_const_value_type>::dot_type
      result;

    Kokkos::parallel_reduce(x.dimension_0(),
        DotFunctor<XVector,YVector>(x,y),
        result);
    return result;
  }
};

#ifdef KOKKOS_HAVE_SERIAL
template<>
struct Dot<const double*,Kokkos::LayoutLeft,Kokkos::Device<Kokkos::Serial,Kokkos::HostSpace>,
           Kokkos::MemoryTraits<Kokkos::Unmanaged>,Kokkos::Impl::ViewDefault,
           const double*,Kokkos::LayoutLeft,Kokkos::Device<Kokkos::Serial,Kokkos::HostSpace>,
           Kokkos::MemoryTraits<Kokkos::Unmanaged>,Kokkos::Impl::ViewDefault>  {
  typedef const double* XT;
  typedef Kokkos::LayoutLeft XL;
  typedef Kokkos::Device<Kokkos::Serial,Kokkos::HostSpace> XD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> XM;
  typedef Kokkos::Impl::ViewDefault XS;
  typedef Kokkos::View<XT,XL,XD,XM,XS> XVector;

  static double
  dot (const XVector& x, const XVector& y);
};
#endif

#ifdef KOKKOS_HAVE_OPENMP
template<>
struct Dot<const double*,Kokkos::LayoutLeft,Kokkos::Device<Kokkos::OpenMP,Kokkos::HostSpace>,
           Kokkos::MemoryTraits<Kokkos::Unmanaged>,Kokkos::Impl::ViewDefault,
           const double*,Kokkos::LayoutLeft,Kokkos::Device<Kokkos::OpenMP,Kokkos::HostSpace>,
           Kokkos::MemoryTraits<Kokkos::Unmanaged>,Kokkos::Impl::ViewDefault>  {
  typedef const double* XT;
  typedef Kokkos::LayoutLeft XL;
  typedef Kokkos::Device<Kokkos::OpenMP,Kokkos::HostSpace> XD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> XM;
  typedef Kokkos::Impl::ViewDefault XS;
  typedef Kokkos::View<XT,XL,XD,XM,XS> XVector;

  static double
  dot (const XVector& x, const XVector& y);
};
#endif

#ifdef KOKKOS_HAVE_PTHREAD
template<>
struct Dot<const double*,Kokkos::LayoutLeft,Kokkos::Device<Kokkos::Threads,Kokkos::HostSpace>,
           Kokkos::MemoryTraits<Kokkos::Unmanaged>,Kokkos::Impl::ViewDefault,
           const double*,Kokkos::LayoutLeft,Kokkos::Device<Kokkos::Threads,Kokkos::HostSpace>,
           Kokkos::MemoryTraits<Kokkos::Unmanaged>,Kokkos::Impl::ViewDefault>  {
  typedef const double* XT;
  typedef Kokkos::LayoutLeft XL;
  typedef Kokkos::Device<Kokkos::Threads,Kokkos::HostSpace> XD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> XM;
  typedef Kokkos::Impl::ViewDefault XS;
  typedef Kokkos::View<XT,XL,XD,XM,XS> XVector;

  static double
  dot (const XVector& x, const XVector& y);
};
#endif

#ifdef KOKKOS_HAVE_CUDA
template<>
struct Dot<const double*,Kokkos::LayoutLeft,Kokkos::Device<Kokkos::Cuda,Kokkos::CudaSpace>,
           Kokkos::MemoryTraits<Kokkos::Unmanaged>,Kokkos::Impl::ViewDefault,
           const double*,Kokkos::LayoutLeft,Kokkos::Device<Kokkos::Cuda,Kokkos::CudaSpace>,
           Kokkos::MemoryTraits<Kokkos::Unmanaged>,Kokkos::Impl::ViewDefault>  {
  typedef const double* XT;
  typedef Kokkos::LayoutLeft XL;
  typedef Kokkos::Device<Kokkos::Cuda,Kokkos::CudaSpace> XD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> XM;
  typedef Kokkos::Impl::ViewDefault XS;
  typedef Kokkos::View<XT,XL,XD,XM,XS> XVector;

  static double
  dot (const XVector& x, const XVector& y);
};

template<>
struct Dot<const double*,Kokkos::LayoutLeft,Kokkos::Device<Kokkos::Cuda,Kokkos::CudaUVMSpace>,
           Kokkos::MemoryTraits<Kokkos::Unmanaged>,Kokkos::Impl::ViewDefault,
           const double*,Kokkos::LayoutLeft,Kokkos::Device<Kokkos::Cuda,Kokkos::CudaUVMSpace>,
           Kokkos::MemoryTraits<Kokkos::Unmanaged>,Kokkos::Impl::ViewDefault>  {
  typedef const double* XT;
  typedef Kokkos::LayoutLeft XL;
  typedef Kokkos::Device<Kokkos::Cuda,Kokkos::CudaUVMSpace> XD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> XM;
  typedef Kokkos::Impl::ViewDefault XS;
  typedef Kokkos::View<XT,XL,XD,XM,XS> XVector;

  static double
  dot (const XVector& x, const XVector& y);
};
#endif


}
}

#endif // KOKKOS_BLAS1_IMPL_HPP_
