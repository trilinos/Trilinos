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
#ifndef KOKKOS_BLAS1_MV_IMPL_FILL_HPP_
#define KOKKOS_BLAS1_MV_IMPL_FILL_HPP_

#include <TpetraKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_InnerProductSpaceTraits.hpp>

namespace KokkosBlas {
namespace Impl {

//
// fill
//

template<class XMV, class SizeType = typename XMV::size_type>
struct MV_FillFunctor {
  typedef typename XMV::execution_space             execution_space;
  typedef SizeType                                        size_type;
  typedef typename XMV::non_const_value_type            xvalue_type;

  const size_type numCols_;
  const xvalue_type val_;
  XMV X_;

  MV_FillFunctor (const XMV& X, const xvalue_type& val) :
    numCols_ (X.dimension_1 ()), val_ (val), X_ (X)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
    for (size_type j = 0; j < numCols_; ++j) {
      X_(i,j) = val_;
    }
  }
};

template<class XV, class SizeType = typename XV::size_type>
struct V_FillFunctor {
  typedef typename XV::execution_space            execution_space;
  typedef SizeType                                      size_type;
  typedef typename XV::non_const_value_type           xvalue_type;

  const xvalue_type val_;
  XV x_;

  V_FillFunctor (const XV& x, const xvalue_type& val) :
    val_ (val), x_ (x)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
    x_(i) = val_;
  }
};

//! Implementation of KokkosBlas::fill for multivectors.
template<class XT, class XL, class XD, class XM, class XS>
struct Fill_MV {
  typedef Kokkos::View<XT,XL,XD,XM,XS> XMV;
  typedef typename XMV::execution_space execution_space;
  typedef typename XMV::size_type size_type;

  static void fill (const XMV& X, const typename XMV::non_const_value_type& val)
  {
    const size_type numRows = X.dimension_0 ();
    const size_type numCols = X.dimension_1 ();

    // The first condition helps avoid overflow with the
    // multiplication in the second condition.
    if (numRows < static_cast<size_type> (INT_MAX) &&
        numRows * numCols < static_cast<size_type> (INT_MAX)) {
      Kokkos::RangePolicy<execution_space, int> policy (0, numRows);
      MV_FillFunctor<XMV, int> op (X, val);
      Kokkos::parallel_for (policy, op);
    }
    else {
      Kokkos::RangePolicy<execution_space, size_type> policy (0, numRows);
      MV_FillFunctor<XMV, size_type> op (X, val);
      Kokkos::parallel_for (policy, op);
    }
  }
};

#ifdef KOKKOS_HAVE_SERIAL
template<>
struct Fill_MV<double**, Kokkos::LayoutLeft, Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>, Kokkos::Impl::ViewDefault> {
  typedef Kokkos::View<double**, Kokkos::LayoutLeft, Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>, Kokkos::Impl::ViewDefault> XMV;
  typedef Kokkos::Serial execution_space;
  typedef XMV::size_type size_type;

  static void fill (const XMV& X, const double& val);
};
#endif // KOKKOS_HAVE_SERIAL

#ifdef KOKKOS_HAVE_OPENMP
template<>
struct Fill_MV<double**, Kokkos::LayoutLeft, Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace>,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>, Kokkos::Impl::ViewDefault> {
  typedef Kokkos::View<double**, Kokkos::LayoutLeft, Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>, Kokkos::Impl::ViewDefault> XMV;
  typedef Kokkos::OpenMP execution_space;
  typedef XMV::size_type size_type;

  static void fill (const XMV& X, const double& val);
};
#endif // KOKKOS_HAVE_OPENMP

#ifdef KOKKOS_HAVE_PTHREAD
template<>
struct Fill_MV<double**, Kokkos::LayoutLeft, Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace>,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>, Kokkos::Impl::ViewDefault> {
  typedef Kokkos::View<double**, Kokkos::LayoutLeft, Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>, Kokkos::Impl::ViewDefault> XMV;
  typedef Kokkos::Threads execution_space;
  typedef XMV::size_type size_type;

  static void fill (const XMV& X, const double& val);
};
#endif // KOKKOS_HAVE_PTHREAD

#ifdef KOKKOS_HAVE_CUDA
template<>
struct Fill_MV<double**, Kokkos::LayoutLeft, Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>, Kokkos::Impl::ViewDefault> {
  typedef Kokkos::View<double**, Kokkos::LayoutLeft, Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>, Kokkos::Impl::ViewDefault> XMV;
  typedef Kokkos::Cuda execution_space;
  typedef XMV::size_type size_type;

  static void fill (const XMV& X, const double& val);
};
#endif // KOKKOS_HAVE_CUDA

#ifdef KOKKOS_HAVE_CUDA
template<>
struct Fill_MV<double**, Kokkos::LayoutLeft, Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace>,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>, Kokkos::Impl::ViewDefault> {
  typedef Kokkos::View<double**, Kokkos::LayoutLeft, Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>, Kokkos::Impl::ViewDefault> XMV;
  typedef Kokkos::Cuda execution_space;
  typedef XMV::size_type size_type;

  static void fill (const XMV& X, const double& val);
};
#endif // KOKKOS_HAVE_CUDA

} // namespace Impl
} // namespace KokkosBlas

#endif // KOKKOS_BLAS1_MV_IMPL_FILL_HPP_
