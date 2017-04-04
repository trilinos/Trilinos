/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
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
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <Kokkos_Blas1_MV_impl_fill.hpp>
#include <climits>


#ifndef KOKKOSKERNELS_KOKKOSBLAS_IMPL_USE_MEMSET
#  define KOKKOSKERNELS_KOKKOSBLAS_IMPL_USE_MEMSET 1
#endif

#ifdef KOKKOSKERNELS_KOKKOSBLAS_IMPL_USE_MEMSET
#  include <cstring> // for memset (see Kokkos::Serial specialization below)
#endif // KOKKOSKERNELS_KOKKOSBLAS_IMPL_USE_MEMSET

namespace KokkosBlas {
namespace Impl {

#ifdef KOKKOSKERNELS_BUILD_EXECUTION_SPACE_SERIAL
#define KOKKOSBLAS_IMPL_MV_EXEC_SPACE Kokkos::Serial
#define KOKKOSBLAS_IMPL_MV_MEM_SPACE Kokkos::HostSpace
#define KOKKOSBLAS_IMPL_MV_SCALAR double

void
Fill<Kokkos::View<KOKKOSBLAS_IMPL_MV_SCALAR**,
                  Kokkos::LayoutLeft,
                  Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >,
     2>::
fill (const XMV& X, const XMV::non_const_value_type& val)
{
  typedef XMV::size_type size_type;
  const size_type numRows = X.dimension_0 ();
  const size_type numCols = X.dimension_1 ();

#ifdef KOKKOSKERNELS_KOKKOSBLAS_IMPL_USE_MEMSET

  // Don't call one of the special cases (memset or 1-D fill) unless
  // the memory in X is contiguous.
  if (X.capacity () == numRows * numCols) {
    if (val == Kokkos::Details::ArithTraits<KOKKOSBLAS_IMPL_MV_SCALAR>::zero ()) {
      // It might not necessarily be true for ALL Scalar types that
      // memset with 0 does the right thing, but it certainly works
      // here.
      memset (X.ptr_on_device (), 0, numRows * numCols * sizeof (KOKKOSBLAS_IMPL_MV_SCALAR));
    }
    else {
      typedef Kokkos::View<KOKKOSBLAS_IMPL_MV_SCALAR*,
        Kokkos::LayoutLeft,
        Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
        Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV1D;

      XV1D X1D (X.ptr_on_device (), X.capacity ());
      Kokkos::deep_copy(X1D, val);

      // mfh 14 Apr 2015: This didn't actually help performance over
      // using ViewFill on the 1-D View.  The key thing is using the
      // right ViewFill specialization: 1-D is faster (.18 s / 26 s,
      // for example).

      // if (numRows < static_cast<size_type> (INT_MAX)) {
      //   V_Fill_Invoke<XV1D, int> (X1D, val);
      // }
      // else {
      //   V_Fill_Invoke<XV1D, size_type> (X1D, val);
      // }
    }
    return;
  }
  //
  // // mfh 14 Apr 2015: Not actually faster for me than the code below
  //
  //   Kokkos::Impl::ViewFill<XMV> (X, val);
  //

#endif // KOKKOSKERNELS_KOKKOSBLAS_IMPL_USE_MEMSET

  // The first condition helps avoid overflow with the
  // multiplication in the second condition.
  if (numRows < static_cast<size_type> (INT_MAX) &&
      numRows * numCols < static_cast<size_type> (INT_MAX)) {
    MV_Fill_Invoke<XMV, int> (X, val);
  }
  else {
    MV_Fill_Invoke<XMV, size_type> (X, val);
  }
}

#undef KOKKOSBLAS_IMPL_MV_EXEC_SPACE
#undef KOKKOSBLAS_IMPL_MV_MEM_SPACE
#undef KOKKOSBLAS_IMPL_MV_SCALAR
#endif // KOKKOSKERNELS_BUILD_EXECUTION_SPACE_SERIAL


#ifdef KOKKOSKERNELS_BUILD_EXECUTION_SPACE_OPENMP

  KOKKOSBLAS_IMPL_MV_FILL_RANK2_DEF( double, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace )

#endif // KOKKOSKERNELS_BUILD_EXECUTION_SPACE_OPENMP


#ifdef KOKKOSKERNELS_BUILD_EXECUTION_SPACE_PTHREAD

  KOKKOSBLAS_IMPL_MV_FILL_RANK2_DEF( double, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace )

#endif // KOKKOSKERNELS_BUILD_EXECUTION_SPACE_PTHREAD


#ifdef KOKKOSKERNELS_BUILD_EXECUTION_SPACE_CUDA

  KOKKOSBLAS_IMPL_MV_FILL_RANK2_DEF( double, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace )

#endif // KOKKOSKERNELS_BUILD_EXECUTION_SPACE_CUDA

} // namespace Impl
} // namespace KokkosBlas

