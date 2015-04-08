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

#include <Kokkos_Blas1_MV_impl.hpp>
#include <climits>

namespace KokkosBlas {
namespace Impl {

#ifdef KOKKOS_HAVE_SERIAL
#define KOKKOSBLAS_IMPL_MV_EXEC_SPACE Kokkos::Serial
#define KOKKOSBLAS_IMPL_MV_MEM_SPACE Kokkos::HostSpace

void
Dot_MV<double*,
       KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
       Kokkos::Impl::ViewDefault,
       const double**,
       Kokkos::LayoutLeft,
       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
       Kokkos::Impl::ViewDefault,
       const double**,
       Kokkos::LayoutLeft,
       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
       Kokkos::Impl::ViewDefault>::
dot (const RV& r, const XMV& X, const XMV& Y)
{
  const size_type numRows = X.dimension_0 ();
  const size_type numVecs = X.dimension_1 ();

  // TODO (mfh 06 Apr 2015) Figure out the optimal max unroll length.
  // It depends on the number of concurrent memory streams that the
  // processor supports, not so much on the SIMD length.

  if (numVecs > 16) {
    MV_Dot_Right_FunctorVector<RV, XMV, YMV> op (r, X, Y);
    Kokkos::parallel_reduce (numRows, op);
  }
  else {
    switch (numVecs) {
    case 16: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 16> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 15: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 15> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 14: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 14> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 13: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 13> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 12: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 12> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 11: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 11> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 10: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 10> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 9: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 9> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 8: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 8> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 7: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 7> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 6: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 6> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 5: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 5> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 4: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 4> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 3: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 3> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 2: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 2> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 1: {
      // RV needs to turn 0-D, and XMV and YMV need to turn 1-D.
      using Kokkos::ALL;
      using Kokkos::subview;
      // RV needs to turn 0-D, and XMV and YMV need to turn 1-D.
      typedef Kokkos::View<RV::value_type, RV::array_layout, RV::device_type, RV::memory_traits, RV::specialize> RV0D;
      typedef Kokkos::View<XMV::value_type*, XMV::array_layout, XMV::device_type, XMV::memory_traits, XMV::specialize> XMV1D;
      typedef Kokkos::View<YMV::value_type*, YMV::array_layout, YMV::device_type, YMV::memory_traits, YMV::specialize> YMV1D;
#if 0
      if (numRows < static_cast<size_type> (INT_MAX) &&
          numRows * numVecs < static_cast<size_type> (INT_MAX)) {
        typedef V_Dot_Functor<RV0D, XMV1D, YMV1D, int> op_type;
        Kokkos::RangePolicy<XMV::execution_space, int> policy (0, numRows);
        op_type op (subview (r, 0), subview (X, ALL (), 0), subview (Y, ALL (), 0));
        Kokkos::parallel_reduce (policy, op);
      }
      else {
        typedef V_Dot_Functor<RV0D, XMV1D, YMV1D, size_type> op_type;
        Kokkos::RangePolicy<XMV::execution_space, size_type> policy (0, numRows);
        op_type op (subview (r, 0), subview (X, ALL (), 0), subview (Y, ALL (), 0));
        Kokkos::parallel_reduce (policy, op);
      }
#else
      RV0D r_0 = subview (r, 0);
      XMV1D X_0 = subview (X, ALL (), 0);
      YMV1D Y_0 = subview (Y, ALL (), 0);

      if (numRows < static_cast<size_type> (INT_MAX) &&
          numRows * numVecs < static_cast<size_type> (INT_MAX)) {
        double sum = 0.0;
        const int lclNumRows = static_cast<int> (numRows);
        for (int i = 0; i < lclNumRows; ++i) {
          sum += X_0(i) * Y_0(i);
        }
        r_0() = sum;
      }
      else {
        double sum = 0.0;
        for (size_type i = 0; i < numRows; ++i) {
          sum += X_0(i) * Y_0(i);
        }
        r_0() = sum;
      }
#endif // 0
      break;
    }
    } // switch
  } // if-else
}

void
Dot_MV<double*,
       KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
       Kokkos::Impl::ViewDefault,
       const double**,
       Kokkos::LayoutLeft,
       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
       Kokkos::Impl::ViewDefault,
       const double**,
       Kokkos::LayoutLeft,
       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
       Kokkos::Impl::ViewDefault>::
dot (const RV& r, const size_t r_col, const XMV& X, const size_t X_col,
     const XMV& Y, const size_t Y_col)
{
  // RV needs to turn 0-D, and XMV and YMV need to turn 1-D.
  using Kokkos::ALL;
  using Kokkos::subview;
  typedef Kokkos::View<RV::value_type, RV::array_layout,
    RV::device_type, RV::memory_traits, RV::specialize> RV0D;
  typedef Kokkos::View<XMV::const_value_type*, XMV::array_layout, XMV::device_type, XMV::memory_traits, XMV::specialize> XMV1D;
  typedef Kokkos::View<YMV::const_value_type*, YMV::array_layout, YMV::device_type, YMV::memory_traits, YMV::specialize> YMV1D;

  const size_type numRows = X.dimension_0 ();
  const size_type numCols = X.dimension_1 ();
  if (numRows < static_cast<size_type> (INT_MAX) &&
      numRows * numCols < static_cast<size_type> (INT_MAX)) {
    typedef V_Dot_Functor<RV0D, XMV1D, YMV1D, int> op_type;
    op_type op (subview (r, r_col), subview (X, ALL (), X_col),
                subview (Y, ALL (), Y_col));
    Kokkos::parallel_reduce (numRows, op);
  }
  else {
    typedef V_Dot_Functor<RV0D, XMV1D, YMV1D, size_type> op_type;
    op_type op (subview (r, r_col), subview (X, ALL (), X_col),
                subview (Y, ALL (), Y_col));
    Kokkos::parallel_reduce (numRows, op);
  }
}


#undef KOKKOSBLAS_IMPL_MV_EXEC_SPACE
#undef KOKKOSBLAS_IMPL_MV_MEM_SPACE
#endif // KOKKOS_HAVE_SERIAL


#ifdef KOKKOS_HAVE_OPENMP
#define KOKKOSBLAS_IMPL_MV_EXEC_SPACE Kokkos::OpenMP
#define KOKKOSBLAS_IMPL_MV_MEM_SPACE Kokkos::HostSpace

void
Dot_MV<double*,
       KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
       Kokkos::Impl::ViewDefault,
       const double**,
       Kokkos::LayoutLeft,
       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
       Kokkos::Impl::ViewDefault,
       const double**,
       Kokkos::LayoutLeft,
       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
       Kokkos::Impl::ViewDefault>::
dot (const RV& r, const XMV& X, const XMV& Y)
{
  const size_type numRows = X.dimension_0 ();
  const size_type numVecs = X.dimension_1 ();

  // TODO (mfh 06 Apr 2015) Figure out the optimal max unroll length.
  // It depends on the number of concurrent memory streams that the
  // processor supports, not so much on the SIMD length.

  if (numVecs > 16) {
    MV_Dot_Right_FunctorVector<RV, XMV, YMV> op (r, X, Y);
    Kokkos::parallel_reduce (numRows, op);
  }
  else {
    switch (numVecs) {
    case 16: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 16> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 15: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 15> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 14: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 14> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 13: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 13> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 12: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 12> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 11: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 11> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 10: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 10> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 9: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 9> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 8: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 8> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 7: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 7> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 6: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 6> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 5: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 5> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 4: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 4> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 3: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 3> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 2: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 2> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 1: {
      // RV needs to turn 0-D, and XMV and YMV need to turn 1-D.
      using Kokkos::ALL;
      using Kokkos::subview;
      // RV needs to turn 0-D, and XMV and YMV need to turn 1-D.
      typedef Kokkos::View<RV::value_type, RV::array_layout, RV::device_type, RV::memory_traits, RV::specialize> RV0D;
      typedef Kokkos::View<XMV::value_type*, XMV::array_layout, XMV::device_type, XMV::memory_traits, XMV::specialize> XMV1D;
      typedef Kokkos::View<YMV::value_type*, YMV::array_layout, YMV::device_type, YMV::memory_traits, YMV::specialize> YMV1D;
      typedef V_Dot_Functor<RV0D, XMV1D, YMV1D> op_type;
      op_type op (subview (r, 0), subview (X, ALL (), 0), subview (Y, ALL (), 0));
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    } // switch
  } // if-else
}

void
Dot_MV<double*,
       KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
       Kokkos::Impl::ViewDefault,
       const double**,
       Kokkos::LayoutLeft,
       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
       Kokkos::Impl::ViewDefault,
       const double**,
       Kokkos::LayoutLeft,
       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
       Kokkos::Impl::ViewDefault>::
dot (const RV& r, const size_t r_col, const XMV& X, const size_t X_col,
     const XMV& Y, const size_t Y_col)
{
  // RV needs to turn 0-D, and XMV and YMV need to turn 1-D.
  using Kokkos::ALL;
  using Kokkos::subview;
  typedef Kokkos::View<RV::value_type, RV::array_layout,
    RV::device_type, RV::memory_traits, RV::specialize> RV0D;
  typedef Kokkos::View<XMV::const_value_type*, XMV::array_layout, XMV::device_type, XMV::memory_traits, XMV::specialize> XMV1D;
  typedef Kokkos::View<YMV::const_value_type*, YMV::array_layout, YMV::device_type, YMV::memory_traits, YMV::specialize> YMV1D;

  const size_type numRows = X.dimension_0 ();
  const size_type numCols = X.dimension_1 ();
  if (numRows < static_cast<size_type> (INT_MAX) &&
      numRows * numCols < static_cast<size_type> (INT_MAX)) {
    typedef V_Dot_Functor<RV0D, XMV1D, YMV1D, int> op_type;
    op_type op (subview (r, r_col), subview (X, ALL (), X_col),
                subview (Y, ALL (), Y_col));
    Kokkos::parallel_reduce (numRows, op);
  }
  else {
    typedef V_Dot_Functor<RV0D, XMV1D, YMV1D, size_type> op_type;
    op_type op (subview (r, r_col), subview (X, ALL (), X_col),
                subview (Y, ALL (), Y_col));
    Kokkos::parallel_reduce (numRows, op);
  }
}

#undef KOKKOSBLAS_IMPL_MV_EXEC_SPACE
#undef KOKKOSBLAS_IMPL_MV_MEM_SPACE
#endif // KOKKOS_HAVE_OPENMP


#ifdef KOKKOS_HAVE_PTHREAD
#define KOKKOSBLAS_IMPL_MV_EXEC_SPACE Kokkos::Threads
#define KOKKOSBLAS_IMPL_MV_MEM_SPACE Kokkos::HostSpace

void
Dot_MV<double*,
       KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
       Kokkos::Impl::ViewDefault,
       const double**,
       Kokkos::LayoutLeft,
       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
       Kokkos::Impl::ViewDefault,
       const double**,
       Kokkos::LayoutLeft,
       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
       Kokkos::Impl::ViewDefault>::
dot (const RV& r, const XMV& X, const XMV& Y)
{
  const size_type numRows = X.dimension_0 ();
  const size_type numVecs = X.dimension_1 ();

  // TODO (mfh 06 Apr 2015) Figure out the optimal max unroll length.
  // It depends on the number of concurrent memory streams that the
  // processor supports, not so much on the SIMD length.

  if (numVecs > 16) {
    MV_Dot_Right_FunctorVector<RV, XMV, YMV> op (r, X, Y);
    Kokkos::parallel_reduce (numRows, op);
  }
  else {
    switch (numVecs) {
    case 16: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 16> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 15: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 15> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 14: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 14> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 13: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 13> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 12: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 12> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 11: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 11> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 10: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 10> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 9: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 9> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 8: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 8> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 7: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 7> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 6: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 6> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 5: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 5> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 4: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 4> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 3: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 3> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 2: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 2> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 1: {
      // RV needs to turn 0-D, and XMV and YMV need to turn 1-D.
      using Kokkos::ALL;
      using Kokkos::subview;
      // RV needs to turn 0-D, and XMV and YMV need to turn 1-D.
      typedef Kokkos::View<RV::value_type, RV::array_layout, RV::device_type, RV::memory_traits, RV::specialize> RV0D;
      typedef Kokkos::View<XMV::value_type*, XMV::array_layout, XMV::device_type, XMV::memory_traits, XMV::specialize> XMV1D;
      typedef Kokkos::View<YMV::value_type*, YMV::array_layout, YMV::device_type, YMV::memory_traits, YMV::specialize> YMV1D;
      typedef V_Dot_Functor<RV0D, XMV1D, YMV1D> op_type;
      op_type op (subview (r, 0), subview (X, ALL (), 0), subview (Y, ALL (), 0));
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    } // switch
  } // if-else
}

void
Dot_MV<double*,
       KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
       Kokkos::Impl::ViewDefault,
       const double**,
       Kokkos::LayoutLeft,
       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
       Kokkos::Impl::ViewDefault,
       const double**,
       Kokkos::LayoutLeft,
       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
       Kokkos::Impl::ViewDefault>::
dot (const RV& r, const size_t r_col, const XMV& X, const size_t X_col,
     const XMV& Y, const size_t Y_col)
{
  // RV needs to turn 0-D, and XMV and YMV need to turn 1-D.
  using Kokkos::ALL;
  using Kokkos::subview;
  typedef Kokkos::View<RV::value_type, RV::array_layout,
    RV::device_type, RV::memory_traits, RV::specialize> RV0D;
  typedef Kokkos::View<XMV::const_value_type*, XMV::array_layout, XMV::device_type, XMV::memory_traits, XMV::specialize> XMV1D;
  typedef Kokkos::View<YMV::const_value_type*, YMV::array_layout, YMV::device_type, YMV::memory_traits, YMV::specialize> YMV1D;

  const size_type numRows = X.dimension_0 ();
  const size_type numCols = X.dimension_1 ();
  if (numRows < static_cast<size_type> (INT_MAX) &&
      numRows * numCols < static_cast<size_type> (INT_MAX)) {
    typedef V_Dot_Functor<RV0D, XMV1D, YMV1D, int> op_type;
    op_type op (subview (r, r_col), subview (X, ALL (), X_col),
                subview (Y, ALL (), Y_col));
    Kokkos::parallel_reduce (numRows, op);
  }
  else {
    typedef V_Dot_Functor<RV0D, XMV1D, YMV1D, size_type> op_type;
    op_type op (subview (r, r_col), subview (X, ALL (), X_col),
                subview (Y, ALL (), Y_col));
    Kokkos::parallel_reduce (numRows, op);
  }
}

#undef KOKKOSBLAS_IMPL_MV_EXEC_SPACE
#undef KOKKOSBLAS_IMPL_MV_MEM_SPACE
#endif // KOKKOS_HAVE_PTHREAD


#ifdef KOKKOS_HAVE_CUDA
#define KOKKOSBLAS_IMPL_MV_EXEC_SPACE Kokkos::Cuda
#define KOKKOSBLAS_IMPL_MV_MEM_SPACE Kokkos::CudaSpace

void
Dot_MV<double*,
       KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
       Kokkos::Impl::ViewDefault,
       const double**,
       Kokkos::LayoutLeft,
       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
       Kokkos::Impl::ViewDefault,
       const double**,
       Kokkos::LayoutLeft,
       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
       Kokkos::Impl::ViewDefault>::
dot (const RV& r, const XMV& X, const XMV& Y)
{
  const size_type numRows = X.dimension_0 ();
  const size_type numVecs = X.dimension_1 ();

  // TODO (mfh 06 Apr 2015) Figure out the optimal max unroll length.
  // It depends on the number of concurrent memory streams that the
  // processor supports, not so much on the SIMD length.

  if (numVecs > 16) {
    MV_Dot_Right_FunctorVector<RV, XMV, YMV> op (r, X, Y);
    Kokkos::parallel_reduce (numRows, op);
  }
  else {
    switch (numVecs) {
    case 16: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 16> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 15: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 15> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 14: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 14> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 13: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 13> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 12: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 12> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 11: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 11> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 10: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 10> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 9: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 9> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 8: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 8> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 7: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 7> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 6: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 6> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 5: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 5> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 4: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 4> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 3: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 3> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 2: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 2> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 1: {
      // RV needs to turn 0-D, and XMV and YMV need to turn 1-D.
      using Kokkos::ALL;
      using Kokkos::subview;
      // RV needs to turn 0-D, and XMV and YMV need to turn 1-D.
      typedef Kokkos::View<RV::value_type, RV::array_layout, RV::device_type, RV::memory_traits, RV::specialize> RV0D;
      typedef Kokkos::View<XMV::value_type*, XMV::array_layout, XMV::device_type, XMV::memory_traits, XMV::specialize> XMV1D;
      typedef Kokkos::View<YMV::value_type*, YMV::array_layout, YMV::device_type, YMV::memory_traits, YMV::specialize> YMV1D;
      typedef V_Dot_Functor<RV0D, XMV1D, YMV1D> op_type;
      op_type op (subview (r, 0), subview (X, ALL (), 0), subview (Y, ALL (), 0));
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    } // switch
  } // if-else
}

void
Dot_MV<double*,
       KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
       Kokkos::Impl::ViewDefault,
       const double**,
       Kokkos::LayoutLeft,
       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
       Kokkos::Impl::ViewDefault,
       const double**,
       Kokkos::LayoutLeft,
       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
       Kokkos::Impl::ViewDefault>::
dot (const RV& r, const size_t r_col, const XMV& X, const size_t X_col,
     const XMV& Y, const size_t Y_col)
{
  // RV needs to turn 0-D, and XMV and YMV need to turn 1-D.
  using Kokkos::ALL;
  using Kokkos::subview;
  typedef Kokkos::View<RV::value_type, RV::array_layout,
    RV::device_type, RV::memory_traits, RV::specialize> RV0D;
  typedef Kokkos::View<XMV::const_value_type*, XMV::array_layout, XMV::device_type, XMV::memory_traits, XMV::specialize> XMV1D;
  typedef Kokkos::View<YMV::const_value_type*, YMV::array_layout, YMV::device_type, YMV::memory_traits, YMV::specialize> YMV1D;

  const size_type numRows = X.dimension_0 ();
  const size_type numCols = X.dimension_1 ();
  if (numRows < static_cast<size_type> (INT_MAX) &&
      numRows * numCols < static_cast<size_type> (INT_MAX)) {
    typedef V_Dot_Functor<RV0D, XMV1D, YMV1D, int> op_type;
    op_type op (subview (r, r_col), subview (X, ALL (), X_col),
                subview (Y, ALL (), Y_col));
    Kokkos::parallel_reduce (numRows, op);
  }
  else {
    typedef V_Dot_Functor<RV0D, XMV1D, YMV1D, size_type> op_type;
    op_type op (subview (r, r_col), subview (X, ALL (), X_col),
                subview (Y, ALL (), Y_col));
    Kokkos::parallel_reduce (numRows, op);
  }
}

#undef KOKKOSBLAS_IMPL_MV_EXEC_SPACE
#undef KOKKOSBLAS_IMPL_MV_MEM_SPACE
#endif // KOKKOS_HAVE_CUDA


#ifdef KOKKOS_HAVE_CUDA
#define KOKKOSBLAS_IMPL_MV_EXEC_SPACE Kokkos::Cuda
#define KOKKOSBLAS_IMPL_MV_MEM_SPACE Kokkos::CudaUVMSpace

void
Dot_MV<double*,
       KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
       Kokkos::Impl::ViewDefault,
       const double**,
       Kokkos::LayoutLeft,
       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
       Kokkos::Impl::ViewDefault,
       const double**,
       Kokkos::LayoutLeft,
       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
       Kokkos::Impl::ViewDefault>::
dot (const RV& r, const XMV& X, const XMV& Y)
{
  const size_type numRows = X.dimension_0 ();
  const size_type numVecs = X.dimension_1 ();

  // TODO (mfh 06 Apr 2015) Figure out the optimal max unroll length.
  // It depends on the number of concurrent memory streams that the
  // processor supports, not so much on the SIMD length.

  if (numVecs > 16) {
    MV_Dot_Right_FunctorVector<RV, XMV, YMV> op (r, X, Y);
    Kokkos::parallel_reduce (numRows, op);
  }
  else {
    switch (numVecs) {
    case 16: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 16> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 15: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 15> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 14: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 14> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 13: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 13> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 12: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 12> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 11: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 11> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 10: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 10> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 9: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 9> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 8: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 8> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 7: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 7> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 6: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 6> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 5: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 5> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 4: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 4> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 3: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 3> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 2: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 2> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    case 1: {
      // RV needs to turn 0-D, and XMV and YMV need to turn 1-D.
      using Kokkos::ALL;
      using Kokkos::subview;
      // RV needs to turn 0-D, and XMV and YMV need to turn 1-D.
      typedef Kokkos::View<RV::value_type, RV::array_layout, RV::device_type, RV::memory_traits, RV::specialize> RV0D;
      typedef Kokkos::View<XMV::value_type*, XMV::array_layout, XMV::device_type, XMV::memory_traits, XMV::specialize> XMV1D;
      typedef Kokkos::View<YMV::value_type*, YMV::array_layout, YMV::device_type, YMV::memory_traits, YMV::specialize> YMV1D;
      typedef V_Dot_Functor<RV0D, XMV1D, YMV1D> op_type;
      op_type op (subview (r, 0), subview (X, ALL (), 0), subview (Y, ALL (), 0));
      Kokkos::parallel_reduce (numRows, op);
      break;
    }
    } // switch
  } // if-else
}

void
Dot_MV<double*,
       KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
       Kokkos::Impl::ViewDefault,
       const double**,
       Kokkos::LayoutLeft,
       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
       Kokkos::Impl::ViewDefault,
       const double**,
       Kokkos::LayoutLeft,
       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
       Kokkos::Impl::ViewDefault>::
dot (const RV& r, const size_t r_col, const XMV& X, const size_t X_col,
     const XMV& Y, const size_t Y_col)
{
  // RV needs to turn 0-D, and XMV and YMV need to turn 1-D.
  using Kokkos::ALL;
  using Kokkos::subview;
  typedef Kokkos::View<RV::value_type, RV::array_layout,
    RV::device_type, RV::memory_traits, RV::specialize> RV0D;
  typedef Kokkos::View<XMV::const_value_type*, XMV::array_layout, XMV::device_type, XMV::memory_traits, XMV::specialize> XMV1D;
  typedef Kokkos::View<YMV::const_value_type*, YMV::array_layout, YMV::device_type, YMV::memory_traits, YMV::specialize> YMV1D;

  const size_type numRows = X.dimension_0 ();
  const size_type numCols = X.dimension_1 ();
  if (numRows < static_cast<size_type> (INT_MAX) &&
      numRows * numCols < static_cast<size_type> (INT_MAX)) {
    typedef V_Dot_Functor<RV0D, XMV1D, YMV1D, int> op_type;
    op_type op (subview (r, r_col), subview (X, ALL (), X_col),
                subview (Y, ALL (), Y_col));
    Kokkos::parallel_reduce (numRows, op);
  }
  else {
    typedef V_Dot_Functor<RV0D, XMV1D, YMV1D, size_type> op_type;
    op_type op (subview (r, r_col), subview (X, ALL (), X_col),
                subview (Y, ALL (), Y_col));
    Kokkos::parallel_reduce (numRows, op);
  }
}

#undef KOKKOSBLAS_IMPL_MV_EXEC_SPACE
#undef KOKKOSBLAS_IMPL_MV_MEM_SPACE
#endif // KOKKOS_HAVE_CUDA

} // namespace Impl
} // namespace KokkosBlas

