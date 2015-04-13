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

#include <Kokkos_Blas1_MV_impl_abs.hpp>
#include <climits>

namespace KokkosBlas {
namespace Impl {

#ifdef KOKKOS_HAVE_SERIAL
#define KOKKOSBLAS_IMPL_MV_EXEC_SPACE Kokkos::Serial
#define KOKKOSBLAS_IMPL_MV_MEM_SPACE Kokkos::HostSpace

void
Abs<Kokkos::View<double**,
                 Kokkos::LayoutLeft,
                 Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                 Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                 Kokkos::Impl::ViewDefault>,
    Kokkos::View<const double**,
                 Kokkos::LayoutLeft,
                 Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                 Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                 Kokkos::Impl::ViewDefault>,
    2>::
abs (const RMV& R, const XMV& X)
{
#ifdef KOKKOS_HAVE_CXX11
  static_assert (Kokkos::Impl::is_view<RMV>::value, "KokkosBlas::Impl::"
                 "Abs<2-D>: RMV is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::"
                 "Abs<2-D>: XMV is not a Kokkos::View.");
  static_assert (RMV::rank == 2, "KokkosBlas::Impl::Abs<2-D>: "
                 "RMV is not rank 2.");
  static_assert (XMV::rank == 2, "KokkosBlas::Impl::Abs<2-D>: "
                 "XMV is not rank 2.");
#endif // KOKKOS_HAVE_CXX11
  const size_type numRows = X.dimension_0 ();
  const size_type numCols = X.dimension_1 ();
  if (numRows < static_cast<size_type> (INT_MAX) &&
      numRows * numCols < static_cast<size_type> (INT_MAX)) {
    typedef int index_type;
    MV_Abs_Generic<RMV, XMV, index_type> (R, X);
  }
  else {
    typedef XMV::size_type index_type;
    MV_Abs_Generic<RMV, XMV, index_type> (R, X);
  }
}

#undef KOKKOSBLAS_IMPL_MV_EXEC_SPACE
#undef KOKKOSBLAS_IMPL_MV_MEM_SPACE
#endif // KOKKOS_HAVE_SERIAL


#ifdef KOKKOS_HAVE_OPENMP
#define KOKKOSBLAS_IMPL_MV_EXEC_SPACE Kokkos::OpenMP
#define KOKKOSBLAS_IMPL_MV_MEM_SPACE Kokkos::HostSpace

void
Abs<Kokkos::View<double**,
                 Kokkos::LayoutLeft,
                 Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                 Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                 Kokkos::Impl::ViewDefault>,
    Kokkos::View<const double**,
                 Kokkos::LayoutLeft,
                 Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                 Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                 Kokkos::Impl::ViewDefault>,
    2>::
abs (const RMV& R, const XMV& X)
{
#ifdef KOKKOS_HAVE_CXX11
  static_assert (Kokkos::Impl::is_view<RMV>::value, "KokkosBlas::Impl::"
                 "Abs<2-D>: RMV is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::"
                 "Abs<2-D>: XMV is not a Kokkos::View.");
  static_assert (RMV::rank == 2, "KokkosBlas::Impl::Abs<2-D>: "
                 "RMV is not rank 2.");
  static_assert (XMV::rank == 2, "KokkosBlas::Impl::Abs<2-D>: "
                 "XMV is not rank 2.");
#endif // KOKKOS_HAVE_CXX11
  const size_type numRows = X.dimension_0 ();
  const size_type numCols = X.dimension_1 ();
  if (numRows < static_cast<size_type> (INT_MAX) &&
      numRows * numCols < static_cast<size_type> (INT_MAX)) {
    typedef int index_type;
    MV_Abs_Generic<RMV, XMV, index_type> (R, X);
  }
  else {
    typedef XMV::size_type index_type;
    MV_Abs_Generic<RMV, XMV, index_type> (R, X);
  }
}

#undef KOKKOSBLAS_IMPL_MV_EXEC_SPACE
#undef KOKKOSBLAS_IMPL_MV_MEM_SPACE
#endif // KOKKOS_HAVE_OPENMP


#ifdef KOKKOS_HAVE_PTHREAD
#define KOKKOSBLAS_IMPL_MV_EXEC_SPACE Kokkos::Threads
#define KOKKOSBLAS_IMPL_MV_MEM_SPACE Kokkos::HostSpace

void
Abs<Kokkos::View<double**,
                 Kokkos::LayoutLeft,
                 Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                 Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                 Kokkos::Impl::ViewDefault>,
    Kokkos::View<const double**,
                 Kokkos::LayoutLeft,
                 Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                 Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                 Kokkos::Impl::ViewDefault>,
    2>::
abs (const RMV& R, const XMV& X)
{
#ifdef KOKKOS_HAVE_CXX11
  static_assert (Kokkos::Impl::is_view<RMV>::value, "KokkosBlas::Impl::"
                 "Abs<2-D>: RMV is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::"
                 "Abs<2-D>: XMV is not a Kokkos::View.");
  static_assert (RMV::rank == 2, "KokkosBlas::Impl::Abs<2-D>: "
                 "RMV is not rank 2.");
  static_assert (XMV::rank == 2, "KokkosBlas::Impl::Abs<2-D>: "
                 "XMV is not rank 2.");
#endif // KOKKOS_HAVE_CXX11
  const size_type numRows = X.dimension_0 ();
  const size_type numCols = X.dimension_1 ();
  if (numRows < static_cast<size_type> (INT_MAX) &&
      numRows * numCols < static_cast<size_type> (INT_MAX)) {
    typedef int index_type;
    MV_Abs_Generic<RMV, XMV, index_type> (R, X);
  }
  else {
    typedef XMV::size_type index_type;
    MV_Abs_Generic<RMV, XMV, index_type> (R, X);
  }
}

#undef KOKKOSBLAS_IMPL_MV_EXEC_SPACE
#undef KOKKOSBLAS_IMPL_MV_MEM_SPACE
#endif // KOKKOS_HAVE_PTHREAD


#ifdef KOKKOS_HAVE_CUDA
#define KOKKOSBLAS_IMPL_MV_EXEC_SPACE Kokkos::Cuda
#define KOKKOSBLAS_IMPL_MV_MEM_SPACE Kokkos::CudaSpace

void
Abs<Kokkos::View<double**,
                 Kokkos::LayoutLeft,
                 Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                 Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                 Kokkos::Impl::ViewDefault>,
    Kokkos::View<const double**,
                 Kokkos::LayoutLeft,
                 Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                 Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                 Kokkos::Impl::ViewDefault>,
    2>::
abs (const RMV& R, const XMV& X)
{
#ifdef KOKKOS_HAVE_CXX11
  static_assert (Kokkos::Impl::is_view<RMV>::value, "KokkosBlas::Impl::"
                 "Abs<2-D>: RMV is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::"
                 "Abs<2-D>: XMV is not a Kokkos::View.");
  static_assert (RMV::rank == 2, "KokkosBlas::Impl::Abs<2-D>: "
                 "RMV is not rank 2.");
  static_assert (XMV::rank == 2, "KokkosBlas::Impl::Abs<2-D>: "
                 "XMV is not rank 2.");
#endif // KOKKOS_HAVE_CXX11
  const size_type numRows = X.dimension_0 ();
  const size_type numCols = X.dimension_1 ();
  if (numRows < static_cast<size_type> (INT_MAX) &&
      numRows * numCols < static_cast<size_type> (INT_MAX)) {
    typedef int index_type;
    MV_Abs_Generic<RMV, XMV, index_type> (R, X);
  }
  else {
    typedef XMV::size_type index_type;
    MV_Abs_Generic<RMV, XMV, index_type> (R, X);
  }
}

#undef KOKKOSBLAS_IMPL_MV_EXEC_SPACE
#undef KOKKOSBLAS_IMPL_MV_MEM_SPACE
#endif // KOKKOS_HAVE_CUDA


#ifdef KOKKOS_HAVE_CUDA
#define KOKKOSBLAS_IMPL_MV_EXEC_SPACE Kokkos::Cuda
#define KOKKOSBLAS_IMPL_MV_MEM_SPACE Kokkos::CudaUVMSpace

void
Abs<Kokkos::View<double**,
                 Kokkos::LayoutLeft,
                 Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                 Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                 Kokkos::Impl::ViewDefault>,
    Kokkos::View<const double**,
                 Kokkos::LayoutLeft,
                 Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                 Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                 Kokkos::Impl::ViewDefault>,
    2>::
abs (const RMV& R, const XMV& X)
{
#ifdef KOKKOS_HAVE_CXX11
  static_assert (Kokkos::Impl::is_view<RMV>::value, "KokkosBlas::Impl::"
                 "Abs<2-D>: RMV is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::"
                 "Abs<2-D>: XMV is not a Kokkos::View.");
  static_assert (RMV::rank == 2, "KokkosBlas::Impl::Abs<2-D>: "
                 "RMV is not rank 2.");
  static_assert (XMV::rank == 2, "KokkosBlas::Impl::Abs<2-D>: "
                 "XMV is not rank 2.");
#endif // KOKKOS_HAVE_CXX11
  const size_type numRows = X.dimension_0 ();
  const size_type numCols = X.dimension_1 ();
  if (numRows < static_cast<size_type> (INT_MAX) &&
      numRows * numCols < static_cast<size_type> (INT_MAX)) {
    typedef int index_type;
    MV_Abs_Generic<RMV, XMV, index_type> (R, X);
  }
  else {
    typedef XMV::size_type index_type;
    MV_Abs_Generic<RMV, XMV, index_type> (R, X);
  }
}

#undef KOKKOSBLAS_IMPL_MV_EXEC_SPACE
#undef KOKKOSBLAS_IMPL_MV_MEM_SPACE
#endif // KOKKOS_HAVE_CUDA

} // namespace Impl
} // namespace KokkosBlas

