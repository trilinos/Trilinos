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

#include <Kokkos_Blas1_MV_impl_nrm1.hpp>

namespace KokkosBlas {
namespace Impl {

#ifdef KOKKOS_HAVE_SERIAL
#define KOKKOSBLAS_IMPL_MV_EXEC_SPACE Kokkos::Serial
#define KOKKOSBLAS_IMPL_MV_MEM_SPACE Kokkos::HostSpace
#define KOKKOSBLAS_IMPL_MV_SCALAR double

void
Nrm1_MV<Kokkos::View<Kokkos::Details::InnerProductSpaceTraits<KOKKOSBLAS_IMPL_MV_SCALAR>::mag_type*,
                     KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
                     Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                     Kokkos::Impl::ViewDefault>,
        Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
                     Kokkos::LayoutLeft,
                     Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                     Kokkos::Impl::ViewDefault>,
        2>::
nrm1 (const RV& r, const XMV& X)
{
  typedef XMV::size_type size_type;
  const size_type numRows = X.dimension_0 ();
  const size_type numCols = X.dimension_1 ();

  // int is generally faster than size_t, but check for overflow first.
  if (numRows < static_cast<size_type> (INT_MAX) &&
      numRows * numCols < static_cast<size_type> (INT_MAX)) {
    MV_Nrm1_Invoke<RV, XMV, int> (r, X);
  }
  else {
    MV_Nrm1_Invoke<RV, XMV, size_type> (r, X);
  }
}

#undef KOKKOSBLAS_IMPL_MV_EXEC_SPACE
#undef KOKKOSBLAS_IMPL_MV_MEM_SPACE
#undef KOKKOSBLAS_IMPL_MV_SCALAR
#endif // KOKKOS_HAVE_SERIAL


#ifdef KOKKOS_HAVE_OPENMP
#define KOKKOSBLAS_IMPL_MV_EXEC_SPACE Kokkos::OpenMP
#define KOKKOSBLAS_IMPL_MV_MEM_SPACE Kokkos::HostSpace
#define KOKKOSBLAS_IMPL_MV_SCALAR double

void
Nrm1_MV<Kokkos::View<Kokkos::Details::InnerProductSpaceTraits<KOKKOSBLAS_IMPL_MV_SCALAR>::mag_type*,
                     KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
                     Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                     Kokkos::Impl::ViewDefault>,
        Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
                     Kokkos::LayoutLeft,
                     Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                     Kokkos::Impl::ViewDefault>,
        2>::
nrm1 (const RV& r, const XMV& X)
{
  typedef XMV::size_type size_type;
  const size_type numRows = X.dimension_0 ();
  const size_type numCols = X.dimension_1 ();

  // int is generally faster than size_t, but check for overflow first.
  if (numRows < static_cast<size_type> (INT_MAX) &&
      numRows * numCols < static_cast<size_type> (INT_MAX)) {
    MV_Nrm1_Invoke<RV, XMV, int> (r, X);
  }
  else {
    MV_Nrm1_Invoke<RV, XMV, size_type> (r, X);
  }
}

#undef KOKKOSBLAS_IMPL_MV_EXEC_SPACE
#undef KOKKOSBLAS_IMPL_MV_MEM_SPACE
#undef KOKKOSBLAS_IMPL_MV_SCALAR
#endif // KOKKOS_HAVE_OPENMP


#ifdef KOKKOS_HAVE_PTHREAD
#define KOKKOSBLAS_IMPL_MV_EXEC_SPACE Kokkos::Threads
#define KOKKOSBLAS_IMPL_MV_MEM_SPACE Kokkos::HostSpace
#define KOKKOSBLAS_IMPL_MV_SCALAR double

void
Nrm1_MV<Kokkos::View<Kokkos::Details::InnerProductSpaceTraits<KOKKOSBLAS_IMPL_MV_SCALAR>::mag_type*,
                     KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
                     Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                     Kokkos::Impl::ViewDefault>,
        Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
                     Kokkos::LayoutLeft,
                     Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                     Kokkos::Impl::ViewDefault>,
        2>::
nrm1 (const RV& r, const XMV& X)
{
  typedef XMV::size_type size_type;
  const size_type numRows = X.dimension_0 ();
  const size_type numCols = X.dimension_1 ();

  // int is generally faster than size_t, but check for overflow first.
  if (numRows < static_cast<size_type> (INT_MAX) &&
      numRows * numCols < static_cast<size_type> (INT_MAX)) {
    MV_Nrm1_Invoke<RV, XMV, int> (r, X);
  }
  else {
    MV_Nrm1_Invoke<RV, XMV, size_type> (r, X);
  }
}

#undef KOKKOSBLAS_IMPL_MV_EXEC_SPACE
#undef KOKKOSBLAS_IMPL_MV_MEM_SPACE
#undef KOKKOSBLAS_IMPL_MV_SCALAR
#endif // KOKKOS_HAVE_PTHREAD


#ifdef KOKKOS_HAVE_CUDA
#define KOKKOSBLAS_IMPL_MV_EXEC_SPACE Kokkos::Cuda
#define KOKKOSBLAS_IMPL_MV_MEM_SPACE Kokkos::CudaSpace
#define KOKKOSBLAS_IMPL_MV_SCALAR double

void
Nrm1_MV<Kokkos::View<Kokkos::Details::InnerProductSpaceTraits<KOKKOSBLAS_IMPL_MV_SCALAR>::mag_type*,
                     KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
                     Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                     Kokkos::Impl::ViewDefault>,
        Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
                     Kokkos::LayoutLeft,
                     Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                     Kokkos::Impl::ViewDefault>,
        2>::
nrm1 (const RV& r, const XMV& X)
{
  typedef XMV::size_type size_type;
  const size_type numRows = X.dimension_0 ();
  const size_type numCols = X.dimension_1 ();

  // int is generally faster than size_t, but check for overflow first.
  if (numRows < static_cast<size_type> (INT_MAX) &&
      numRows * numCols < static_cast<size_type> (INT_MAX)) {
    MV_Nrm1_Invoke<RV, XMV, int> (r, X);
  }
  else {
    MV_Nrm1_Invoke<RV, XMV, size_type> (r, X);
  }
}

#undef KOKKOSBLAS_IMPL_MV_EXEC_SPACE
#undef KOKKOSBLAS_IMPL_MV_MEM_SPACE
#undef KOKKOSBLAS_IMPL_MV_SCALAR
#endif // KOKKOS_HAVE_CUDA


#ifdef KOKKOS_HAVE_CUDA
#define KOKKOSBLAS_IMPL_MV_EXEC_SPACE Kokkos::Cuda
#define KOKKOSBLAS_IMPL_MV_MEM_SPACE Kokkos::CudaUVMSpace
#define KOKKOSBLAS_IMPL_MV_SCALAR double

void
Nrm1_MV<Kokkos::View<Kokkos::Details::InnerProductSpaceTraits<KOKKOSBLAS_IMPL_MV_SCALAR>::mag_type*,
                     KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
                     Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                     Kokkos::Impl::ViewDefault>,
        Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
                     Kokkos::LayoutLeft,
                     Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                     Kokkos::Impl::ViewDefault>,
        2>::
nrm1 (const RV& r, const XMV& X)
{
  typedef XMV::size_type size_type;
  const size_type numRows = X.dimension_0 ();
  const size_type numCols = X.dimension_1 ();

  // int is generally faster than size_t, but check for overflow first.
  if (numRows < static_cast<size_type> (INT_MAX) &&
      numRows * numCols < static_cast<size_type> (INT_MAX)) {
    MV_Nrm1_Invoke<RV, XMV, int> (r, X);
  }
  else {
    MV_Nrm1_Invoke<RV, XMV, size_type> (r, X);
  }
}

#undef KOKKOSBLAS_IMPL_MV_EXEC_SPACE
#undef KOKKOSBLAS_IMPL_MV_MEM_SPACE
#undef KOKKOSBLAS_IMPL_MV_SCALAR
#endif // KOKKOS_HAVE_CUDA

} // namespace Impl
} // namespace KokkosBlas

