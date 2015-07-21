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

#include <Kokkos_Sparse.hpp>

namespace KokkosSparse {
namespace Impl {

#ifdef KOKKOS_HAVE_OPENMP

KOKKOSSPARSE_IMPL_SPMV_DEF( double, int, size_t, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace )

#define KOKKOSSPARSE_IMPL_MV_EXEC_SPACE Kokkos::OpenMP
#define KOKKOSSPARSE_IMPL_MV_MEM_SPACE Kokkos::HostSpace
#define KOKKOSSPARSE_IMPL_MV_DEVICE_TYPE Kokkos::Device<KOKKOSSPARSE_IMPL_MV_EXEC_SPACE,KOKKOSSPARSE_IMPL_MV_MEM_SPACE>
#define KOKKOSSPARSE_IMPL_MV_SCALAR double

template<>
struct SPMV<const KOKKOSSPARSE_IMPL_MV_SCALAR,
            const size_t,
            KOKKOSSPARSE_IMPL_MV_DEVICE_TYPE,
            Kokkos::MemoryTraits<Kokkos::Unmanaged>,
            size_t,
            const KOKKOSSPARSE_IMPL_MV_SCALAR*,
            Kokkos::LayoutLeft,
            KOKKOSSPARSE_IMPL_MV_DEVICE_TYPE,
            Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess>,
            Kokkos::Impl::ViewDefault,
            KOKKOSSPARSE_IMPL_MV_SCALAR*,
            Kokkos::LayoutLeft,
            KOKKOSSPARSE_IMPL_MV_DEVICE_TYPE,
            Kokkos::MemoryTraits<Kokkos::Unmanaged>,
            Kokkos::Impl::ViewDefault>
{
  typedef CrsMatrix<const KOKKOSSPARSE_IMPL_MV_SCALAR,
                    const size_t,
                    KOKKOSSPARSE_IMPL_MV_DEVICE_TYPE,
                    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                    size_t> AMatrix;
  typedef Kokkos::View<const KOKKOSSPARSE_IMPL_MV_SCALAR*,
                       Kokkos::LayoutLeft,
                       KOKKOSSPARSE_IMPL_MV_DEVICE_TYPE,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess>,
                       Kokkos::Impl::ViewDefault> XVector;
  typedef Kokkos::View<KOKKOSSPARSE_IMPL_MV_SCALAR*,
                       Kokkos::LayoutLeft,
                       KOKKOSSPARSE_IMPL_MV_DEVICE_TYPE,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> YVector;
  typedef typename YVector::non_const_value_type Scalar;

  static void
  spmv (const char mode[], const Scalar& alpha, const AMatrix& A,
        const XVector& x, const Scalar& beta, const YVector& y)
  {
    if (alpha == Kokkos::Details::ArithTraits<Scalar>::zero ()) {
      spmv_alpha<AMatrix,XVector,YVector,0> (mode, alpha, A, x, beta, y);
      return;
    }
    if (alpha == Kokkos::Details::ArithTraits<Scalar>::one ()) {
      spmv_alpha<AMatrix,XVector,YVector,1> (mode, alpha, A, x, beta, y);
      return;
    }
    if (alpha == -Kokkos::Details::ArithTraits<Scalar>::one ()) {
      spmv_alpha<AMatrix,XVector,YVector,-1> (mode, alpha, A, x, beta, y);
      return;
    }
    spmv_alpha<AMatrix,XVector,YVector,2> (mode, alpha, A, x, beta, y);
  }

};
#endif // KOKKOS_HAVE_OPENMP

} // namespace Impl
} // namespace KokkosSparse

