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
#ifndef KOKKOSBLAS3_GEMV_HPP_
#define KOKKOSBLAS3_GEMV_HPP_

/// \file KokkosBlas3_gemm.hpp

#include <KokkosKernels_Macros.hpp>
#include <KokkosBlas3_gemm_spec.hpp>
#include <KokkosKernels_helpers.hpp>
#include <sstream>
#include <type_traits>

namespace KokkosBlas {

/// \brief Dense matrix-matrix multiply: C = beta*C + alpha*op(A)*op(B).
///
/// \tparam AViewType Input matrix, as a 2-D Kokkos::View
/// \tparam BViewType Input matrix, as a 2-D Kokkos::View
/// \tparam CViewType Output matrix, as a nonconst 2-D Kokkos::View
///
/// \param transA [in] "N" for non-transpose, "T" for transpose, "C"
///   for conjugate transpose.  All characters after the first are
///   ignored.  This works just like the BLAS routines.
/// \param transB [in] "N" for non-transpose, "T" for transpose, "C"
///   for conjugate transpose.  All characters after the first are
///   ignored.  This works just like the BLAS routines.
/// \param alpha [in] Input coefficient of A*x
/// \param A [in] Input matrix, as a 2-D Kokkos::View
/// \param B [in] Input matrix, as a 2-D Kokkos::View
/// \param beta [in] Input coefficient of C
/// \param C [in/out] Output vector, as a nonconst 2-D Kokkos::View
template<class AViewType,
         class BViewType,
         class CViewType>
void
gemm (const char transA[],
      const char transB[],
      typename AViewType::const_value_type& alpha,
      const AViewType& A,
      const BViewType& B,
      typename CViewType::const_value_type& beta,
      const CViewType& C)
{

  #if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
  static_assert (Kokkos::Impl::is_view<AViewType>::value,
                 "AViewType must be a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<BViewType>::value,
                 "BViewType must be a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<CViewType>::value,
                 "CViewType must be a Kokkos::View.");
  static_assert (static_cast<int> (AViewType::rank) == 2,
                 "AViewType must have rank 2.");
  static_assert (static_cast<int> (BViewType::rank) == 2,
                 "BViewType must have rank 2.");
  static_assert (static_cast<int> (CViewType::rank) == 2,
                 "CViewType must have rank 2.");

  // Check validity of transpose argument
  bool valid_transA = (transA[0] == 'N') || (transA[0] == 'n') ||
                      (transA[0] == 'T') || (transA[0] == 't') ||
                      (transA[0] == 'C') || (transA[0] == 'c');
  bool valid_transB = (transB[0] == 'N') || (transB[0] == 'n') ||
                      (transB[0] == 'T') || (transB[0] == 't') ||
                      (transB[0] == 'C') || (transB[0] == 'c');
  if(!(valid_transA && valid_transB)) {
    std::ostringstream os;
    os << "KokkosBlas::gemm: transA[0] = '" << transA[0] << " transB[0] = '" << transB[0] << "'. " <<
      "Valid values include 'N' or 'n' (No transpose), 'T' or 't' (Transpose), "
      "and 'C' or 'c' (Conjugate transpose).";
    Kokkos::Impl::throw_runtime_exception (os.str ());
  }

  // Check compatibility of dimensions at run time.
  bool A_t = !(transA[0] == 'N' || transA[0] == 'n');
  bool B_t = !(transB[0] == 'N' || transB[0] == 'n');
  int64_t A0 = A.extent(0);
  int64_t A1 = A.extent(1);
  int64_t B0 = B.extent(0);
  int64_t B1 = B.extent(1);
  int64_t C0 = C.extent(0);
  int64_t C1 = C.extent(1);

  if ( ((A_t?A1:A0) != C0) ||
       ((B_t?B0:B1) != C1) ||
       ((A_t?A0:A1) != (B_t?B1:B0)) ) {
      std::ostringstream os;
      os << "KokkosBlas::gemm: Dimensions of A, B, and C do not match: "
         << "transA: " << transA[0] << " transB: " << transB[0]
         << " A: " << A.extent(0) << " x " << A.extent(1)
         << " B: " << B.extent(0) << " x " << B.extent(1)
         << " C: " << C.extent(0) << " x " << C.extent(1);
      Kokkos::Impl::throw_runtime_exception (os.str ());
    }
  #endif // KOKKOSKERNELS_DEBUG_LEVEL > 0

  // Return if degenerated matrices are provided
  if((A.extent(0) == 0) || (A.extent(1) == 0) || (C.extent(1) == 0))
    return;

  // Minimize the number of Impl::GEMV instantiations, by
  // standardizing on particular View specializations for its template
  // parameters.
  typedef Kokkos::View<typename AViewType::const_value_type**,
    typename AViewType::array_layout,
    typename AViewType::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged> > AVT;
  typedef Kokkos::View<typename BViewType::const_value_type**,
    typename BViewType::array_layout,
    typename BViewType::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged> > BVT;
  typedef Kokkos::View<typename CViewType::non_const_value_type**,
    typename CViewType::array_layout,
    typename CViewType::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged> > CVT;
  typedef Impl::GEMM<AVT, BVT, CVT> impl_type;
  impl_type::gemm (transA, transB, alpha, A, B, beta, C);
}

} // namespace KokkosBlas

#endif // KOKKOS_BLAS3_MV_HPP_
