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
#ifndef KOKKOSBLAS2_GEMV_HPP_
#define KOKKOSBLAS2_GEMV_HPP_

/// \file Kokkos_Blas2_MV.hpp
/// \brief BLAS 2 kernels specifically optimized for typical
///   Tpetra::MultiVector use cases.

#include <KokkosBlas2_gemv_spec.hpp>
#include <KokkosKernels_helpers.hpp>
#include <sstream>
#include <type_traits> // requires C++11, but so does Kokkos

namespace KokkosBlas {

/// \brief Dense matrix-vector multiply: y = beta*y + alpha*A*x.
///
/// \tparam AViewType Input matrix, as a 2-D Kokkos::View
/// \tparam XViewType Input vector, as a 1-D Kokkos::View
/// \tparam YViewType Output vector, as a nonconst 1-D Kokkos::View
/// \tparam AlphaCoeffType Type of input coefficient alpha
/// \tparam BetaCoeffType Type of input coefficient beta
///
/// \param trans [in] "N" for non-transpose, "T" for transpose, "C"
///   for conjugate transpose.  All characters after the first are
///   ignored.  This works just like the BLAS routines.
/// \param alpha [in] Input coefficient of A*x
/// \param A [in] Input matrix, as a 2-D Kokkos::View
/// \param x [in] Input vector, as a 1-D Kokkos::View
/// \param beta [in] Input coefficient of y
/// \param y [in/out] Output vector, as a nonconst 1-D Kokkos::View
template<class AViewType,
         class XViewType,
         class YViewType>
void
gemv (const char trans[],
      typename AViewType::const_value_type& alpha,
      const AViewType& A,
      const XViewType& x,
      typename YViewType::const_value_type& beta,
      const YViewType& y)
{
  static_assert (Kokkos::Impl::is_view<AViewType>::value,
                 "AViewType must be a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<XViewType>::value,
                 "XViewType must be a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<YViewType>::value,
                 "YViewType must be a Kokkos::View.");
  static_assert (static_cast<int> (AViewType::rank) == 2,
                 "AViewType must have rank 2.");
  static_assert (static_cast<int> (XViewType::rank) == 1,
                 "XViewType must have rank 1.");
  static_assert (static_cast<int> (YViewType::rank) == 1,
                 "YViewType must have rank 1.");

  // Check compatibility of dimensions at run time.
  if (trans[0] == 'N' || trans[0] == 'n') {
    if (A.extent(0) != y.extent(0) || A.extent(1) != x.extent(0)) {
      std::ostringstream os;
      os << "KokkosBlas::gemv: Dimensions of A, x, and y do not match: "
         << "A: " << A.extent(0) << " x " << A.extent(1)
         << ", x: " << x.extent(0) << ", y: " << y.extent(0);
      Kokkos::Impl::throw_runtime_exception (os.str ());
    }
  }
  else if (trans[0] == 'T' || trans[0] == 't' ||
           trans[0] == 'C' || trans[0] == 'c' ||
           trans[0] == 'H' || trans[0] == 'h') {
    if (A.extent(1) != y.extent(0) || A.extent(0) != x.extent(0)) {
      std::ostringstream os;
      os << "KokkosBlas::dot: Dimensions of A, x, and y do not match: "
         << "A: " << A.extent(0) << " x " << A.extent(1)
         << ", x: " << x.extent(0) << ", y: " << y.extent(0);
      Kokkos::Impl::throw_runtime_exception (os.str ());
    }
  }
  else {
    std::ostringstream os;
    os << "KokkosBlas::gemv: trans[0] = '" << trans[0] << "'.  Valid values "
      "include 'N' (No transpose), 'T' (Transpose), and 'C' (Conjugate "
      "transpose).";
    Kokkos::Impl::throw_runtime_exception (os.str ());
  }

  // Minimize the number of Impl::GEMV instantiations, by
  // standardizing on particular View specializations for its template
  // parameters.
  typedef Kokkos::View<typename AViewType::const_value_type**,
    typename AViewType::array_layout,
    typename AViewType::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged> > AVT;
  typedef Kokkos::View<typename XViewType::const_value_type*,
    typename KokkosKernels::Impl::GetUnifiedLayout<XViewType>::array_layout,
    typename XViewType::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged> > XVT;
  typedef Kokkos::View<typename YViewType::non_const_value_type*,
    typename KokkosKernels::Impl::GetUnifiedLayout<YViewType>::array_layout,
    typename YViewType::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged> > YVT;
  typedef Impl::GEMV<AVT, XVT, YVT> impl_type;
  impl_type::gemv (trans, alpha, A, x, beta, y);
}

} // namespace KokkosBlas

#endif // KOKKOS_BLAS2_MV_HPP_
