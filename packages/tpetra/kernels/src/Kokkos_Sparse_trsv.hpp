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

#ifndef KOKKOS_SPARSE_TRSV_HPP_
#define KOKKOS_SPARSE_TRSV_HPP_

#include <Kokkos_Sparse_impl_trsm.hpp>
#ifdef KOKKOS_HAVE_CXX11
#include <type_traits>
#endif // KOKKOS_HAVE_CXX11

namespace KokkosSparse {

/// \brief Solve the triangular sparse linear system Op(A) x = b.
///
/// \tparam AMatrix KokkosSparse::CrsMatrix specialization.
/// \tparam BMV The type of the input (right-hand side) (multi)vector.
/// \tparam XMV The type of the output (left-hand side) (multi)vector.
///
/// \param uplo [in] "U" (for upper triangular) or "L" (for lower
///   triangular).
/// \param trans [in] "C" (for conjugate transpose), "T" (for
///   transpose), or "N" (for no transpose).
/// \param diag [in] "U" (for implicit unit diagonal) or "N" (for
///   not).
/// \param A [in] The input matrix A; must be upper triangular or
///   lower triangular.
/// \param b [in] The input (right-hand side) (multi)vector.
/// \param x [in] The output (left-hand side) (multi)vector.
template <class AMatrix, class BMV, class XMV>
void
trsv (const char uplo[],
      const char trans[],
      const char diag[],
      const AMatrix& A,
      const BMV& b,
      const XMV& x)
{
  // FIXME (mfh 23 Apr 2015) Need to implement rank-1 version of this function.
  static_assert (BMV::rank == 2, "KokkosBlas::trsv: Rank-1 version of this "
                 "function has not yet been implemented.");

  static_assert (Kokkos::Impl::is_view<BMV>::value,
                 "KokkosBlas::trsv: b is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<XMV>::value,
                 "KokkosBlas::trsv: x is not a Kokkos::View.");
  static_assert ((int) BMV::rank == (int) XMV::rank,
                 "KokkosBlas::trsv: The ranks of b and x do not match.");
  static_assert (BMV::rank == 1 || BMV::rank == 2,
                 "KokkosBlas::trsv: b and x must both either have rank 1, or rank 2.");
  static_assert (Kokkos::Impl::is_same<typename XMV::value_type,
                 typename XMV::non_const_value_type>::value,
                 "KokkosBlas::trsv: The output x must be nonconst.");

  if (uplo[0] != 'U' && uplo[0] != 'u' && uplo[0] != 'L' && uplo[0] != 'l') {
    std::ostringstream os;
    os << "Invalid uplo[0] = \'" << uplo << "\'";
    Kokkos::Impl::throw_runtime_exception (os.str ());
  }
  if (trans[0] != 'C' && trans[0] != 'c' && trans[0] != 'T' && trans[0] != 't' && trans[0] != 'N' && trans[0] != 'n') {
    std::ostringstream os;
    os << "Invalid trans[0] = \'" << trans << "\'";
    Kokkos::Impl::throw_runtime_exception (os.str ());
  }
  if (diag[0] != 'U' && diag[0] != 'u' && diag[0] != 'N' && diag[0] != 'n') {
    std::ostringstream os;
    os << "Invalid diag[0] = \'" << diag << "\'";
    Kokkos::Impl::throw_runtime_exception (os.str ());
  }

  typedef typename BMV::size_type size_type;
  const size_type numRows = static_cast<size_type> (A.numRows ());
  const size_type numCols = static_cast<size_type> (A.numCols ());

  const bool transpose = ! trans[0] == 'N' && ! trans[0] == 'n';
  if (! transpose && (numCols != x.dimension_0 () || numRows != b.dimension_0 ())) {
    std::ostringstream os;
    os << "Dimensions do not match (non-transpose case).  "
       << "A is " << numRows << " x " << numCols
       << ", x is " << x.dimension_0 () << " x " << x.dimension_1 ()
       << ", and b is " << b.dimension_0 () << " x " << b.dimension_1 ();
    Kokkos::Impl::throw_runtime_exception (os.str ());
  }
  if (transpose && (numRows != x.dimension_0 () || numCols != b.dimension_0 ())) {
    std::ostringstream os;
    os << "Dimensions do not match (transpose or conjugate transpose case).  "
       << "A is " << numRows << " x " << numCols
       << ", x is " << x.dimension_0 () << " x " << x.dimension_1 ()
       << ", and b is " << b.dimension_0 () << " x " << b.dimension_1 ();
    Kokkos::Impl::throw_runtime_exception (os.str ());
  }

  typedef KokkosSparse::CrsMatrix<typename AMatrix::const_value_type,
    typename AMatrix::const_ordinal_type,
    typename AMatrix::device_type,
    typename AMatrix::memory_traits,
    typename AMatrix::const_size_type> AMatrix_Internal;
  AMatrix_Internal A_i = A;

  using KokkosSparse::Impl::Sequential::Trsv;
  Trsv<AMatrix_Internal, XMV, BMV>::trsv (uplo, trans, diag, A, b, x);
}

} // namespace KokkosSparse

#endif // KOKKOS_SPARSE_TRSV_HPP_

