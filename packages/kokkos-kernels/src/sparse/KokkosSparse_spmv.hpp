/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
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
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
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




#ifndef KOKKOSSPARSE_SPMV_HPP_
#define KOKKOSSPARSE_SPMV_HPP_

#include "KokkosKernels_helpers.hpp"
#include "KokkosKernels_Controls.hpp"
#include "KokkosSparse_spmv_spec.hpp"
#include "KokkosSparse_spmv_struct_spec.hpp"
#include <type_traits>
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosBlas1_scal.hpp"
#include "KokkosKernels_Utils.hpp"

namespace KokkosSparse {

namespace {
  struct RANK_ONE{};
  struct RANK_TWO{};
}

template <class AlphaType, class AMatrix, class XVector, class BetaType, class YVector>
void
spmv (KokkosKernels::Experimental::Controls controls,
      const char mode[],
      const AlphaType& alpha,
      const AMatrix& A,
      const XVector& x,
      const BetaType& beta,
      const YVector& y,
      const RANK_ONE)
{
  // Make sure that both x and y have the same rank.
  static_assert (static_cast<int> (XVector::rank) ==
                 static_cast<int> (YVector::rank),
    "KokkosSparse::spmv: Vector ranks do not match.");
  // Make sure that x (and therefore y) is rank 1.
  static_assert (static_cast<int> (XVector::rank) == 1,
    "KokkosSparse::spmv: Both Vector inputs must have rank 1 "
    "in order to call this specialization of spmv.");
  // Make sure that y is non-const.
  static_assert (std::is_same<typename YVector::value_type,
                   typename YVector::non_const_value_type>::value,
    "KokkosSparse::spmv: Output Vector must be non-const.");

  // Check compatibility of dimensions at run time.
  if ((mode[0] == NoTranspose[0]) || (mode[0] == Conjugate[0])) {
    if ((x.extent(1) != y.extent(1)) ||
        (static_cast<size_t> (A.numCols ()) > static_cast<size_t> (x.extent(0))) ||
        (static_cast<size_t> (A.numRows ()) > static_cast<size_t> (y.extent(0)))) {
      std::ostringstream os;
      os << "KokkosSparse::spmv: Dimensions do not match: "
         << ", A: " << A.numRows () << " x " << A.numCols()
         << ", x: " << x.extent(0) << " x " << x.extent(1)
         << ", y: " << y.extent(0) << " x " << y.extent(1)
         ;
      Kokkos::Impl::throw_runtime_exception (os.str ());
    }
  }
  else {
    if ((x.extent(1) != y.extent(1)) ||
        (static_cast<size_t> (A.numCols ()) > static_cast<size_t> (y.extent(0))) ||
        (static_cast<size_t> (A.numRows ()) > static_cast<size_t> (x.extent(0)))) {
      std::ostringstream os;
      os << "KokkosSparse::spmv: Dimensions do not match (transpose): "
         << ", A: " << A.numRows () << " x " << A.numCols()
         << ", x: " << x.extent(0) << " x " << x.extent(1)
         << ", y: " << y.extent(0) << " x " << y.extent(1)
         ;
      Kokkos::Impl::throw_runtime_exception (os.str ());
    }
  }

  typedef KokkosSparse::CrsMatrix<
              typename AMatrix::const_value_type,
              typename AMatrix::const_ordinal_type,
              typename AMatrix::device_type,
              Kokkos::MemoryTraits<Kokkos::Unmanaged>,
              typename AMatrix::const_size_type>          AMatrix_Internal;

  typedef Kokkos::View<
            typename XVector::const_value_type*,
            typename KokkosKernels::Impl::GetUnifiedLayout<XVector>::array_layout,
            typename XVector::device_type,
            Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> > XVector_Internal;

  typedef Kokkos::View<
            typename YVector::non_const_value_type*,
            typename KokkosKernels::Impl::GetUnifiedLayout<YVector>::array_layout,
            typename YVector::device_type,
            Kokkos::MemoryTraits<Kokkos::Unmanaged> > YVector_Internal;

  AMatrix_Internal A_i = A;
  XVector_Internal x_i = x;
  YVector_Internal y_i = y;

  if(alpha == Kokkos::ArithTraits<AlphaType>::zero() ||
      A_i.numRows() == 0 || A_i.numCols() == 0 || A_i.nnz() == 0)
  {
    //This is required to maintain semantics of KokkosKernels native SpMV:
    //if y contains NaN but beta = 0, the result y should be filled with 0.
    //For example, this is useful for passing in uninitialized y and beta=0.
    if(beta == Kokkos::ArithTraits<BetaType>::zero())
      Kokkos::deep_copy(y_i, Kokkos::ArithTraits<BetaType>::zero());
    else
      KokkosBlas::scal(y_i, beta, y_i);
    return;
  }
  return Impl::SPMV<
    typename AMatrix_Internal::value_type,
    typename AMatrix_Internal::ordinal_type,
    typename AMatrix_Internal::device_type,
    typename AMatrix_Internal::memory_traits,
    typename AMatrix_Internal::size_type,
    typename XVector_Internal::value_type*,
    typename XVector_Internal::array_layout,
    typename XVector_Internal::device_type,
    typename XVector_Internal::memory_traits,
    typename YVector_Internal::value_type*,
    typename YVector_Internal::array_layout,
    typename YVector_Internal::device_type,
    typename YVector_Internal::memory_traits>
      ::spmv (controls, mode, alpha, A_i, x_i, beta, y_i);
}


template<class AlphaType, class AMatrix, class XVector, class BetaType, class YVector ,
         class XLayout = typename XVector::array_layout>
struct SPMV2D1D {
  static bool spmv2d1d (const char mode[],
        const AlphaType& alpha,
        const AMatrix& A,
        const XVector& x,
        const BetaType& beta,
        const YVector& y);
};


template<class AlphaType, class AMatrix, class XVector, class BetaType, class YVector>
struct SPMV2D1D<AlphaType, AMatrix, XVector, BetaType, YVector, Kokkos::LayoutStride>{
  static bool spmv2d1d (const char mode[],
        const AlphaType& alpha,
        const AMatrix& A,
        const XVector& x,
        const BetaType& beta,
        const YVector& y)
  {
#if defined (KOKKOSKERNELS_INST_LAYOUTSTRIDE) || !defined(KOKKOSKERNELS_ETI_ONLY)
    spmv (mode, alpha, A, x, beta, y);
    return true;
#else
    return false;
#endif
  }
};

template<class AlphaType, class AMatrix, class XVector, class BetaType, class YVector>
struct SPMV2D1D<AlphaType, AMatrix, XVector, BetaType, YVector, Kokkos::LayoutLeft>{
  static bool spmv2d1d (const char mode[],
        const AlphaType& alpha,
        const AMatrix& A,
        const XVector& x,
        const BetaType& beta,
        const YVector& y)
  {
#if defined (KOKKOSKERNELS_INST_LAYOUTLEFT) || !defined(KOKKOSKERNELS_ETI_ONLY)
    spmv (mode, alpha, A, x, beta, y);
    return true;
#else
    return false;
#endif
  }
};


template<class AlphaType, class AMatrix, class XVector, class BetaType, class YVector>
struct SPMV2D1D<AlphaType, AMatrix, XVector, BetaType, YVector, Kokkos::LayoutRight>{
  static bool spmv2d1d (const char mode[],
        const AlphaType& alpha,
        const AMatrix& A,
        const XVector& x,
        const BetaType& beta,
        const YVector& y)
  {
#if defined (KOKKOSKERNELS_INST_LAYOUTLEFT) || !defined(KOKKOSKERNELS_ETI_ONLY)
    spmv (mode, alpha, A, x, beta, y);
    return true;
#else
    return false;
#endif
  }
};

template<class AlphaType, class AMatrix, class XVector, class BetaType, class YVector>
void
spmv (KokkosKernels::Experimental::Controls controls,
      const char mode[],
      const AlphaType& alpha,
      const AMatrix& A,
      const XVector& x,
      const BetaType& beta,
      const YVector& y,
      const RANK_TWO)
{
  // Make sure that both x and y have the same rank.
  static_assert (static_cast<int> (XVector::rank) ==
                 static_cast<int> (YVector::rank),
    "KokkosBlas::spmv: Vector ranks do not match.");
  // Make sure that y is non-const.
  static_assert (std::is_same<typename YVector::value_type,
                   typename YVector::non_const_value_type>::value,
    "KokkosBlas::spmv: Output Vector must be non-const.");

  // Check compatibility of dimensions at run time.
  if ((mode[0] == NoTranspose[0]) || (mode[0] == Conjugate[0])) {
    if ((x.extent(1) != y.extent(1)) ||
        (static_cast<size_t> (A.numCols ()) > static_cast<size_t> (x.extent(0))) ||
        (static_cast<size_t> (A.numRows ()) > static_cast<size_t> (y.extent(0)))) {
      std::ostringstream os;
      os << "KokkosBlas::spmv: Dimensions do not match: "
         << ", A: " << A.numRows () << " x " << A.numCols()
         << ", x: " << x.extent(0) << " x " << x.extent(1)
         << ", y: " << y.extent(0) << " x " << y.extent(1);
      Kokkos::Impl::throw_runtime_exception (os.str ());
    }
  } else {
    if ((x.extent(1) != y.extent(1)) ||
        (static_cast<size_t> (A.numCols ()) > static_cast<size_t> (y.extent(0))) ||
        (static_cast<size_t> (A.numRows ()) > static_cast<size_t> (x.extent(0)))) {
      std::ostringstream os;
      os << "KokkosBlas::spmv: Dimensions do not match (transpose): "
         << ", A: " << A.numRows () << " x " << A.numCols()
         << ", x: " << x.extent(0) << " x " << x.extent(1)
         << ", y: " << y.extent(0) << " x " << y.extent(1);
      Kokkos::Impl::throw_runtime_exception (os.str ());
    }
  }

  typedef KokkosSparse::CrsMatrix<
        typename AMatrix::const_value_type,
        typename AMatrix::const_ordinal_type,
        typename AMatrix::device_type,
        Kokkos::MemoryTraits<Kokkos::Unmanaged>,
        typename AMatrix::const_size_type>              AMatrix_Internal;

  AMatrix_Internal A_i = A;

  // Call single-vector version if appropriate
  if (x.extent(1) == 1) {
    typedef Kokkos::View<typename XVector::const_value_type*,
      typename Kokkos::Impl::if_c<std::is_same<typename YVector::array_layout, Kokkos::LayoutLeft>::value,
                                  Kokkos::LayoutLeft, Kokkos::LayoutStride>::type,
      typename XVector::device_type,
      Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> > XVector_SubInternal;
    typedef Kokkos::View<typename YVector::non_const_value_type*,
      typename Kokkos::Impl::if_c<std::is_same<typename YVector::array_layout,Kokkos::LayoutLeft>::value,
                                  Kokkos::LayoutLeft,Kokkos::LayoutStride>::type,
      typename YVector::device_type,
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > YVector_SubInternal;

    XVector_SubInternal x_i = Kokkos::subview (x, Kokkos::ALL (), 0);
    YVector_SubInternal y_i = Kokkos::subview (y, Kokkos::ALL (), 0);

    //spmv (mode, alpha, A, x_i, beta, y_i);
    using impl_type = SPMV2D1D<AlphaType, AMatrix_Internal,
      XVector_SubInternal, BetaType, YVector_SubInternal,
      typename XVector_SubInternal::array_layout>;
    if (impl_type::spmv2d1d (mode, alpha, A, x_i, beta, y_i)) {
      return;
    }
  }
  {
    typedef Kokkos::View<
            typename XVector::const_value_type**,
            typename XVector::array_layout,
            typename XVector::device_type,
            Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> > XVector_Internal;

    typedef Kokkos::View<
              typename YVector::non_const_value_type**,
              typename YVector::array_layout,
              typename YVector::device_type,
              Kokkos::MemoryTraits<Kokkos::Unmanaged> > YVector_Internal;

    XVector_Internal x_i = x;
    YVector_Internal y_i = y;

    return Impl::SPMV_MV<typename AMatrix_Internal::value_type,
                         typename AMatrix_Internal::ordinal_type,
                         typename AMatrix_Internal::device_type,
                         typename AMatrix_Internal::memory_traits,
                         typename AMatrix_Internal::size_type,
                         typename XVector_Internal::value_type**,
                         typename XVector_Internal::array_layout,
                         typename XVector_Internal::device_type,
                         typename XVector_Internal::memory_traits,
                         typename YVector_Internal::value_type**,
                         typename YVector_Internal::array_layout,
                         typename YVector_Internal::device_type,
                         typename YVector_Internal::memory_traits>::spmv_mv (mode, alpha, A_i, x_i, beta, y_i);
  }
}

/// \brief Public interface to local sparse matrix-vector multiply.
///
/// Compute y = beta*y + alpha*Op(A)*x, where x and y are either both
/// rank 1 (single vectors) or rank 2 (multivectors) Kokkos::View
/// instances, A is a KokkosSparse::CrsMatrix, and Op(A) is determined
/// by \c mode.  If beta == 0, ignore and overwrite the initial
/// entries of y; if alpha == 0, ignore the entries of A and x.
///
/// \param controls [in] kokkos-kernels control structure
/// \param mode [in] "N" for no transpose, "T" for transpose, or "C"
///   for conjugate transpose.
/// \param alpha [in] Scalar multiplier for the matrix A.
/// \param A [in] The sparse matrix; KokkosSparse::CrsMatrix instance.
/// \param x [in] Either a single vector (rank-1 Kokkos::View) or
///   multivector (rank-2 Kokkos::View).
/// \param beta [in] Scalar multiplier for the (multi)vector y.
/// \param y [in/out] Either a single vector (rank-1 Kokkos::View) or
///   multivector (rank-2 Kokkos::View).  It must have the same number
///   of columns as x.
template <class AlphaType, class AMatrix, class XVector, class BetaType, class YVector>
void spmv(KokkosKernels::Experimental::Controls controls,
	  const char mode[],
	  const AlphaType& alpha,
	  const AMatrix& A,
	  const XVector& x,
	  const BetaType& beta,
	  const YVector& y) {
  using RANK_SPECIALISE =
    typename std::conditional<static_cast<int> (XVector::rank) == 2,
                              RANK_TWO, RANK_ONE>::type;
  spmv (controls, mode, alpha, A, x, beta, y, RANK_SPECIALISE ());
}

// Overload for backward compatibility and also just simpler
// interface for users that are happy with the kernel default settings
template <class AlphaType, class AMatrix, class XVector, class BetaType, class YVector>
void spmv(const char mode[],
	  const AlphaType& alpha,
	  const AMatrix& A,
	  const XVector& x,
	  const BetaType& beta,
	  const YVector& y) {

  KokkosKernels::Experimental::Controls controls;
  spmv(controls, mode, alpha, A, x, beta, y);

}

  namespace Experimental {

    template <class AlphaType, class AMatrix, class XVector, class BetaType, class YVector>
    void
    spmv_struct (const char mode[],
                 const int stencil_type,
                 const Kokkos::View<typename AMatrix::non_const_ordinal_type*, Kokkos::HostSpace>& structure,
                 const AlphaType& alpha,
                 const AMatrix& A,
                 const XVector& x,
                 const BetaType& beta,
                 const YVector& y,
                 const RANK_ONE)
    {
      // Make sure that both x and y have the same rank.
      static_assert ((int) XVector::rank == (int) YVector::rank,
                     "KokkosSparse::spmv_struct: Vector ranks do not match.");
      // Make sure that x (and therefore y) is rank 1.
      static_assert ((int) XVector::rank == 1,
                     "KokkosSparse::spmv_struct: Both Vector inputs must have rank 1 in "
                     "order to call this specialization of spmv.");
      // Make sure that y is non-const.
      static_assert (std::is_same<typename YVector::value_type,
                     typename YVector::non_const_value_type>::value,
                     "KokkosSparse::spmv_struct: Output Vector must be non-const.");

      // Check compatibility of dimensions at run time.
      if ((mode[0] == NoTranspose[0]) || (mode[0] == Conjugate[0])) {
        if ((x.extent(1) != y.extent(1)) ||
            (static_cast<size_t> (A.numCols ()) > static_cast<size_t> (x.extent(0))) ||
            (static_cast<size_t> (A.numRows ()) > static_cast<size_t> (y.extent(0)))) {
          std::ostringstream os;
          os << "KokkosSparse::spmv_struct: Dimensions do not match: "
             << ", A: " << A.numRows () << " x " << A.numCols()
             << ", x: " << x.extent(0) << " x " << x.extent(1)
             << ", y: " << y.extent(0) << " x " << y.extent(1)
            ;

          Kokkos::Impl::throw_runtime_exception (os.str ());
        }
      }
      else {
        if ((x.extent(1) != y.extent(1)) ||
            (static_cast<size_t> (A.numCols ()) > static_cast<size_t> (y.extent(0))) ||
            (static_cast<size_t> (A.numRows ()) > static_cast<size_t> (x.extent(0)))) {
          std::ostringstream os;
          os << "KokkosSparse::spmv_struct: Dimensions do not match (transpose): "
             << ", A: " << A.numRows () << " x " << A.numCols()
             << ", x: " << x.extent(0) << " x " << x.extent(1)
             << ", y: " << y.extent(0) << " x " << y.extent(1)
            ;

          Kokkos::Impl::throw_runtime_exception (os.str ());
        }
      }


      typedef KokkosSparse::CrsMatrix<
        typename AMatrix::const_value_type,
                                            typename AMatrix::const_ordinal_type,
                                            typename AMatrix::device_type,
                                            Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                                            typename AMatrix::const_size_type>          AMatrix_Internal;

      typedef Kokkos::View<
        typename XVector::const_value_type*,
        typename KokkosKernels::Impl::GetUnifiedLayout<XVector>::array_layout,
        typename XVector::device_type,
        Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> > XVector_Internal;

      typedef Kokkos::View<
        typename YVector::non_const_value_type*,
        typename KokkosKernels::Impl::GetUnifiedLayout<YVector>::array_layout,
        typename YVector::device_type,
        Kokkos::MemoryTraits<Kokkos::Unmanaged> > YVector_Internal;

      AMatrix_Internal A_i = A;
      XVector_Internal x_i = x;
      YVector_Internal y_i = y;

      return Impl::SPMV_STRUCT<
        typename AMatrix_Internal::value_type,
        typename AMatrix_Internal::ordinal_type,
        typename AMatrix_Internal::device_type,
        typename AMatrix_Internal::memory_traits,
        typename AMatrix_Internal::size_type,
        typename XVector_Internal::value_type*,
        typename XVector_Internal::array_layout,
        typename XVector_Internal::device_type,
        typename XVector_Internal::memory_traits,
        typename YVector_Internal::value_type*,
        typename YVector_Internal::array_layout,
        typename YVector_Internal::device_type,
        typename YVector_Internal::memory_traits>::spmv_struct (mode, stencil_type, structure,
                                                                alpha, A_i, x_i, beta, y_i);
    }


    template<class AlphaType, class AMatrix, class XVector, class BetaType, class YVector ,
             class XLayout = typename XVector::array_layout>
    struct SPMV2D1D_STRUCT{
      static bool spmv2d1d_struct (const char mode[],
                                   const int stencil_type,
                                   const Kokkos::View<typename AMatrix::non_const_ordinal_type*, Kokkos::HostSpace>& structure,
                                   const AlphaType& alpha,
                                   const AMatrix& A,
                                   const XVector& x,
                                   const BetaType& beta,
                                   const YVector& y);
    };


    template<class AlphaType, class AMatrix, class XVector, class BetaType, class YVector>
    struct SPMV2D1D_STRUCT<AlphaType, AMatrix, XVector, BetaType, YVector, Kokkos::LayoutStride>{
      static bool spmv2d1d_struct (const char mode[],
                                   const int stencil_type,
                                   const Kokkos::View<typename AMatrix::non_const_ordinal_type*, Kokkos::HostSpace>& structure,
                                   const AlphaType& alpha,
                                   const AMatrix& A,
                                   const XVector& x,
                                   const BetaType& beta,
                                   const YVector& y){
#if defined (KOKKOSKERNELS_INST_LAYOUTSTRIDE) || !defined(KOKKOSKERNELS_ETI_ONLY)
        spmv_struct (mode, stencil_type, structure, alpha, A, x, beta, y, RANK_ONE());
        return true;
#else
        return false;
#endif
      }
    };

    template<class AlphaType, class AMatrix, class XVector, class BetaType, class YVector>
    struct SPMV2D1D_STRUCT<AlphaType, AMatrix, XVector, BetaType, YVector, Kokkos::LayoutLeft>{
      static bool spmv2d1d_struct (const char mode[],
                                   const int stencil_type,
                                   const Kokkos::View<typename AMatrix::non_const_ordinal_type*, Kokkos::HostSpace>& structure,
                                   const AlphaType& alpha,
                                   const AMatrix& A,
                                   const XVector& x,
                                   const BetaType& beta,
                                   const YVector& y){
#if defined (KOKKOSKERNELS_INST_LAYOUTLEFT) || !defined(KOKKOSKERNELS_ETI_ONLY)
        spmv_struct (mode, stencil_type, structure, alpha, A, x, beta, y, RANK_ONE());
        return true;
#else
        return false;
#endif
      }
    };


    template<class AlphaType, class AMatrix, class XVector, class BetaType, class YVector>
    struct SPMV2D1D_STRUCT<AlphaType, AMatrix, XVector, BetaType, YVector, Kokkos::LayoutRight>{
      static bool spmv2d1d_struct (const char mode[],
                                   const int stencil_type,
                                   const Kokkos::View<typename AMatrix::non_const_ordinal_type*, Kokkos::HostSpace>& structure,
                                   const AlphaType& alpha,
                                   const AMatrix& A,
                                   const XVector& x,
                                   const BetaType& beta,
                                   const YVector& y){
#if defined (KOKKOSKERNELS_INST_LAYOUTLEFT) || !defined(KOKKOSKERNELS_ETI_ONLY)
        spmv_struct (mode, stencil_type, structure, alpha, A, x, beta, y, RANK_ONE());
        return true;
#else
        return false;
#endif
      }
    };

    template<class AlphaType, class AMatrix, class XVector, class BetaType, class YVector>
    void
    spmv_struct (const char mode[],
                 const int stencil_type,
                 const Kokkos::View<typename AMatrix::non_const_ordinal_type*, Kokkos::HostSpace>& structure,
                 const AlphaType& alpha,
                 const AMatrix& A,
                 const XVector& x,
                 const BetaType& beta,
                 const YVector& y,
                 const RANK_TWO)
    {
      // Make sure that both x and y have the same rank.
      static_assert (XVector::rank == YVector::rank,
                     "KokkosBlas::spmv: Vector ranks do not match.");
      // Make sure that y is non-const.
      static_assert (std::is_same<typename YVector::value_type,
                     typename YVector::non_const_value_type>::value,
                     "KokkosBlas::spmv: Output Vector must be non-const.");

      // Check compatibility of dimensions at run time.
      if ((mode[0] == NoTranspose[0]) || (mode[0] == Conjugate[0])) {
        if ((x.extent(1) != y.extent(1)) ||
            (static_cast<size_t> (A.numCols ()) > static_cast<size_t> (x.extent(0))) ||
            (static_cast<size_t> (A.numRows ()) > static_cast<size_t> (y.extent(0)))) {
          std::ostringstream os;
          os << "KokkosBlas::spmv: Dimensions do not match: "
             << ", A: " << A.numRows () << " x " << A.numCols()
             << ", x: " << x.extent(0) << " x " << x.extent(1)
             << ", y: " << y.extent(0) << " x " << y.extent(1);
          Kokkos::Impl::throw_runtime_exception (os.str ());
        }
      } else {
        if ((x.extent(1) != y.extent(1)) ||
            (static_cast<size_t> (A.numCols ()) > static_cast<size_t> (y.extent(0))) ||
            (static_cast<size_t> (A.numRows ()) > static_cast<size_t> (x.extent(0)))) {
          std::ostringstream os;
          os << "KokkosBlas::spmv: Dimensions do not match (transpose): "
             << ", A: " << A.numRows () << " x " << A.numCols()
             << ", x: " << x.extent(0) << " x " << x.extent(1)
             << ", y: " << y.extent(0) << " x " << y.extent(1);
          Kokkos::Impl::throw_runtime_exception (os.str ());
        }
      }

      typedef KokkosSparse::CrsMatrix<
        typename AMatrix::const_value_type,
                                            typename AMatrix::const_ordinal_type,
                                            typename AMatrix::device_type,
                                            Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                                            typename AMatrix::const_size_type>              AMatrix_Internal;

      AMatrix_Internal A_i = A;

      // Call single-vector version if appropriate
      if (x.extent(1) == 1) {
        typedef Kokkos::View<typename XVector::const_value_type*,
                             typename Kokkos::Impl::if_c<std::is_same<typename YVector::array_layout, Kokkos::LayoutLeft>::value,
                                                         Kokkos::LayoutLeft, Kokkos::LayoutStride>::type,
                             typename XVector::device_type,
                             Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> > XVector_SubInternal;
        typedef Kokkos::View<typename YVector::non_const_value_type*,
                             typename Kokkos::Impl::if_c<std::is_same<typename YVector::array_layout,Kokkos::LayoutLeft>::value,
                                                         Kokkos::LayoutLeft,Kokkos::LayoutStride>::type,
                             typename YVector::device_type,
                             Kokkos::MemoryTraits<Kokkos::Unmanaged> > YVector_SubInternal;

        XVector_SubInternal x_i = Kokkos::subview (x, Kokkos::ALL (), 0);
        YVector_SubInternal y_i = Kokkos::subview (y, Kokkos::ALL (), 0);



        //spmv_struct (mode, alpha, A, x_i, beta, y_i);
        if (SPMV2D1D_STRUCT  <AlphaType, AMatrix_Internal, XVector_SubInternal,
            BetaType, YVector_SubInternal, typename XVector_SubInternal::array_layout>::spmv2d1d_struct(mode, stencil_type, structure, alpha, A, x_i, beta, y_i)) {
          return;
        }
      }

      // Call true rank 2 vector implementation
      {
        typedef Kokkos::View<
          typename XVector::const_value_type**,
          typename XVector::array_layout,
          typename XVector::device_type,
          Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> > XVector_Internal;

        typedef Kokkos::View<
          typename YVector::non_const_value_type**,
          typename YVector::array_layout,
          typename YVector::device_type,
          Kokkos::MemoryTraits<Kokkos::Unmanaged> > YVector_Internal;

        XVector_Internal x_i = x;
        YVector_Internal y_i = y;

        return Impl::SPMV_MV<typename AMatrix_Internal::value_type,
                             typename AMatrix_Internal::ordinal_type,
                             typename AMatrix_Internal::device_type,
                             typename AMatrix_Internal::memory_traits,
                             typename AMatrix_Internal::size_type,
                             typename XVector_Internal::value_type**,
                             typename XVector_Internal::array_layout,
                             typename XVector_Internal::device_type,
                             typename XVector_Internal::memory_traits,
                             typename YVector_Internal::value_type**,
                             typename YVector_Internal::array_layout,
                             typename YVector_Internal::device_type,
                             typename YVector_Internal::memory_traits>::spmv_mv (mode, alpha, A_i, x_i, beta, y_i);
      }
    }

    /// \brief Public interface to structured local sparse matrix-vector multiply.
    ///
    /// Compute y = beta*y + alpha*Op(A)*x, where x and y are either both
    /// rank 1 (single vectors) or rank 2 (multivectors) Kokkos::View
    /// instances, A is a KokkosSparse::CrsMatrix, and Op(A) is determined
    /// by \c mode.  If beta == 0, ignore and overwrite the initial
    /// entries of y; if alpha == 0, ignore the entries of A and x.
    ///
    /// \param mode [in] "N" for no transpose, "T" for transpose, or "C"
    ///   for conjugate transpose.
    /// \param structure [in] this 1D view stores the # rows in each dimension (i,j,k)
    /// \param alpha [in] Scalar multiplier for the matrix A.
    /// \param A [in] The sparse matrix; KokkosSparse::CrsMatrix instance.
    /// \param x [in] Either a single vector (rank-1 Kokkos::View) or
    ///   multivector (rank-2 Kokkos::View).
    /// \param beta [in] Scalar multiplier for the (multi)vector y.
    /// \param y [in/out] Either a single vector (rank-1 Kokkos::View) or
    ///   multivector (rank-2 Kokkos::View).  It must have the same number
    ///   of columns as x.
    template <class AlphaType, class AMatrix, class XVector, class BetaType, class YVector>
    void
    spmv_struct (const char mode[],
                 const int stencil_type,
                 const Kokkos::View<typename AMatrix::non_const_ordinal_type*, Kokkos::HostSpace>& structure,
                 const AlphaType& alpha,
                 const AMatrix& A,
                 const XVector& x,
                 const BetaType& beta,
                 const YVector& y) {
      typedef typename Kokkos::Impl::if_c<XVector::rank == 2, RANK_TWO, RANK_ONE>::type RANK_SPECIALISE;
      spmv_struct (mode, stencil_type, structure, alpha, A, x, beta, y, RANK_SPECIALISE ());
    }

  } // namespace Experimental
} // namespace KokkosSparse

#endif

