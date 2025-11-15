// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_CLASSICALDROPPING_HPP
#define MUELU_CLASSICALDROPPING_HPP

#include "MueLu_DroppingCommon.hpp"
#include "Kokkos_Core.hpp"
#if KOKKOS_VERSION >= 40799
#include "KokkosKernels_ArithTraits.hpp"
#else
#include "Kokkos_ArithTraits.hpp"
#endif
#include "Xpetra_Matrix.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu::ClassicalDropping {

/*!
  @class DropFunctor
  @brief Classical dropping criterion

  Depending on the value of measure evaluates the dropping criterion

  SmoothedAggregationMeasure

  \f[
  \frac{|A_{ij}|^2}{|A_{ii}| |A_{jj}|} \le \theta^2
  \f]

  SignedRugeStuebenMeasure

  \f[
  \frac{-\operatorname{Re}A_{ij}}{| max_j -A_{ij}|} \le \theta
  \f]

  SignedSmoothedAggregationMeasure

  \f[
  \frac{-\operatorname{sign}(A_{ij}) |A_{ij}|^2}{|A_{ii}| |A_{jj}|} \le \theta^2
  \f]
*/
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, Misc::StrengthMeasure measure>
class DropFunctor {
 public:
  using matrix_type        = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using diag_vec_type      = Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using local_matrix_type  = typename matrix_type::local_matrix_device_type;
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using memory_space       = typename local_matrix_type::memory_space;
  using diag_view_type     = typename Kokkos::DualView<const scalar_type*, Kokkos::LayoutStride, typename Node::device_type, Kokkos::MemoryUnmanaged>::t_dev;

  using results_view = Kokkos::View<DecisionType*, memory_space>;

#if KOKKOS_VERSION >= 40799
  using ATS = KokkosKernels::ArithTraits<scalar_type>;
#else
  using ATS  = Kokkos::ArithTraits<scalar_type>;
#endif
  using magnitudeType = typename ATS::magnitudeType;
#if KOKKOS_VERSION >= 40799
  using mATS = KokkosKernels::ArithTraits<magnitudeType>;
#else
  using mATS = Kokkos::ArithTraits<magnitudeType>;
#endif
  using boundary_nodes_view = Kokkos::View<const bool*, memory_space>;

 private:
  local_matrix_type A;
  Teuchos::RCP<diag_vec_type> diagVec;
  diag_view_type diag;  // corresponds to overlapped diagonal
  magnitudeType eps;
  results_view results;

 public:
  DropFunctor(matrix_type& A_, magnitudeType threshold, results_view& results_)
    : A(A_.getLocalMatrixDevice())
    , eps(threshold)
    , results(results_) {
    // Construct ghosted matrix diagonal
    if constexpr ((measure == Misc::SmoothedAggregationMeasure) || (measure == Misc::SignedSmoothedAggregationMeasure)) {
      diagVec        = Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetMatrixOverlappedDiagonal(A_);
      auto lclDiag2d = diagVec->getLocalViewDevice(Tpetra::Access::ReadOnly);
      diag           = Kokkos::subview(lclDiag2d, Kokkos::ALL(), 0);
    } else if constexpr (measure == Misc::SignedRugeStuebenMeasure) {
      diagVec    = Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetMatrixMaxMinusOffDiagonal(A_);
      auto lcl2d = diagVec->getLocalViewDevice(Tpetra::Access::ReadOnly);
      diag       = Kokkos::subview(lcl2d, Kokkos::ALL(), 0);
    }
  }

  KOKKOS_FORCEINLINE_FUNCTION
  void operator()(const local_ordinal_type rlid) const {
    auto row      = A.rowConst(rlid);
    size_t offset = A.graph.row_map(rlid);

#ifdef MUELU_COALESCE_DROP_DEBUG
    {
      Kokkos::printf("SoC:        ");
      for (local_ordinal_type k = 0; k < row.length; ++k) {
        auto clid = row.colidx(k);

        auto val = row.value(k);

        if constexpr (measure == Misc::SmoothedAggregationMeasure) {
          auto aiiajj = ATS::magnitude(diag(rlid)) * ATS::magnitude(diag(clid));  // |a_ii|*|a_jj|
          auto aij2   = ATS::magnitude(val) * ATS::magnitude(val);                // |a_ij|^2

          Kokkos::printf("%5f ", ATS::sqrt(aij2 / aiiajj));
        } else if constexpr (measure == Misc::SignedRugeStuebenMeasure) {
          auto neg_aij     = -ATS::real(val);
          auto max_neg_aik = eps * ATS::real(diag(rlid));

          Kokkos::printf("%5f ", neg_aij / max_neg_aik);
        } else if constexpr (measure == Misc::SignedSmoothedAggregationMeasure) {
          auto aiiajj               = ATS::magnitude(diag(rlid)) * ATS::magnitude(diag(clid));  // |a_ii|*|a_jj|
          const bool is_nonpositive = ATS::real(val) <= mATS::zero();
          magnitudeType aij2        = ATS::magnitude(val) * ATS::magnitude(val);  // |a_ij|^2
          // + |a_ij|^2, if a_ij < 0, - |a_ij|^2 if a_ij >=0
          if (!is_nonpositive)
            aij2 = -aij2;
          Kokkos::printf("%5f ", ATS::sqrt(aij2 / aiiajj));
        }
      }
      Kokkos::printf("\n");
    }
#endif

    for (local_ordinal_type k = 0; k < row.length; ++k) {
      auto clid = row.colidx(k);

      auto val = row.value(k);

      if constexpr (measure == Misc::SmoothedAggregationMeasure) {
        auto aiiajj = ATS::magnitude(diag(rlid)) * ATS::magnitude(diag(clid));  // |a_ii|*|a_jj|
        auto aij2   = ATS::magnitude(val) * ATS::magnitude(val);                // |a_ij|^2

        results(offset + k) = Kokkos::max((aij2 <= eps * eps * aiiajj) ? DROP : KEEP,
                                          results(offset + k));
      } else if constexpr (measure == Misc::SignedRugeStuebenMeasure) {
        auto neg_aij        = -ATS::real(val);
        auto max_neg_aik    = eps * ATS::real(diag(rlid));
        results(offset + k) = Kokkos::max((neg_aij < max_neg_aik) ? DROP : KEEP,
                                          results(offset + k));
      } else if constexpr (measure == Misc::SignedSmoothedAggregationMeasure) {
        auto aiiajj               = ATS::magnitude(diag(rlid)) * ATS::magnitude(diag(clid));  // |a_ii|*|a_jj|
        const bool is_nonpositive = ATS::real(val) <= mATS::zero();
        magnitudeType aij2        = ATS::magnitude(val) * ATS::magnitude(val);  // |a_ij|^2
        // + |a_ij|^2, if a_ij < 0, - |a_ij|^2 if a_ij >=0
        if (!is_nonpositive)
          aij2 = -aij2;
        results(offset + k) = Kokkos::max((aij2 <= eps * eps * aiiajj) ? DROP : KEEP,
                                          results(offset + k));
      }
    }
  }
};

// helper function to allow partial template deduction
template <Misc::StrengthMeasure measure, class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
auto make_drop_functor(Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A_,
                       typename DropFunctor<Scalar, LocalOrdinal, GlobalOrdinal, Node, measure>::magnitudeType threshold,
                       typename DropFunctor<Scalar, LocalOrdinal, GlobalOrdinal, Node, measure>::results_view& results_) {
  auto functor = DropFunctor<Scalar, LocalOrdinal, GlobalOrdinal, Node, measure>(A_, threshold, results_);
  return functor;
}

}  // namespace MueLu::ClassicalDropping

#endif
