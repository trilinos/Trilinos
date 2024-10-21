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
#include "Kokkos_ArithTraits.hpp"
#include "Xpetra_Matrix.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu::ClassicalDropping {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class AbsDropFunctor {
 private:
  using matrix_type        = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using diag_vec_type      = Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using local_matrix_type  = typename matrix_type::local_matrix_type;
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using memory_space       = typename local_matrix_type::memory_space;
  using diag_view_type     = typename Kokkos::DualView<const scalar_type*, Kokkos::LayoutStride, typename Node::device_type, Kokkos::MemoryUnmanaged>::t_dev;

  using results_view = Kokkos::View<DecisionType*, memory_space>;

  using ATS                 = Kokkos::ArithTraits<scalar_type>;
  using magnitudeType       = typename ATS::magnitudeType;
  using boundary_nodes_view = Kokkos::View<const bool*, memory_space>;

  local_matrix_type A;
  Teuchos::RCP<diag_vec_type> diagVec;
  diag_view_type diag;  // corresponds to overlapped diagonal
  magnitudeType eps;
  results_view results;

 public:
  AbsDropFunctor(matrix_type& A_, magnitudeType threshold, results_view& results_)
    : A(A_.getLocalMatrixDevice())
    , eps(threshold)
    , results(results_) {
    diagVec        = Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetMatrixOverlappedDiagonal(A_);
    auto lclDiag2d = diagVec->getDeviceLocalView(Xpetra::Access::ReadOnly);
    diag           = Kokkos::subview(lclDiag2d, Kokkos::ALL(), 0);
  }

  KOKKOS_FORCEINLINE_FUNCTION
  void operator()(const local_ordinal_type rlid) const {
    auto row      = A.rowConst(rlid);
    size_t offset = A.graph.row_map(rlid);
    for (local_ordinal_type k = 0; k < row.length; ++k) {
      auto clid = row.colidx(k);

      auto val    = row.value(k);
      auto aiiajj = ATS::magnitude(diag(rlid)) * ATS::magnitude(diag(clid));  // |a_ii|*|a_jj|
      auto aij2   = ATS::magnitude(val) * ATS::magnitude(val);                // |a_ij|^2

      results(offset + k) = Kokkos::max((aij2 <= eps * eps * aiiajj) ? DROP : KEEP,
                                        results(offset + k));
    }
  }
};

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class SignedClassicalRSDropFunctor {
 private:
  using matrix_type        = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using local_matrix_type  = typename matrix_type::local_matrix_type;
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using memory_space       = typename local_matrix_type::memory_space;

  using results_view = Kokkos::View<DecisionType*, memory_space>;

  using ATS                 = Kokkos::ArithTraits<scalar_type>;
  using magnitudeType       = typename ATS::magnitudeType;
  using boundary_nodes_view = Kokkos::View<const bool*, memory_space>;

  using diag_vec_type  = Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using diag_view_type = typename Kokkos::DualView<const scalar_type*, Kokkos::LayoutStride, typename Node::device_type, Kokkos::MemoryUnmanaged>::t_dev;

  local_matrix_type A;
  Teuchos::RCP<diag_vec_type> diagVec;
  diag_view_type diag;  // corresponds to overlapped diagonal
  magnitudeType eps;
  results_view results;

 public:
  SignedClassicalRSDropFunctor(matrix_type& A_, magnitudeType threshold, results_view& results_)
    : A(A_.getLocalMatrixDevice())
    , eps(threshold)
    , results(results_) {
    diagVec        = Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetMatrixMaxMinusOffDiagonal(A_);
    auto lclDiag2d = diagVec->getDeviceLocalView(Xpetra::Access::ReadOnly);
    diag           = Kokkos::subview(lclDiag2d, Kokkos::ALL(), 0);
  }

  KOKKOS_FORCEINLINE_FUNCTION
  void operator()(const local_ordinal_type rlid) const {
    auto row      = A.rowConst(rlid);
    size_t offset = A.graph.row_map(rlid);
    for (local_ordinal_type k = 0; k < row.length; ++k) {
      auto val            = row.value(k);
      auto neg_aij        = -ATS::real(val);
      auto max_neg_aik    = eps * ATS::real(diag(rlid));
      results(offset + k) = Kokkos::max((neg_aij <= max_neg_aik) ? DROP : KEEP,
                                        results(offset + k));
    }
  }
};

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class SignedClassicalSADropFunctor {
 private:
  using matrix_type        = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using diag_vec_type      = Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using local_matrix_type  = typename matrix_type::local_matrix_type;
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using memory_space       = typename local_matrix_type::memory_space;
  using diag_view_type     = typename Kokkos::DualView<const scalar_type*, Kokkos::LayoutStride, typename Node::device_type, Kokkos::MemoryUnmanaged>::t_dev;

  using results_view = Kokkos::View<DecisionType*, memory_space>;

  using ATS                 = Kokkos::ArithTraits<scalar_type>;
  using magnitudeType       = typename ATS::magnitudeType;
  using mATS                = Kokkos::ArithTraits<magnitudeType>;
  using boundary_nodes_view = Kokkos::View<const bool*, memory_space>;

  local_matrix_type A;
  Teuchos::RCP<diag_vec_type> diagVec;
  diag_view_type diag;  // corresponds to overlapped diagonal
  magnitudeType eps;
  results_view results;

 public:
  SignedClassicalSADropFunctor(matrix_type& A_, magnitudeType threshold, results_view& results_)
    : A(A_.getLocalMatrixDevice())
    , eps(threshold)
    , results(results_) {
    // Construct ghosted matrix diagonal
    diagVec        = Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetMatrixOverlappedDiagonal(A_);
    auto lclDiag2d = diagVec->getDeviceLocalView(Xpetra::Access::ReadOnly);
    diag           = Kokkos::subview(lclDiag2d, Kokkos::ALL(), 0);
  }

  KOKKOS_FORCEINLINE_FUNCTION
  void operator()(const local_ordinal_type rlid) const {
    auto row      = A.rowConst(rlid);
    size_t offset = A.graph.row_map(rlid);
    for (local_ordinal_type k = 0; k < row.length; ++k) {
      auto clid = row.colidx(k);

      auto val                  = row.value(k);
      auto aiiajj               = ATS::magnitude(diag(rlid)) * ATS::magnitude(diag(clid));  // |a_ii|*|a_jj|
      const bool is_nonpositive = ATS::real(val) <= mATS::zero();
      magnitudeType aij2        = ATS::magnitude(val) * ATS::magnitude(val);  // |a_ij|^2
      // + |a_ij|^2, if a_ij < 0, - |a_ij|^2 if a_ij >=0
      if (is_nonpositive)
        aij2 = -aij2;
      results(offset + k) = Kokkos::max((aij2 <= eps * eps * aiiajj) ? DROP : KEEP,
                                        results(offset + k));
    }
  }
};

}  // namespace MueLu::ClassicalDropping

#endif
