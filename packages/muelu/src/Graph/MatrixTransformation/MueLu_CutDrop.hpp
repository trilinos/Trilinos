// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_CUTDROP_HPP
#define MUELU_CUTDROP_HPP

#include "Kokkos_Core.hpp"
#include "Kokkos_ArithTraits.hpp"
#include "MueLu_DroppingCommon.hpp"
#include "MueLu_Utilities.hpp"
#include "Xpetra_Matrix.hpp"
#include "Xpetra_MultiVector.hpp"
#include "MueLu_DistanceLaplacianDropping.hpp"

namespace MueLu::CutDrop {

enum decisionAlgoType { defaultAlgo,
                        unscaled_cut,
                        scaled_cut,
                        scaled_cut_symmetric };

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class UnscaledComparison {
 public:
  using matrix_type = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  using local_matrix_type  = typename matrix_type::local_matrix_type;
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using memory_space       = typename local_matrix_type::memory_space;
  using results_view       = Kokkos::View<DecisionType*, memory_space>;

  local_matrix_type A;
  results_view results;

 private:
  using ATS           = Kokkos::ArithTraits<scalar_type>;
  using magnitudeType = typename ATS::magnitudeType;

 public:
  UnscaledComparison(matrix_type& A_, results_view& results_)
    : A(A_.getLocalMatrixDevice())
    , results(results_) {}

  template <class local_matrix_type2>
  struct Comparator {
   private:
    using scalar_type        = typename local_matrix_type2::value_type;
    using local_ordinal_type = typename local_matrix_type2::ordinal_type;
    using memory_space       = typename local_matrix_type2::memory_space;
    using results_view       = Kokkos::View<DecisionType*, memory_space>;

    using ATS           = Kokkos::ArithTraits<scalar_type>;
    using magnitudeType = typename ATS::magnitudeType;

    const local_matrix_type2 A;
    const local_ordinal_type offset;
    const results_view results;

   public:
    KOKKOS_INLINE_FUNCTION
    Comparator(const local_matrix_type2& A_, local_ordinal_type rlid_, const results_view& results_)
      : A(A_)
      , offset(A_.graph.row_map(rlid_))
      , results(results_) {}

    KOKKOS_INLINE_FUNCTION
    magnitudeType get_value(size_t x) const {
      return ATS::magnitude(A.values(offset + x) * A.values(offset + x));
    }

    KOKKOS_INLINE_FUNCTION
    bool operator()(size_t x, size_t y) const {
      if (results(offset + x) != UNDECIDED) {
        if (results(offset + y) != UNDECIDED) {
          // does not matter
          return (x < y);
        } else {
          // sort undecided to the right
          return true;
        }
      } else {
        if (results(offset + y) != UNDECIDED) {
          // sort undecided to the right
          return false;
        } else {
          return get_value(x) > get_value(y);
        }
      }
    }
  };

  using comparator_type = Comparator<local_matrix_type>;

  KOKKOS_INLINE_FUNCTION
  comparator_type getComparator(local_ordinal_type rlid) const {
    return comparator_type(A, rlid, results);
  }
};

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class ScaledComparison {
 public:
  using matrix_type        = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using local_matrix_type  = typename matrix_type::local_matrix_type;
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using memory_space       = typename local_matrix_type::memory_space;
  using diag_vec_type      = Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using diag_view_type     = typename Kokkos::DualView<const scalar_type*, Kokkos::LayoutStride, typename Node::device_type, Kokkos::MemoryUnmanaged>::t_dev;
  using results_view       = Kokkos::View<DecisionType*, memory_space>;

  local_matrix_type A;
  results_view results;

 private:
  using ATS           = Kokkos::ArithTraits<scalar_type>;
  using magnitudeType = typename ATS::magnitudeType;

  Teuchos::RCP<diag_vec_type> diagVec;
  diag_view_type diag;

 public:
  ScaledComparison(matrix_type& A_, results_view& results_)
    : A(A_.getLocalMatrixDevice())
    , results(results_) {
    diagVec        = Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetMatrixOverlappedDiagonal(A_);
    auto lclDiag2d = diagVec->getDeviceLocalView(Xpetra::Access::ReadOnly);
    diag           = Kokkos::subview(lclDiag2d, Kokkos::ALL(), 0);
  }

  template <class local_matrix_type2, class diag_view_type2>
  struct Comparator {
   private:
    using scalar_type        = typename local_matrix_type2::value_type;
    using local_ordinal_type = typename local_matrix_type2::ordinal_type;
    using memory_space       = typename local_matrix_type2::memory_space;
    using results_view       = Kokkos::View<DecisionType*, memory_space>;

    using ATS           = Kokkos::ArithTraits<scalar_type>;
    using magnitudeType = typename ATS::magnitudeType;

    const local_matrix_type2 A;
    const diag_view_type2 diag;
    const local_ordinal_type rlid;
    const local_ordinal_type offset;
    const results_view results;

   public:
    KOKKOS_INLINE_FUNCTION
    Comparator(const local_matrix_type2& A_, const diag_view_type2& diag_, const local_ordinal_type rlid_, const results_view& results_)
      : A(A_)
      , diag(diag_)
      , rlid(rlid_)
      , offset(A_.graph.row_map(rlid_))
      , results(results_) {}

    KOKKOS_INLINE_FUNCTION
    magnitudeType get_value(size_t x) const {
      auto x_aij    = ATS::magnitude(A.values(offset + x) * A.values(offset + x));
      auto x_aiiajj = ATS::magnitude(diag(rlid) * diag(A.graph.entries(offset + x)));
      return (x_aij / x_aiiajj);
    }

    KOKKOS_INLINE_FUNCTION
    bool operator()(size_t x, size_t y) const {
      if (results(offset + x) != UNDECIDED) {
        if (results(offset + y) != UNDECIDED) {
          // does not matter
          return (x < y);
        } else {
          // sort undecided to the right
          return true;
        }
      } else {
        if (results(offset + y) != UNDECIDED) {
          // sort undecided to the right
          return false;
        } else {
          return get_value(x) > get_value(y);
        }
      }
    }
  };

  using comparator_type = Comparator<local_matrix_type, diag_view_type>;

  KOKKOS_INLINE_FUNCTION
  comparator_type getComparator(local_ordinal_type rlid) const {
    return comparator_type(A, diag, rlid, results);
  }
};

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class DistanceFunctorType>
class UnscaledDistanceLaplacianComparison {
 public:
  using matrix_type        = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using local_matrix_type  = typename matrix_type::local_matrix_type;
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using memory_space       = typename local_matrix_type::memory_space;
  using diag_vec_type      = Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using diag_view_type     = typename Kokkos::DualView<const scalar_type*, Kokkos::LayoutStride, typename Node::device_type, Kokkos::MemoryUnmanaged>::t_dev;
  using results_view       = Kokkos::View<DecisionType*, memory_space>;

  local_matrix_type A;
  results_view results;

 private:
  using ATS           = Kokkos::ArithTraits<scalar_type>;
  using magnitudeType = typename ATS::magnitudeType;

  Teuchos::RCP<diag_vec_type> diagVec;
  diag_view_type diag;
  DistanceFunctorType dist2;

 public:
  UnscaledDistanceLaplacianComparison(matrix_type& A_, DistanceFunctorType& dist2_, results_view& results_)
    : A(A_.getLocalMatrixDevice())
    , results(results_)
    , dist2(dist2_) {
    // Construct ghosted distance Laplacian diagonal
    diagVec        = DistanceLaplacian::getDiagonal(A_, dist2);
    auto lclDiag2d = diagVec->getDeviceLocalView(Xpetra::Access::ReadOnly);
    diag           = Kokkos::subview(lclDiag2d, Kokkos::ALL(), 0);
  }

  template <class local_matrix_type2, class DistanceFunctorType2, class diag_view_type2>
  struct Comparator {
   private:
    using scalar_type        = typename local_matrix_type2::value_type;
    using local_ordinal_type = typename local_matrix_type2::ordinal_type;
    using memory_space       = typename local_matrix_type2::memory_space;
    using results_view       = Kokkos::View<DecisionType*, memory_space>;

    using ATS           = Kokkos::ArithTraits<scalar_type>;
    using magnitudeType = typename ATS::magnitudeType;

    const local_matrix_type2 A;
    const diag_view_type2 diag;
    const DistanceFunctorType2* dist2;
    const local_ordinal_type rlid;
    const local_ordinal_type offset;
    const results_view results;

    const scalar_type one = ATS::one();

   public:
    KOKKOS_INLINE_FUNCTION
    Comparator(const local_matrix_type2& A_, const diag_view_type2& diag_, const DistanceFunctorType2* dist2_, local_ordinal_type rlid_, const results_view& results_)
      : A(A_)
      , diag(diag_)
      , dist2(dist2_)
      , rlid(rlid_)
      , offset(A_.graph.row_map(rlid_))
      , results(results_) {}

    KOKKOS_INLINE_FUNCTION
    magnitudeType get_value(size_t x) const {
      auto clid = A.graph.entries(offset + x);
      scalar_type val;
      if (rlid != clid) {
        val = one / dist2->distance2(rlid, clid);
      } else {
        val = diag(rlid);
      }
      auto aij2 = ATS::magnitude(val) * ATS::magnitude(val);  // |a_ij|^2
      return aij2;
    }

    KOKKOS_INLINE_FUNCTION
    bool operator()(size_t x, size_t y) const {
      if (results(offset + x) != UNDECIDED) {
        if (results(offset + y) != UNDECIDED) {
          // does not matter
          return (x < y);
        } else {
          // sort undecided to the right
          return true;
        }
      } else {
        if (results(offset + y) != UNDECIDED) {
          // sort undecided to the right
          return false;
        } else {
          return get_value(x) > get_value(y);
        }
      }
    }
  };

  using comparator_type = Comparator<local_matrix_type, DistanceFunctorType, diag_view_type>;

  KOKKOS_INLINE_FUNCTION
  comparator_type getComparator(local_ordinal_type rlid) const {
    return comparator_type(A, diag, &dist2, rlid, results);
  }
};

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class DistanceFunctorType>
class ScaledDistanceLaplacianComparison {
 public:
  using matrix_type        = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using local_matrix_type  = typename matrix_type::local_matrix_type;
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using memory_space       = typename local_matrix_type::memory_space;
  using diag_vec_type      = Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using diag_view_type     = typename Kokkos::DualView<const scalar_type*, Kokkos::LayoutStride, typename Node::device_type, Kokkos::MemoryUnmanaged>::t_dev;
  using results_view       = Kokkos::View<DecisionType*, memory_space>;

  local_matrix_type A;
  results_view results;

 private:
  using ATS           = Kokkos::ArithTraits<scalar_type>;
  using magnitudeType = typename ATS::magnitudeType;

  Teuchos::RCP<diag_vec_type> diagVec;
  diag_view_type diag;
  DistanceFunctorType dist2;

 public:
  ScaledDistanceLaplacianComparison(matrix_type& A_, DistanceFunctorType& dist2_, results_view& results_)
    : A(A_.getLocalMatrixDevice())
    , results(results_)
    , dist2(dist2_) {
    // Construct ghosted distance Laplacian diagonal
    diagVec        = DistanceLaplacian::getDiagonal(A_, dist2);
    auto lclDiag2d = diagVec->getDeviceLocalView(Xpetra::Access::ReadOnly);
    diag           = Kokkos::subview(lclDiag2d, Kokkos::ALL(), 0);
  }

  template <class local_matrix_type2, class DistanceFunctorType2, class diag_view_type2>
  struct Comparator {
   private:
    using scalar_type        = typename local_matrix_type2::value_type;
    using local_ordinal_type = typename local_matrix_type2::ordinal_type;
    using memory_space       = typename local_matrix_type2::memory_space;
    using results_view       = Kokkos::View<DecisionType*, memory_space>;

    using ATS           = Kokkos::ArithTraits<scalar_type>;
    using magnitudeType = typename ATS::magnitudeType;

    const local_matrix_type2 A;
    const diag_view_type2 diag;
    const DistanceFunctorType2* dist2;
    const local_ordinal_type rlid;
    const local_ordinal_type offset;
    const results_view results;

    const scalar_type one = ATS::one();

   public:
    KOKKOS_INLINE_FUNCTION
    Comparator(const local_matrix_type2& A_, const diag_view_type2& diag_, const DistanceFunctorType2* dist2_, local_ordinal_type rlid_, const results_view& results_)
      : A(A_)
      , diag(diag_)
      , dist2(dist2_)
      , rlid(rlid_)
      , offset(A_.graph.row_map(rlid_))
      , results(results_) {}

    KOKKOS_INLINE_FUNCTION
    magnitudeType get_value(size_t x) const {
      auto clid = A.graph.entries(offset + x);
      scalar_type val;
      if (rlid != clid) {
        val = one / dist2->distance2(rlid, clid);
      } else {
        val = diag(rlid);
      }
      auto aiiajj = ATS::magnitude(diag(rlid)) * ATS::magnitude(diag(clid));  // |a_ii|*|a_jj|
      auto aij2   = ATS::magnitude(val) * ATS::magnitude(val);                // |a_ij|^2
      return (aij2 / aiiajj);
    }

    KOKKOS_INLINE_FUNCTION
    bool operator()(size_t x, size_t y) const {
      if (results(offset + x) != UNDECIDED) {
        if (results(offset + y) != UNDECIDED) {
          // does not matter
          return (x < y);
        } else {
          // sort undecided to the right
          return true;
        }
      } else {
        if (results(offset + y) != UNDECIDED) {
          // sort undecided to the right
          return false;
        } else {
          return get_value(x) > get_value(y);
        }
      }
    }
  };

  using comparator_type = Comparator<local_matrix_type, DistanceFunctorType, diag_view_type>;

  KOKKOS_INLINE_FUNCTION
  comparator_type getComparator(local_ordinal_type rlid) const {
    return comparator_type(A, diag, &dist2, rlid, results);
  }
};

template <class view_type, class comparator_type>
KOKKOS_INLINE_FUNCTION void serialHeapSort(view_type& v, comparator_type comparator) {
  auto N       = v.extent(0);
  size_t start = N / 2;
  size_t end   = N;
  while (end > 1) {
    if (start > 0)
      start = start - 1;
    else {
      end       = end - 1;
      auto temp = v(0);
      v(0)      = v(end);
      v(end)    = temp;
    }
    size_t root = start;
    while (2 * root + 1 < end) {
      size_t child = 2 * root + 1;
      if ((child + 1 < end) and (comparator(v(child), v(child + 1))))
        ++child;

      if (comparator(v(root), v(child))) {
        auto temp = v(root);
        v(root)   = v(child);
        v(child)  = temp;
        root      = child;
      } else
        break;
    }
  }
}

template <class comparison_type>
class CutDropFunctor {
 private:
  using local_matrix_type  = typename comparison_type::local_matrix_type;
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using memory_space       = typename local_matrix_type::memory_space;
  using results_view       = Kokkos::View<DecisionType*, memory_space>;

  using ATS                 = Kokkos::ArithTraits<scalar_type>;
  using magnitudeType       = typename ATS::magnitudeType;
  using boundary_nodes_view = Kokkos::View<const bool*, memory_space>;

  local_matrix_type A;
  comparison_type comparison;
  magnitudeType eps;
  results_view results;
  Kokkos::View<local_ordinal_type*, memory_space> index;

 public:
  CutDropFunctor(comparison_type& comparison_, magnitudeType threshold)
    : A(comparison_.A)
    , comparison(comparison_)
    , eps(threshold)
    , results(comparison_.results) {
    index = Kokkos::View<local_ordinal_type*, memory_space>("indices", A.nnz());
  }

  KOKKOS_FORCEINLINE_FUNCTION
  void operator()(const local_ordinal_type& rlid) const {
    auto row   = A.rowConst(rlid);
    size_t nnz = row.length;

    auto drop_view       = Kokkos::subview(results, Kokkos::make_pair(A.graph.row_map(rlid), A.graph.row_map(rlid + 1)));
    auto row_permutation = Kokkos::subview(index, Kokkos::make_pair(A.graph.row_map(rlid), A.graph.row_map(rlid + 1)));

    auto comparator = comparison.getComparator(rlid);

    for (size_t i = 0; i < nnz; ++i) {
      row_permutation(i) = i;
    }
    serialHeapSort(row_permutation, comparator);

    size_t keepStart = 0;
    size_t dropStart = nnz;
    // find index where dropping starts
    for (size_t i = 1; i < nnz; ++i) {
      auto const& x = row_permutation(i - 1);
      auto const& y = row_permutation(i);
      if ((drop_view(x) != UNDECIDED) && (drop_view(y) == UNDECIDED))
        keepStart = i;
      if ((drop_view(x) != UNDECIDED) || (drop_view(y) != UNDECIDED))
        continue;
      magnitudeType x_aij = comparator.get_value(x);
      magnitudeType y_aij = comparator.get_value(y);
      if (eps * eps * x_aij > y_aij) {
        if (i < dropStart) {
          dropStart = i;
        }
      }
    }

    // drop everything to the right of where values stop passing threshold
    for (size_t i = keepStart; i < nnz; ++i) {
      drop_view(row_permutation(i)) = Kokkos::max(dropStart <= i ? DROP : KEEP, drop_view(row_permutation(i)));
    }
  }
};

}  // namespace MueLu::CutDrop

#endif
