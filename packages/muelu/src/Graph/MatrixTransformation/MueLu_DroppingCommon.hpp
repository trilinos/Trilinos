// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_DROPPINGCOMMON_HPP
#define MUELU_DROPPINGCOMMON_HPP

#include "Kokkos_Core.hpp"
#include "Kokkos_ArithTraits.hpp"
#include "Xpetra_Access.hpp"
#include "Xpetra_Matrix.hpp"

namespace MueLu {

/*! Possible decision for a single entry.
  Once we are done with dropping, we should have no UNDECIDED entries left.
  Normally, both DROP and BOUNDARY entries will be dropped, but we distinguish them in case we want to keep boundaries.
 */
enum DecisionType : char {
  UNDECIDED = 0,  // no decision has been taken yet, used for initialization
  KEEP      = 1,  // keeep the entry
  DROP      = 2,  // drop it
  BOUNDARY  = 3   // entry is a boundary
};

namespace Misc {

template <class local_ordinal_type>
class NoOpFunctor {
 public:
  NoOpFunctor() {}

  KOKKOS_FORCEINLINE_FUNCTION
  void operator()(local_ordinal_type rlid) const {
  }
};

/*!
  @class PointwiseDropBoundaryFunctor
  @brief Functor that drops boundary nodes for a blockSize == 1 problem.
*/
template <class local_matrix_type>
class PointwiseDropBoundaryFunctor {
 private:
  using scalar_type         = typename local_matrix_type::value_type;
  using local_ordinal_type  = typename local_matrix_type::ordinal_type;
  using memory_space        = typename local_matrix_type::memory_space;
  using results_view        = Kokkos::View<DecisionType*, memory_space>;
  using boundary_nodes_view = Kokkos::View<const bool*, memory_space>;

  local_matrix_type A;
  boundary_nodes_view boundaryNodes;
  results_view results;

 public:
  PointwiseDropBoundaryFunctor(local_matrix_type& A_, boundary_nodes_view boundaryNodes_, results_view& results_)
    : A(A_)
    , boundaryNodes(boundaryNodes_)
    , results(results_) {}

  KOKKOS_FORCEINLINE_FUNCTION
  void operator()(local_ordinal_type rlid) const {
    auto row                 = A.rowConst(rlid);
    const size_t offset      = A.graph.row_map(rlid);
    const bool isBoundaryRow = boundaryNodes(rlid);
    if (isBoundaryRow) {
      for (local_ordinal_type k = 0; k < row.length; ++k) {
        auto clid           = row.colidx(k);
        results(offset + k) = Kokkos::max(rlid == clid ? KEEP : DROP,
                                          results(offset + k));
      }
    }
  }
};

/*!
  @class VectorDropBoundaryFunctor
  @brief Functor that drops boundary nodes for a blockSize > 1 problem.
*/
template <class local_matrix_type>
class VectorDropBoundaryFunctor {
 private:
  using scalar_type             = typename local_matrix_type::value_type;
  using local_ordinal_type      = typename local_matrix_type::ordinal_type;
  using memory_space            = typename local_matrix_type::memory_space;
  using results_view            = Kokkos::View<DecisionType*, memory_space>;
  using boundary_nodes_view     = Kokkos::View<const bool*, memory_space>;
  using block_indices_view_type = Kokkos::View<local_ordinal_type*, memory_space>;

  local_matrix_type A;
  block_indices_view_type point_to_block;
  boundary_nodes_view boundaryNodes;
  results_view results;

 public:
  VectorDropBoundaryFunctor(local_matrix_type& A_, block_indices_view_type point_to_block_, boundary_nodes_view boundaryNodes_, results_view& results_)
    : A(A_)
    , point_to_block(point_to_block_)
    , boundaryNodes(boundaryNodes_)
    , results(results_) {}

  KOKKOS_FORCEINLINE_FUNCTION
  void operator()(local_ordinal_type rlid) const {
    auto row                 = A.rowConst(rlid);
    const size_t offset      = A.graph.row_map(rlid);
    const bool isBoundaryRow = boundaryNodes(point_to_block(rlid));
    if (isBoundaryRow) {
      for (local_ordinal_type k = 0; k < row.length; ++k) {
        auto clid           = row.colidx(k);
        results(offset + k) = Kokkos::max(rlid == clid ? KEEP : DROP,
                                          results(offset + k));
      }
    }
  }
};

/*!
@class KeepDiagonalFunctor
@brief Functor that marks diagonal as kept, unless the are already marked as boundary.
*/
template <class local_matrix_type>
class KeepDiagonalFunctor {
 private:
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using memory_space       = typename local_matrix_type::memory_space;
  using results_view       = Kokkos::View<DecisionType*, memory_space>;

  local_matrix_type A;
  results_view results;

 public:
  KeepDiagonalFunctor(local_matrix_type& A_, results_view& results_)
    : A(A_)
    , results(results_) {}

  KOKKOS_FORCEINLINE_FUNCTION
  void operator()(local_ordinal_type rlid) const {
    auto row            = A.rowConst(rlid);
    const size_t offset = A.graph.row_map(rlid);
    for (local_ordinal_type k = 0; k < row.length; ++k) {
      auto clid = row.colidx(k);
      if ((rlid == clid) && (results(offset + k) != BOUNDARY)) {
        results(offset + k) = KEEP;
        break;
      }
    }
  }
};

/*!
@class DropOffRankFunctor
@brief Functor that drops off-rank entries
*/
template <class local_matrix_type>
class DropOffRankFunctor {
 private:
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using memory_space       = typename local_matrix_type::memory_space;
  using results_view       = Kokkos::View<DecisionType*, memory_space>;

  local_matrix_type A;
  results_view results;

 public:
  DropOffRankFunctor(local_matrix_type& A_, results_view& results_)
    : A(A_)
    , results(results_) {}

  KOKKOS_FORCEINLINE_FUNCTION
  void operator()(local_ordinal_type rlid) const {
    auto row            = A.rowConst(rlid);
    const size_t offset = A.graph.row_map(rlid);
    for (local_ordinal_type k = 0; k < row.length; ++k) {
      auto clid = row.colidx(k);
      if (clid >= A.numRows()) {
        results(offset + k) = Kokkos::max(DROP, results(offset + k));
      }
    }
  }
};

/*!
@class MarkSingletonFunctor
@brief Functor that marks singletons (all off-diagonal entries in a row are dropped) as boundary.
*/
template <class local_matrix_type>
class MarkSingletonFunctor {
 private:
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using memory_space       = typename local_matrix_type::memory_space;
  using results_view       = Kokkos::View<DecisionType*, memory_space>;

  using boundary_nodes_view = Kokkos::View<bool*, memory_space>;

  local_matrix_type A;
  boundary_nodes_view boundaryNodes;
  results_view results;

 public:
  MarkSingletonFunctor(local_matrix_type& A_, boundary_nodes_view boundaryNodes_, results_view& results_)
    : A(A_)
    , boundaryNodes(boundaryNodes_)
    , results(results_) {}

  KOKKOS_FORCEINLINE_FUNCTION
  void operator()(local_ordinal_type rlid) const {
    auto row            = A.rowConst(rlid);
    const size_t offset = A.graph.row_map(rlid);
    for (local_ordinal_type k = 0; k < row.length; ++k) {
      auto clid = row.colidx(k);
      if ((results(offset + k) == KEEP) && (rlid != clid))
        return;
    }
    boundaryNodes(rlid) = true;
    for (local_ordinal_type k = 0; k < row.length; ++k) {
      auto clid = row.colidx(k);
      if (rlid == clid)
        results(offset + k) = KEEP;
      else
        results(offset + k) = BOUNDARY;
    }
  }
};

/*!
@class MarkSingletonVectorFunctor
@brief Functor that marks singletons (all off-diagonal entries in a row are dropped) as boundary.
*/
template <class local_matrix_type>
class MarkSingletonVectorFunctor {
 private:
  using scalar_type             = typename local_matrix_type::value_type;
  using local_ordinal_type      = typename local_matrix_type::ordinal_type;
  using memory_space            = typename local_matrix_type::memory_space;
  using results_view            = Kokkos::View<DecisionType*, memory_space>;
  using block_indices_view_type = Kokkos::View<local_ordinal_type*, memory_space>;

  using boundary_nodes_view = Kokkos::View<bool*, memory_space>;

  local_matrix_type A;
  block_indices_view_type point_to_block;
  boundary_nodes_view boundaryNodes;
  results_view results;

 public:
  MarkSingletonVectorFunctor(local_matrix_type& A_, block_indices_view_type point_to_block_, boundary_nodes_view boundaryNodes_, results_view& results_)
    : A(A_)
    , point_to_block(point_to_block_)
    , boundaryNodes(boundaryNodes_)
    , results(results_) {}

  KOKKOS_FORCEINLINE_FUNCTION
  void operator()(local_ordinal_type rlid) const {
    auto row            = A.rowConst(rlid);
    const size_t offset = A.graph.row_map(rlid);
    for (local_ordinal_type k = 0; k < row.length; ++k) {
      auto clid = row.colidx(k);
      if ((results(offset + k) == KEEP) && (rlid != clid))
        return;
    }
    auto brlid           = point_to_block(rlid);
    boundaryNodes(brlid) = true;
    for (local_ordinal_type k = 0; k < row.length; ++k) {
      auto clid = row.colidx(k);
      if (rlid == clid)
        results(offset + k) = KEEP;
      else
        results(offset + k) = BOUNDARY;
    }
  }
};

/*!
@class BlockDiagonalizeFunctor
@brief Functor that drops all entries that are not on the block diagonal.
*/
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class BlockDiagonalizeFunctor {
 private:
  using matrix_type       = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using local_matrix_type = typename matrix_type::local_matrix_type;

  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using memory_space       = typename local_matrix_type::memory_space;
  using results_view       = Kokkos::View<DecisionType*, memory_space>;

  using block_indices_type            = Xpetra::MultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>;
  using local_block_indices_view_type = typename block_indices_type::dual_view_type_const::t_dev;

  local_matrix_type A;
  local_block_indices_view_type point_to_block;
  local_block_indices_view_type ghosted_point_to_block;
  results_view results;

 public:
  BlockDiagonalizeFunctor(matrix_type& A_, block_indices_type& point_to_block_, block_indices_type& ghosted_point_to_block_, results_view& results_)
    : A(A_.getLocalMatrixDevice())
    , point_to_block(point_to_block_.getDeviceLocalView(Xpetra::Access::ReadOnly))
    , ghosted_point_to_block(ghosted_point_to_block_.getDeviceLocalView(Xpetra::Access::ReadOnly))
    , results(results_) {}

  KOKKOS_FORCEINLINE_FUNCTION
  void operator()(local_ordinal_type rlid) const {
    auto row            = A.rowConst(rlid);
    const size_t offset = A.graph.row_map(rlid);
    for (local_ordinal_type k = 0; k < row.length; ++k) {
      auto clid = row.colidx(k);
      if (point_to_block(rlid, 0) == ghosted_point_to_block(clid, 0)) {
        results(offset + k) = Kokkos::max(KEEP, results(offset + k));
      } else {
        results(offset + k) = Kokkos::max(DROP, results(offset + k));
      }
    }
  }
};

/*!
@class DebugFunctor
@brief Functor that checks that all entries have been marked.
*/
template <class local_matrix_type>
class DebugFunctor {
 private:
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using memory_space       = typename local_matrix_type::memory_space;
  using results_view       = Kokkos::View<DecisionType*, memory_space>;

  using boundary_nodes_view = Kokkos::View<bool*, memory_space>;

  local_matrix_type A;
  results_view results;

 public:
  DebugFunctor(local_matrix_type& A_, results_view& results_)
    : A(A_)
    , results(results_) {}

  KOKKOS_FORCEINLINE_FUNCTION
  void operator()(local_ordinal_type rlid) const {
    auto row            = A.rowConst(rlid);
    const size_t offset = A.graph.row_map(rlid);
    for (local_ordinal_type k = 0; k < row.length; ++k) {
      if (results(offset + k) == UNDECIDED) {
        Kokkos::printf("No dropping decision was taken for entry (%d, %d)\n", rlid, row.colidx(k));
        assert(false);
      }
    }
  }
};

/*!
@class SymmetrizeFunctor
@brief Functor that symmetrizes the dropping decisions.
*/
template <class local_matrix_type>
class SymmetrizeFunctor {
 private:
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using memory_space       = typename local_matrix_type::memory_space;
  using results_view       = Kokkos::View<DecisionType*, memory_space>;

  local_matrix_type A;
  results_view results;

 public:
  SymmetrizeFunctor(local_matrix_type& A_, results_view& results_)
    : A(A_)
    , results(results_) {}

  KOKKOS_FORCEINLINE_FUNCTION
  void operator()(local_ordinal_type rlid) const {
    auto row            = A.rowConst(rlid);
    const size_t offset = A.graph.row_map(rlid);
    for (local_ordinal_type k = 0; k < row.length; ++k) {
      if (results(offset + k) == KEEP) {
        auto clid = row.colidx(k);
        if (clid >= A.numRows())
          continue;
        auto row2            = A.rowConst(clid);
        const size_t offset2 = A.graph.row_map(clid);
        for (local_ordinal_type k2 = 0; k2 < row2.length; ++k2) {
          auto clid2 = row2.colidx(k2);
          if (clid2 == rlid) {
            if (results(offset2 + k2) == DROP)
              results(offset2 + k2) = KEEP;
            break;
          }
        }
      }
    }
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

}  // namespace Misc

}  // namespace MueLu

#endif
