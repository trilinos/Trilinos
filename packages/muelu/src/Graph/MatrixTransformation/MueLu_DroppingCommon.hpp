#ifndef MUELU_DROPPINGCOMMON_HPP
#define MUELU_DROPPINGCOMMON_HPP

#include "Kokkos_Core.hpp"
#include "Kokkos_ArithTraits.hpp"

namespace MueLu {

enum DecisionType {
  UNDECIDED = 0,
  KEEP      = 1,
  DROP      = 2,
  BOUNDARY  = 3
};

namespace Misc {

template <class local_matrix_type>
class DropBoundaryFunctor {
 private:
  using scalar_type         = typename local_matrix_type::value_type;
  using local_ordinal_type  = typename local_matrix_type::ordinal_type;
  using memory_space        = typename local_matrix_type::memory_space;
  using results_view        = Kokkos::View<DecisionType*, memory_space>;
  using boundary_nodes_view = Kokkos::View<bool*, memory_space>;

  local_matrix_type A;
  boundary_nodes_view boundaryNodes;
  results_view results;

 public:
  DropBoundaryFunctor(local_matrix_type& A_, boundary_nodes_view boundaryNodes_, results_view& results_)
    : A(A_)
    , boundaryNodes(boundaryNodes_)
    , results(results_) {}

  KOKKOS_FORCEINLINE_FUNCTION
  bool operator()(local_ordinal_type rlid) const {
    auto row                 = A.rowConst(rlid);
    const size_t offset      = A.graph.row_map(rlid);
    const bool isBoundaryRow = boundaryNodes(rlid);
    if (isBoundaryRow) {
      for (local_ordinal_type k = 0; k < row.length; ++k) {
        results(offset + k) = KEEP;
      }
    }
    return isBoundaryRow;
  }
};

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
  bool operator()(local_ordinal_type rlid) const {
    auto row            = A.rowConst(rlid);
    const size_t offset = A.graph.row_map(rlid);
    for (local_ordinal_type k = 0; k < row.length; ++k) {
      auto clid = row.colidx(k);
      if (rlid == clid) {
        results(offset + k) = KEEP;
        break;
      }
    }
    return false;
  }
};

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
  bool operator()(local_ordinal_type rlid) const {
    auto row            = A.rowConst(rlid);
    const size_t offset = A.graph.row_map(rlid);
    for (local_ordinal_type k = 0; k < row.length; ++k) {
      auto clid = row.colidx(k);
      if ((results(offset + k) == KEEP) && (rlid != clid))
        return false;
    }
    boundaryNodes(rlid) = true;
    for (local_ordinal_type k = 0; k < row.length; ++k) {
      auto clid = row.colidx(k);
      if (rlid == clid)
        results(offset + k) = KEEP;
      else
        results(offset + k) = BOUNDARY;
    }
    return true;
  }
};

}  // namespace Misc

namespace MatrixConstruction {

template <class local_matrix_type, class FunctorsType>
class PointwiseCountingFunctor {
 private:
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using memory_space       = typename local_matrix_type::memory_space;
  using results_view       = Kokkos::View<DecisionType*, memory_space>;

  using rowptr_type = typename local_matrix_type::row_map_type::non_const_type;

  local_matrix_type A;
  results_view results;
  rowptr_type rowptr;
  FunctorsType functors;

 public:
  PointwiseCountingFunctor(local_matrix_type& A_, results_view& results_, rowptr_type& rowptr_, FunctorsType& functors_)
    : A(A_)
    , results(results_)
    , rowptr(rowptr_)
    , functors(functors_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const local_ordinal_type rlid, local_ordinal_type& nnz, const bool& final) const {
    bool doneWithRow = false;
    if constexpr (std::tuple_size<FunctorsType>::value >= 1) {
      auto functor = std::get<0>(functors);
      doneWithRow  = functor(rlid);
    }
    if constexpr (std::tuple_size<FunctorsType>::value >= 2) {
      if (!doneWithRow) {
        auto functor = std::get<1>(functors);
        doneWithRow  = functor(rlid);
      }
    }
    if constexpr (std::tuple_size<FunctorsType>::value >= 3) {
      if (!doneWithRow) {
        auto functor = std::get<2>(functors);
        doneWithRow  = functor(rlid);
      }
    }
    if constexpr (std::tuple_size<FunctorsType>::value >= 4) {
      if (!doneWithRow) {
        auto functor = std::get<3>(functors);
        doneWithRow  = functor(rlid);
      }
    }
    if constexpr (std::tuple_size<FunctorsType>::value >= 5) {
      if (!doneWithRow) {
        auto functor = std::get<4>(functors);
        doneWithRow  = functor(rlid);
      }
    }
    if constexpr (std::tuple_size<FunctorsType>::value >= 6) {
      if (!doneWithRow) {
        auto functor = std::get<5>(functors);
        doneWithRow  = functor(rlid);
      }
    }
    if constexpr (std::tuple_size<FunctorsType>::value >= 7) {
      if (!doneWithRow) {
        auto functor = std::get<6>(functors);
        doneWithRow  = functor(rlid);
      }
    }

    size_t start = A.graph.row_map(rlid);
    size_t end   = A.graph.row_map(rlid + 1);
    for (size_t i = start; i < end; ++i) {
      if (results(i) == KEEP) {
        ++nnz;
      }
    }
    if (final)
      rowptr(rlid + 1) = nnz;
  }
};

template <class local_matrix_type, bool lumping>
class PointwiseFillFunctor {
 private:
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using memory_space       = typename local_matrix_type::memory_space;
  using results_view       = Kokkos::View<DecisionType*, memory_space>;
  using ATS                = Kokkos::ArithTraits<scalar_type>;

  local_matrix_type A;
  results_view results;
  local_matrix_type filteredA;
  const scalar_type zero = ATS::zero();

 public:
  PointwiseFillFunctor(local_matrix_type& A_, results_view& results_, local_matrix_type& filteredA_)
    : A(A_)
    , results(results_)
    , filteredA(filteredA_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const local_ordinal_type rlid) const {
    auto rowA                     = A.row(rlid);
    size_t K                      = A.graph.row_map(rlid);
    auto rowFilteredA             = filteredA.row(rlid);
    local_ordinal_type j          = 0;
    scalar_type diagCorrection    = zero;
    local_ordinal_type diagOffset = -1;
    for (local_ordinal_type k = 0; k < rowA.length; ++k) {
      if constexpr (lumping) {
        local_ordinal_type clid = rowA.colidx(k);
        if (rlid == clid) {
          diagOffset = j;
        }
      }
      if (results(K + k) == KEEP) {
        rowFilteredA.colidx(j) = rowA.colidx(k);
        rowFilteredA.value(j)  = rowA.value(k);
        ++j;
      } else if constexpr (lumping) {
        diagCorrection += rowA.value(k);
      }
    }
    if constexpr (lumping) {
      rowA.value(diagOffset) += diagCorrection;
    }
  }
};

}  // namespace MatrixConstruction

}  // namespace MueLu

#endif
