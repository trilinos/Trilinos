#ifndef MUELU_DROPPINGCOMMON_HPP
#define MUELU_DROPPINGCOMMON_HPP

#include "Kokkos_Core.hpp"
#include "Kokkos_ArithTraits.hpp"

// For demangling function names
#include <cxxabi.h>

namespace MueLu {

enum DecisionType {
  UNDECIDED = 0,  // no decision has been taken yet, used for initialization
  KEEP      = 1,  // keeep the entry
  DROP      = 2,  // drop it
  BOUNDARY  = 3   // entry is a boundary
};

namespace Misc {

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
  bool operator()(local_ordinal_type rlid) const {
    auto row                 = A.rowConst(rlid);
    const size_t offset      = A.graph.row_map(rlid);
    const bool isBoundaryRow = boundaryNodes(rlid);
    if (isBoundaryRow) {
      for (local_ordinal_type k = 0; k < row.length; ++k) {
        auto clid           = row.colidx(k);
        results(offset + k) = std::max(rlid == clid ? KEEP : DROP,
                                       results(offset + k));
      }
    }
    return false;
  }
};

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
  bool operator()(local_ordinal_type rlid) const {
    auto row                 = A.rowConst(rlid);
    const size_t offset      = A.graph.row_map(rlid);
    const bool isBoundaryRow = boundaryNodes(point_to_block(rlid));
    if (isBoundaryRow) {
      for (local_ordinal_type k = 0; k < row.length; ++k) {
        auto clid           = row.colidx(k);
        results(offset + k) = std::max(rlid == clid ? KEEP : DROP,
                                       results(offset + k));
      }
    }
    return false;
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
      if ((rlid == clid) && (results(offset + k) != BOUNDARY)) {
        results(offset + k) = KEEP;
        break;
      }
    }
    return false;
  }
};

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
  bool operator()(local_ordinal_type rlid) const {
    auto row            = A.rowConst(rlid);
    const size_t offset = A.graph.row_map(rlid);
    for (local_ordinal_type k = 0; k < row.length; ++k) {
      auto clid = row.colidx(k);
      if (clid >= A.numRows()) {
        results(offset + k) = std::max(DROP, results(offset + k));
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
  bool operator()(local_ordinal_type rlid) const {
    auto row            = A.rowConst(rlid);
    const size_t offset = A.graph.row_map(rlid);
    for (local_ordinal_type k = 0; k < row.length; ++k) {
      auto clid = row.colidx(k);
      if ((results(offset + k) == KEEP) && (rlid != clid))
        return false;
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
    return true;
  }
};

template <class local_matrix_type, class block_indices_view_type>
class BlockDiagonalizeFunctor {
 private:
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using memory_space       = typename local_matrix_type::memory_space;
  using results_view       = Kokkos::View<DecisionType*, memory_space>;

  local_matrix_type A;
  block_indices_view_type point_to_block;
  block_indices_view_type ghosted_point_to_block;
  results_view results;

 public:
  BlockDiagonalizeFunctor(local_matrix_type& A_, block_indices_view_type point_to_block_, block_indices_view_type ghosted_point_to_block_, results_view& results_)
    : A(A_)
    , point_to_block(point_to_block_)
    , ghosted_point_to_block(ghosted_point_to_block_)
    , results(results_) {}

  KOKKOS_FORCEINLINE_FUNCTION
  bool operator()(local_ordinal_type rlid) const {
    auto row            = A.rowConst(rlid);
    const size_t offset = A.graph.row_map(rlid);
    for (local_ordinal_type k = 0; k < row.length; ++k) {
      auto clid = row.colidx(k);
      if (point_to_block(rlid, 0) == ghosted_point_to_block(clid, 0)) {
        results(offset + k) = std::max(KEEP, results(offset + k));
      } else {
        results(offset + k) = std::max(DROP, results(offset + k));
      }
    }
    return false;
  }
};

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
  bool operator()(local_ordinal_type rlid) const {
    auto row            = A.rowConst(rlid);
    const size_t offset = A.graph.row_map(rlid);
    for (local_ordinal_type k = 0; k < row.length; ++k) {
      if (results(offset + k) == UNDECIDED)
        assert(false);
    }
    return false;
  }
};

}  // namespace Misc

namespace MatrixConstruction {

template <class local_matrix_type, class FunctorsType, bool debug = false>
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

    if constexpr (debug) {
      Kokkos::printf("\nStarting on row %d\n", rlid);

      auto row = A.rowConst(rlid);

      Kokkos::printf("indices:    ");
      for (local_ordinal_type k = 0; k < row.length; ++k) {
        auto clid = row.colidx(k);
        Kokkos::printf("%5d ", clid);
      }
      Kokkos::printf("\n");

      Kokkos::printf("values:     ");
      for (local_ordinal_type k = 0; k < row.length; ++k) {
        auto val = row.value(k);
        Kokkos::printf("%5f ", val);
      }
      Kokkos::printf("\n");
    }

    if constexpr (std::tuple_size<FunctorsType>::value >= 1) {
      auto functor = std::get<0>(functors);
      doneWithRow  = functor(rlid);
      if constexpr (debug) {
        std::string functorName    = typeid(functor).name();
        int status                 = 0;
        char* demangledFunctorName = 0;
        demangledFunctorName       = abi::__cxa_demangle(functorName.c_str(), 0, 0, &status);
        Kokkos::printf("%s\n", demangledFunctorName);

        auto row            = A.rowConst(rlid);
        const size_t offset = A.graph.row_map(rlid);

        Kokkos::printf("decisions:  ");
        for (local_ordinal_type k = 0; k < row.length; ++k) {
          Kokkos::printf("%5d ", results(offset + k));
        }
        Kokkos::printf("\n");
      }
    }
    if constexpr (std::tuple_size<FunctorsType>::value >= 2) {
      if (!doneWithRow) {
        auto functor = std::get<1>(functors);
        doneWithRow  = functor(rlid);
        if constexpr (debug) {
          std::string functorName    = typeid(functor).name();
          int status                 = 0;
          char* demangledFunctorName = 0;
          demangledFunctorName       = abi::__cxa_demangle(functorName.c_str(), 0, 0, &status);
          Kokkos::printf("%s\n", demangledFunctorName);

          auto row            = A.rowConst(rlid);
          const size_t offset = A.graph.row_map(rlid);

          Kokkos::printf("decisions:  ");
          for (local_ordinal_type k = 0; k < row.length; ++k) {
            Kokkos::printf("%5d ", results(offset + k));
          }
          Kokkos::printf("\n");
        }
      }
    }
    if constexpr (std::tuple_size<FunctorsType>::value >= 3) {
      if (!doneWithRow) {
        auto functor = std::get<2>(functors);
        doneWithRow  = functor(rlid);
        if constexpr (debug) {
          std::string functorName    = typeid(functor).name();
          int status                 = 0;
          char* demangledFunctorName = 0;
          demangledFunctorName       = abi::__cxa_demangle(functorName.c_str(), 0, 0, &status);
          Kokkos::printf("%s\n", demangledFunctorName);

          auto row            = A.rowConst(rlid);
          const size_t offset = A.graph.row_map(rlid);

          Kokkos::printf("decisions:  ");
          for (local_ordinal_type k = 0; k < row.length; ++k) {
            Kokkos::printf("%5d ", results(offset + k));
          }
          Kokkos::printf("\n");
        }
      }
    }
    if constexpr (std::tuple_size<FunctorsType>::value >= 4) {
      if (!doneWithRow) {
        auto functor = std::get<3>(functors);
        doneWithRow  = functor(rlid);
        if constexpr (debug) {
          std::string functorName    = typeid(functor).name();
          int status                 = 0;
          char* demangledFunctorName = 0;
          demangledFunctorName       = abi::__cxa_demangle(functorName.c_str(), 0, 0, &status);
          Kokkos::printf("%s\n", demangledFunctorName);

          auto row            = A.rowConst(rlid);
          const size_t offset = A.graph.row_map(rlid);

          Kokkos::printf("decisions:  ");
          for (local_ordinal_type k = 0; k < row.length; ++k) {
            Kokkos::printf("%5d ", results(offset + k));
          }
          Kokkos::printf("\n");
        }
      }
    }
    if constexpr (std::tuple_size<FunctorsType>::value >= 5) {
      if (!doneWithRow) {
        auto functor = std::get<4>(functors);
        doneWithRow  = functor(rlid);
        if constexpr (debug) {
          std::string functorName    = typeid(functor).name();
          int status                 = 0;
          char* demangledFunctorName = 0;
          demangledFunctorName       = abi::__cxa_demangle(functorName.c_str(), 0, 0, &status);
          Kokkos::printf("%s\n", demangledFunctorName);

          auto row            = A.rowConst(rlid);
          const size_t offset = A.graph.row_map(rlid);

          Kokkos::printf("decisions:  ");
          for (local_ordinal_type k = 0; k < row.length; ++k) {
            Kokkos::printf("%5d ", results(offset + k));
          }
          Kokkos::printf("\n");
        }
      }
    }
    if constexpr (std::tuple_size<FunctorsType>::value >= 6) {
      if (!doneWithRow) {
        auto functor = std::get<5>(functors);
        doneWithRow  = functor(rlid);
        if constexpr (debug) {
          std::string functorName    = typeid(functor).name();
          int status                 = 0;
          char* demangledFunctorName = 0;
          demangledFunctorName       = abi::__cxa_demangle(functorName.c_str(), 0, 0, &status);
          Kokkos::printf("%s\n", demangledFunctorName);

          auto row            = A.rowConst(rlid);
          const size_t offset = A.graph.row_map(rlid);

          Kokkos::printf("decisions:  ");
          for (local_ordinal_type k = 0; k < row.length; ++k) {
            Kokkos::printf("%5d ", results(offset + k));
          }
          Kokkos::printf("\n");
        }
      }
    }
    if constexpr (std::tuple_size<FunctorsType>::value >= 7) {
      if (!doneWithRow) {
        auto functor = std::get<6>(functors);
        doneWithRow  = functor(rlid);
        if constexpr (debug) {
          std::string functorName    = typeid(functor).name();
          int status                 = 0;
          char* demangledFunctorName = 0;
          demangledFunctorName       = abi::__cxa_demangle(functorName.c_str(), 0, 0, &status);
          Kokkos::printf("%s\n", demangledFunctorName);

          auto row            = A.rowConst(rlid);
          const size_t offset = A.graph.row_map(rlid);

          Kokkos::printf("decisions:  ");
          for (local_ordinal_type k = 0; k < row.length; ++k) {
            Kokkos::printf("%5d ", results(offset + k));
          }
          Kokkos::printf("\n");
        }
      }
    }

    if constexpr (debug) {
      Kokkos::printf("Done with row %d\n", rlid);
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

template <class local_matrix_type, class local_graph_type, bool lumping>
class PointwiseFillReuseFunctor {
 private:
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using memory_space       = typename local_matrix_type::memory_space;
  using results_view       = Kokkos::View<DecisionType*, memory_space>;
  using ATS                = Kokkos::ArithTraits<scalar_type>;

  local_matrix_type A;
  results_view results;
  local_matrix_type filteredA;
  local_graph_type graph;
  const scalar_type zero = ATS::zero();

 public:
  PointwiseFillReuseFunctor(local_matrix_type& A_, results_view& results_, local_matrix_type& filteredA_, local_graph_type& graph_)
    : A(A_)
    , results(results_)
    , filteredA(filteredA_)
    , graph(graph_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const local_ordinal_type rlid) const {
    auto rowA                       = A.row(rlid);
    size_t K                        = A.graph.row_map(rlid);
    auto rowFilteredA               = filteredA.row(rlid);
    local_ordinal_type j            = 0;
    local_ordinal_type jj           = 0;
    local_ordinal_type graph_offset = graph.row_map(rlid);
    scalar_type diagCorrection      = zero;
    local_ordinal_type diagOffset   = -1;
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
        graph.entries(graph_offset + jj) = rowA.colidx(k);
        ++jj;
      } else if constexpr (lumping) {
        diagCorrection += rowA.value(k);
        rowFilteredA.colidx(j) = rowA.colidx(k);
        rowFilteredA.value(j)  = zero;
        ++j;
      } else {
        rowFilteredA.colidx(j) = rowA.colidx(k);
        rowFilteredA.value(j)  = zero;
        ++j;
      }
    }
    if constexpr (lumping) {
      rowFilteredA.value(diagOffset) += diagCorrection;
    }
  }
};

template <class local_matrix_type, bool lumping>
class PointwiseFillNoReuseFunctor {
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
  PointwiseFillNoReuseFunctor(local_matrix_type& A_, results_view& results_, local_matrix_type& filteredA_)
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
      rowFilteredA.value(diagOffset) += diagCorrection;
    }
  }
};

template <class local_matrix_type, class FunctorsType, bool debug = false>
class VectorCountingFunctor {
 private:
  using scalar_type             = typename local_matrix_type::value_type;
  using local_ordinal_type      = typename local_matrix_type::ordinal_type;
  using memory_space            = typename local_matrix_type::memory_space;
  using results_view            = Kokkos::View<DecisionType*, memory_space>;
  using block_indices_view_type = Kokkos::View<local_ordinal_type*, memory_space>;

  using rowptr_type = typename local_matrix_type::row_map_type::non_const_type;
  using ATS         = Kokkos::ArithTraits<local_ordinal_type>;

  local_matrix_type A;
  local_ordinal_type blockSize;
  block_indices_view_type ghosted_point_to_block;
  results_view results;
  rowptr_type filtered_rowptr;
  rowptr_type graph_rowptr;
  FunctorsType functors;

 public:
  VectorCountingFunctor(local_matrix_type& A_, local_ordinal_type blockSize_, block_indices_view_type ghosted_point_to_block_, results_view& results_, rowptr_type& filtered_rowptr_, rowptr_type& graph_rowptr_, FunctorsType& functors_)
    : A(A_)
    , blockSize(blockSize_)
    , ghosted_point_to_block(ghosted_point_to_block_)
    , results(results_)
    , filtered_rowptr(filtered_rowptr_)
    , graph_rowptr(graph_rowptr_)
    , functors(functors_) {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const local_ordinal_type brlid, std::tuple<local_ordinal_type, local_ordinal_type>& nnz, const bool& final) const {
    auto nnz_filtered = &std::get<0>(nnz);
    auto nnz_graph    = &std::get<1>(nnz);

    if constexpr (debug) {
      Kokkos::printf("\nStarting on block row %d\n", brlid);
    }

    for (local_ordinal_type rlid = blockSize * brlid; rlid < blockSize * (brlid + 1); ++rlid) {
      if constexpr (debug) {
        Kokkos::printf("\nStarting on row %d\n", rlid);

        auto row = A.rowConst(rlid);

        Kokkos::printf("indices:    ");
        for (local_ordinal_type k = 0; k < row.length; ++k) {
          auto clid = row.colidx(k);
          Kokkos::printf("%5d ", clid);
        }
        Kokkos::printf("\n");

        Kokkos::printf("values:     ");
        for (local_ordinal_type k = 0; k < row.length; ++k) {
          auto val = row.value(k);
          Kokkos::printf("%5f ", val);
        }
        Kokkos::printf("\n");
      }

      bool doneWithRow = false;
      if constexpr (std::tuple_size<FunctorsType>::value >= 1) {
        auto functor = std::get<0>(functors);
        doneWithRow  = functor(rlid);
        if constexpr (debug) {
          std::string functorName    = typeid(functor).name();
          int status                 = 0;
          char* demangledFunctorName = 0;
          demangledFunctorName       = abi::__cxa_demangle(functorName.c_str(), 0, 0, &status);
          Kokkos::printf("%s\n", demangledFunctorName);

          auto row            = A.rowConst(rlid);
          const size_t offset = A.graph.row_map(rlid);

          Kokkos::printf("decisions:  ");
          for (local_ordinal_type k = 0; k < row.length; ++k) {
            Kokkos::printf("%5d ", results(offset + k));
          }
          Kokkos::printf("\n");
        }
      }
      if constexpr (std::tuple_size<FunctorsType>::value >= 2) {
        if (!doneWithRow) {
          auto functor = std::get<1>(functors);
          doneWithRow  = functor(rlid);
          if constexpr (debug) {
            std::string functorName    = typeid(functor).name();
            int status                 = 0;
            char* demangledFunctorName = 0;
            demangledFunctorName       = abi::__cxa_demangle(functorName.c_str(), 0, 0, &status);
            Kokkos::printf("%s\n", demangledFunctorName);

            auto row            = A.rowConst(rlid);
            const size_t offset = A.graph.row_map(rlid);

            Kokkos::printf("decisions:  ");
            for (local_ordinal_type k = 0; k < row.length; ++k) {
              Kokkos::printf("%5d ", results(offset + k));
            }
            Kokkos::printf("\n");
          }
        }
      }
      if constexpr (std::tuple_size<FunctorsType>::value >= 3) {
        if (!doneWithRow) {
          auto functor = std::get<2>(functors);
          doneWithRow  = functor(rlid);
          if constexpr (debug) {
            std::string functorName    = typeid(functor).name();
            int status                 = 0;
            char* demangledFunctorName = 0;
            demangledFunctorName       = abi::__cxa_demangle(functorName.c_str(), 0, 0, &status);
            Kokkos::printf("%s\n", demangledFunctorName);

            auto row            = A.rowConst(rlid);
            const size_t offset = A.graph.row_map(rlid);

            Kokkos::printf("decisions:  ");
            for (local_ordinal_type k = 0; k < row.length; ++k) {
              Kokkos::printf("%5d ", results(offset + k));
            }
            Kokkos::printf("\n");
          }
        }
      }
      if constexpr (std::tuple_size<FunctorsType>::value >= 4) {
        if (!doneWithRow) {
          auto functor = std::get<3>(functors);
          doneWithRow  = functor(rlid);
          if constexpr (debug) {
            std::string functorName    = typeid(functor).name();
            int status                 = 0;
            char* demangledFunctorName = 0;
            demangledFunctorName       = abi::__cxa_demangle(functorName.c_str(), 0, 0, &status);
            Kokkos::printf("%s\n", demangledFunctorName);

            auto row            = A.rowConst(rlid);
            const size_t offset = A.graph.row_map(rlid);

            Kokkos::printf("decisions:  ");
            for (local_ordinal_type k = 0; k < row.length; ++k) {
              Kokkos::printf("%5d ", results(offset + k));
            }
            Kokkos::printf("\n");
          }
        }
      }
      if constexpr (std::tuple_size<FunctorsType>::value >= 5) {
        if (!doneWithRow) {
          auto functor = std::get<4>(functors);
          doneWithRow  = functor(rlid);
          if constexpr (debug) {
            std::string functorName    = typeid(functor).name();
            int status                 = 0;
            char* demangledFunctorName = 0;
            demangledFunctorName       = abi::__cxa_demangle(functorName.c_str(), 0, 0, &status);
            Kokkos::printf("%s\n", demangledFunctorName);

            auto row            = A.rowConst(rlid);
            const size_t offset = A.graph.row_map(rlid);

            Kokkos::printf("decisions:  ");
            for (local_ordinal_type k = 0; k < row.length; ++k) {
              Kokkos::printf("%5d ", results(offset + k));
            }
            Kokkos::printf("\n");
          }
        }
      }
      if constexpr (std::tuple_size<FunctorsType>::value >= 6) {
        if (!doneWithRow) {
          auto functor = std::get<5>(functors);
          doneWithRow  = functor(rlid);
          if constexpr (debug) {
            std::string functorName    = typeid(functor).name();
            int status                 = 0;
            char* demangledFunctorName = 0;
            demangledFunctorName       = abi::__cxa_demangle(functorName.c_str(), 0, 0, &status);
            Kokkos::printf("%s\n", demangledFunctorName);

            auto row            = A.rowConst(rlid);
            const size_t offset = A.graph.row_map(rlid);

            Kokkos::printf("decisions:  ");
            for (local_ordinal_type k = 0; k < row.length; ++k) {
              Kokkos::printf("%5d ", results(offset + k));
            }
            Kokkos::printf("\n");
          }
        }
      }
      if constexpr (std::tuple_size<FunctorsType>::value >= 7) {
        if (!doneWithRow) {
          auto functor = std::get<6>(functors);
          doneWithRow  = functor(rlid);
          if constexpr (debug) {
            std::string functorName    = typeid(functor).name();
            int status                 = 0;
            char* demangledFunctorName = 0;
            demangledFunctorName       = abi::__cxa_demangle(functorName.c_str(), 0, 0, &status);
            Kokkos::printf("%s\n", demangledFunctorName);

            auto row            = A.rowConst(rlid);
            const size_t offset = A.graph.row_map(rlid);

            Kokkos::printf("decisions:  ");
            for (local_ordinal_type k = 0; k < row.length; ++k) {
              Kokkos::printf("%5d ", results(offset + k));
            }
            Kokkos::printf("\n");
          }
        }
      }

      if constexpr (debug) {
        Kokkos::printf("Done with row %d\n", rlid);
      }

      size_t start = A.graph.row_map(rlid);
      size_t end   = A.graph.row_map(rlid + 1);
      for (size_t i = start; i < end; ++i) {
        if (results(i) == KEEP) {
          ++(*nnz_filtered);
        }
      }
      if (final)
        filtered_rowptr(rlid + 1) = *nnz_filtered;
    }

    if constexpr (debug) {
      Kokkos::printf("Done with block row %d\nGraph indices ", brlid);
    }

    local_ordinal_type* nextIndices = new local_ordinal_type[blockSize];
    for (local_ordinal_type block_index = 0; block_index < blockSize; ++block_index) {
      nextIndices[block_index] = 0;
    }
    local_ordinal_type prev_bclid = -1;
    while (true) {
      local_ordinal_type min_block_index = -1;
      local_ordinal_type min_clid        = ATS::max();
      local_ordinal_type min_offset      = -1;
      for (local_ordinal_type block_index = 0; block_index < blockSize; ++block_index) {
        auto rlid   = blockSize * brlid + block_index;
        auto offset = A.graph.row_map(rlid) + nextIndices[block_index];
        if (offset == A.graph.row_map(rlid + 1))
          continue;
        auto clid = A.graph.entries(offset);
        if (clid < min_clid) {
          min_block_index = block_index;
          min_clid        = clid;
          min_offset      = offset;
        }
      }
      if (min_block_index == -1)
        break;
      ++nextIndices[min_block_index];
      auto bclid = ghosted_point_to_block(min_clid);
      if (prev_bclid < bclid) {
        if (results(min_offset) == KEEP) {
          ++(*nnz_graph);
          if constexpr (debug) {
            Kokkos::printf("%5d ", bclid);
          }
          prev_bclid = bclid;
        }
      }
    }
    if constexpr (debug) {
      Kokkos::printf("\n");
    }
    if (final)
      graph_rowptr(brlid + 1) = *nnz_graph;
  }
};

template <class local_matrix_type, bool lumping, bool reuse>
class VectorFillFunctor {
 private:
  using scalar_type             = typename local_matrix_type::value_type;
  using local_ordinal_type      = typename local_matrix_type::ordinal_type;
  using local_graph_type        = typename local_matrix_type::staticcrsgraph_type;
  using memory_space            = typename local_matrix_type::memory_space;
  using results_view            = Kokkos::View<DecisionType*, memory_space>;
  using ATS                     = Kokkos::ArithTraits<scalar_type>;
  using OTS                     = Kokkos::ArithTraits<local_ordinal_type>;
  using block_indices_view_type = Kokkos::View<local_ordinal_type*, memory_space>;

  local_matrix_type A;
  local_ordinal_type blockSize;
  block_indices_view_type ghosted_point_to_block;
  results_view results;
  local_matrix_type filteredA;
  local_graph_type graph;
  const scalar_type zero = ATS::zero();

 public:
  VectorFillFunctor(local_matrix_type& A_, local_ordinal_type blockSize_, block_indices_view_type ghosted_point_to_block_, results_view& results_, local_matrix_type& filteredA_, local_graph_type& graph_)
    : A(A_)
    , blockSize(blockSize_)
    , ghosted_point_to_block(ghosted_point_to_block_)
    , results(results_)
    , filteredA(filteredA_)
    , graph(graph_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const local_ordinal_type brlid) const {
    for (local_ordinal_type rlid = blockSize * brlid; rlid < blockSize * (brlid + 1); ++rlid) {
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
          if constexpr (reuse) {
            rowFilteredA.colidx(j) = rowA.colidx(k);
            rowFilteredA.value(j)  = zero;
            ++j;
          }
        } else if constexpr (reuse) {
          rowFilteredA.colidx(j) = rowA.colidx(k);
          rowFilteredA.value(j)  = zero;
          ++j;
        }
      }
      if constexpr (lumping) {
        rowFilteredA.value(diagOffset) += diagCorrection;
      }
    }

    local_ordinal_type* nextIndices = new local_ordinal_type[blockSize];
    for (local_ordinal_type block_index = 0; block_index < blockSize; ++block_index) {
      nextIndices[block_index] = 0;
    }
    local_ordinal_type prev_bclid = -1;

    local_ordinal_type j = graph.row_map(brlid);
    while (true) {
      local_ordinal_type min_block_index = -1;
      local_ordinal_type min_clid        = OTS::max();
      local_ordinal_type min_offset      = -1;
      for (local_ordinal_type block_index = 0; block_index < blockSize; ++block_index) {
        auto rlid   = blockSize * brlid + block_index;
        auto offset = A.graph.row_map(rlid) + nextIndices[block_index];
        if (offset == A.graph.row_map(rlid + 1))
          continue;
        auto clid = A.graph.entries(offset);
        if (clid < min_clid) {
          min_block_index = block_index;
          min_clid        = clid;
          min_offset      = offset;
        }
      }
      if (min_block_index == -1)
        break;
      ++nextIndices[min_block_index];
      auto bclid = ghosted_point_to_block(min_clid);
      if (prev_bclid < bclid) {
        if (results(min_offset) == KEEP) {
          graph.entries(j) = bclid;
          ++j;
          prev_bclid = bclid;
        }
      }
    }
  }
};

}  // namespace MatrixConstruction

}  // namespace MueLu

#endif
