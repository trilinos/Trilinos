#ifndef MUELU_BOUNDARYDETECTION_HPP
#define MUELU_BOUNDARYDETECTION_HPP

#include <cstddef>
#include <tuple>
#include "Kokkos_Core.hpp"
#include "Kokkos_ArithTraits.hpp"
#include "MueLu_LWGraph_kokkos.hpp"
#include "MueLu_Utilities.hpp"
#include "Teuchos_RCP.hpp"
#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_CrsGraph.hpp"
#include "Xpetra_MultiVector.hpp"

namespace MueLu::BoundaryDetection {

// These functors all assume that the boundaryNodes view has been initialized to false.

// Marks rows as Dirichlet based on value threshold and number of off-diagonal entries.
template <class local_matrix_type>
class PointDirichletFunctor {
 private:
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using memory_space       = typename local_matrix_type::memory_space;

  using ATS                 = Kokkos::ArithTraits<scalar_type>;
  using magnitudeType       = typename ATS::magnitudeType;
  using boundary_nodes_view = Kokkos::View<bool*, memory_space>;

  local_matrix_type A;
  boundary_nodes_view boundaryNodes;
  magnitudeType dirichletThreshold;
  local_ordinal_type dirichletNonzeroThreshold;

 public:
  PointDirichletFunctor(local_matrix_type& A_, boundary_nodes_view boundaryNodes_, magnitudeType dirichletThreshold_, local_ordinal_type dirichletNonzeroThreshold_)
    : A(A_)
    , boundaryNodes(boundaryNodes_)
    , dirichletThreshold(dirichletThreshold_)
    , dirichletNonzeroThreshold(dirichletNonzeroThreshold_) {}

  KOKKOS_FORCEINLINE_FUNCTION
  void operator()(const local_ordinal_type rlid) const {
    auto row               = A.row(rlid);
    local_ordinal_type nnz = 0;
    for (local_ordinal_type k = 0; k < row.length; ++k) {
      auto clid = row.colidx(k);
      auto val  = row.value(k);
      if ((rlid != static_cast<local_ordinal_type>(clid)) && (ATS::magnitude(val) > dirichletThreshold)) {
        ++nnz;
        if (nnz == dirichletNonzeroThreshold) {
          return;
        }
      }
    }
    boundaryNodes(rlid) = true;
  }
};

// Marks rows as Dirichlet based on abs(rowsum) and abs(diag).
template <class local_matrix_type>
class RowSumFunctor {
 private:
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using memory_space       = typename local_matrix_type::memory_space;

  using ATS                 = Kokkos::ArithTraits<scalar_type>;
  using magnitudeType       = typename ATS::magnitudeType;
  using magATS              = Kokkos::ArithTraits<magnitudeType>;
  using boundary_nodes_view = Kokkos::View<bool*, memory_space>;

  local_matrix_type A;
  boundary_nodes_view boundaryNodes;
  magnitudeType rowSumTol;

 public:
  RowSumFunctor(local_matrix_type& A_, boundary_nodes_view boundaryNodes_, magnitudeType rowSumTol_)
    : A(A_)
    , boundaryNodes(boundaryNodes_)
    , rowSumTol(rowSumTol_) {}

  KOKKOS_FORCEINLINE_FUNCTION
  void operator()(const local_ordinal_type rlid) const {
    scalar_type rowsum  = ATS::zero();
    scalar_type diagval = ATS::zero();
    auto row            = A.row(rlid);
    for (local_ordinal_type k = 0; k < row.length; ++k) {
      auto clid = row.colidx(k);
      auto val  = row.value(k);
      if (rlid == static_cast<local_ordinal_type>(clid))
        diagval = val;
      rowsum += val;
    }
    if (ATS::magnitude(rowsum) > ATS::magnitude(diagval) * rowSumTol) {
      boundaryNodes(rlid) = true;
    }
  }
};

template <class local_matrix_type, class FunctorsType>
class BoundaryFunctor {
 private:
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using memory_space       = typename local_matrix_type::memory_space;

  local_matrix_type A;
  FunctorsType functors;

 public:
  BoundaryFunctor(local_matrix_type& A_, FunctorsType& functors_)
    : A(A_)
    , functors(functors_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const local_ordinal_type rlid) const {
    if constexpr (std::tuple_size<FunctorsType>::value >= 1) {
      auto functor = std::get<0>(functors);
      functor(rlid);
    }
    if constexpr (std::tuple_size<FunctorsType>::value >= 2) {
      auto functor = std::get<1>(functors);
      functor(rlid);
    }
    if constexpr (std::tuple_size<FunctorsType>::value >= 3) {
      auto functor = std::get<2>(functors);
      functor(rlid);
    }
    if constexpr (std::tuple_size<FunctorsType>::value >= 4) {
      auto functor = std::get<3>(functors);
      functor(rlid);
    }
  }
};

// Marks rows as Dirichlet based on value threshold and number of off-diagonal entries.
// Marks blocks as Dirichlet when one row is Dirichlet (useGreedyDirichlet==true) or when all rows are Dirichlet (useGreedyDirichlet==false).
template <class local_matrix_type, bool useGreedyDirichlet>
class VectorDirichletFunctor {
 private:
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using memory_space       = typename local_matrix_type::memory_space;

  using ATS                 = Kokkos::ArithTraits<scalar_type>;
  using magnitudeType       = typename ATS::magnitudeType;
  using boundary_nodes_view = Kokkos::View<bool*, memory_space>;

  local_matrix_type A;
  local_ordinal_type blockSize;
  boundary_nodes_view boundaryNodes;
  magnitudeType dirichletThreshold;
  local_ordinal_type dirichletNonzeroThreshold;

 public:
  VectorDirichletFunctor(local_matrix_type& A_, local_ordinal_type blockSize_, boundary_nodes_view boundaryNodes_, magnitudeType dirichletThreshold_, local_ordinal_type dirichletNonzeroThreshold_)
    : A(A_)
    , blockSize(blockSize_)
    , boundaryNodes(boundaryNodes_)
    , dirichletThreshold(dirichletThreshold_)
    , dirichletNonzeroThreshold(dirichletNonzeroThreshold_) {}

  KOKKOS_FORCEINLINE_FUNCTION
  void operator()(const local_ordinal_type rblid) const {
    for (local_ordinal_type rlid = rblid * blockSize; rlid < (rblid + 1) * blockSize; ++rlid) {
      auto row               = A.row(rlid);
      local_ordinal_type nnz = 0;
      bool rowIsDirichlet    = true;
      for (local_ordinal_type k = 0; k < row.length; ++k) {
        auto clid = row.colidx(k);
        auto val  = row.value(k);
        if ((rlid != static_cast<local_ordinal_type>(clid)) && (ATS::magnitude(val) > dirichletThreshold)) {
          ++nnz;
          if (nnz == dirichletNonzeroThreshold) {
            rowIsDirichlet = false;
            break;
          }
        }
      }
      if constexpr (useGreedyDirichlet) {
        if (rowIsDirichlet) {
          boundaryNodes(rblid) = true;
          return;
        }
      } else {
        if (!rowIsDirichlet) {
          boundaryNodes(rblid) = false;
          return;
        }
      }
    }
    if constexpr (useGreedyDirichlet)
      boundaryNodes(rblid) = false;
    else
      boundaryNodes(rblid) = true;
  }
};

}  // namespace MueLu::BoundaryDetection

#endif
