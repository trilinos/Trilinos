// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_BOUNDARYDETECTION_HPP
#define MUELU_BOUNDARYDETECTION_HPP

#include <cstddef>
#include <type_traits>
#include "Kokkos_Core.hpp"
#if KOKKOS_VERSION >= 40799
#include "KokkosKernels_ArithTraits.hpp"
#else
#include "Kokkos_ArithTraits.hpp"
#endif
#include "MueLu_LWGraph_kokkos.hpp"
#include "MueLu_Utilities.hpp"
#include "Teuchos_RCP.hpp"
#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_CrsGraph.hpp"
#include "Xpetra_MultiVector.hpp"

namespace MueLu::BoundaryDetection {

/*!
  @class PointDirichletFunctor
  @brief Functor for marking nodes as Dirichlet.

  A row is marked as Dirichlet boundary if fewer than dirichletNonzeroThreshold entries are larger in absolute value than dirichletThreshold.
  It is assumed that boundaryNodes was initialized to false.
*/
template <class local_matrix_type>
class PointDirichletFunctor {
 private:
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using memory_space       = typename local_matrix_type::memory_space;

#if KOKKOS_VERSION >= 40799
  using ATS = KokkosKernels::ArithTraits<scalar_type>;
#else
  using ATS    = Kokkos::ArithTraits<scalar_type>;
#endif
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
    auto row               = A.rowConst(rlid);
    local_ordinal_type nnz = 0;
    for (local_ordinal_type k = 0; k < row.length; ++k) {
      local_ordinal_type clid = row.colidx(k);
      scalar_type val         = row.value(k);
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

/*!
  @class VectorDirichletFunctor
  @brief Functor for marking nodes as Dirichlet in a block operator.

  Assumes a single fixed block size specified by blockSize.
  Marks blocks as Dirichlet when one row is Dirichlet (useGreedyDirichlet==true) or when all rows are Dirichlet (useGreedyDirichlet==false).
  A row is marked as Dirichlet boundary if fewer than dirichletNonzeroThreshold entries are larger in absolute value than dirichletThreshold.
  It is assumed that boundaryNodes was initialized to false.
*/
template <class local_matrix_type, bool useGreedyDirichlet>
class VectorDirichletFunctor {
 private:
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using memory_space       = typename local_matrix_type::memory_space;

#if KOKKOS_VERSION >= 40799
  using ATS = KokkosKernels::ArithTraits<scalar_type>;
#else
  using ATS    = Kokkos::ArithTraits<scalar_type>;
#endif
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
      auto row               = A.rowConst(rlid);
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
          return;
        }
      }
    }
    if constexpr (!useGreedyDirichlet)
      boundaryNodes(rblid) = true;
  }
};

/*!
  @class RowSumFunctor
  @brief Functor for marking nodes as Dirichlet based on rowsum.

  A row is marked as Dirichlet boundary if the sum of off-diagonal values is smaller in absolute value than the diagonal multiplied by the threshold rowSumTol.
  It is assumed that boundaryNodes was initialized to false.
*/
template <class local_matrix_type>
class RowSumFunctor {
 private:
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using memory_space       = typename local_matrix_type::memory_space;

#if KOKKOS_VERSION >= 40799
  using ATS = KokkosKernels::ArithTraits<scalar_type>;
#else
  using ATS    = Kokkos::ArithTraits<scalar_type>;
#endif
  using magnitudeType = typename ATS::magnitudeType;
#if KOKKOS_VERSION >= 40799
  using magATS = KokkosKernels::ArithTraits<magnitudeType>;
#else
  using magATS = Kokkos::ArithTraits<magnitudeType>;
#endif
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
    auto row            = A.rowConst(rlid);
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

/*!
  @class BoundaryFunctor
  @brief Functor that serially applies sub-functors to rows.
*/
template <class local_matrix_type, class Functor, class... RemainingFunctors>
class BoundaryFunctor {
 private:
  using local_ordinal_type = typename local_matrix_type::ordinal_type;

  Functor functor;
  BoundaryFunctor<local_matrix_type, RemainingFunctors...> remainingFunctors;

 public:
  BoundaryFunctor(local_matrix_type& A_, Functor& functor_, RemainingFunctors&... remainingFunctors_)
    : functor(functor_)
    , remainingFunctors(A_, remainingFunctors_...) {}

  KOKKOS_FUNCTION void operator()(const local_ordinal_type rlid) const {
    functor(rlid);
    remainingFunctors(rlid);
  }
};

template <class local_matrix_type, class Functor>
class BoundaryFunctor<local_matrix_type, Functor> {
 private:
  using local_ordinal_type = typename local_matrix_type::ordinal_type;

  local_matrix_type A;
  Functor functor;

 public:
  BoundaryFunctor(local_matrix_type& A_, Functor& functor_)
    : A(A_)
    , functor(functor_) {}

  KOKKOS_FUNCTION void operator()(const local_ordinal_type rlid) const {
    functor(rlid);
  }
};

}  // namespace MueLu::BoundaryDetection

#endif
