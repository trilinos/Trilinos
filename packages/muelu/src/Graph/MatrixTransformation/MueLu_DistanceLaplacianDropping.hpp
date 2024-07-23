// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_DISTANCELAPLACIANDROPPING_HPP
#define MUELU_DISTANCELAPLACIANDROPPING_HPP

#include "MueLu_DroppingCommon.hpp"
#include "Kokkos_Core.hpp"
#include "Kokkos_ArithTraits.hpp"
#include "Teuchos_RCP.hpp"
#include "Xpetra_Matrix.hpp"
#include "Xpetra_MultiVector.hpp"
#include "Xpetra_MultiVectorFactory.hpp"

namespace MueLu::DistanceLaplacian {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class DistanceFunctor {
 private:
  using matrix_type        = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using local_matrix_type  = typename matrix_type::local_matrix_type;
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = LocalOrdinal;
  using ATS                = Kokkos::ArithTraits<scalar_type>;
  using magnitudeType      = typename ATS::magnitudeType;
  using magATS             = Kokkos::ArithTraits<magnitudeType>;
  using coords_type        = Xpetra::MultiVector<magnitudeType, LocalOrdinal, GlobalOrdinal, Node>;
  using local_coords_type  = typename coords_type::dual_view_type_const::t_dev;

  Teuchos::RCP<coords_type> coordsMV;
  Teuchos::RCP<coords_type> ghostedCoordsMV;

  local_coords_type coords;
  local_coords_type ghostedCoords;

 public:
  DistanceFunctor(matrix_type& A, Teuchos::RCP<coords_type>& coords_) {
    coordsMV      = coords_;
    auto importer = A.getCrsGraph()->getImporter();
    if (!importer.is_null()) {
      ghostedCoordsMV = Xpetra::MultiVectorFactory<magnitudeType, LocalOrdinal, GlobalOrdinal, Node>::Build(importer->getTargetMap(), coordsMV->getNumVectors());
      ghostedCoordsMV->doImport(*coordsMV, *importer, Xpetra::INSERT);
      coords        = coordsMV->getDeviceLocalView(Xpetra::Access::ReadOnly);
      ghostedCoords = ghostedCoordsMV->getDeviceLocalView(Xpetra::Access::ReadOnly);
    } else {
      coords        = coordsMV->getDeviceLocalView(Xpetra::Access::ReadOnly);
      ghostedCoords = coords;
    }
  }

  KOKKOS_FORCEINLINE_FUNCTION
  magnitudeType distance2(const local_ordinal_type row, const local_ordinal_type col) const {
    magnitudeType d = magATS::zero();
    magnitudeType s;
    for (size_t j = 0; j < coords.extent(1); ++j) {
      s = coords(row, j) - ghostedCoords(col, j);
      d += s * s;
    }
    return d;
  }
};

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class DistanceFunctorType>
Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
getDiagonal(Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
            DistanceFunctorType& distFunctor) {
  using scalar_type         = Scalar;
  using local_ordinal_type  = LocalOrdinal;
  using global_ordinal_type = GlobalOrdinal;
  using node_type           = Node;
  using ATS                 = Kokkos::ArithTraits<scalar_type>;
  using impl_scalar_type    = typename ATS::val_type;
  using implATS             = Kokkos::ArithTraits<impl_scalar_type>;
  using magnitudeType       = typename implATS::magnitudeType;
  using execution_space     = typename Node::execution_space;
  using range_type          = Kokkos::RangePolicy<LocalOrdinal, execution_space>;

  auto diag = Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(A.getRowMap(), 1);
  {
    auto lclA    = A.getLocalMatrixDevice();
    auto lclDiag = diag->getDeviceLocalView(Xpetra::Access::OverwriteAll);

    Kokkos::parallel_for(
        "MueLu:CoalesceDropF:Build:scalar_filter:laplacian_diag",
        range_type(0, lclA.numRows()),
        KOKKOS_LAMBDA(const local_ordinal_type& row) {
          auto rowView = lclA.rowConst(row);
          auto length  = rowView.length;

          magnitudeType d;
          impl_scalar_type d2 = implATS::zero();
          for (local_ordinal_type colID = 0; colID < length; colID++) {
            auto col = rowView.colidx(colID);
            if (row != col) {
              d = distFunctor.distance2(row, col);
              d2 += implATS::one() / d;
            }
          }
          lclDiag(row, 0) = d2;
        });
  }
  auto importer = A.getCrsGraph()->getImporter();
  if (!importer.is_null()) {
    auto ghostedDiag = Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(A.getColMap(), 1);
    ghostedDiag->doImport(*diag, *importer, Xpetra::INSERT);
    return ghostedDiag;
  } else {
    return diag;
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class DistanceFunctorType>
class DropFunctor {
 private:
  using matrix_type        = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using local_matrix_type  = typename matrix_type::local_matrix_type;
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using memory_space       = typename local_matrix_type::memory_space;
  using diag_vec_type      = Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using diag_view_type     = typename Kokkos::DualView<const scalar_type*, Kokkos::LayoutStride, typename Node::device_type, Kokkos::MemoryUnmanaged>::t_dev;

  using results_view = Kokkos::View<DecisionType*, memory_space>;

  using ATS                 = Kokkos::ArithTraits<scalar_type>;
  using magnitudeType       = typename ATS::magnitudeType;
  using boundary_nodes_view = Kokkos::View<const bool*, memory_space>;

  local_matrix_type A;
  magnitudeType eps;
  Teuchos::RCP<diag_vec_type> diagVec;
  diag_view_type diag;  // corresponds to overlapped diagonal
  DistanceFunctorType dist2;
  results_view results;
  const scalar_type one = ATS::one();

 public:
  DropFunctor(matrix_type& A_, magnitudeType threshold, DistanceFunctorType& dist2_, results_view& results_)
    : A(A_.getLocalMatrixDevice())
    , eps(threshold)
    , dist2(dist2_)
    , results(results_) {
    diagVec        = getDiagonal(A_, dist2);
    auto lclDiag2d = diagVec->getDeviceLocalView(Xpetra::Access::ReadOnly);
    diag           = Kokkos::subview(lclDiag2d, Kokkos::ALL(), 0);
  }

  KOKKOS_FORCEINLINE_FUNCTION
  void operator()(local_ordinal_type rlid) const {
    auto row            = A.rowConst(rlid);
    const size_t offset = A.graph.row_map(rlid);
    for (local_ordinal_type k = 0; k < row.length; ++k) {
      auto clid = row.colidx(k);

      scalar_type val;
      if (rlid != clid) {
        val = one / dist2.distance2(rlid, clid);
      } else {
        val = diag(rlid);
      }
      auto aiiajj = ATS::magnitude(diag(rlid)) * ATS::magnitude(diag(clid));  // |a_ii|*|a_jj|
      auto aij2   = ATS::magnitude(val) * ATS::magnitude(val);                // |a_ij|^2

      results(offset + k) = Kokkos::max((aij2 <= eps * eps * aiiajj) ? DROP : KEEP,
                                        results(offset + k));
    }
  }
};

}  // namespace MueLu::DistanceLaplacian

#endif
