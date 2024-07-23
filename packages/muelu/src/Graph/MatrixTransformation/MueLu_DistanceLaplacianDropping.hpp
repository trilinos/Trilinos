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

template <class CoordsType, class local_ordinal_type>
class DistanceFunctor {
 private:
  using scalar_type   = typename CoordsType::non_const_value_type;
  using ATS           = Kokkos::ArithTraits<scalar_type>;
  using magnitudeType = typename ATS::magnitudeType;

  CoordsType coords;
  CoordsType ghostedCoords;

 public:
  DistanceFunctor(CoordsType coords_, CoordsType ghostedCoords_)
    : coords(coords_)
    , ghostedCoords(ghostedCoords_) {}

  KOKKOS_INLINE_FUNCTION
  magnitudeType distance2(const local_ordinal_type row, const local_ordinal_type col) const {
    scalar_type d = ATS::zero();
    scalar_type s;
    for (size_t j = 0; j < coords.extent(1); j++) {
      s = coords(row, j) - ghostedCoords(col, j);
      d += s * s;
    }
    return ATS::magnitude(d);
  }
};

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class DistanceFunctorType>
Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
getDiagonal(const Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& A,
            DistanceFunctorType& distFunctor) {
  using scalar_type         = Scalar;
  using local_ordinal_type  = LocalOrdinal;
  using global_ordinal_type = GlobalOrdinal;
  using node_type           = Node;
  using ATS                 = Kokkos::ArithTraits<scalar_type>;
  using execution_space     = typename Node::execution_space;
  using range_type          = Kokkos::RangePolicy<LocalOrdinal, execution_space>;

  auto diag = Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(A->getRowMap(), 1);
  {
    auto lclA    = A->getLocalMatrixDevice();
    auto lclDiag = diag->getDeviceLocalView(Xpetra::Access::OverwriteAll);

    Kokkos::parallel_for(
        "MueLu:CoalesceDropF:Build:scalar_filter:laplacian_diag",
        range_type(0, lclA.numRows()),
        KOKKOS_LAMBDA(const local_ordinal_type row) {
          auto rowView = lclA.rowConst(row);
          auto length  = rowView.length;

          scalar_type d = ATS::zero();
          for (local_ordinal_type colID = 0; colID < length; colID++) {
            auto col = rowView.colidx(colID);
            if (row != col)
              d += ATS::one() / distFunctor.distance2(row, col);
          }
          lclDiag(row, 0) = d;
        });
  }
  auto importer = A->getCrsGraph()->getImporter();
  if (!importer.is_null()) {
    auto ghostedDiag = Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(A->getColMap(), 1);
    ghostedDiag->doImport(*diag, *importer, Xpetra::INSERT);
    return ghostedDiag;
  } else {
    return diag;
  }
}

template <class local_matrix_type, class DistanceFunctorType, class diag_view_type>
class DropFunctor {
 private:
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using memory_space       = typename local_matrix_type::memory_space;
  using results_view       = Kokkos::View<DecisionType*, memory_space>;

  using ATS                 = Kokkos::ArithTraits<scalar_type>;
  using magnitudeType       = typename ATS::magnitudeType;
  using boundary_nodes_view = Kokkos::View<const bool*, memory_space>;

  local_matrix_type A;
  magnitudeType eps;
  diag_view_type diag;  // corresponds to overlapped diagonal
  DistanceFunctorType dist2;
  results_view results;
  const scalar_type one = ATS::one();

 public:
  DropFunctor(local_matrix_type& A_, magnitudeType threshold, diag_view_type diag_, DistanceFunctorType& dist2_, results_view& results_)
    : A(A_)
    , eps(threshold)
    , diag(diag_)
    , dist2(dist2_)
    , results(results_) {}

  KOKKOS_FORCEINLINE_FUNCTION
  bool operator()(local_ordinal_type rlid) {
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

      results(offset + k) = std::max((aij2 <= eps * eps * aiiajj) ? DROP : KEEP,
                                     results(offset + k));
    }
    return false;
  }
};

}  // namespace MueLu::DistanceLaplacian

#endif
