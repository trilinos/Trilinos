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

#include "Kokkos_Macros.hpp"
#include "KokkosBatched_LU_Decl.hpp"
#include "KokkosBatched_LU_Serial_Impl.hpp"
#include "KokkosBatched_Trsv_Decl.hpp"
#include "KokkosBatched_Trsv_Serial_Impl.hpp"
#include "MueLu_DroppingCommon.hpp"
#include "Kokkos_Core.hpp"
#include "Kokkos_ArithTraits.hpp"
#include "Teuchos_RCP.hpp"
#include "Xpetra_Matrix.hpp"
#include "Xpetra_MultiVector.hpp"
#include "Xpetra_MultiVectorFactory.hpp"

namespace MueLu::DistanceLaplacian {

/*!
@class DistanceFunctor
@brief Computes the unscaled distance Laplacian.
*/
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class UnweightedDistanceFunctor {
 private:
  using matrix_type        = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using local_matrix_type  = typename matrix_type::local_matrix_type;
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = LocalOrdinal;
  using ATS                = Kokkos::ArithTraits<scalar_type>;
  using impl_scalar_type   = typename ATS::val_type;
  using implATS            = Kokkos::ArithTraits<impl_scalar_type>;
  using magnitudeType      = typename implATS::magnitudeType;
  using magATS             = Kokkos::ArithTraits<magnitudeType>;
  using coords_type        = Xpetra::MultiVector<magnitudeType, LocalOrdinal, GlobalOrdinal, Node>;
  using local_coords_type  = typename coords_type::dual_view_type_const::t_dev;

  Teuchos::RCP<coords_type> coordsMV;
  Teuchos::RCP<coords_type> ghostedCoordsMV;

  local_coords_type coords;
  local_coords_type ghostedCoords;

 public:
  UnweightedDistanceFunctor(matrix_type& A, Teuchos::RCP<coords_type>& coords_) {
    coordsMV      = coords_;
    auto importer = A.getCrsGraph()->getImporter();
    if (!importer.is_null()) {
      ghostedCoordsMV = Xpetra::MultiVectorFactory<magnitudeType, LocalOrdinal, GlobalOrdinal, Node>::Build(importer->getTargetMap(), coordsMV->getNumVectors());
      ghostedCoordsMV->doImport(*coordsMV, *importer, Xpetra::INSERT);
      coords        = coordsMV->getLocalViewDevice(Xpetra::Access::ReadOnly);
      ghostedCoords = ghostedCoordsMV->getLocalViewDevice(Xpetra::Access::ReadOnly);
    } else {
      coords        = coordsMV->getLocalViewDevice(Xpetra::Access::ReadOnly);
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

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class ScalarMaterialDistanceFunctor {
 private:
  using matrix_type         = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using local_matrix_type   = typename matrix_type::local_matrix_type;
  using scalar_type         = typename local_matrix_type::value_type;
  using local_ordinal_type  = LocalOrdinal;
  using ATS                 = Kokkos::ArithTraits<scalar_type>;
  using impl_scalar_type    = typename ATS::val_type;
  using implATS             = Kokkos::ArithTraits<impl_scalar_type>;
  using magnitudeType       = typename implATS::magnitudeType;
  using magATS              = Kokkos::ArithTraits<magnitudeType>;
  using coords_type         = Xpetra::MultiVector<magnitudeType, LocalOrdinal, GlobalOrdinal, Node>;
  using local_coords_type   = typename coords_type::dual_view_type_const::t_dev;
  using material_type       = Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using local_material_type = typename material_type::dual_view_type_const::t_dev;

  Teuchos::RCP<coords_type> coordsMV;
  Teuchos::RCP<coords_type> ghostedCoordsMV;

  local_coords_type coords;
  local_coords_type ghostedCoords;

  Teuchos::RCP<material_type> materialMV;
  Teuchos::RCP<material_type> ghostedMaterialMV;

  local_material_type material;
  local_material_type ghostedMaterial;

 public:
  ScalarMaterialDistanceFunctor(matrix_type& A, Teuchos::RCP<coords_type>& coords_, Teuchos::RCP<material_type>& material_) {
    coordsMV      = coords_;
    materialMV    = material_;
    auto importer = A.getCrsGraph()->getImporter();
    if (!importer.is_null()) {
      ghostedCoordsMV = Xpetra::MultiVectorFactory<magnitudeType, LocalOrdinal, GlobalOrdinal, Node>::Build(importer->getTargetMap(), coordsMV->getNumVectors());
      ghostedCoordsMV->doImport(*coordsMV, *importer, Xpetra::INSERT);
      coords        = coordsMV->getLocalViewDevice(Xpetra::Access::ReadOnly);
      ghostedCoords = ghostedCoordsMV->getLocalViewDevice(Xpetra::Access::ReadOnly);

      ghostedMaterialMV = Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(importer->getTargetMap(), materialMV->getNumVectors());
      ghostedMaterialMV->doImport(*materialMV, *importer, Xpetra::INSERT);
      material        = materialMV->getLocalViewDevice(Xpetra::Access::ReadOnly);
      ghostedMaterial = ghostedMaterialMV->getLocalViewDevice(Xpetra::Access::ReadOnly);
    } else {
      coords        = coordsMV->getLocalViewDevice(Xpetra::Access::ReadOnly);
      ghostedCoords = coords;

      material        = materialMV->getLocalViewDevice(Xpetra::Access::ReadOnly);
      ghostedMaterial = material;
    }
  }

  KOKKOS_INLINE_FUNCTION
  magnitudeType distance2(const local_ordinal_type row, const local_ordinal_type col) const {
    // || x_row - x_col ||_S^2
    // where
    // S = 1/material(row) * Identity
    magnitudeType d     = magATS::zero();
    magnitudeType d_row = magATS::zero();
    magnitudeType d_col = magATS::zero();
    magnitudeType s;

    for (size_t j = 0; j < coords.extent(1); ++j) {
      s = coords(row, j) - ghostedCoords(col, j);
      d += s * s;
    }

    d_row = d / implATS::magnitude(ghostedMaterial(row, 0));
    d_col = d / implATS::magnitude(ghostedMaterial(col, 0));

    return Kokkos::max(d_row, d_col);
  }
};

template <class local_ordinal_type, class material_vector_type, class material_matrix_type>
class TensorInversion {
 private:
  material_vector_type material_vector;
  material_matrix_type material_matrix;

 public:
  TensorInversion(material_vector_type& material_vector_, material_matrix_type& material_matrix_)
    : material_vector(material_vector_)
    , material_matrix(material_matrix_) {}

  KOKKOS_FORCEINLINE_FUNCTION
  void operator()(local_ordinal_type i) const {
    for (size_t j = 0; j < material_matrix.extent(1); ++j) {
      for (size_t k = 0; k < material_matrix.extent(2); ++k) {
        material_matrix(i, j, k) = material_vector(i, j * material_matrix.extent(1) + k);
      }
    }
    auto matrix_material = Kokkos::subview(material_matrix, i, Kokkos::ALL(), Kokkos::ALL());
    KokkosBatched::SerialLU<KokkosBatched::Algo::LU::Unblocked>::invoke(matrix_material);
  }
};

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class TensorMaterialDistanceFunctor {
 private:
  using matrix_type        = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using local_matrix_type  = typename matrix_type::local_matrix_type;
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = LocalOrdinal;
  using ATS                = Kokkos::ArithTraits<scalar_type>;
  using impl_scalar_type   = typename ATS::val_type;
  using implATS            = Kokkos::ArithTraits<impl_scalar_type>;
  using magnitudeType      = typename implATS::magnitudeType;
  using magATS             = Kokkos::ArithTraits<magnitudeType>;
  using coords_type        = Xpetra::MultiVector<magnitudeType, LocalOrdinal, GlobalOrdinal, Node>;
  using local_coords_type  = typename coords_type::dual_view_type_const::t_dev;
  using material_type      = Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using memory_space       = typename local_matrix_type::memory_space;

  using local_material_type = Kokkos::View<impl_scalar_type***, memory_space>;
  using local_dist_type     = Kokkos::View<impl_scalar_type**, memory_space>;

  Teuchos::RCP<coords_type> coordsMV;
  Teuchos::RCP<coords_type> ghostedCoordsMV;

  local_coords_type coords;
  local_coords_type ghostedCoords;

  local_material_type material;

  local_dist_type lcl_dist;

  const scalar_type one = ATS::one();

 public:
  TensorMaterialDistanceFunctor(matrix_type& A, Teuchos::RCP<coords_type>& coords_, Teuchos::RCP<material_type>& material_) {
    coordsMV = coords_;

    auto importer = A.getCrsGraph()->getImporter();
    if (!importer.is_null()) {
      ghostedCoordsMV = Xpetra::MultiVectorFactory<magnitudeType, LocalOrdinal, GlobalOrdinal, Node>::Build(importer->getTargetMap(), coordsMV->getNumVectors());
      ghostedCoordsMV->doImport(*coordsMV, *importer, Xpetra::INSERT);
      coords        = coordsMV->getLocalViewDevice(Xpetra::Access::ReadOnly);
      ghostedCoords = ghostedCoordsMV->getLocalViewDevice(Xpetra::Access::ReadOnly);
    } else {
      coords        = coordsMV->getLocalViewDevice(Xpetra::Access::ReadOnly);
      ghostedCoords = coords;
    }

    {
      Teuchos::RCP<material_type> ghostedMaterial;
      if (!importer.is_null()) {
        ghostedMaterial = Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(importer->getTargetMap(), material_->getNumVectors());
        ghostedMaterial->doImport(*material_, *importer, Xpetra::INSERT);
      } else {
        ghostedMaterial = material_;
      }

      using execution_space = typename Node::execution_space;
      using range_type      = Kokkos::RangePolicy<LocalOrdinal, execution_space>;

      local_ordinal_type dim = std::sqrt(material_->getNumVectors());
      auto lclMaterial       = ghostedMaterial->getLocalViewDevice(Xpetra::Access::ReadOnly);
      material               = local_material_type("material", lclMaterial.extent(0), dim, dim);
      lcl_dist               = local_dist_type("material", lclMaterial.extent(0), dim);
      TensorInversion<local_ordinal_type, typename material_type::dual_view_type::t_dev_const_um, local_material_type> functor(lclMaterial, material);
      Kokkos::parallel_for("MueLu:TensorMaterialDistanceFunctor::inversion", range_type(0, lclMaterial.extent(0)), functor);
    }
  }

  KOKKOS_INLINE_FUNCTION
  magnitudeType distance2(const local_ordinal_type row, const local_ordinal_type col) const {
    // || x_row - x_col ||_S^2
    // where
    // S = inv(material(col))

    // row material
    impl_scalar_type d_row = implATS::zero();
    {
      auto matrix_row_material = Kokkos::subview(material, row, Kokkos::ALL(), Kokkos::ALL());
      auto dist                = Kokkos::subview(lcl_dist, row, Kokkos::ALL());

      for (size_t j = 0; j < coords.extent(1); ++j) {
        dist(j) = coords(row, j) - ghostedCoords(col, j);
      }

      KokkosBatched::SerialTrsv<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::NoTranspose, KokkosBatched::Diag::Unit, KokkosBatched::Algo::Trsv::Unblocked>::invoke(one, matrix_row_material, dist);
      KokkosBatched::SerialTrsv<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::NoTranspose, KokkosBatched::Diag::NonUnit, KokkosBatched::Algo::Trsv::Unblocked>::invoke(one, matrix_row_material, dist);

      for (size_t j = 0; j < coords.extent(1); ++j) {
        d_row += dist(j) * (coords(row, j) - ghostedCoords(col, j));
      }
    }

    // column material
    impl_scalar_type d_col = implATS::zero();
    {
      auto matrix_col_material = Kokkos::subview(material, col, Kokkos::ALL(), Kokkos::ALL());
      auto dist                = Kokkos::subview(lcl_dist, row, Kokkos::ALL());

      for (size_t j = 0; j < coords.extent(1); ++j) {
        dist(j) = coords(row, j) - ghostedCoords(col, j);
      }

      KokkosBatched::SerialTrsv<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::NoTranspose, KokkosBatched::Diag::Unit, KokkosBatched::Algo::Trsv::Unblocked>::invoke(one, matrix_col_material, dist);
      KokkosBatched::SerialTrsv<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::NoTranspose, KokkosBatched::Diag::NonUnit, KokkosBatched::Algo::Trsv::Unblocked>::invoke(one, matrix_col_material, dist);

      for (size_t j = 0; j < coords.extent(1); ++j) {
        d_col += dist(j) * (coords(row, j) - ghostedCoords(col, j));
      }
    }

    return Kokkos::max(implATS::magnitude(d_row), implATS::magnitude(d_col));
  }
};

/*!
Method to compute ghosted distance Laplacian diagonal.
*/
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
    auto lclDiag = diag->getLocalViewDevice(Xpetra::Access::OverwriteAll);

    Kokkos::parallel_for(
        "MueLu:CoalesceDropF:Build:scalar_filter:laplacian_diag",
        range_type(0, lclA.numRows()),
        KOKKOS_LAMBDA(const local_ordinal_type& row) {
          auto rowView = lclA.rowConst(row);
          auto length  = rowView.length;

          magnitudeType d;
          impl_scalar_type d2  = implATS::zero();
          bool haveAddedToDiag = false;
          for (local_ordinal_type colID = 0; colID < length; colID++) {
            auto col = rowView.colidx(colID);
            if (row != col) {
              d = distFunctor.distance2(row, col);
              d2 += implATS::one() / d;
              haveAddedToDiag = true;
            }
          }

          // Deal with the situation where boundary conditions have only been enforced on rows, but not on columns.
          // We enforce dropping of these entries by assigning a very large number to the diagonal entries corresponding to BCs.
          lclDiag(row, 0) = !haveAddedToDiag ? implATS::squareroot(implATS::rmax()) : d2;
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

/*!
@class DropFunctor
@brief Drops entries the unscaled distance Laplacian.

Evaluates the dropping criterion
\f[
\frac{|d_{ij}|^2}{|d_{ii}| |d_{jj}|} \le \theta^2
\f]
where \f$d_{ij}\f$ is a distance metric.
*/
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
    auto lclDiag2d = diagVec->getLocalViewDevice(Xpetra::Access::ReadOnly);
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

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class DistanceFunctorType>
class VectorDropFunctor {
 private:
  using matrix_type             = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using local_matrix_type       = typename matrix_type::local_matrix_type;
  using scalar_type             = typename local_matrix_type::value_type;
  using local_ordinal_type      = typename local_matrix_type::ordinal_type;
  using memory_space            = typename local_matrix_type::memory_space;
  using block_indices_view_type = Kokkos::View<local_ordinal_type*, memory_space>;
  using diag_vec_type           = Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using diag_view_type          = typename Kokkos::DualView<const scalar_type*, Kokkos::LayoutStride, typename Node::device_type, Kokkos::MemoryUnmanaged>::t_dev;

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
  block_indices_view_type point_to_block;
  block_indices_view_type ghosted_point_to_block;
  const scalar_type one = ATS::one();

 public:
  VectorDropFunctor(matrix_type& A_, matrix_type& mergedA_, magnitudeType threshold, DistanceFunctorType& dist2_, results_view& results_, block_indices_view_type point_to_block_, block_indices_view_type ghosted_point_to_block_)
    : A(A_.getLocalMatrixDevice())
    , eps(threshold)
    , dist2(dist2_)
    , results(results_)
    , point_to_block(point_to_block_)
    , ghosted_point_to_block(ghosted_point_to_block_) {
    diagVec        = getDiagonal(mergedA_, dist2);
    auto lclDiag2d = diagVec->getLocalViewDevice(Xpetra::Access::ReadOnly);
    diag           = Kokkos::subview(lclDiag2d, Kokkos::ALL(), 0);
  }

  KOKKOS_FORCEINLINE_FUNCTION
  void operator()(local_ordinal_type rlid) const {
    auto brlid          = point_to_block(rlid);
    auto row            = A.rowConst(rlid);
    const size_t offset = A.graph.row_map(rlid);
    for (local_ordinal_type k = 0; k < row.length; ++k) {
      auto clid  = row.colidx(k);
      auto bclid = ghosted_point_to_block(clid);

      scalar_type val;
      if (brlid != bclid) {
        val = one / dist2.distance2(brlid, bclid);
      } else {
        val = diag(brlid);
      }
      auto aiiajj = ATS::magnitude(diag(brlid)) * ATS::magnitude(diag(bclid));  // |a_ii|*|a_jj|
      auto aij2   = ATS::magnitude(val) * ATS::magnitude(val);                  // |a_ij|^2

      results(offset + k) = Kokkos::max((aij2 <= eps * eps * aiiajj) ? DROP : KEEP,
                                        results(offset + k));
    }
  }
};

}  // namespace MueLu::DistanceLaplacian

#endif
