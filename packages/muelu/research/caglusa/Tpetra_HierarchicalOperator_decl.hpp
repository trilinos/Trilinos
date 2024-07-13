// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_HIERARCHICALOPERATOR_DECL_HPP
#define TPETRA_HIERARCHICALOPERATOR_DECL_HPP

#include <Teuchos_RCP.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_RowMatrix.hpp>
#include <Tpetra_OperatorWithDiagonal.hpp>
#include <Tpetra_BlockedMap_decl.hpp>
#include <Tpetra_BlockedMatrix_decl.hpp>
#include "Teuchos_VerbosityLevel.hpp"

namespace Tpetra {

/*

A container class for hierarchical matrices of different types.
In particular, both H- and H2-matrices are supported.

The unknowns of the kernel approximations are collected in the clusterMap.
For H-matrices, this is just a concatenation.
For H2-matrices, the map also contains the intermediate clusters that might be needed in upward/downward pass.


H = nearField
    + basisMatrix *
      ((I+transferMatrices[K-1]^T) * ... * (I+transferMatrices[0]^T)) *
      kernelApproximations *
      ((I+transferMatrices[0]) * ... * (I+transferMatrices[K-1])) *
      basisMatrix^T

nearField and basisMatrix are standard (point) CRS matrices.
kernelApproximations and transferMatrices[.] are blocked CRS matrices

I is the identity matrix and is not explicitely saved.

Maps:
map (standard): domain and range of H;
               domain, range, row of nearField
clusterMap (blocked map): domain, range, row, column of transferMatrices;
                         domain, range, row of kernelApproximations
ghosted_clusterMap (blocked map): column of kernelApproximations


For H-matrices:
K = 0, i.e. there are no transfer matrices

For H2-matrices:
upward and downward pass in the cluster hierarchy are encoded in transfer matrices
*/

template <class Scalar        = Tpetra::Operator<>::scalar_type,
          class LocalOrdinal  = typename Tpetra::Operator<Scalar>::local_ordinal_type,
          class GlobalOrdinal = typename Tpetra::Operator<Scalar, LocalOrdinal>::global_ordinal_type,
          class Node          = typename Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
class HierarchicalOperator : public Tpetra::OperatorWithDiagonal<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
 public:
  using matrix_type = Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using mv_type     = Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using map_type    = Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;

  //! The RowMatrix representing the base class of CrsMatrix
  using row_matrix_type = RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  using impl_scalar_type = typename row_matrix_type::impl_scalar_type;
  using mag_type         = typename Kokkos::ArithTraits<impl_scalar_type>::mag_type;

  using local_inds_device_view_type =
      typename row_matrix_type::local_inds_device_view_type;
  using local_inds_host_view_type =
      typename row_matrix_type::local_inds_host_view_type;
  using nonconst_local_inds_host_view_type =
      typename row_matrix_type::nonconst_local_inds_host_view_type;

  using global_inds_device_view_type =
      typename row_matrix_type::global_inds_device_view_type;
  using global_inds_host_view_type =
      typename row_matrix_type::global_inds_host_view_type;
  using nonconst_global_inds_host_view_type =
      typename row_matrix_type::nonconst_global_inds_host_view_type;

  using values_device_view_type =
      typename row_matrix_type::values_device_view_type;
  using values_host_view_type =
      typename row_matrix_type::values_host_view_type;
  using nonconst_values_host_view_type =
      typename row_matrix_type::nonconst_values_host_view_type;

  using blocked_matrix_type = BlockedMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using blocked_map_type    = BlockedMap<LocalOrdinal, GlobalOrdinal, Node>;

  //! @name Constructor/Destructor
  //@{

  //! Constructor
  HierarchicalOperator(const Teuchos::RCP<matrix_type>& nearField,
                       const Teuchos::RCP<blocked_matrix_type>& kernelApproximations,
                       const Teuchos::RCP<matrix_type>& basisMatrix,
                       std::vector<Teuchos::RCP<blocked_matrix_type> >& transferMatrices,
                       const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

  //! Returns the Tpetra::Map object associated with the domain of this operator.
  Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > getDomainMap() const {
    return nearField_->getDomainMap();
  }

  //! Returns the Tpetra::Map object associated with the range of this operator.
  Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > getRangeMap() const {
    return nearField_->getRangeMap();
  }

  //! Returns in Y the result of a Tpetra::Operator applied to a Tpetra::MultiVector X.
  /*!
    \param[in]  X - Tpetra::MultiVector of dimension NumVectors to multiply with matrix.
    \param[out] Y -Tpetra::MultiVector of dimension NumVectors containing result.
  */
  void apply(const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
             Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Y,
             Teuchos::ETransp mode = Teuchos::NO_TRANS,
             Scalar alpha          = Teuchos::ScalarTraits<Scalar>::one(),
             Scalar beta           = Teuchos::ScalarTraits<Scalar>::zero()) const;

  Teuchos::RCP<HierarchicalOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > restrict(const Teuchos::RCP<matrix_type>& P);

  Teuchos::RCP<matrix_type> toMatrix();

  double getCompression() {
    size_t nnz = (nearField_->getGlobalNumEntries() +
                  kernelApproximations_->pointA_->getGlobalNumEntries() +
                  basisMatrix_->getGlobalNumEntries());
    for (size_t i = 0; i < transferMatrices_.size(); i++)
      nnz += transferMatrices_[i]->pointA_->getGlobalNumEntries();
    return Teuchos::as<double>(nnz) / (getDomainMap()->getGlobalNumElements() * getDomainMap()->getGlobalNumElements());
  }

  size_t numClusterPairsOnLevelGlobal(size_t level) const {
    using vec_type = typename Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    using Teuchos::RCP;
    const Scalar ONE  = Teuchos::ScalarTraits<Scalar>::one();
    const Scalar ZERO = Teuchos::ScalarTraits<Scalar>::zero();
    if (level == 0) {
      // RCP<vec_type> tempV         = Teuchos::rcp(new vec_type(kernelApproximations_->blockMap_->blockMap_, false));
      // RCP<vec_type> tempV2        = Teuchos::rcp(new vec_type(kernelApproximations_->blockMap_->blockMap_, false));

      // tempV->putScalar(ONE);
      // transferMatrices_[0]->blockA_->apply(*tempV, *tempV2, Teuchos::TRANS);
      // tempV->putScalar(ZERO);
      // kernelApproximations_->blockA_->apply(*tempV2, *tempV);
      // return Teuchos::as<size_t>(tempV->dot(*tempV2));
      return 0;
    } else if (level <= transferMatrices_.size()) {
      RCP<vec_type> tempV  = Teuchos::rcp(new vec_type(kernelApproximations_->blockMap_->blockMap_, false));
      RCP<vec_type> tempV2 = Teuchos::rcp(new vec_type(kernelApproximations_->blockMap_->blockMap_, false));

      tempV->putScalar(ONE);
      transferMatrices_[level - 1]->blockA_->apply(*tempV, *tempV2, Teuchos::TRANS);
      tempV->putScalar(ZERO);
      kernelApproximations_->blockA_->apply(*tempV2, *tempV);
      return Teuchos::as<size_t>(tempV->dot(*tempV2));
    } else {
      TEUCHOS_ASSERT(false);
    }
  }

  size_t numClustersOnLevelGlobal(size_t level) const {
    if (level > 0)
      return Teuchos::as<size_t>(transferMatrices_[level - 1]->blockA_->getGlobalNumEntries());
    else {
      const double treeCoarseningFactor = params_->get<double>("treeCoarseningFactor");
      return Teuchos::as<size_t>(transferMatrices_[0]->blockA_->getGlobalNumEntries() / treeCoarseningFactor);
    }
  }

  size_t numClustersOnLevelLocal(size_t level) const {
    if (level > 0)
      return Teuchos::as<size_t>(transferMatrices_[level - 1]->blockA_->getLocalNumEntries());
    else {
      const double treeCoarseningFactor = params_->get<double>("treeCoarseningFactor");
      return Teuchos::as<size_t>(transferMatrices_[0]->blockA_->getLocalNumEntries() / treeCoarseningFactor);
    }
  }

  size_t numLevels() const {
    return transferMatrices_.size() + 1;
  }

  void setDebugOutput(const bool debugOutput) {
    debugOutput_ = debugOutput;
  }

  Teuchos::RCP<matrix_type> nearFieldMatrix() {
    return nearField_;
  }

  // Fake RowMatrix interface
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > getRowMap() const {
    return nearField_->getRowMap();
  }

  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > getColMap() const {
    return nearField_->getColMap();
  }

  Teuchos::RCP<const Teuchos::Comm<int> > getComm() const {
    return nearField_->getDomainMap()->getComm();
  }

  Teuchos::RCP<const RowGraph<LocalOrdinal, GlobalOrdinal, Node> > getGraph() const {
    return nearField_->getCrsGraph();
  }

  global_size_t getGlobalNumRows() const {
    return nearField_->getGlobalNumRows();
  }

  global_size_t getGlobalNumCols() const {
    return nearField_->getGlobalNumCols();
  }

  size_t getLocalNumRows() const {
    return nearField_->getLocalNumRows();
  }

  size_t getLocalNumCols() const {
    return nearField_->getLocalNumCols();
  }

  GlobalOrdinal getIndexBase() const {
    return nearField_->getIndexBase();
  }

  global_size_t getGlobalNumEntries() const {
    return nearField_->getGlobalNumEntries();
  }

  size_t getLocalNumEntries() const {
    return nearField_->getLocalNumEntries();
  }

  void getLocalDiagCopy(Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& diag) const {
    nearField_->getLocalDiagCopy(diag);
  }

  void describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const {
    describe(out, verbLevel, true);
  }

  void describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel, const bool printHeader) const {
    using std::endl;
    using std::setw;
    const size_t numRows         = nearField_->getRowMap()->getGlobalNumElements();
    const size_t nnzNearField    = nearField_->getGlobalNumEntries();
    const double nnzNearPerRow   = Teuchos::as<double>(nnzNearField) / numRows;
    const size_t nnzKernelApprox = kernelApproximations_->pointA_->getGlobalNumEntries();
    const size_t numClusterPairs = kernelApproximations_->blockA_->getGlobalNumEntries();
    const size_t nnzBasis        = basisMatrix_->getGlobalNumEntries();
    size_t numTransfers          = transferMatrices_.size();
    size_t nnzTransfer           = 0;
    for (size_t i = 0; i < transferMatrices_.size(); i++)
      nnzTransfer += transferMatrices_[i]->pointA_->getGlobalNumEntries();
    const size_t nnzTotal       = nnzNearField + nnzKernelApprox + nnzBasis + nnzTransfer;
    const double nnzTotalPerRow = Teuchos::as<double>(nnzTotal) / numRows;
    std::ostringstream oss;
    oss << std::left;
    if (printHeader)
      oss << setw(9) << "rows" << setw(12)
          << "nnz(near)" << setw(14)
          << "nnz(near)/row" << setw(12)
          << "nnz(basis)" << setw(15)
          << "#cluster pairs" << setw(12)
          << "nnz(kernel)" << setw(14)
          << "#transfers" << setw(14)
          << "nnz(transfer)" << setw(12)
          << "nnz(total)" << setw(14)
          << "nnz(total)/row" << endl;
    oss << setw(9) << numRows << setw(12)
        << nnzNearField << setw(14)
        << nnzNearPerRow << setw(12)
        << nnzBasis << setw(15)
        << numClusterPairs << setw(12)
        << nnzKernelApprox << setw(14)
        << numTransfers << setw(14)
        << nnzTransfer << setw(12)
        << nnzTotal << setw(14)
        << nnzTotalPerRow << endl;

    if (Teuchos::includesVerbLevel(verbLevel, Teuchos::VERB_EXTREME)) {
      oss << setw(9) << "level" << setw(9) << "numClusters" << setw(9) << "numClusterPairs" << endl;
      for (size_t lvl = 0; lvl < numLevels(); ++lvl) {
        size_t numClustersLvl_gbl = numClustersOnLevelGlobal(lvl);
        size_t numClusterPairsLvl = numClusterPairsOnLevelGlobal(lvl);
        // size_t numClustersLvl_lcl = numClustersOnLevelLocal(lvl);

        oss << setw(9) << lvl << setw(9) << numClustersLvl_gbl << setw(9) << numClusterPairsLvl << endl;
      }
    }

    out << oss.str();
  }

  bool hasFarField() const {
    return kernelApproximations_->blockA_->getGlobalNumEntries() > 0;
  }

  bool hasTransferMatrices() const {
    return !transferMatrices_.empty();
  }

  bool denserThanDenseMatrix() const {
    const size_t numRows      = nearField_->getRowMap()->getGlobalNumElements();
    const size_t nnzNearField = nearField_->getGlobalNumEntries();
    // const double nnzNearPerRow = Teuchos::as<double>(nnzNearField)/numRows;
    const size_t nnzKernelApprox = kernelApproximations_->pointA_->getGlobalNumEntries();
    // const size_t numClusterPairs = kernelApproximations_->blockA_->getGlobalNumEntries();
    const size_t nnzBasis = basisMatrix_->getGlobalNumEntries();
    size_t nnzTransfer    = 0;
    for (size_t i = 0; i < transferMatrices_.size(); i++)
      nnzTransfer += transferMatrices_[i]->pointA_->getGlobalNumEntries();
    const size_t nnzTotal       = nnzNearField + nnzKernelApprox + nnzBasis + nnzTransfer;
    const double nnzTotalPerRow = Teuchos::as<double>(nnzTotal) / numRows;

    return (nnzTotalPerRow >= numRows);
  }

 private:
  void allocateMemory(size_t numVectors) const;

  void applyWithTransposes(const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
                           Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Y,
                           Teuchos::ETransp mode = Teuchos::NO_TRANS,
                           Scalar alpha          = Teuchos::ScalarTraits<Scalar>::one(),
                           Scalar beta           = Teuchos::ScalarTraits<Scalar>::zero()) const;

  void applyWithoutTransposes(const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
                              Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Y,
                              Teuchos::ETransp mode = Teuchos::NO_TRANS,
                              Scalar alpha          = Teuchos::ScalarTraits<Scalar>::one(),
                              Scalar beta           = Teuchos::ScalarTraits<Scalar>::zero()) const;

  bool canApplyWithoutTransposes_;
  std::string coarseningCriterion_;
  bool debugOutput_;

  Teuchos::RCP<matrix_type> nearField_;
  Teuchos::RCP<blocked_matrix_type> kernelApproximations_;
  Teuchos::RCP<matrix_type> basisMatrix_;
  Teuchos::RCP<matrix_type> basisMatrixT_;
  std::vector<Teuchos::RCP<blocked_matrix_type> > transferMatrices_;
  std::vector<Teuchos::RCP<blocked_matrix_type> > transferMatricesT_;
  Teuchos::RCP<const map_type> clusterCoeffMap_;
  mutable Teuchos::RCP<mv_type> coefficients_, coefficients2_;
  mutable Teuchos::RCP<mv_type> X_colmap_, coefficients_colmap_;

  Teuchos::RCP<Teuchos::ParameterList> params_;
};
}  // namespace Tpetra

#endif  // TPETRA_HIERARCHICALOPERATOR_DECL_HPP
