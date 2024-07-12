// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_HIERARCHICALOPERATOR_DEF_HPP
#define TPETRA_HIERARCHICALOPERATOR_DEF_HPP

#include <TpetraExt_MatrixMatrix.hpp>
#include <Tpetra_RowMatrixTransposer.hpp>
#include <Tpetra_HierarchicalOperator_decl.hpp>

namespace Tpetra {

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
removeSmallEntries(Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& A,
                   typename Teuchos::ScalarTraits<Scalar>::magnitudeType tol) {
  using crs_matrix   = Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using row_ptr_type = typename crs_matrix::local_graph_device_type::row_map_type::non_const_type;
  using col_idx_type = typename crs_matrix::local_graph_device_type::entries_type::non_const_type;
  using vals_type    = typename crs_matrix::local_matrix_device_type::values_type;

  using ATS      = Kokkos::ArithTraits<Scalar>;
  using impl_SC  = typename ATS::val_type;
  using impl_ATS = Kokkos::ArithTraits<impl_SC>;

  auto lclA = A->getLocalMatrixDevice();

  auto rowptr = row_ptr_type("rowptr", lclA.numRows() + 1);

  LocalOrdinal nnz;
  Kokkos::parallel_scan(
      "removeSmallEntries::rowptr",
      Kokkos::RangePolicy<LocalOrdinal>(0, lclA.numRows()),
      KOKKOS_LAMBDA(const LocalOrdinal rlid, LocalOrdinal& partial_nnz, bool is_final) {
        auto row = lclA.row(rlid);
        for (LocalOrdinal k = 0; k < row.length; ++k) {
          if (impl_ATS::magnitude(row.value(k)) > tol) {
            partial_nnz += 1;
          }
        }
        if (is_final)
          rowptr(rlid + 1) = partial_nnz;
      },
      nnz);

  auto idx  = col_idx_type("idx", nnz);
  auto vals = vals_type("vals", nnz);

  Kokkos::parallel_for(
      "removeSmallEntries::indicesValues",
      Kokkos::RangePolicy<LocalOrdinal>(0, lclA.numRows()),
      KOKKOS_LAMBDA(const LocalOrdinal rlid) {
        auto row = lclA.row(rlid);
        auto I   = rowptr(rlid);
        for (LocalOrdinal k = 0; k < row.length; ++k) {
          if (impl_ATS::magnitude(row.value(k)) > tol) {
            idx(I)  = row.colidx(k);
            vals(I) = row.value(k);
            I += 1;
          }
        }
      });

  Kokkos::fence();

  auto newA = Teuchos::rcp(new crs_matrix(A->getRowMap(), A->getColMap(), rowptr, idx, vals));
  newA->fillComplete(A->getDomainMap(),
                     A->getRangeMap());
  return newA;
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
constructSubMatrix(const Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& A,
                   std::pair<LocalOrdinal, LocalOrdinal> lid_range) {
  using crs_matrix        = Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using local_graph_type  = typename crs_matrix::local_graph_device_type;
  using local_matrix_type = typename crs_matrix::local_matrix_device_type;
  using row_ptr_type      = typename local_graph_type::row_map_type::non_const_type;
  using col_idx_type      = typename local_graph_type::entries_type::non_const_type;
  using vals_type         = typename local_matrix_type::values_type;

  auto lclA = A->getLocalMatrixDevice();

  auto old_rowptr = lclA.graph.row_map;
  auto rowptr     = row_ptr_type("rowptr", lclA.numRows() + 1);

  LocalOrdinal nnz;
  Kokkos::parallel_scan(
      "constructSubMatrix::rowptr",
      Kokkos::RangePolicy<LocalOrdinal>(0, lclA.numRows()),
      KOKKOS_LAMBDA(const LocalOrdinal rlid, LocalOrdinal& partial_nnz, bool is_final) {
        if ((lid_range.first <= rlid) && (rlid < lid_range.second)) {
          partial_nnz += old_rowptr(rlid + 1) - old_rowptr(rlid);
        }
        if (is_final)
          rowptr(rlid + 1) = partial_nnz;
      },
      nnz);

  typename local_graph_type::size_type start;
  typename local_graph_type::size_type end;
  Kokkos::deep_copy(start, Kokkos::subview(old_rowptr, lid_range.first));
  Kokkos::deep_copy(end, Kokkos::subview(old_rowptr, lid_range.second));
  auto vals_range = Kokkos::make_pair(start, end);

  TEUCHOS_ASSERT_EQUALITY(nnz, static_cast<LocalOrdinal>(end - start));

  auto idx  = Kokkos::subview(lclA.graph.entries, vals_range);
  auto vals = Kokkos::subview(lclA.values, vals_range);

  auto newA = Teuchos::rcp(new crs_matrix(A->getRowMap(), A->getColMap(), rowptr, idx, vals));
  newA->fillComplete(A->getDomainMap(),
                     A->getRangeMap());
  return newA;
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
Teuchos::RCP<Tpetra::BlockedMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
transpose(Teuchos::RCP<Tpetra::BlockedMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& A) {
  TEUCHOS_ASSERT(A->ghosted_blockMap_.is_null());

  Teuchos::RCP<Teuchos::ParameterList> transposeParams = rcp(new Teuchos::ParameterList);

  Tpetra::RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node> transposerPoint(A->pointA_);
  auto pointAT = transposerPoint.createTranspose(transposeParams);

  Tpetra::RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node> transposerBlock(A->blockA_);
  auto blockAT = transposerBlock.createTranspose(transposeParams);

  auto AT = Teuchos::rcp(new Tpetra::BlockedMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(pointAT, blockAT, A->blockMap_));

  return AT;
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
buildIdentityMatrix(Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> >& map) {
  using matrix_type                               = Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  Teuchos::RCP<matrix_type> identity              = Teuchos::rcp(new matrix_type(map, 1));
  Teuchos::ArrayView<const GlobalOrdinal> gblRows = map->getLocalElementList();
  for (auto it = gblRows.begin(); it != gblRows.end(); ++it) {
    Teuchos::Array<GlobalOrdinal> col(1, *it);
    Teuchos::Array<Scalar> val(1, Teuchos::ScalarTraits<Scalar>::one());
    identity->insertGlobalValues(*it, col(), val());
  }
  identity->fillComplete();
  return identity;
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
HierarchicalOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    HierarchicalOperator(const Teuchos::RCP<matrix_type>& nearField,
                         const Teuchos::RCP<blocked_matrix_type>& kernelApproximations,
                         const Teuchos::RCP<matrix_type>& basisMatrix,
                         std::vector<Teuchos::RCP<blocked_matrix_type> >& transferMatrices,
                         const Teuchos::RCP<Teuchos::ParameterList>& params)
  : nearField_(nearField)
  , kernelApproximations_(kernelApproximations)
  , basisMatrix_(basisMatrix)
  , transferMatrices_(transferMatrices)
  , params_(params) {
  auto map         = nearField_->getDomainMap();
  clusterCoeffMap_ = basisMatrix_->getDomainMap();

  bool setupTransposes;
  bool doDebugChecks;
  std::string sendTypeNearField;
  std::string sendTypeBasisMatrix;
  std::string sendTypeKernelApproximations;

  Teuchos::ParameterList defaultParams("Default params");
  defaultParams.set("setupTransposes", true);
  defaultParams.set("doDebugChecks", true);
  defaultParams.set("Send type nearField", "Isend");
  defaultParams.set("Send type basisMatrix", "Isend");
  defaultParams.set("Send type kernelApproximations", "Alltoall");
  defaultParams.set("Coarsening criterion", "transferLevels");
  defaultParams.set("debugOutput", false);
  defaultParams.set("keepTransfers", -1);
  defaultParams.set("treeCoarseningFactor", 2.0);
  defaultParams.set("leftOverFactor", 1.0);
  defaultParams.set("batchSize", 50);
  defaultParams.set("numBatches", 36);
  if (params_.is_null())
    params_ = Teuchos::rcp(new Teuchos::ParameterList(""));
  params_->validateParametersAndSetDefaults(defaultParams);

  setupTransposes              = params_->get<bool>("setupTransposes");
  doDebugChecks                = params_->get<bool>("doDebugChecks");
  sendTypeNearField            = params_->get<std::string>("Send type nearField");
  sendTypeBasisMatrix          = params_->get<std::string>("Send type basisMatrix");
  sendTypeKernelApproximations = params_->get<std::string>("Send type kernelApproximations");
  coarseningCriterion_         = params_->get<std::string>("Coarsening criterion");
  TEUCHOS_ASSERT((coarseningCriterion_ == "numClusters") || (coarseningCriterion_ == "equivalentDense") || (coarseningCriterion_ == "transferLevels"));
  setDebugOutput(params_->get<bool>("debugOutput"));

  auto comm = getComm();
  if (debugOutput_ && (comm->getRank() == 0))
    std::cout << *params_ << std::endl;

  if (doDebugChecks) {
    // near field matrix lives on map and is nonlocal
    TEUCHOS_ASSERT(map->isSameAs(*nearField_->getRangeMap()));
    TEUCHOS_ASSERT(map->isSameAs(*nearField_->getRowMap()));

    // basis matrix is entirely local and maps from clusterCoeffMap_ to map.
    TEUCHOS_ASSERT(map->isSameAs(*basisMatrix->getRangeMap()));
    TEUCHOS_ASSERT(map->isSameAs(*basisMatrix->getRowMap()));
    // TEUCHOS_ASSERT(clusterCoeffMap_->isSameAs(*basisMatrix->getDomainMap()));

    // kernel approximations live on clusterCoeffMap and are nonlocal
    TEUCHOS_ASSERT(clusterCoeffMap_->isSameAs(*kernelApproximations_->pointA_->getDomainMap()));
    TEUCHOS_ASSERT(clusterCoeffMap_->isSameAs(*kernelApproximations_->pointA_->getRangeMap()));
    TEUCHOS_ASSERT(clusterCoeffMap_->isSameAs(*kernelApproximations_->pointA_->getRowMap()));

    for (size_t i = 0; i < transferMatrices_.size(); i++) {
      // transfer matrices are entirely local, block diagonal on clusterCoeffMap
      TEUCHOS_ASSERT(clusterCoeffMap_->isSameAs(*transferMatrices_[i]->pointA_->getDomainMap()));
      TEUCHOS_ASSERT(clusterCoeffMap_->isSameAs(*transferMatrices_[i]->pointA_->getColMap()));
      TEUCHOS_ASSERT(clusterCoeffMap_->isSameAs(*transferMatrices_[i]->pointA_->getRowMap()));
      TEUCHOS_ASSERT(clusterCoeffMap_->isSameAs(*transferMatrices_[i]->pointA_->getRangeMap()));
    }
  }

  // Set the send types
  Teuchos::RCP<Teuchos::ParameterList> distParams = rcp(new Teuchos::ParameterList());
  {
    distParams->set("Send type", sendTypeNearField);
    Teuchos::RCP<const Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > nearFieldImporter = nearField_->getGraph()->getImporter();
    nearFieldImporter->getDistributor().setParameterList(distParams);
    auto revDistor = nearFieldImporter->getDistributor().getReverse(false);
    if (!revDistor.is_null())
      revDistor->setParameterList(distParams);
  }

  {
    distParams->set("Send type", sendTypeBasisMatrix);
    Teuchos::RCP<const Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > basisMatrixImporter = basisMatrix_->getGraph()->getImporter();
    if (!basisMatrixImporter.is_null()) {
      basisMatrixImporter->getDistributor().setParameterList(distParams);
      auto revDistor = basisMatrixImporter->getDistributor().getReverse(false);
      if (!revDistor.is_null())
        revDistor->setParameterList(distParams);
    }
  }

  {
    distParams->set("Send type", sendTypeKernelApproximations);
    Teuchos::RCP<const Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > kernelApproximationsImporter = kernelApproximations_->pointA_->getGraph()->getImporter();
    kernelApproximationsImporter->getDistributor().setParameterList(distParams);
    auto revDistor = kernelApproximationsImporter->getDistributor().getReverse(false);
    if (!revDistor.is_null())
      revDistor->setParameterList(distParams);
  }

  if (setupTransposes) {
    Teuchos::RCP<Teuchos::ParameterList> transposeParams = rcp(new Teuchos::ParameterList);

    Tpetra::RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node> transposerBasisMatrix(basisMatrix_);
    basisMatrixT_ = transposerBasisMatrix.createTranspose(transposeParams);

    for (size_t i = 0; i < transferMatrices_.size(); i++) {
      transferMatricesT_.push_back(transpose(transferMatrices_[i]));
    }

    canApplyWithoutTransposes_ = true;
  } else
    canApplyWithoutTransposes_ = false;

  // Allocate memory for apply with vectors
  allocateMemory(1);
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void HierarchicalOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    apply(const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
          Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Y,
          Teuchos::ETransp mode,
          Scalar alpha,
          Scalar beta) const {
  if (canApplyWithoutTransposes_)
    applyWithoutTransposes(X, Y, mode, alpha, beta);
  else
    applyWithTransposes(X, Y, mode, alpha, beta);
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void HierarchicalOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    applyWithTransposes(const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
                        Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Y,
                        Teuchos::ETransp mode,
                        Scalar alpha,
                        Scalar beta) const {
  using Teuchos::RCP;
  const Scalar one  = Teuchos::ScalarTraits<Scalar>::one();
  const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
  bool flip         = true;

  allocateMemory(X.getNumVectors());

  // upward pass
  {
    Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("upward pass")));

    basisMatrix_->apply(X, *coefficients_, Teuchos::TRANS);

    for (int i = Teuchos::as<int>(transferMatrices_.size()) - 1; i >= 0; i--)
      if (flip) {
        coefficients2_->assign(*coefficients_);
        transferMatrices_[i]->localApply(*coefficients_, *coefficients2_, Teuchos::NO_TRANS, one, one);
        flip = false;
      } else {
        coefficients_->assign(*coefficients2_);
        transferMatrices_[i]->localApply(*coefficients2_, *coefficients_, Teuchos::NO_TRANS, one, one);
        flip = true;
      }
  }

  // far field interactions - part 1
  {
    Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("far field 1")));

    RCP<const Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > kernelApproximationsImporter = kernelApproximations_->pointA_->getGraph()->getImporter();
    if (flip) {
      if (mode == Teuchos::NO_TRANS) {
        coefficients_colmap_->beginImport(*coefficients_, *kernelApproximationsImporter, INSERT);
      } else if (mode == Teuchos::TRANS) {
        kernelApproximations_->localApply(*coefficients_, *coefficients_colmap_, mode, alpha);
        coefficients2_->putScalar(zero);
        coefficients2_->beginExport(*coefficients_colmap_, *kernelApproximationsImporter, ADD_ASSIGN);
      }
    } else {
      if (mode == Teuchos::NO_TRANS) {
        coefficients_colmap_->beginImport(*coefficients2_, *kernelApproximationsImporter, INSERT);
      } else if (mode == Teuchos::TRANS) {
        kernelApproximations_->localApply(*coefficients2_, *coefficients_colmap_, mode, alpha);
        coefficients_->putScalar(zero);
        coefficients_->beginExport(*coefficients_colmap_, *kernelApproximationsImporter, ADD_ASSIGN);
      }
    }
  }

  // near field - part 1
  RCP<const Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > nearFieldImporter = nearField_->getGraph()->getImporter();
  {
    Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("near field 1")));
    if (mode == Teuchos::NO_TRANS) {
      X_colmap_->beginImport(X, *nearFieldImporter, INSERT);
    } else if (mode == Teuchos::TRANS) {
      nearField_->localApply(X, *X_colmap_, mode, alpha, zero);
      Y.scale(beta);
      Y.beginExport(*X_colmap_, *nearFieldImporter, ADD_ASSIGN);
    }
  }

  // far field interactions - part 2
  {
    Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("far field 2")));

    RCP<const Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > kernelApproximationsImporter = kernelApproximations_->pointA_->getGraph()->getImporter();
    if (flip) {
      if (mode == Teuchos::NO_TRANS) {
        coefficients_colmap_->endImport(*coefficients_, *kernelApproximationsImporter, INSERT);
        kernelApproximations_->localApply(*coefficients_colmap_, *coefficients2_, mode, alpha);
      } else if (mode == Teuchos::TRANS) {
        coefficients2_->endExport(*coefficients_colmap_, *kernelApproximationsImporter, ADD_ASSIGN);
      }
    } else {
      if (mode == Teuchos::NO_TRANS) {
        coefficients_colmap_->endImport(*coefficients2_, *kernelApproximationsImporter, INSERT);
        kernelApproximations_->localApply(*coefficients_colmap_, *coefficients_, mode, alpha);
      } else if (mode == Teuchos::TRANS) {
        coefficients_->endExport(*coefficients_colmap_, *kernelApproximationsImporter, ADD_ASSIGN);
      }
    }
  }

  // near field - part 2
  {
    Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("near field 2")));

    if (mode == Teuchos::NO_TRANS) {
      X_colmap_->endImport(X, *nearFieldImporter, INSERT);
      nearField_->localApply(*X_colmap_, Y, mode, alpha, beta);
    } else if (mode == Teuchos::TRANS) {
      Y.endExport(*X_colmap_, *nearFieldImporter, ADD_ASSIGN);
    }
  }

  // downward pass
  {
    Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("downward pass")));

    for (size_t i = 0; i < transferMatrices_.size(); i++)
      if (flip) {
        coefficients_->assign(*coefficients2_);
        transferMatrices_[i]->localApply(*coefficients2_, *coefficients_, Teuchos::TRANS, one, one);
        flip = false;
      } else {
        coefficients2_->assign(*coefficients_);
        transferMatrices_[i]->localApply(*coefficients_, *coefficients2_, Teuchos::TRANS, one, one);
        flip = true;
      }
    if (flip)
      basisMatrix_->apply(*coefficients2_, Y, Teuchos::NO_TRANS, one, one);
    else
      basisMatrix_->apply(*coefficients_, Y, Teuchos::NO_TRANS, one, one);
  }
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void HierarchicalOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    applyWithoutTransposes(const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
                           Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Y,
                           Teuchos::ETransp mode,
                           Scalar alpha,
                           Scalar beta) const {
  using Teuchos::RCP;
  const Scalar one  = Teuchos::ScalarTraits<Scalar>::one();
  const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
  bool flip         = true;

  allocateMemory(X.getNumVectors());

  // upward pass
  {
    Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("upward pass")));

    basisMatrixT_->apply(X, *coefficients_, Teuchos::NO_TRANS);

    for (int i = Teuchos::as<int>(transferMatrices_.size()) - 1; i >= 0; i--)
      if (flip) {
        coefficients2_->assign(*coefficients_);
        transferMatrices_[i]->localApply(*coefficients_, *coefficients2_, Teuchos::NO_TRANS, one, one);
        flip = false;
      } else {
        coefficients_->assign(*coefficients2_);
        transferMatrices_[i]->localApply(*coefficients2_, *coefficients_, Teuchos::NO_TRANS, one, one);
        flip = true;
      }
  }

  // far field interactions - part 1
  {
    Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("far field 1")));

    RCP<const Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > kernelApproximationsImporter = kernelApproximations_->pointA_->getGraph()->getImporter();
    if (flip) {
      if (mode == Teuchos::NO_TRANS) {
        coefficients_colmap_->beginImport(*coefficients_, *kernelApproximationsImporter, INSERT);
      } else if (mode == Teuchos::TRANS) {
        kernelApproximations_->localApply(*coefficients_, *coefficients_colmap_, mode, alpha);
        coefficients2_->putScalar(zero);
        coefficients2_->beginExport(*coefficients_colmap_, *kernelApproximationsImporter, ADD_ASSIGN);
      }
    } else {
      if (mode == Teuchos::NO_TRANS) {
        coefficients_colmap_->beginImport(*coefficients2_, *kernelApproximationsImporter, INSERT);
      } else if (mode == Teuchos::TRANS) {
        kernelApproximations_->localApply(*coefficients2_, *coefficients_colmap_, mode, alpha);
        coefficients_->putScalar(zero);
        coefficients_->beginExport(*coefficients_colmap_, *kernelApproximationsImporter, ADD_ASSIGN);
      }
    }
  }

  // near field - part 1
  RCP<const Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > nearFieldImporter = nearField_->getGraph()->getImporter();
  {
    Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("near field 1")));
    if (mode == Teuchos::NO_TRANS) {
      X_colmap_->beginImport(X, *nearFieldImporter, INSERT);
    } else if (mode == Teuchos::TRANS) {
      nearField_->localApply(X, *X_colmap_, mode, alpha, zero);
      Y.scale(beta);
      Y.beginExport(*X_colmap_, *nearFieldImporter, ADD_ASSIGN);
    }
  }

  // far field interactions - part 2
  {
    Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("far field 2")));

    RCP<const Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > kernelApproximationsImporter = kernelApproximations_->pointA_->getGraph()->getImporter();
    if (flip) {
      if (mode == Teuchos::NO_TRANS) {
        coefficients_colmap_->endImport(*coefficients_, *kernelApproximationsImporter, INSERT);
        kernelApproximations_->localApply(*coefficients_colmap_, *coefficients2_, mode, alpha);
      } else if (mode == Teuchos::TRANS) {
        coefficients2_->endExport(*coefficients_colmap_, *kernelApproximationsImporter, ADD_ASSIGN);
      }
    } else {
      if (mode == Teuchos::NO_TRANS) {
        coefficients_colmap_->endImport(*coefficients2_, *kernelApproximationsImporter, INSERT);
        kernelApproximations_->localApply(*coefficients_colmap_, *coefficients_, mode, alpha);
      } else if (mode == Teuchos::TRANS) {
        coefficients_->endExport(*coefficients_colmap_, *kernelApproximationsImporter, ADD_ASSIGN);
      }
    }
  }

  // near field - part 2
  {
    Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("near field 2")));

    if (mode == Teuchos::NO_TRANS) {
      X_colmap_->endImport(X, *nearFieldImporter, INSERT);
      nearField_->localApply(*X_colmap_, Y, mode, alpha, beta);
    } else if (mode == Teuchos::TRANS) {
      Y.endExport(*X_colmap_, *nearFieldImporter, ADD_ASSIGN);
    }
  }

  // downward pass
  {
    Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("downward pass")));

    for (size_t i = 0; i < transferMatrices_.size(); i++)
      if (flip) {
        coefficients_->assign(*coefficients2_);
        transferMatricesT_[i]->localApply(*coefficients2_, *coefficients_, Teuchos::NO_TRANS, one, one);
        flip = false;
      } else {
        coefficients2_->assign(*coefficients_);
        transferMatricesT_[i]->localApply(*coefficients_, *coefficients2_, Teuchos::NO_TRANS, one, one);
        flip = true;
      }
    if (flip)
      basisMatrix_->apply(*coefficients2_, Y, Teuchos::NO_TRANS, one, one);
    else
      basisMatrix_->apply(*coefficients_, Y, Teuchos::NO_TRANS, one, one);
  }
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
Teuchos::RCP<HierarchicalOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    HierarchicalOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
        restrict(const Teuchos::RCP<matrix_type>& P) {
  Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("Galerkin product")));

  // H_c = P^T * H * P
  using lo_vec_type   = typename blocked_map_type::lo_vec_type;
  using vec_type      = typename Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using MagnitudeType = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
  using Teuchos::RCP;
  using Teuchos::rcp;
  const Scalar ONE  = Teuchos::ScalarTraits<Scalar>::one();
  const Scalar ZERO = Teuchos::ScalarTraits<Scalar>::zero();
  const Scalar HALF = ONE / (ONE + ONE);

  RCP<Teuchos::ParameterList> coarseParams = rcp(new Teuchos::ParameterList(*params_));

  // newBasisMatrix = P^T * basisMatrix
  RCP<matrix_type> newBasisMatrix = rcp(new matrix_type(P->getDomainMap(), clusterCoeffMap_, 0));
  MatrixMatrix::Multiply(*P, true, *basisMatrix_, false, *newBasisMatrix);

  //
  auto clusterSizes         = kernelApproximations_->blockMap_->blockSizes_;
  auto ghosted_clusterMap   = kernelApproximations_->blockA_->getColMap();
  auto ghosted_clusterSizes = kernelApproximations_->ghosted_blockMap_->blockSizes_;

  // Get number of unknowns associated with each cluster via new basisMatrix.
  // numUnknownsPerCluster = \prod transfer_k * graph(newBasisMatrix)^T * ones
  RCP<vec_type> numUnknownsPerCluster;
  RCP<vec_type> ghosted_numUnknownsPerCluster;
  if ((coarseningCriterion_ == "equivalentDense") ||
      (coarseningCriterion_ == "numClusters")) {
    {
      numUnknownsPerCluster          = rcp(new vec_type(kernelApproximations_->blockA_->getRowMap(), false));
      auto lcl_clusterSizes          = clusterSizes->getLocalViewHost(Tpetra::Access::ReadOnly);
      auto lcl_numUnknownsPerCluster = numUnknownsPerCluster->getLocalViewHost(Tpetra::Access::OverwriteAll);
      // Compute the transpose of the newBasisMatrix.
      RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > newBasisMatrixT;
      Tpetra::RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node> transposer(newBasisMatrix);
      RCP<Teuchos::ParameterList> transposeParams = rcp(new Teuchos::ParameterList);
      newBasisMatrixT                             = transposer.createTranspose(transposeParams);

      // TODO: parallel_for
      auto rowptr               = newBasisMatrixT->getLocalRowPtrsHost();
      LocalOrdinal clusterStart = 0;
      LocalOrdinal clusterEnd   = 0;
      for (LocalOrdinal cluster = 0; cluster < lcl_clusterSizes.extent_int(0); ++cluster) {
        clusterStart = clusterEnd;
        clusterEnd += lcl_clusterSizes(cluster, 0);
        LocalOrdinal maxEntries = 0;
        for (LocalOrdinal row = clusterStart; row < clusterEnd; ++row) {
          LocalOrdinal numEntriesPerRow = rowptr(row + 1) - rowptr(row);
          maxEntries                    = std::max(maxEntries, numEntriesPerRow);
        }
        lcl_numUnknownsPerCluster(cluster, 0) = maxEntries;
      }
      TEUCHOS_ASSERT_EQUALITY(clusterEnd + 1, rowptr.extent_int(0));
    }

    // sum from child nodes to parents via transfer operators
    for (int i = Teuchos::as<int>(transferMatrices_.size()) - 1; i >= 0; i--)
      transferMatrices_[i]->blockA_->apply(*numUnknownsPerCluster, *numUnknownsPerCluster, Teuchos::NO_TRANS, ONE, ONE);

    // get ghosted numUnknownsPerCluster
    ghosted_numUnknownsPerCluster = rcp(new vec_type(ghosted_clusterMap, false));
    auto import                   = kernelApproximations_->blockA_->getCrsGraph()->getImporter();
    ghosted_numUnknownsPerCluster->doImport(*numUnknownsPerCluster, *import, Tpetra::INSERT);
  }

  // coarse cluster pair graph
  RCP<matrix_type> newKernelBlockGraph;
  // point entries of cluster pairs that should be moved to the near field
  RCP<matrix_type> diffKernelApprox;
  // coarse point matrix of cluster pairs
  Teuchos::RCP<matrix_type> newKernelApprox;

  // Determine which cluster pairs should be moved to the near field.
  // We are constructing the coarse block matrix newKernelBlockGraph
  // and the point matrix diffKernelApprox.
  {
    typename vec_type::dual_view_type::t_host::const_type lcl_numUnknownsPerCluster;
    typename vec_type::dual_view_type::t_host::const_type lcl_ghosted_numUnknownsPerCluster;
    auto lcl_offsets         = Kokkos::create_mirror_view(kernelApproximations_->blockMap_->offsets_);
    auto lcl_ghosted_offsets = Kokkos::create_mirror_view(kernelApproximations_->ghosted_blockMap_->offsets_);
    Kokkos::deep_copy(lcl_offsets, kernelApproximations_->blockMap_->offsets_);
    Kokkos::deep_copy(lcl_ghosted_offsets, kernelApproximations_->ghosted_blockMap_->offsets_);
    auto lcl_clusterSizes         = clusterSizes->getLocalViewHost(Tpetra::Access::ReadOnly);
    auto lcl_ghosted_clusterSizes = ghosted_clusterSizes->getLocalViewHost(Tpetra::Access::ReadOnly);

    if ((coarseningCriterion_ == "equivalentDense") ||
        (coarseningCriterion_ == "numClusters")) {
      lcl_numUnknownsPerCluster         = numUnknownsPerCluster->getLocalViewHost(Tpetra::Access::ReadOnly);
      lcl_ghosted_numUnknownsPerCluster = ghosted_numUnknownsPerCluster->getLocalViewHost(Tpetra::Access::ReadOnly);
    }

    // Criterion: "numClusters"
    // Compute all cluster pair sizes, sort them, and pick cut-off
    // so that the number of cluster pairs decreases propotionally
    // to the number of unknowns.
    size_t tgt_clusterPairSize = 0;
    if (coarseningCriterion_ == "numClusters") {
      auto lcl_BlockGraph = kernelApproximations_->blockA_->getLocalMatrixHost();
      std::vector<size_t> clusterPairSizes;
      for (LocalOrdinal brlid = 0; brlid < lcl_BlockGraph.numRows(); ++brlid) {
        auto brow = lcl_BlockGraph.row(brlid);
        for (LocalOrdinal k = 0; k < brow.length; ++k) {
          // Entries of the block matrix for kernelApproximations
          // decide whether the cluster pair is present and only take
          // values 1 or 0.
          if (brow.value(k) > HALF) {
            LocalOrdinal bclid = brow.colidx(k);
            clusterPairSizes.push_back(lcl_numUnknownsPerCluster(brlid, 0) * lcl_ghosted_numUnknownsPerCluster(bclid, 0));
          }
        }
      }
      std::sort(clusterPairSizes.begin(), clusterPairSizes.end());
      double coarseningRate = Teuchos::as<double>(P->getGlobalNumCols()) / Teuchos::as<double>(P->getGlobalNumRows());
      tgt_clusterPairSize   = clusterPairSizes[Teuchos::as<size_t>(static_cast<double>(clusterPairSizes.size()) * (1 - coarseningRate))];
      // std::cout << "HERE " << clusterPairSizes[0] << " " << tgt_clusterPairSize << " " << clusterPairSizes[clusterPairSizes.size()-1] << std::endl;
    }

    // Criterion: "transferLevels"
    // Drop cluster pairs by level in the tree.
    auto comm = getComm();
    std::set<LocalOrdinal> blidsToDrop;
    LocalOrdinal minDrop = kernelApproximations_->blockMap_->blockMap_->getLocalNumElements() + 1;
    LocalOrdinal maxDrop = -1;
    if (coarseningCriterion_ == "transferLevels") {
      double coarseningRate       = Teuchos::as<double>(P->getGlobalNumCols()) / Teuchos::as<double>(P->getGlobalNumRows());
      size_t droppedClusterPairs  = 0;
      size_t totalNumClusterPairs = kernelApproximations_->blockA_->getGlobalNumEntries();
      RCP<vec_type> tempV         = Teuchos::rcp(new vec_type(kernelApproximations_->blockMap_->blockMap_, false));
      RCP<vec_type> tempV2        = Teuchos::rcp(new vec_type(kernelApproximations_->blockMap_->blockMap_, false));
      // keepTransfers == transferMatrices_.size(): keep all transfers
      // keepTransfers == 0: keep no transfers
      int keepTransfers                 = params_->get<int>("keepTransfers", -1);
      const double treeCoarseningFactor = params_->get<double>("treeCoarseningFactor");
      if (keepTransfers == -1) {
        double leftOverFactor = params_->get<double>("leftOverFactor");
        keepTransfers         = transferMatrices_.size();
        double temp           = (1.0 / coarseningRate) * leftOverFactor;
        while (temp >= 1.0) {
          --keepTransfers;
          temp /= treeCoarseningFactor;
        }
        keepTransfers = std::max(keepTransfers, 0);
        coarseParams->set("leftOverFactor", temp);
      }
      coarseParams->set("keepTransfers", -1);
      TEUCHOS_ASSERT((0 <= keepTransfers) && (keepTransfers <= Teuchos::as<int>(transferMatrices_.size())));

      size_t droppedTransfers = 0;
      for (int k = Teuchos::as<int>(transferMatrices_.size()) - 1; k >= 0; --k) {
        size_t numTreeEdgesBetweenLevels = transferMatrices_[k]->blockA_->getGlobalNumEntries();
        size_t numClusters_k1            = numTreeEdgesBetweenLevels;
        auto numClusters_k               = Teuchos::as<size_t>(static_cast<double>(numTreeEdgesBetweenLevels) / treeCoarseningFactor);

        if (debugOutput_ && (comm->getRank() == 0))
          std::cout << "transfer " << k << " between levels " << k + 1 << " and " << k << " maps " << numClusters_k1 << " to " << numClusters_k << " clusters" << std::endl;

        tempV->putScalar(ONE);
        transferMatrices_[k]->blockA_->apply(*tempV, *tempV2, Teuchos::TRANS);
        tempV->putScalar(ZERO);
        kernelApproximations_->blockA_->apply(*tempV2, *tempV);
        Scalar numClusterPairs = tempV->dot(*tempV2);
        if (debugOutput_ && (comm->getRank() == 0))
          std::cout << "cluster pairs on level " << k + 1 << ": " << numClusterPairs << std::endl;

        bool doDrop;
        if (keepTransfers >= 0) {
          doDrop = (keepTransfers <= k);
        } else {
          doDrop = (droppedClusterPairs + numClusterPairs < (1.0 - coarseningRate) * static_cast<double>(totalNumClusterPairs));
        }
        if (doDrop) {
          auto lcl_transfer       = transferMatrices_[k]->blockA_->getLocalMatrixHost();
          auto lcl_transfer_graph = lcl_transfer.graph;
          for (LocalOrdinal j = 0; j < lcl_transfer_graph.entries.extent_int(0); j++) {
            blidsToDrop.insert(lcl_transfer_graph.entries(j));
            minDrop = std::min(minDrop, lcl_transfer_graph.entries(j));
            maxDrop = std::max(maxDrop, lcl_transfer_graph.entries(j));
          }

          droppedTransfers += 1;
          droppedClusterPairs += numClusterPairs;
        }
      }
      if (debugOutput_ && (comm->getRank() == 0))
        std::cout << "Dropped " << droppedTransfers << " transfers of " << transferMatrices_.size() << ", dropped cluster pairs: " << droppedClusterPairs << std::endl;
    }

    bool dropContiguousRange = (static_cast<size_t>(maxDrop + 1 - minDrop) == blidsToDrop.size());

    // number of cluster pairs dropped
    int dropped = 0;
    // number of cluster pairs we kept
    int kept = 0;
    // number of cluster pairs that were no longer present
    int ignored = 0;

    if (dropContiguousRange) {
      newKernelBlockGraph = constructSubMatrix(kernelApproximations_->blockA_, {0, minDrop});
      dropped             = maxDrop + 1 - minDrop;
      kept                = minDrop;
      LocalOrdinal minLid = lcl_offsets(minDrop);
      LocalOrdinal maxLid = lcl_offsets(maxDrop + 1);
      diffKernelApprox    = constructSubMatrix(kernelApproximations_->pointA_, {minLid, maxLid});
      newKernelApprox     = constructSubMatrix(kernelApproximations_->pointA_, {0, minLid});

    } else {
      newKernelBlockGraph = rcp(new matrix_type(kernelApproximations_->blockA_->getCrsGraph()));
      newKernelBlockGraph->resumeFill();
      diffKernelApprox = rcp(new matrix_type(kernelApproximations_->pointA_->getCrsGraph()));

      // loop over cluster pairs
      // TODO: parallel_for
      auto lcl_BlockGraph       = kernelApproximations_->blockA_->getLocalMatrixHost();
      auto lcl_newBlockGraph    = newKernelBlockGraph->getLocalMatrixHost();
      auto lcl_KernelApprox     = kernelApproximations_->pointA_->getLocalMatrixHost();
      auto lcl_diffKernelApprox = diffKernelApprox->getLocalMatrixHost();
      for (LocalOrdinal brlid = 0; brlid < lcl_BlockGraph.numRows(); ++brlid) {
        size_t brsize = lcl_clusterSizes(brlid, 0);
        auto brow     = lcl_BlockGraph.row(brlid);
        auto new_brow = lcl_newBlockGraph.row(brlid);
        for (LocalOrdinal k = 0; k < brow.length; ++k) {
          // Entries of the block matrix for kernelApproximations
          // decide whether the cluster pair is present and only take
          // values 1 or 0.
          if (brow.value(k) > HALF) {
            LocalOrdinal bclid = brow.colidx(k);
            size_t bcsize      = lcl_ghosted_clusterSizes(bclid, 0);

            // criterium for removing a cluster pair from the far field
            bool removeCluster = false;
            if (coarseningCriterion_ == "equivalentDense") {
              // Size of the sparse cluster approximation >= size of dense equivalent
              removeCluster = (brsize * bcsize >= lcl_numUnknownsPerCluster(brlid, 0) * lcl_ghosted_numUnknownsPerCluster(bclid, 0));
            } else if (coarseningCriterion_ == "numClusters") {
              removeCluster = (lcl_numUnknownsPerCluster(brlid, 0) * lcl_ghosted_numUnknownsPerCluster(bclid, 0) < tgt_clusterPairSize);
            } else if (coarseningCriterion_ == "transferLevels") {
              removeCluster = ((blidsToDrop.find(brlid) != blidsToDrop.end()) ||
                               (blidsToDrop.find(bclid) != blidsToDrop.end()));
            }
            if (removeCluster) {
              // we are dropping the cluster pair from the far field
              ++dropped;
              new_brow.value(k) = ZERO;

              // loop over the point matrix and add the entries to diffKernelApprox
              const LocalOrdinal row_start = lcl_offsets(brlid);
              const LocalOrdinal row_end   = lcl_offsets(brlid + 1);
              const LocalOrdinal col_start = lcl_ghosted_offsets(bclid);
              const LocalOrdinal col_end   = lcl_ghosted_offsets(bclid + 1);
              TEUCHOS_ASSERT_EQUALITY(Teuchos::as<size_t>(row_end - row_start), brsize);
              TEUCHOS_ASSERT_EQUALITY(Teuchos::as<size_t>(col_end - col_start), bcsize);
              for (LocalOrdinal rlid = row_start; rlid < row_end; ++rlid) {
                auto diff_row  = lcl_diffKernelApprox.row(rlid);
                auto row       = lcl_KernelApprox.row(rlid);
                size_t removed = 0;
                for (LocalOrdinal n = 0; n < row.length; ++n) {
                  if ((col_start <= row.colidx(n)) && (col_end > row.colidx(n))) {
                    diff_row.value(n) = row.value(n);
                    ++removed;
                  }
                }
                if (removed != bcsize) {
                  std::ostringstream oss;
                  oss << "brlid " << brlid << " row " << rlid << std::endl;
                  oss << "col_start " << col_start << " col_end " << col_end << std::endl;
                  for (LocalOrdinal n = 0; n < row.length; ++n) {
                    oss << row.colidx(n) << " " << row.value(n) << std::endl;
                  }
                  std::cout << oss.str();
                }
                TEUCHOS_ASSERT_EQUALITY(removed, bcsize);
              }
            } else {
              // We are keeping the cluster pair.
              ++kept;
              new_brow.value(k) = brow.value(k);
            }
          } else {
            // The cluster pair has already been dropped on the fine level.
            ++ignored;
            new_brow.value(k) = brow.value(k);
          }
        }
      }

      newKernelBlockGraph->fillComplete(kernelApproximations_->blockA_->getDomainMap(),
                                        kernelApproximations_->blockA_->getRangeMap());
      newKernelBlockGraph = removeSmallEntries(newKernelBlockGraph, Teuchos::ScalarTraits<MagnitudeType>::eps());
      diffKernelApprox->fillComplete(clusterCoeffMap_,
                                     clusterCoeffMap_);

      {
        Teuchos::RCP<matrix_type> temp = MatrixMatrix::add(ONE, false, *kernelApproximations_->pointA_, -ONE, false, *diffKernelApprox);
        newKernelApprox                = removeSmallEntries(temp, Teuchos::ScalarTraits<MagnitudeType>::eps());
      }
    }
    if (debugOutput_) {
      // number of cluster pairs dropped
      int gbl_dropped = 0;
      // number of cluster pairs we kept
      int gbl_kept = 0;
      // number of cluster pairs that were no longer present
      int gbl_ignored = 0;
      Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &dropped, &gbl_dropped);
      Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &kept, &gbl_kept);
      Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &ignored, &gbl_ignored);
      if (comm->getRank() == 0)
        std::cout << "dropped " << gbl_dropped << " cluster pairs, kept " << gbl_kept << " cluster pairs, ignored " << gbl_ignored << " cluster pairs" << std::endl;
    }
  }

  // construct identity on clusterCoeffMap_
  Teuchos::RCP<matrix_type> identity = buildIdentityMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(clusterCoeffMap_);

  Teuchos::RCP<blocked_matrix_type> newBlockedKernelApproximation = rcp(new blocked_matrix_type(newKernelApprox, newKernelBlockGraph, kernelApproximations_->blockMap_, kernelApproximations_->ghosted_blockMap_));

  // select subset of transfer matrices for coarse operator
  std::vector<Teuchos::RCP<blocked_matrix_type> > newTransferMatrices;
  {
    Teuchos::TimeMonitor tM_basis_transfer(*Teuchos::TimeMonitor::getNewTimer(std::string("Coarse basis and transfers")));

    auto comm = getComm();

    RCP<vec_type> v_temp          = rcp(new vec_type(newKernelBlockGraph->getDomainMap()));
    RCP<vec_type> clusterUseCount = rcp(new vec_type(newKernelBlockGraph->getDomainMap()));
    v_temp->putScalar(ONE);
    clusterUseCount->putScalar(ZERO);
    newKernelBlockGraph->apply(*v_temp, *clusterUseCount, Teuchos::NO_TRANS);
    newKernelBlockGraph->apply(*v_temp, *clusterUseCount, Teuchos::TRANS, ONE, ONE);

    for (int i = Teuchos::as<int>(transferMatrices_.size()) - 1; i >= 0; i--) {
      // We drop a transfer operator T_i when
      //  sum(T_i * clusterUseCount) == 0
      // Since we need to use Scalar, we instead check for < 0.5
      transferMatrices_[i]->blockA_->localApply(*clusterUseCount, *v_temp, Teuchos::NO_TRANS);
      Scalar gbl_use_count = v_temp->norm1();
      // if (comm->getRank() == 0)
      //   std::cout << "Transfer " << i << " count " << gbl_use_count << std::endl;

      if (gbl_use_count < HALF) {
        // We do not keep the i-th transfer for the coarse operator.
        // newBasisMatrix := newBasisMatrix * (I+transferMatrices_[i])^T
        Teuchos::RCP<matrix_type> temp2 = MatrixMatrix::add(ONE, false, *identity, ONE, false, *transferMatrices_[i]->pointA_);
        RCP<matrix_type> temp           = rcp(new matrix_type(newBasisMatrix->getRowMap(), clusterCoeffMap_, 0));
        MatrixMatrix::Multiply(*newBasisMatrix, false, *temp2, true, *temp);
        newBasisMatrix = temp;
      } else {
        // We keep the i-th transfer for the coarse operator.
        newTransferMatrices.insert(newTransferMatrices.begin(), transferMatrices_[i]);
      }
    }
  }

  // Coarse near field
  RCP<matrix_type> newNearField;
  {
    Teuchos::TimeMonitor tM_near(*Teuchos::TimeMonitor::getNewTimer(std::string("Coarse near field")));

    // diffKernelApprox := (identity + newTransferMatrices[K-1])^T * ... * (identity + newTransferMatrices[0])^T
    //                      * diffKernelApprox
    //                      * (identity + newTransferMatrices[0]) * ... * (identity + newTransferMatrices[K-1])
    {
      Teuchos::TimeMonitor tM_near_1(*Teuchos::TimeMonitor::getNewTimer(std::string("Coarse near field 1")));
      for (int i = 0; i < Teuchos::as<int>(newTransferMatrices.size()); ++i) {
        // diffKernelApprox := (I + newTransferMatrices[i])^T * diffKernelApprox * (I + newTransferMatrices[i])
        Teuchos::RCP<matrix_type> temp  = MatrixMatrix::add(ONE, false, *identity, ONE, false, *newTransferMatrices[i]->pointA_);
        Teuchos::RCP<matrix_type> temp2 = rcp(new matrix_type(clusterCoeffMap_, 0));
        MatrixMatrix::Multiply(*temp, true, *diffKernelApprox, false, *temp2);
        MatrixMatrix::Multiply(*temp2, false, *temp, false, *diffKernelApprox);
      }
    }

    // diffKernelApprox := (newBasisMatrix * diffKernelApprox) * newBasisMatrix^T
    {
      Teuchos::TimeMonitor tM_near2(*Teuchos::TimeMonitor::getNewTimer(std::string("Coarse near field 2")));

      Teuchos::RCP<matrix_type> temp = rcp(new matrix_type(newBasisMatrix->getRowMap(), 0));
      {
        Teuchos::TimeMonitor tM_near2a(*Teuchos::TimeMonitor::getNewTimer(std::string("Coarse near field 2a")));
        MatrixMatrix::Multiply(*newBasisMatrix, false, *diffKernelApprox, false, *temp);
      }

      {
        Teuchos::TimeMonitor tM_near2a(*Teuchos::TimeMonitor::getNewTimer(std::string("Coarse near field 2b")));
        diffKernelApprox = rcp(new matrix_type(newBasisMatrix->getRowMap(), 0));
        MatrixMatrix::Multiply(*temp, false, *newBasisMatrix, true, *diffKernelApprox);
      }
    }

    // newNearField = (P^T * nearField * P) + diffKernelApprox
    {
      Teuchos::TimeMonitor tM_near3(*Teuchos::TimeMonitor::getNewTimer(std::string("Coarse near field 3")));
      RCP<matrix_type> temp = rcp(new matrix_type(nearField_->getRowMap(), 0));
      MatrixMatrix::Multiply(*nearField_, false, *P, false, *temp);
      RCP<matrix_type> temp2 = rcp(new matrix_type(P->getDomainMap(), 0));
      MatrixMatrix::Multiply(*P, true, *temp, false, *temp2);
      newNearField     = MatrixMatrix::add(ONE, false, *temp2, ONE, false, *diffKernelApprox);
      diffKernelApprox = Teuchos::null;
      // newNearField = removeSmallEntries(newNearField, Teuchos::ScalarTraits<MagnitudeType>::eps());
    }
  }

  return Teuchos::rcp(new HierarchicalOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>(newNearField,
                                                                                          newBlockedKernelApproximation,
                                                                                          newBasisMatrix,
                                                                                          newTransferMatrices,
                                                                                          coarseParams));
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
HierarchicalOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    toMatrix() {
  Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("Conversion from H-matrix to CSR")));

  using Teuchos::RCP;
  using Teuchos::rcp;

  const Scalar ONE = Teuchos::ScalarTraits<Scalar>::one();

  if (hasFarField()) {
    RCP<matrix_type> kernelApproximations;

    if (hasTransferMatrices()) {
      kernelApproximations = rcp(new matrix_type(*kernelApproximations_->pointA_));

      // construct identity on clusterCoeffMap_
      Teuchos::RCP<matrix_type> identity = buildIdentityMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(clusterCoeffMap_);

      Teuchos::TimeMonitor tM_near_1(*Teuchos::TimeMonitor::getNewTimer(std::string("Densify far field 1")));
      for (int i = 0; i < Teuchos::as<int>(transferMatrices_.size()); ++i) {
        // kernelApproximations := (I + newTransferMatrices[i])^T * kernelApproximations * (I + newTransferMatrices[i])
        Teuchos::RCP<matrix_type> temp  = MatrixMatrix::add(ONE, false, *identity, ONE, false, *transferMatrices_[i]->pointA_);
        Teuchos::RCP<matrix_type> temp2 = rcp(new matrix_type(clusterCoeffMap_, 0));
        MatrixMatrix::Multiply(*temp, true, *kernelApproximations, false, *temp2);
        Teuchos::RCP<matrix_type> temp3 = rcp(new matrix_type(clusterCoeffMap_, 0));
        MatrixMatrix::Multiply(*temp2, false, *temp, false, *temp3);
        kernelApproximations = temp3;
      }
    } else {
      kernelApproximations = kernelApproximations_->pointA_;
    }

    // farField = (basisMatrix_ * kernelApproximations) * basisMatrix_^T
    RCP<matrix_type> farField;
    {
      Teuchos::TimeMonitor tM_near_2(*Teuchos::TimeMonitor::getNewTimer(std::string("Densify far field 2")));
      int rank       = getComm()->getRank();
      int size       = getComm()->getSize();
      int numBatches = params_->get<int>("numBatches");
      if (numBatches < 0) {
        int batchSize = params_->get<int>("batchSize");
        numBatches    = size / batchSize;
      }
      numBatches = std::max(numBatches, 1);
      for (int batchNo = 0; batchNo < numBatches; ++batchNo) {
        RCP<matrix_type> kernelApproximationsSlice;
        {
          Teuchos::TimeMonitor tM_near_2a(*Teuchos::TimeMonitor::getNewTimer(std::string("Densify far field 2 0")));

          using local_matrix_type = typename matrix_type::local_matrix_device_type;
          local_matrix_type lclKernelApproximationsSlice;

          if (rank % numBatches == batchNo) {
            lclKernelApproximationsSlice = kernelApproximations->getLocalMatrixDevice();

          } else {
            using local_graph_type = typename matrix_type::local_graph_device_type;
            using row_ptr_type     = typename local_graph_type::row_map_type::non_const_type;
            using col_idx_type     = typename local_graph_type::entries_type::non_const_type;
            using vals_type        = typename local_matrix_type::values_type;

            auto lclKernelApproximations = kernelApproximations->getLocalMatrixDevice();
            auto rowptr                  = row_ptr_type("rowptr", lclKernelApproximations.numRows() + 1);
            auto idx                     = col_idx_type("colidx", 0);
            auto vals                    = vals_type("vals", 0);

            auto lclKernelApproximationsSliceGraph = local_graph_type(idx, rowptr);
            lclKernelApproximationsSlice           = local_matrix_type("slice", lclKernelApproximations.numCols(), vals, lclKernelApproximationsSliceGraph);
          }
          kernelApproximationsSlice = rcp(new matrix_type(lclKernelApproximationsSlice,
                                                          kernelApproximations->getRowMap(),
                                                          kernelApproximations->getColMap(),
                                                          kernelApproximations->getDomainMap(),
                                                          kernelApproximations->getRangeMap()));
        }
        {
          RCP<matrix_type> temp;
          {
            Teuchos::TimeMonitor tM_near_2a(*Teuchos::TimeMonitor::getNewTimer(std::string("Densify far field 2a")));
            temp = rcp(new matrix_type(basisMatrix_->getRowMap(), 0));
            MatrixMatrix::Multiply(*basisMatrix_, false, *kernelApproximationsSlice, false, *temp);
          }

          {
            Teuchos::TimeMonitor tM_near_2b(*Teuchos::TimeMonitor::getNewTimer(std::string("Densify far field 2b")));
            auto temp2 = rcp(new matrix_type(basisMatrix_->getRowMap(), 0));
            MatrixMatrix::Multiply(*temp, false, *basisMatrix_, true, *temp2);
            temp = Teuchos::null;
            if (batchNo == 0)
              farField = temp2;
            else
              farField = MatrixMatrix::add(ONE, false, *farField, ONE, false, *temp2);
          }
        }
      }
    }

    // nearField_ + farField
    RCP<matrix_type> dense;
    {
      Teuchos::TimeMonitor tM_near_3(*Teuchos::TimeMonitor::getNewTimer(std::string("Densify far field 3")));
      dense = MatrixMatrix::add(ONE, false, *nearField_, ONE, false, *farField);
    }
    return dense;

  } else

    return nearField_;
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void HierarchicalOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    allocateMemory(size_t numVectors) const {
  if (coefficients_.is_null() || coefficients_->getNumVectors() != numVectors) {
    coefficients_        = Teuchos::rcp(new mv_type(clusterCoeffMap_, numVectors));
    coefficients2_       = Teuchos::rcp(new mv_type(clusterCoeffMap_, numVectors));
    X_colmap_            = Teuchos::rcp(new mv_type(nearField_->getColMap(), numVectors));
    coefficients_colmap_ = Teuchos::rcp(new mv_type(kernelApproximations_->pointA_->getColMap(), numVectors));
  }
}
}  // namespace Tpetra

#endif  // TPETRA_HIERARCHICALOPERATOR_DEF_HPP
