#ifndef TPETRA_HIERARCHICALOPERATOR_DEF_HPP
#define TPETRA_HIERARCHICALOPERATOR_DEF_HPP

#include <TpetraExt_MatrixMatrix.hpp>
#include <Tpetra_RowMatrixTransposer.hpp>

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

  typedef Kokkos::ArithTraits<Scalar> ATS;
  using impl_SC  = typename ATS::val_type;
  using impl_ATS = Kokkos::ArithTraits<impl_SC>;

  auto lclA = A->getLocalMatrixDevice();

  auto rowptr = row_ptr_type("rowptr", lclA.numRows() + 1);

  Kokkos::parallel_for(
      "removeSmallEntries::rowptr1",
      Kokkos::RangePolicy<LocalOrdinal>(0, lclA.numRows()),
      KOKKOS_LAMBDA(const LocalOrdinal rlid) {
        auto row = lclA.row(rlid);
        for (LocalOrdinal k = 0; k < row.length; ++k) {
          if (impl_ATS::magnitude(row.value(k)) > tol) {
            rowptr(rlid + 1) += 1;
          }
        }
      });
  LocalOrdinal nnz;
  Kokkos::parallel_scan(
      "removeSmallEntries::rowptr2",
      Kokkos::RangePolicy<LocalOrdinal>(0, lclA.numRows()),
      KOKKOS_LAMBDA(const LocalOrdinal rlid, LocalOrdinal& partial_nnz, bool is_final) {
        partial_nnz += rowptr(rlid + 1);
        if (is_final)
          rowptr(rlid + 1) = partial_nnz;
      },
      nnz);

  // auto nnz = rowptr(lclA.numRows());

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
  debugOutput_ = params_->get<bool>("debugOutput");

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
  RCP<matrix_type> newKernelBlockGraph = rcp(new matrix_type(kernelApproximations_->blockA_->getCrsGraph()));
  newKernelBlockGraph->resumeFill();
  // point entries of cluster pairs that should be moved to the near field
  RCP<matrix_type> diffKernelApprox = rcp(new matrix_type(kernelApproximations_->pointA_->getCrsGraph()));

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
      tgt_clusterPairSize   = clusterPairSizes[Teuchos::as<size_t>(clusterPairSizes.size() * (1 - coarseningRate))];
      // std::cout << "HERE " << clusterPairSizes[0] << " " << tgt_clusterPairSize << " " << clusterPairSizes[clusterPairSizes.size()-1] << std::endl;
    }

    // Criterion: "transferLevels"
    // Drop cluster pairs by level in the tree.
    auto comm = getComm();
    std::set<LocalOrdinal> blidsToDrop;
    if (coarseningCriterion_ == "transferLevels") {
      double coarseningRate       = Teuchos::as<double>(P->getGlobalNumCols()) / Teuchos::as<double>(P->getGlobalNumRows());
      size_t droppedClusterPairs  = 0;
      size_t totalNumClusterPairs = kernelApproximations_->blockA_->getGlobalNumEntries();
      RCP<vec_type> tempV         = Teuchos::rcp(new vec_type(kernelApproximations_->blockMap_->blockMap_, false));
      RCP<vec_type> tempV2        = Teuchos::rcp(new vec_type(kernelApproximations_->blockMap_->blockMap_, false));
      int keepTransfers           = params_->get<int>("keepTransfers", -1);
      if (keepTransfers == -1) {
        double leftOverFactor             = params_->get<double>("leftOverFactor");
        keepTransfers                     = transferMatrices_.size();
        double temp                       = (1.0 / coarseningRate) * leftOverFactor;
        const double treeCoarseningFactor = params_->get<double>("treeCoarseningFactor");
        while (temp >= 2.0) {
          --keepTransfers;
          temp /= treeCoarseningFactor;
        }
        keepTransfers = std::max(keepTransfers, 0);
        params_->set("leftOverFactor", temp);
      }

      for (int k = Teuchos::as<int>(transferMatrices_.size()) - 1; k >= 0; --k) {
        size_t clustersInLevel = transferMatrices_[k]->blockA_->getGlobalNumEntries();

        if (debugOutput_ && (comm->getRank() == 0))
          std::cout << "level " << k << " clustersInLevel " << clustersInLevel << std::endl;

        tempV->putScalar(ONE);
        transferMatrices_[k]->blockA_->apply(*tempV, *tempV2, Teuchos::TRANS);

        size_t numClusters = tempV2->norm1();
        if (debugOutput_ && (comm->getRank() == 0))
          std::cout << "numClusters " << numClusters << std::endl;
        tempV->putScalar(ZERO);
        kernelApproximations_->blockA_->apply(*tempV2, *tempV);

        Scalar numClusterPairs = tempV->dot(*tempV2);
        if (debugOutput_ && (comm->getRank() == 0))
          std::cout << "numClusterPairs " << numClusterPairs << std::endl;

        bool doDrop;
        if (keepTransfers >= 0) {
          doDrop = (keepTransfers <= k);
        } else {
          doDrop = (droppedClusterPairs + numClusterPairs < (1.0 - coarseningRate) * totalNumClusterPairs);
        }
        if (doDrop) {
          auto lcl_transfer       = transferMatrices_[k]->blockA_->getLocalMatrixHost();
          auto lcl_transfer_graph = lcl_transfer.graph;
          for (LocalOrdinal j = 0; j < lcl_transfer_graph.entries.extent_int(0); j++)
            blidsToDrop.insert(lcl_transfer_graph.entries(j));

          droppedClusterPairs += numClusterPairs;
        } else {
          if (debugOutput_ && (comm->getRank() == 0))
            std::cout << "Dropped " << transferMatrices_.size() - 1 - k << " transfers of " << transferMatrices_.size() << " dropped cp: " << droppedClusterPairs << std::endl;
          break;
        }
      }
    }

    // number of cluster pairs dropped
    int dropped = 0;
    // number of cluster pairs we kept
    int kept = 0;
    // number of cluster pairs that were no longer present
    int ignored = 0;
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
        std::cout << "dropped " << gbl_dropped << " kept " << gbl_kept << " ignored " << gbl_ignored << std::endl;
    }
  }

  newKernelBlockGraph->fillComplete(kernelApproximations_->blockA_->getDomainMap(),
                                    kernelApproximations_->blockA_->getRangeMap());
  newKernelBlockGraph = removeSmallEntries(newKernelBlockGraph, Teuchos::ScalarTraits<MagnitudeType>::eps());
  diffKernelApprox->fillComplete(clusterCoeffMap_,
                                 clusterCoeffMap_);

  // coarse point matrix of cluster pairs
  Teuchos::RCP<matrix_type> newKernelApprox;
  {
    Teuchos::RCP<matrix_type> temp = MatrixMatrix::add(ONE, false, *kernelApproximations_->pointA_, -ONE, false, *diffKernelApprox);
    newKernelApprox                = removeSmallEntries(temp, Teuchos::ScalarTraits<MagnitudeType>::eps());
  }

  // construct identity on clusterCoeffMap_
  Teuchos::RCP<matrix_type> identity = buildIdentityMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(clusterCoeffMap_);

  Teuchos::RCP<blocked_matrix_type> newBlockedKernelApproximation = rcp(new blocked_matrix_type(newKernelApprox, newKernelBlockGraph, kernelApproximations_->blockMap_, kernelApproximations_->ghosted_blockMap_));

  // select subset of transfer matrices for coarse operator
  std::vector<Teuchos::RCP<blocked_matrix_type> > newTransferMatrices;
  {
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
    // transfer = newBasisMatrix * (identity + newTransferMatrices[K-1]^T) * ... * (identity + newTransferMatrices[0])^T
    Teuchos::RCP<matrix_type> transfer = rcp(new matrix_type(*newBasisMatrix));
    for (int i = Teuchos::as<int>(newTransferMatrices.size()) - 1; i >= 0; i--) {
      Teuchos::RCP<matrix_type> temp  = MatrixMatrix::add(ONE, false, *identity, ONE, false, *newTransferMatrices[i]->pointA_);
      Teuchos::RCP<matrix_type> temp2 = rcp(new matrix_type(newBasisMatrix->getRowMap(), 0));
      MatrixMatrix::Multiply(*transfer, false, *temp, true, *temp2);
      transfer = temp2;
    }

    // diffFarField = transfer * diffKernelApprox * transfer^T
    RCP<matrix_type> diffFarField;
    {
      Teuchos::RCP<matrix_type> temp = rcp(new matrix_type(newBasisMatrix->getRowMap(), 0));
      MatrixMatrix::Multiply(*transfer, false, *diffKernelApprox, false, *temp);
      diffFarField = rcp(new matrix_type(newBasisMatrix->getRowMap(), 0));
      MatrixMatrix::Multiply(*temp, false, *transfer, true, *diffFarField);
    }

    // newNearField = P^T * nearField * P + diffFarField
    {
      RCP<matrix_type> temp = rcp(new matrix_type(nearField_->getRowMap(), 0));
      MatrixMatrix::Multiply(*nearField_, false, *P, false, *temp);
      RCP<matrix_type> temp2 = rcp(new matrix_type(P->getDomainMap(), 0));
      MatrixMatrix::Multiply(*P, true, *temp, false, *temp2);
      newNearField = MatrixMatrix::add(ONE, false, *temp2, ONE, false, *diffFarField);
      newNearField = removeSmallEntries(newNearField, Teuchos::ScalarTraits<MagnitudeType>::eps());
    }
  }

  return Teuchos::rcp(new HierarchicalOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>(newNearField,
                                                                                          newBlockedKernelApproximation,
                                                                                          newBasisMatrix,
                                                                                          newTransferMatrices,
                                                                                          params_));
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
    // transfer = basisMatrix_ * (identity + transferMatrices_[K-1]) * ... * (identity + transferMatrices_[0])
    RCP<matrix_type> transfer = rcp(new matrix_type(*basisMatrix_));

    if (hasTransferMatrices()) {
      // construct identity on clusterCoeffMap_
      Teuchos::RCP<matrix_type> identity = buildIdentityMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(clusterCoeffMap_);

      for (int i = Teuchos::as<int>(transferMatrices_.size()) - 1; i >= 0; i--) {
        RCP<matrix_type> temp  = MatrixMatrix::add(ONE, false, *identity, ONE, false, *transferMatrices_[i]->pointA_);
        RCP<matrix_type> temp2 = rcp(new matrix_type(basisMatrix_->getRowMap(), 0));
        MatrixMatrix::Multiply(*transfer, false, *temp, true, *temp2);
        transfer = temp2;
      }
    }

    // farField = transfer * kernelApproximations_ * transfer^T
    RCP<matrix_type> temp = rcp(new matrix_type(basisMatrix_->getRowMap(), 0));
    MatrixMatrix::Multiply(*transfer, false, *kernelApproximations_->pointA_, false, *temp);
    RCP<matrix_type> farField = rcp(new matrix_type(basisMatrix_->getRowMap(), 0));
    MatrixMatrix::Multiply(*temp, false, *transfer, true, *farField);

    // nearField_ + farField
    return MatrixMatrix::add(ONE, false, *nearField_, ONE, false, *farField);

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
