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

    using crs_matrix = Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    using row_ptr_type = typename crs_matrix::local_graph_device_type::row_map_type::non_const_type;
    using col_idx_type = typename crs_matrix::local_graph_device_type::entries_type::non_const_type;
    using vals_type = typename crs_matrix::local_matrix_device_type::values_type;

    typedef Kokkos::ArithTraits<Scalar> ATS;
    using impl_SC = typename ATS::val_type;
    using impl_ATS = Kokkos::ArithTraits<impl_SC>;

    auto lclA = A->getLocalMatrixDevice();

    auto rowptr = row_ptr_type("rowptr", lclA.numRows()+1);

    Kokkos::parallel_for("removeSmallEntries::rowptr1",
                         Kokkos::RangePolicy<LocalOrdinal>(0, lclA.numRows()),
                         KOKKOS_LAMBDA(const LocalOrdinal rlid) {
                           auto row = lclA.row(rlid);
                           for (LocalOrdinal k = 0; k<row.length; ++k) {
                             if (impl_ATS::magnitude(row.value(k)) > tol) {
                               rowptr(rlid+1) += 1;
                             }
                           }
                         });
    LocalOrdinal nnz;
    Kokkos::parallel_scan("removeSmallEntries::rowptr2",
                          Kokkos::RangePolicy<LocalOrdinal>(0, lclA.numRows()),
                          KOKKOS_LAMBDA(const LocalOrdinal rlid, LocalOrdinal& partial_nnz, bool is_final) {

                            partial_nnz += rowptr(rlid+1);
                            if (is_final)
                              rowptr(rlid+1) = partial_nnz;

                          }, nnz);

    // auto nnz = rowptr(lclA.numRows());

    auto idx = col_idx_type("idx", nnz);
    auto vals = vals_type("vals", nnz);

    Kokkos::parallel_for("removeSmallEntries::indicesValues",
                         Kokkos::RangePolicy<LocalOrdinal>(0, lclA.numRows()),
                         KOKKOS_LAMBDA(const LocalOrdinal rlid) {
                           auto row = lclA.row(rlid);
                           auto I = rowptr(rlid);
                           for (LocalOrdinal k = 0; k<row.length; ++k) {
                             if (impl_ATS::magnitude(row.value(k)) > tol) {
                               idx(I) = row.colidx(k);
                               vals(I) = row.value(k);
                               I += 1;
                             }
                           }
                         });

    auto newA =  Teuchos::rcp(new crs_matrix(A->getRowMap(), A->getColMap(), rowptr, idx, vals));
    newA->fillComplete(A->getDomainMap(),
                       A->getRangeMap());
    return newA;
  }


  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class Node>
  HierarchicalOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  HierarchicalOperator(const Teuchos::RCP<matrix_type>& nearField,
                       const Teuchos::RCP<blocked_matrix_type>& kernelApproximations,
                       const Teuchos::RCP<matrix_type>& basisMatrix,
                       std::vector<Teuchos::RCP<blocked_matrix_type> >& transferMatrices)
      :
      nearField_(nearField),
      kernelApproximations_(kernelApproximations),
      basisMatrix_(basisMatrix),
      transferMatrices_(transferMatrices)
    {
      auto map = nearField_->getDomainMap();
      clusterCoeffMap_ = basisMatrix_->getDomainMap();

      const bool doDebugChecks = true;

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

        for (size_t i = 0; i<transferMatrices_.size(); i++) {
          // transfer matrices are entirely local, block diagonal on clusterCoeffMap
          TEUCHOS_ASSERT(clusterCoeffMap_->isSameAs(*transferMatrices_[i]->pointA_->getDomainMap()));
          TEUCHOS_ASSERT(clusterCoeffMap_->isSameAs(*transferMatrices_[i]->pointA_->getColMap()));
          TEUCHOS_ASSERT(clusterCoeffMap_->isSameAs(*transferMatrices_[i]->pointA_->getRowMap()));
          TEUCHOS_ASSERT(clusterCoeffMap_->isSameAs(*transferMatrices_[i]->pointA_->getRangeMap()));
        }
      }

      // Set the two importers to Isend
      Teuchos::RCP<Teuchos::ParameterList> distParams = rcp(new Teuchos::ParameterList());
      distParams->set("Send type", "Isend");
      {
        Teuchos::RCP<const Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> > nearFieldImporter = nearField_->getGraph()->getImporter();
        nearFieldImporter->getDistributor().setParameterList(distParams);
        auto revDistor = nearFieldImporter->getDistributor().getReverse(false);
        if (!revDistor.is_null())
          revDistor->setParameterList(distParams);
      }

      {
        Teuchos::RCP<const Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> > kernelApproximationsImporter = kernelApproximations_->pointA_->getGraph()->getImporter();
        kernelApproximationsImporter->getDistributor().setParameterList(distParams);
        auto revDistor = kernelApproximationsImporter->getDistributor().getReverse(false);
        if (!revDistor.is_null())
          revDistor->setParameterList(distParams);
      }

      // Allocate memory for apply with vectors
      allocateMemory(1);
    }


  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class Node>
  void
  HierarchicalOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  apply(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
        Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
        Teuchos::ETransp mode,
        Scalar alpha,
        Scalar beta) const {
    using Teuchos::RCP;
    const Scalar one = Teuchos::ScalarTraits<Scalar>::one();
    const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
    bool flip = true;

    allocateMemory(X.getNumVectors());

    // near field - part 1
    RCP<const Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> > nearFieldImporter = nearField_->getGraph()->getImporter();
    {
      if (mode == Teuchos::NO_TRANS) {
        X_colmap_->beginImport(X, *nearFieldImporter, INSERT);
      } else if (mode == Teuchos::TRANS) {
        nearField_->localApply(X, *X_colmap_, mode, alpha, zero);
        Y.scale (beta);
        Y.beginExport(*X_colmap_, *nearFieldImporter, ADD_ASSIGN);
      }
    }

    // upward pass
    {
      Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("upward pass")));

      basisMatrix_->apply(X, *coefficients_, Teuchos::TRANS);

      for (int i = Teuchos::as<int>(transferMatrices_.size())-1; i>=0; i--)
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

      RCP<const Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> > kernelApproximationsImporter = kernelApproximations_->pointA_->getGraph()->getImporter();
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

    // far field interactions - part 2
    {
      Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("far field 2")));

      RCP<const Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> > kernelApproximationsImporter = kernelApproximations_->pointA_->getGraph()->getImporter();
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

    // downward pass
    {
      Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("downward pass")));

      for (size_t i = 0; i<transferMatrices_.size(); i++)
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
  Teuchos::RCP<HierarchicalOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  HierarchicalOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  restrict(const Teuchos::RCP<matrix_type>& P) {
    // H_c = P^T * H * P
    using lo_vec_type = typename blocked_map_type::lo_vec_type;
    using vec_type = typename Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
    using MagnitudeType = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
    using Teuchos::RCP;
    using Teuchos::rcp;

    // P^T * basisMatrix
    RCP<matrix_type> newBasisMatrix = rcp(new matrix_type(P->getDomainMap(), clusterCoeffMap_, 0));
    MatrixMatrix::Multiply(*P, true, *basisMatrix_, false, *newBasisMatrix);

    //
    auto clusterMap = kernelApproximations_->blockA_->getRowMap();
    auto clusterSizes = kernelApproximations_->blockMap_->blockSizes_;
    auto ghosted_clusterMap = kernelApproximations_->blockA_->getColMap();
    auto ghosted_clusterSizes = kernelApproximations_->ghosted_blockMap_->blockSizes_;
    auto lcl_clusterSizes = clusterSizes->getLocalViewHost(Tpetra::Access::ReadOnly);
    auto lcl_ghosted_clusterSizes = ghosted_clusterSizes->getLocalViewHost(Tpetra::Access::ReadOnly);

    // get number of unknowns associated with cluster via new basisMatrix
    RCP<vec_type> numUnknownsPerCluster = rcp(new vec_type(clusterMap, false));
    {
      auto lcl_numUnknownsPerCluster = numUnknownsPerCluster->getLocalViewHost(Tpetra::Access::OverwriteAll);
      // Compute the transpose of the newBasisMatrix.
      RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > newBasisMatrixT;
      Tpetra::RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node> transposer(newBasisMatrix);
      RCP<Teuchos::ParameterList> transposeParams = rcp(new Teuchos::ParameterList);
      newBasisMatrixT = transposer.createTranspose(transposeParams);

      auto rowptr = newBasisMatrixT->getLocalRowPtrsHost();
      LocalOrdinal clusterStart = 0;
      LocalOrdinal clusterEnd = 0;
      for (LocalOrdinal cluster = 0; cluster < lcl_clusterSizes.extent_int(0); ++cluster) {
        clusterStart = clusterEnd;
        clusterEnd += lcl_clusterSizes(cluster, 0);
        LocalOrdinal maxEntries = 0;
        for (LocalOrdinal row = clusterStart; row < clusterEnd; ++row) {
          LocalOrdinal numEntriesPerRow = rowptr(row+1)-rowptr(row);
          maxEntries = std::max(maxEntries, numEntriesPerRow);
        }
        lcl_numUnknownsPerCluster(cluster, 0) = maxEntries;
      }
      TEUCHOS_ASSERT_EQUALITY(clusterEnd+1, rowptr.extent_int(0));
    }

    const Scalar one = Teuchos::ScalarTraits<Scalar>::one();
    for (int i = Teuchos::as<int>(transferMatrices_.size())-1; i>=0; i--)
      transferMatrices_[i]->blockA_->apply(*numUnknownsPerCluster, *numUnknownsPerCluster, Teuchos::NO_TRANS, one, one);

    // get ghosted numUnknownsPerCluster
    RCP<vec_type> ghosted_numUnknownsPerCluster = rcp(new vec_type(ghosted_clusterMap, false));
    auto import = kernelApproximations_->blockA_->getCrsGraph()->getImporter();
    ghosted_numUnknownsPerCluster->doImport(*numUnknownsPerCluster, *import, Tpetra::INSERT);

    // coarse cluster pair graph
    RCP<matrix_type> newKernelBlockGraph = rcp(new matrix_type(kernelApproximations_->blockA_->getCrsGraph()));
    newKernelBlockGraph->resumeFill();
    // point entries of cluster pairs that should be moved to the near field
    RCP<matrix_type> diffKernelApprox = rcp(new matrix_type(kernelApproximations_->pointA_->getCrsGraph()));

    {
      auto lcl_BlockGraph = kernelApproximations_->blockA_->getLocalMatrixHost();
      auto lcl_newBlockGraph = newKernelBlockGraph->getLocalMatrixHost();
      auto lcl_KernelApprox = kernelApproximations_->pointA_->getLocalMatrixHost();
      auto lcl_diffKernelApprox = diffKernelApprox->getLocalMatrixHost();
      auto lcl_numUnknownsPerCluster = numUnknownsPerCluster->getLocalViewHost(Tpetra::Access::ReadOnly);
      auto lcl_ghosted_numUnknownsPerCluster = ghosted_numUnknownsPerCluster->getLocalViewHost(Tpetra::Access::ReadOnly);
      auto lcl_offsets = Kokkos::create_mirror_view(kernelApproximations_->blockMap_->offsets_);
      auto lcl_ghosted_offsets = Kokkos::create_mirror_view(kernelApproximations_->ghosted_blockMap_->offsets_);
      Kokkos::deep_copy(lcl_offsets, kernelApproximations_->blockMap_->offsets_);
      Kokkos::deep_copy(lcl_ghosted_offsets, kernelApproximations_->ghosted_blockMap_->offsets_);

      int dropped = 0;
      int kept = 0;
      int ignored = 0;
      // loop over cluster pairs
      for (LocalOrdinal brlid = 0; brlid < lcl_BlockGraph.numRows(); ++brlid) {
        size_t brsize = lcl_clusterSizes(brlid, 0);
        auto brow = lcl_BlockGraph.row(brlid);
        auto new_brow = lcl_newBlockGraph.row(brlid);
        for (LocalOrdinal k = 0; k < brow.length; ++k) {
          if (brow.value(k) > 0.5) {
            LocalOrdinal bclid = brow.colidx(k);
            size_t bcsize = lcl_ghosted_clusterSizes(bclid, 0);

            if (brsize * bcsize >= lcl_numUnknownsPerCluster(brlid, 0) * lcl_ghosted_numUnknownsPerCluster(bclid, 0)) {
              ++dropped;
              new_brow.value(k) = 0;

              const LocalOrdinal row_start = lcl_offsets(brlid);
              const LocalOrdinal row_end = lcl_offsets(brlid+1);
              const LocalOrdinal col_start = lcl_ghosted_offsets(bclid);
              const LocalOrdinal col_end = lcl_ghosted_offsets(bclid+1);
              TEUCHOS_ASSERT_EQUALITY(Teuchos::as<size_t>(row_end-row_start), brsize);
              TEUCHOS_ASSERT_EQUALITY(Teuchos::as<size_t>(col_end-col_start), bcsize);
              // loop over rows of kernelApproximations in pointwise indexing
              for (LocalOrdinal rlid = row_start; rlid <  row_end; ++rlid) {
                auto diff_row = lcl_diffKernelApprox.row(rlid);
                auto row = lcl_KernelApprox.row(rlid);
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
                  oss << "col_start "<< col_start << " col_end " << col_end << std::endl;
                  for (LocalOrdinal n = 0; n < row.length; ++n) {
                    oss << row.colidx(n) << " " << row.value(n) << std::endl;
                  }
                  std::cout << oss.str();
                }
                TEUCHOS_ASSERT_EQUALITY(removed, bcsize);
              }

            } else {
              ++kept;
              new_brow.value(k) = brow.value(k);
            }
          } else {
            ++ignored;
            new_brow.value(k) = brow.value(k);
          }
        }
      }

      // std::cout << "dropped " << dropped << " kept " << kept << " ignored " << ignored << std::endl;
    }
    newKernelBlockGraph->fillComplete(kernelApproximations_->blockA_->getDomainMap(),
                                      kernelApproximations_->blockA_->getRangeMap());
    newKernelBlockGraph = removeSmallEntries(newKernelBlockGraph, Teuchos::ScalarTraits<MagnitudeType>::eps());
    diffKernelApprox->fillComplete(clusterCoeffMap_,
                                   clusterCoeffMap_);

    std::vector<Teuchos::RCP<blocked_matrix_type> > newTransferMatrices;
    {
      const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
      Teuchos::RCP<mv_type> mv_temp = rcp(new mv_type(newKernelBlockGraph->getDomainMap(), 1));
      mv_temp->putScalar(one);
      RCP<mv_type> clusterUseCount = rcp(new mv_type(newKernelBlockGraph->getDomainMap(), 1));
      clusterUseCount->putScalar(zero);
      newKernelBlockGraph->apply(*mv_temp, *clusterUseCount, Teuchos::NO_TRANS);
      newKernelBlockGraph->apply(*mv_temp, *clusterUseCount, Teuchos::TRANS, one, one);

      std::vector<int> keepTransfers;
      auto comm = getComm();
      // std::ostringstream oss;
      for (int i = Teuchos::as<int>(transferMatrices_.size())-1; i>=0; i--) {
        transferMatrices_[i]->blockA_->localApply(*clusterUseCount, *mv_temp, Teuchos::TRANS);
        auto lcl_counts = mv_temp->getLocalViewHost(Tpetra::Access::ReadOnly);
        Scalar lcl_count = zero, gbl_count = zero;
        for (LocalOrdinal n = 0; n < lcl_counts.extent_int(0); ++n) {
          lcl_count += lcl_counts(n, 0);
        }
        Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &lcl_count, &gbl_count);
        // oss << "Transfer " << i << " count " << gbl_count << std::endl;

        if (gbl_count < 0.5) {
          RCP<matrix_type> temp = rcp(new matrix_type(P->getDomainMap(), clusterCoeffMap_, 0));
          MatrixMatrix::Multiply(*newBasisMatrix, false, *transferMatrices_[i]->pointA_, true, *temp);
          newBasisMatrix = temp;
        } else
          keepTransfers.push_back(i);
      }
      // std::cout << oss.str();

      for (auto it = keepTransfers.begin(); it != keepTransfers.end(); ++it) {
        newTransferMatrices.insert(newTransferMatrices.begin(), transferMatrices_[*it]);
      }
    }

    // coarse point matrix of cluster pairs
    Teuchos::RCP<matrix_type> newKernelApprox;
    {
      Teuchos::RCP<matrix_type> temp = MatrixMatrix::add(one, false, *kernelApproximations_->pointA_, -one, false, *diffKernelApprox);
      newKernelApprox = removeSmallEntries(temp, Teuchos::ScalarTraits<MagnitudeType>::eps());
    }

    // static int lvlNo = 1;
    // Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >::writeSparseFile("kernel"+std::to_string(lvlNo), *newKernelApprox);
    // Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >::writeSparseFile("diffKernel"+std::to_string(lvlNo), *diffKernelApprox);
    // Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >::writeSparseFile("kernelGraph"+std::to_string(lvlNo), *newKernelBlockGraph);
    // Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >::writeDenseFile("numUnknownsPerCluster"+std::to_string(lvlNo), numUnknownsPerCluster);
    // ++lvlNo;

    Teuchos::RCP<blocked_matrix_type> newBlockedKernelApproximation = rcp(new blocked_matrix_type(newKernelApprox, newKernelBlockGraph, kernelApproximations_->blockMap_, kernelApproximations_->ghosted_blockMap_));

    // construct identity on clusterCoeffMap_
    Teuchos::RCP<matrix_type> identity = rcp(new matrix_type(clusterCoeffMap_, 1));
    Teuchos::ArrayView<const GlobalOrdinal> gblRows = clusterCoeffMap_->getLocalElementList ();
    for (auto it = gblRows.begin (); it != gblRows.end (); ++it) {
      Teuchos::Array<GlobalOrdinal> col (1, *it);
      Teuchos::Array<Scalar> val (1, one);
      identity->insertGlobalValues (*it, col (), val ());
    }
    identity->fillComplete ();

    // transfer = basisMatrix_ * (identity + transferMatrices_[0]) * ... * (identity + transferMatrices_[n-1])
    Teuchos::RCP<matrix_type> transfer = rcp(new matrix_type(*newBasisMatrix));
    for (size_t i = 0; i<transferMatrices_.size(); i++) {
      Teuchos::RCP<matrix_type> temp = MatrixMatrix::add(one, false, *identity, one, false, *transferMatrices_[i]->pointA_);
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

    // P^T * nearField * P
    RCP<matrix_type> newNearField;
    {
      RCP<matrix_type> temp = rcp(new matrix_type(nearField_->getRowMap(), 0));
      MatrixMatrix::Multiply(*nearField_, false, *P, false, *temp);
      RCP<matrix_type> temp2 = rcp(new matrix_type(P->getDomainMap(), 0));
      MatrixMatrix::Multiply(*P, true, *temp, false, *temp2);
      newNearField = MatrixMatrix::add(one, false, *temp2, one, false, *diffFarField);
    }

    return Teuchos::rcp(new HierarchicalOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node>(newNearField, newBlockedKernelApproximation, newBasisMatrix, newTransferMatrices));
  }


  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class Node>
  Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  HierarchicalOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  toMatrix() {
    using Teuchos::RCP;
    using Teuchos::rcp;

    const Scalar one = Teuchos::ScalarTraits<Scalar>::one();

    // construct identity on clusterCoeffMap_
    RCP<matrix_type> identity = rcp(new matrix_type(clusterCoeffMap_, 1));
    Teuchos::ArrayView<const GlobalOrdinal> gblRows = clusterCoeffMap_->getLocalElementList ();
    for (auto it = gblRows.begin (); it != gblRows.end (); ++it) {
      Teuchos::Array<GlobalOrdinal> col (1, *it);
      Teuchos::Array<Scalar> val (1, one);
      identity->insertGlobalValues (*it, col (), val ());
    }
    identity->fillComplete ();

    // transfer = basisMatrix_ * (identity + transferMatrices_[0]) * ... * (identity + transferMatrices_[n-1])
    RCP<matrix_type> transfer = rcp(new matrix_type(*basisMatrix_));
    for (size_t i = 0; i<transferMatrices_.size(); i++) {
      RCP<matrix_type> temp = MatrixMatrix::add(one, false, *identity, one, false, *transferMatrices_[i]->pointA_);
      RCP<matrix_type> temp2 = rcp(new matrix_type(basisMatrix_->getRowMap(), 0));
      MatrixMatrix::Multiply(*transfer, false, *temp, true, *temp2);
      transfer = temp2;
    }

    // farField = transfer * kernelApproximations_ * transfer^T
    RCP<matrix_type> temp = rcp(new matrix_type(basisMatrix_->getRowMap(), 0));
    MatrixMatrix::Multiply(*transfer, false, *kernelApproximations_->pointA_, false, *temp);
    RCP<matrix_type> farField = rcp(new matrix_type(basisMatrix_->getRowMap(), 0));
    MatrixMatrix::Multiply(*temp, false, *transfer, true, *farField);

    // nearField_ + farField
    return MatrixMatrix::add(one, false, *nearField_, one, false, *farField);

  }


  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class Node>
  void
  HierarchicalOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  allocateMemory(size_t numVectors) const {
    if (coefficients_.is_null() || coefficients_->getNumVectors() != numVectors) {
      coefficients_  = Teuchos::rcp(new mv_type(clusterCoeffMap_, numVectors));
      coefficients2_ = Teuchos::rcp(new mv_type(clusterCoeffMap_, numVectors));
      X_colmap_ = Teuchos::rcp(new mv_type(nearField_->getColMap(), numVectors));
      coefficients_colmap_  = Teuchos::rcp(new mv_type(kernelApproximations_->pointA_->getColMap(), numVectors));
    }
  }
}

#endif // TPETRA_HIERARCHICALOPERATOR_DEF_HPP
