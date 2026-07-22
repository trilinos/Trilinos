// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifndef AMESOS2_TPETRACRSMATRIX_MATRIXADAPTER_DEF_HPP
#define AMESOS2_TPETRACRSMATRIX_MATRIXADAPTER_DEF_HPP

#include "Amesos2_TpetraCrsMatrix_MatrixAdapter_decl.hpp"
#include "Amesos2_TpetraRowMatrix_AbstractMatrixAdapter_def.hpp"
#include "Amesos2_MatrixAdapter_def.hpp"

namespace Amesos2 {

  template <typename Scalar,
            typename LocalOrdinal,
            typename GlobalOrdinal,
            typename Node>
  ConcreteMatrixAdapter<
    Tpetra::CrsMatrix<Scalar,
                      LocalOrdinal,
                      GlobalOrdinal,
                      Node>
    >::ConcreteMatrixAdapter(Teuchos::RCP<matrix_t> m)
      : AbstractConcreteMatrixAdapter<Tpetra::RowMatrix<Scalar,
                                                        LocalOrdinal,
                                                        GlobalOrdinal,
                                                        Node>,
                                      Tpetra::CrsMatrix<Scalar,
                                                        LocalOrdinal,
                                                        GlobalOrdinal,
                                                        Node> >(m) // with implicit cast
    {}

  template <typename Scalar,
            typename LocalOrdinal,
            typename GlobalOrdinal,
            typename Node>
  Teuchos::RCP<const MatrixAdapter<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > >
  ConcreteMatrixAdapter<
    Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>
    >::get_impl(const Teuchos::Ptr<const map_t> map, EDistribution distribution) const
    {
      using Teuchos::RCP;
      using Teuchos::rcp;
      using Teuchos::rcpFromPtr;
      typedef Tpetra::Import<local_ordinal_t, global_ordinal_t, node_t> import_t;

      RCP<import_t> importer =
        rcp (new import_t (this->getRowMap(), rcpFromPtr (map)));

      RCP<matrix_t> t_mat;

      t_mat = Tpetra::importAndFillCompleteCrsMatrix<matrix_t>( (this->mat_), *importer ); // DomainMap, RangeMap, params inputs default to Teuchos::null

      // Case for non-contiguous GIDs
      if ( distribution == CONTIGUOUS_AND_ROOTED ) {

        auto myRank = map->getComm()->getRank();

        auto local_matrix = t_mat->getLocalMatrixDevice();
        const size_t global_num_contiguous_entries = t_mat->getGlobalNumRows();
        const size_t local_num_contiguous_entries = (myRank == 0) ? t_mat->getGlobalNumRows() : 0;

        //create maps
        //typedef Tpetra::Map< local_ordinal_t, global_ordinal_t, node_t>  contiguous_map_type;
        RCP<const map_t> contiguousRowMap = rcp( new map_t(global_num_contiguous_entries, local_num_contiguous_entries, 0, (t_mat->getComm() ) ) );
        RCP<const map_t> contiguousColMap = rcp( new map_t(global_num_contiguous_entries, local_num_contiguous_entries, 0, (t_mat->getComm() ) ) );
        RCP<const map_t> contiguousDomainMap = rcp( new map_t(global_num_contiguous_entries, local_num_contiguous_entries, 0, (t_mat->getComm() ) ) );
        RCP<const map_t> contiguousRangeMap  = rcp( new map_t(global_num_contiguous_entries, local_num_contiguous_entries, 0, (t_mat->getComm() ) ) );

        RCP<matrix_t> contiguous_t_mat = rcp( new matrix_t(contiguousRowMap, contiguousColMap, local_matrix) );
        contiguous_t_mat->resumeFill();
        contiguous_t_mat->expertStaticFillComplete(contiguousDomainMap, contiguousRangeMap);

        return rcp (new ConcreteMatrixAdapter<matrix_t> (contiguous_t_mat));
      } //end if not contiguous

      return rcp (new ConcreteMatrixAdapter<matrix_t> (t_mat));
    }



  template <typename Scalar,
            typename LocalOrdinal,
            typename GlobalOrdinal,
            typename Node>
  Teuchos::RCP<const MatrixAdapter<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > >
  ConcreteMatrixAdapter<
    Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>
    >::reindex_impl(Teuchos::RCP<const map_t> &contigRowMap,
                    Teuchos::RCP<const map_t> &contigColMap,
                    const EPhase current_phase) const
    {
      typedef Tpetra::Map< local_ordinal_t, global_ordinal_t, node_t> contiguous_map_type;
      using Teuchos::RCP;
      using Teuchos::rcp;
#ifdef HAVE_AMESOS2_TIMERS
      auto reindexTimer = Teuchos::TimeMonitor::getNewTimer("Time to re-index matrix gids");
      Teuchos::TimeMonitor ReindexTimer(*reindexTimer);
#endif

      auto rowMap = this->mat_->getRowMap();
      auto colMap = this->mat_->getColMap();
      auto rowComm = rowMap->getComm();
      auto colComm = colMap->getComm();

      GlobalOrdinal indexBase = rowMap->getIndexBase();
      GlobalOrdinal numDoFs = this->mat_->getGlobalNumRows();
      LocalOrdinal nRows = this->mat_->getLocalNumRows();
      LocalOrdinal nCols = colMap->getLocalNumElements();

      RCP<matrix_t> contiguous_t_mat;
      // check when to recompute contigRowMap & contigColMap
      if(current_phase == PREORDERING || current_phase == SYMBFACT) {
        auto tmpMap = rcp (new contiguous_map_type (numDoFs, nRows, indexBase, rowComm));
        global_ordinal_t frow = tmpMap->getMinGlobalIndex();

        // Create new GID list for RowMap
        Kokkos::View<global_ordinal_t*, HostExecSpaceType> rowIndexList ("rowIndexList", nRows);
        for (local_ordinal_t k = 0; k < nRows; k++) {
          rowIndexList(k) = frow+k; // based on index-base of rowMap
        }
        // Create new GID list for ColMap
        Kokkos::View<global_ordinal_t*, HostExecSpaceType> colIndexList ("colIndexList", nCols);
        // initialize to catch col GIDs that are not in row GIDs
        // they will be all assigned to (n+1)th columns
        for (local_ordinal_t k = 0; k < nCols; k++) {
          colIndexList(k) = numDoFs+indexBase;
        }
        typedef Tpetra::MultiVector<global_ordinal_t,
                                    local_ordinal_t,
                                    global_ordinal_t,
                                    node_t> gid_mv_t;
        Teuchos::ArrayView<const global_ordinal_t> rowIndexArray(rowIndexList.data(), nRows);
        Teuchos::ArrayView<const global_ordinal_t> colIndexArray(colIndexList.data(), nCols);
        gid_mv_t row_mv (rowMap, rowIndexArray, nRows, 1);
        gid_mv_t col_mv (colMap, colIndexArray, nCols, 1);
        typedef Tpetra::Import<local_ordinal_t, global_ordinal_t, node_t> import_t;
        RCP<import_t> importer = rcp (new import_t (rowMap, colMap));
        col_mv.doImport (row_mv, *importer, Tpetra::INSERT);
        {
          // col_mv is imported from rowIndexList, which is based on index-base of rowMap
          auto col_view = col_mv.getLocalViewHost(Tpetra::Access::ReadOnly);
          for(int i=0; i<nCols; i++) colIndexList(i) = col_view(i,0);
        }
        // Create new Row & Col Maps (both based on indexBase of rowMap)
        contigRowMap = rcp (new contiguous_map_type (numDoFs, rowIndexList.data(), nRows, indexBase, rowComm));
        contigColMap = rcp (new contiguous_map_type (numDoFs, colIndexList.data(), nCols, indexBase, colComm));

        // Create contiguous Matrix
        auto lclMatrix = this->mat_->getLocalMatrixDevice();
        contiguous_t_mat = rcp( new matrix_t(contigRowMap, contigColMap, lclMatrix));
      } else {
        // Build Matrix with contiguous Maps
        auto lclMatrix = this->mat_->getLocalMatrixDevice();
        auto importer  = this->mat_->getCrsGraph()->getImporter();
        auto exporter  = this->mat_->getCrsGraph()->getExporter();
        contiguous_t_mat = rcp( new matrix_t(lclMatrix, contigRowMap, contigColMap, contigRowMap, contigColMap, importer,exporter));
      }

      // Return new matrix adapter
      return rcp (new ConcreteMatrixAdapter<matrix_t> (contiguous_t_mat));
    }


  template <typename Scalar,
            typename LocalOrdinal,
            typename GlobalOrdinal,
            typename Node>
  template<typename KV_S, typename KV_GO, typename KV_GS, typename host_ordinal_type_array, typename host_scalar_type_array>
  LocalOrdinal
  ConcreteMatrixAdapter<
    Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>
    >::gather_impl(KV_S& nzvals, KV_GO& indices, KV_GS& pointers,
                   host_ordinal_type_array &perm_g2l,
                   host_ordinal_type_array &recvCountRows, host_ordinal_type_array &recvDisplRows,
                   host_ordinal_type_array &recvCounts, host_ordinal_type_array &recvDispls,
                   host_ordinal_type_array &transpose_map, host_scalar_type_array &nzvals_t,
                   bool column_major, EPhase current_phase) const
    {
      LocalOrdinal ret = -1;
      {
#ifdef HAVE_AMESOS2_TIMERS
        Teuchos::RCP< Teuchos::Time > gatherTime = Teuchos::TimeMonitor::getNewCounter ("Amesos2::gather");
        Teuchos::TimeMonitor GatherTimer(*gatherTime);
#endif
        auto rowMap = this->mat_->getRowMap();
        auto colMap = this->mat_->getColMap();
        auto comm = rowMap->getComm();
        auto nRanks = comm->getSize();
        auto myRank = comm->getRank();

        global_ordinal_t nRows = this->mat_->getGlobalNumRows();
        auto lclMatrix = this->mat_->getLocalMatrixDevice();

        // check when to recompute communication patterns
        if(current_phase == PREORDERING || current_phase == SYMBFACT) {
          // grab rowptr and colind on host
          auto lclRowptr_d = lclMatrix.graph.row_map;
          auto lclColind_d = lclMatrix.graph.entries;
          auto lclRowptr = Kokkos::create_mirror_view(lclRowptr_d);
          auto lclColind = Kokkos::create_mirror_view(lclColind_d);
          Kokkos::deep_copy(lclRowptr, lclRowptr_d);
          Kokkos::deep_copy(lclColind, lclColind_d);

          // index-bases
          global_ordinal_t rowIndexBase = rowMap->getIndexBase();
          global_ordinal_t colIndexBase = colMap->getIndexBase();
          // map from global to local
          host_ordinal_type_array  perm_l2g; // map from GIDs
          // true uses 'original' contiguous row inds 0:(n-1) (no need to perm sol or rhs), 
          // false uses GIDs given by map (need to swap sol & rhs, but use matrix perm, e.g., fill-reducing)
          bool swap_cols = false; //true;
          // workspace to transpose
          KV_GS pointers_t;
          KV_GO indices_t;
          // communication counts & displacements
          LocalOrdinal myNRows = this->mat_->getLocalNumRows();
          Kokkos::resize(recvCounts, nRanks);
          Kokkos::resize(recvDispls, nRanks+1);
          bool need_to_perm = false;
          {
#ifdef HAVE_AMESOS2_TIMERS
            Teuchos::RCP< Teuchos::Time > gatherTime_ = Teuchos::TimeMonitor::getNewCounter ("Amesos2::gather(rowptr)");
            Teuchos::TimeMonitor GatherTimer_(*gatherTime_);
#endif
            Teuchos::gather<int, LocalOrdinal> (&myNRows, 1, recvCounts.data(), 1, 0, *comm);
            if (myRank == 0) {
              Kokkos::resize(recvDispls, nRanks+1);
              recvDispls(0) = 0;
              for (int p = 1; p <= nRanks; p++) {
                recvDispls(p) = recvDispls(p-1) + recvCounts(p-1);
              }
              if (recvDispls(nRanks) != nRows) {
                TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Amesos2_TpetraCrsMatrix_MatrixAdapter::gather_impl : mismatch between gathered(local nrows) and global nrows.");
              }
            } else {
              for (int p = 0; p < nRanks; p++) {
                recvCounts(p) = 0;
                recvDispls(p) = 0;
              }
              recvDispls(nRanks) = 0;
            }
            // gether g2l perm (convert to 0-base)
            {
              host_ordinal_type_array lclMap;
              Kokkos::resize(lclMap, myNRows);
              if (myRank == 0) {
                Kokkos::resize(perm_g2l, nRows);
                Kokkos::resize(perm_l2g, nRows);
              } else {
                Kokkos::resize(perm_g2l, 1);
              }
              for (int i=0; i < myNRows; i++) {
                lclMap(i) = rowMap->getGlobalElement(i);
              }
              Teuchos::gatherv<int, LocalOrdinal> (lclMap.data(), myNRows, perm_g2l.data(),
                                                   recvCounts.data(), recvDispls.data(),
                                                   0, *comm);
              if (myRank == 0) {
                for (int i=0; i < nRows; i++) {
                  perm_g2l(i) -= rowIndexBase; // map to GIDs
                  perm_l2g(perm_g2l(i)) = i;   // map from GIDs
                  if (i != perm_g2l(i)) need_to_perm = true;
                }
              }
            }
            // gether rowptr
            // - making sure same ordinal type
            KV_GS lclRowptr_ ("localRowptr_", lclRowptr.extent(0));
            for (int i = 0; i < int(lclRowptr.extent(0)); i++) lclRowptr_(i) = lclRowptr(i);
            if (myRank == 0 && (column_major || need_to_perm)) {
              Kokkos::resize(pointers_t, nRows+1);
            } else if (myRank != 0) {
              Kokkos::resize(pointers_t, 2);
            }
            LocalOrdinal sendIdx = (myNRows > 0 ? 1 : 0); // To skip sending the first rowptr entry (note: 0, if local matrix is empty)
            LocalOrdinal *pointers_ = ((myRank != 0 || (column_major || need_to_perm)) ? pointers_t.data() : pointers.data());
            Teuchos::gatherv<int, LocalOrdinal> (&lclRowptr_(sendIdx), myNRows, &pointers_[1],
                                                 recvCounts.data(), recvDispls.data(),
                                                 0, *comm);

            // save row counts/displs
            Kokkos::resize(recvCountRows, nRanks);
            Kokkos::resize(recvDisplRows, nRanks+1);
            Kokkos::deep_copy(recvCountRows, recvCounts);
            Kokkos::deep_copy(recvDisplRows, recvDispls);
            if (myRank == 0) {
              // shift to global pointers
              pointers_[0] = 0;
              recvCounts(0) = pointers_[recvDispls(1)];
              LocalOrdinal displs = recvCounts(0);
              for (int p = 1; p < nRanks; p++) {
                // skip "Empty" submatrix (no rows)
                //  recvCounts(p) is zero, while pointers_[recvDispls(p+1)] now contains nnz from p-1
                if (recvDispls(p+1) > recvDispls(p)) {
                  // save recvCounts from pth MPI
                  recvCounts(p) = pointers_[recvDispls(p+1)];
                  // shift pointers for pth MPI to global
                  for (int i = 1+recvDispls(p); i <= recvDispls(p+1); i++) {
                    pointers_[i] += displs;
                  }
                  displs += recvCounts(p);
                }
              }
              ret = pointers_[nRows];
            }
          }
          {
#ifdef HAVE_AMESOS2_TIMERS
            Teuchos::RCP< Teuchos::Time > gatherTime_ = Teuchos::TimeMonitor::getNewCounter ("Amesos2::gather(colind)");
            Teuchos::TimeMonitor GatherTimer_(*gatherTime_);
#endif
            // gather colinds
            if (myRank == 0) {
              recvDispls(0) = 0;
              for (int p = 0; p < nRanks; p++) {
                recvDispls(p+1) = recvDispls(p) + recvCounts(p);
              }
            }
            // -- convert to global colids & ** convert to base-zero **
            KV_GO lclColind_ ("localColind_", lclColind.extent(0));
            for (int i = 0; i < int(lclColind.extent(0)); i++) {
              lclColind_(i) = (colMap->getGlobalElement((lclColind(i))) - colIndexBase);
            }
            if (column_major || need_to_perm) {
              Kokkos::resize(indices_t, indices.extent(0));
              Teuchos::gatherv<int, LocalOrdinal> (lclColind_.data(), lclColind_.extent(0), indices_t.data(),
                                                   recvCounts.data(), recvDispls.data(),
                                                   0, *comm);
            } else {
              Teuchos::gatherv<int, LocalOrdinal> (lclColind_.data(), lclColind_.extent(0), indices.data(),
                                                   recvCounts.data(), recvDispls.data(),
                                                   0, *comm);
            }
          }
          if (myRank == 0) {
#ifdef HAVE_AMESOS2_TIMERS
            Teuchos::RCP< Teuchos::Time > gatherTime_ = Teuchos::TimeMonitor::getNewCounter ("Amesos2::gather(transpose index)");
            Teuchos::TimeMonitor GatherTimer_(*gatherTime_);
#endif
            if (swap_cols && need_to_perm) {
              // convert col GIDs to 0:(n-1)
              for (int i=0; i<ret; i++) {
                if (column_major || need_to_perm) {
                  indices_t(i) = perm_l2g(indices_t(i));
                } else {
                  indices(i) = perm_l2g(indices(i));
                }
              }
            }
            // (note: column idexes are now in base-0)
            if (column_major) {
              // Map to transpose
              Kokkos::resize(transpose_map, ret);
              // Transopose to convert to CSC
              for (int i=0; i<=nRows; i++) {
                pointers(i) = 0;
              }
              for (int k=0; k<ret; k++) {
                int col = indices_t(k);
                if (col < nRows-1) {
                  pointers(col+2) ++;
                }
              }
              for (int i=1; i < nRows; i++) {
                pointers(i+1) += pointers(i);
              }
              // nonzeroes in each column are sorted in ascending row
              for (int row=0; row<nRows; row++) {
                // if swap cols, then just use the original 0:(n-1), otherwise GIDs
                int i = (swap_cols ? row : perm_l2g(row));
                for (int k=pointers_t(i); k<pointers_t(i+1); k++) {
                  int col = indices_t(k);
                  if (col < nRows) {
                    transpose_map(k) = pointers(1+col);
                    indices(pointers(1+col)) = row;
                    pointers(1+col) ++;
                  } else {
                    // extra columns
                    transpose_map(k) = -1;
                  }
                }
              }
            } else if (need_to_perm) {
              // Map to permute
              Kokkos::resize(transpose_map, ret);
              for (int i=0; i<nRows; i++) {
                int row = perm_g2l(i);
                pointers(row+2) = pointers_t(i+1)-pointers_t(i);
              }
              for (int i=1; i < nRows; i++) {
                pointers(i+1) += pointers(i);
              }
              for (int i=0; i<nRows; i++) {
                int row = perm_g2l(i);
                for (int k=pointers_t(i); k<pointers_t(i+1); k++) {
                  int col = indices_t(k);
                  if (col < nRows) {
                    transpose_map(k) = pointers(1+row);
                    indices(pointers(1+row)) = col;
                    pointers(1+row) ++;
                  } else {
                    transpose_map(k) = -1;
                  }
                }
              }
            } else {
              Kokkos::resize(transpose_map, 0);
            }
          }
          if (!need_to_perm || swap_cols) {
            // no need to perm rhs or sol
            Kokkos::resize(perm_g2l, 0);
          }
        }
        //if(current_phase == NUMFACT) // Numerical values may be used in symbolic (e.g, MWM)
        {
          {
#ifdef HAVE_AMESOS2_TIMERS
            Teuchos::RCP< Teuchos::Time > gatherTime_ = Teuchos::TimeMonitor::getNewCounter ("Amesos2::gather(nzvals)");
            Teuchos::TimeMonitor GatherTimer_(*gatherTime_);
#endif
            // grab numerical values on host
            auto lclNzvals_d = lclMatrix.values;
            auto lclNzvals = Kokkos::create_mirror_view(lclNzvals_d);;
            Kokkos::deep_copy(lclNzvals, lclNzvals_d);

            // gather nzvals
            if (transpose_map.extent(0) > 0) {
              Kokkos::resize(nzvals_t, nzvals.extent(0));
              Teuchos::gatherv<int, Scalar> (reinterpret_cast<Scalar*> (lclNzvals.data()), lclNzvals.extent(0),
                                             reinterpret_cast<Scalar*> (nzvals_t.data()), recvCounts.data(), recvDispls.data(),
                                             0, *comm);
            } else {
              Teuchos::gatherv<int, Scalar> (reinterpret_cast<Scalar*> (lclNzvals.data()), lclNzvals.extent(0),
                                             reinterpret_cast<Scalar*> (nzvals.data()), recvCounts.data(), recvDispls.data(),
                                             0, *comm);
            }
          }
          if (myRank == 0) {
            // Insert Numerical values to transopose matrix
            ret = pointers(nRows);
#ifdef HAVE_AMESOS2_TIMERS
            Teuchos::RCP< Teuchos::Time > gatherTime_ = Teuchos::TimeMonitor::getNewCounter ("Amesos2::gather(transpose values)");
            Teuchos::TimeMonitor GatherTimer_(*gatherTime_);
#endif
            if (transpose_map.extent(0) > 0) {
              //for (int k=0; k<ret; k++) {
              //  if (transpose_map(k) >= 0) {
              //    nzvals(transpose_map(k)) = nzvals_t(k);
              //  }
              //}
              Kokkos::parallel_for("Amesos2::TpetraCrsMatrixAdapter::gather", Kokkos::RangePolicy<HostExecSpaceType>(0, ret),
                KOKKOS_LAMBDA(const int k) { if (transpose_map(k) >= 0) nzvals(transpose_map(k)) = nzvals_t(k); });
            }
          }
        }
        // broadcast return value
        Teuchos::broadcast<int, LocalOrdinal>(*comm, 0, 1, &ret);
      }
      return ret;
    }


  template <typename Scalar,
            typename LocalOrdinal,
            typename GlobalOrdinal,
            typename Node>
  void
  ConcreteMatrixAdapter<
    Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>
    >::describe (Teuchos::FancyOStream& os,
                 const Teuchos::EVerbosityLevel verbLevel) const
    {
      this->mat_->describe(os, verbLevel);
    }
} // end namespace Amesos2

#endif  // AMESOS2_TPETRACRSMATRIX_MATRIXADAPTER_DEF_HPP
