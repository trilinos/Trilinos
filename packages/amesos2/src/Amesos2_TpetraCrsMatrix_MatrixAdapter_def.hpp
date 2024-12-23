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
                    Teuchos::RCP<const map_t> &contigColMap) const
    {
      typedef Tpetra::Map< local_ordinal_t, global_ordinal_t, node_t> contiguous_map_type;
      auto rowMap = this->mat_->getRowMap();
      auto colMap = this->mat_->getColMap();
      auto rowComm = rowMap->getComm();
      auto colComm = colMap->getComm();

#ifdef HAVE_AMESOS2_TIMERS
      auto reindexTimer = Teuchos::TimeMonitor::getNewTimer("Time to re-index matrix gids");
      Teuchos::TimeMonitor ReindexTimer(*reindexTimer);
#endif

      global_ordinal_t indexBase = rowMap->getIndexBase();
      global_ordinal_t numDoFs = this->mat_->getGlobalNumRows();
      local_ordinal_t nRows = this->mat_->getLocalNumRows();
      local_ordinal_t nCols = colMap->getLocalNumElements();

      RCP<matrix_t> contiguous_t_mat;
      // if-checks when to recompute contigRowMap & contigColMap
      // TODO: this is currentlly based on the global matrix dimesions
      if ((contigRowMap.is_null() || contigColMap.is_null()) ||
          (contigRowMap->getGlobalNumElements() != numDoFs || contigColMap->getGlobalNumElements() != numDoFs)) {
        auto tmpMap = rcp (new contiguous_map_type (numDoFs, nRows, indexBase, rowComm));
        global_ordinal_t frow = tmpMap->getMinGlobalIndex();

        // Create new GID list for RowMap
        typedef Kokkos::DefaultHostExecutionSpace HostExecSpaceType;
        Kokkos::View<global_ordinal_t*, HostExecSpaceType> rowIndexList ("indexList", nRows);
        for (local_ordinal_t k = 0; k < nRows; k++) {
          rowIndexList(k) = frow+k;
        }
        // Create new GID list for ColMap
        Kokkos::View<global_ordinal_t*, HostExecSpaceType> colIndexList ("indexList", nCols);
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
          auto col_view = col_mv.getLocalViewHost(Tpetra::Access::ReadOnly);
          for(int i=0; i<nCols; i++) colIndexList(i) = col_view(i,0);
        }
        // Create new Row & Col Maps
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
      return rcp (new ConcreteMatrixAdapter<matrix_t> (contiguous_t_mat));
    }


  template <typename Scalar,
            typename LocalOrdinal,
            typename GlobalOrdinal,
            typename Node>
  template<typename KV_S, typename KV_GO, typename KV_GS>
  LocalOrdinal
  ConcreteMatrixAdapter<
    Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>
    >::gather_impl(KV_S& nzvals, KV_GO& indices, KV_GS& pointers, bool column_major, EPhase current_phase) const
    {
      LocalOrdinal ret = 0;
      {
#ifdef HAVE_AMESOS2_TIMERS
        Teuchos::RCP< Teuchos::Time > gatherTime = Teuchos::TimeMonitor::getNewCounter ("Amesos2::gather");
        Teuchos::TimeMonitor GatherTimer(*gatherTime);
#endif
        auto lclMatrix = this->mat_->getLocalMatrixHost();
        auto lclRowptr = lclMatrix.graph.row_map;
        auto lclColind = lclMatrix.graph.entries;
        auto lclNzvals = lclMatrix.values;


        auto rowMap = this->mat_->getRowMap();
        auto colMap = this->mat_->getColMap();
        auto comm = rowMap->getComm();
        auto nRanks = comm->getSize();
        auto myRank = comm->getRank();

        global_ordinal_t nRows = this->mat_->getGlobalNumRows();
        if(current_phase == SYMBFACT) {
          // workspace to transpose
          KV_GS pointers_t;
          KV_GO indices_t;
          // gether rowptr
          LocalOrdinal myNRows = this->mat_->getLocalNumRows();
          Kokkos::resize(recvCounts, nRanks);
          Kokkos::resize(recvDispls, nRanks+1);
          {
#ifdef HAVE_AMESOS2_TIMERS
            Teuchos::RCP< Teuchos::Time > gatherTime = Teuchos::TimeMonitor::getNewCounter ("Amesos2::gather(rowptr)");
            Teuchos::TimeMonitor GatherTimer(*gatherTime);
#endif
            Teuchos::gather<int, LocalOrdinal> (&myNRows, 1, recvCounts.data(), 1, 0, *comm);
            if (myRank == 0) {
              Kokkos::resize(recvDispls, nRanks+1);
              recvDispls(0) = 0;
              for (int p = 1; p <= nRanks; p++) {
                recvDispls(p) = recvDispls(p-1) + recvCounts(p-1);
              }
              fflush(stdout);
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
            // - making sure same ordinal type
            KV_GS lclRowptr_ ("localRowptr_", lclRowptr.extent(0));
            for (int i = 0; i < int(lclRowptr.extent(0)); i++) lclRowptr_(i) = lclRowptr(i);
            if (myRank == 0 && column_major) {
              Kokkos::resize(pointers_t, nRows+1);
            } else if (myRank != 0) {
              Kokkos::resize(pointers_t, 2);
            }
            LocalOrdinal *pointers_ = (myRank != 0 || column_major ? pointers_t.data() : pointers.data());
            Teuchos::gatherv<int, LocalOrdinal> (&lclRowptr_(1), myNRows, &pointers_[1], 
                                                 recvCounts.data(), recvDispls.data(),
                                                 0, *comm);

            if (myRank == 0) {
              // shift to global pointers
              pointers_[0] = 0;
              recvCounts(0) = pointers_[recvDispls(1)];
              LocalOrdinal displs = recvCounts(0);
              for (int p = 1; p < nRanks; p++) {
                // save recvCounts from pth MPI
                recvCounts(p) = pointers_[recvDispls(p+1)];
                // shift pointers for pth MPI to global
                for (int i = 1+recvDispls(p); i <= recvDispls(p+1); i++) {
                  pointers_[i] += displs;
                }
                displs += recvCounts(p);
              }
              ret = pointers_[nRows];
            }
          }
          {
#ifdef HAVE_AMESOS2_TIMERS
            Teuchos::RCP< Teuchos::Time > gatherTime = Teuchos::TimeMonitor::getNewCounter ("Amesos2::gather(colind)");
            Teuchos::TimeMonitor GatherTimer(*gatherTime);
#endif
            // gather colinds & nzvals
            if (myRank == 0) {
              recvDispls(0) = 0;
              for (int p = 0; p < nRanks; p++) {
                recvDispls(p+1) = recvDispls(p) + recvCounts(p);
              }
            }
            // -- convert to global colids
            KV_GO lclColind_ ("localColind_", lclColind.extent(0));
            for (int i = 0; i < int(lclColind.extent(0)); i++) lclColind_(i) = colMap->getGlobalElement((lclColind(i)));
            if (column_major) {
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
          if (myRank == 0 && column_major) {
#ifdef HAVE_AMESOS2_TIMERS
            Teuchos::RCP< Teuchos::Time > gatherTime = Teuchos::TimeMonitor::getNewCounter ("Amesos2::gather(transpose index)");
            Teuchos::TimeMonitor GatherTimer(*gatherTime);
#endif
            // Map to transpose
            Kokkos::resize(this->transpose_map_, ret);
            // Transopose to convert to CSC
            for (int i=0; i<=nRows; i++) {
              pointers(i) = 0;
            }
            for (int k=0; k<ret; k++) {
              if (indices_t(k) < nRows-1) {
                pointers(indices_t(k)+2) ++;
              }
            }
            for (int i=1; i < nRows; i++) {
              pointers(i+1) += pointers(i);
            }
            for (int i=0; i<nRows; i++) {
              for (int k=pointers_t(i); k<pointers_t(i+1); k++) {
                this->transpose_map_(k) = pointers(1+indices_t(k));
                indices(pointers(1+indices_t(k))) = i;
                pointers(1+indices_t(k)) ++;
              }
            }
          }
        }
        //if(current_phase == NUMFACT) // Numerical values may be used in symbolic (e.g, MWM)
        {
          // workspace to transpose
          {
#ifdef HAVE_AMESOS2_TIMERS
            Teuchos::RCP< Teuchos::Time > gatherTime = Teuchos::TimeMonitor::getNewCounter ("Amesos2::gather(nzvals)");
            Teuchos::TimeMonitor GatherTimer(*gatherTime);
#endif
            if (column_major) {
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
          if (myRank == 0 && column_major) {
            // Insert Numerical values to transopose matrix
            ret = pointers(nRows);
#ifdef HAVE_AMESOS2_TIMERS
            Teuchos::RCP< Teuchos::Time > gatherTime = Teuchos::TimeMonitor::getNewCounter ("Amesos2::gather(transpose values)");
            Teuchos::TimeMonitor GatherTimer(*gatherTime);
#endif
            for (int k=0; k<ret; k++) {
              nzvals(this->transpose_map_(k)) = nzvals_t(k);
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
