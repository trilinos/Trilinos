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
    >::get_impl(const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > map, EDistribution distribution) const
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
        typedef Tpetra::Map< local_ordinal_t, global_ordinal_t, node_t>  contiguous_map_type;
        RCP<const contiguous_map_type> contiguousRowMap = rcp( new contiguous_map_type(global_num_contiguous_entries, local_num_contiguous_entries, 0, (t_mat->getComm() ) ) );
        RCP<const contiguous_map_type> contiguousColMap = rcp( new contiguous_map_type(global_num_contiguous_entries, local_num_contiguous_entries, 0, (t_mat->getComm() ) ) );
        RCP<const contiguous_map_type> contiguousDomainMap = rcp( new contiguous_map_type(global_num_contiguous_entries, local_num_contiguous_entries, 0, (t_mat->getComm() ) ) );
        RCP<const contiguous_map_type> contiguousRangeMap = rcp( new contiguous_map_type(global_num_contiguous_entries, local_num_contiguous_entries, 0, (t_mat->getComm() ) ) );

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
    >::reindex_impl() const
    {
      typedef Kokkos::DefaultHostExecutionSpace HostExecSpaceType;
      typedef Tpetra::Map< local_ordinal_t, global_ordinal_t, node_t>  contiguous_map_type;
      auto rowMap = this->mat_->getRowMap();
      auto colMap = this->mat_->getColMap();
      auto rowComm = rowMap->getComm();
      auto colComm = colMap->getComm();

      global_ordinal_t indexBase = rowMap->getIndexBase();
      global_ordinal_t numDoFs = this->mat_->getGlobalNumRows();
      local_ordinal_t nRows = this->mat_->getLocalNumRows();
      local_ordinal_t nCols = colMap->getLocalNumElements();

      auto tmpMap = rcp (new contiguous_map_type (numDoFs, nRows, indexBase, rowComm));
      global_ordinal_t frow = tmpMap->getMinGlobalIndex();

      // Create new GID list for RowMap
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
      Teuchos::RCP<const contiguous_map_type> newRowMap = rcp (new contiguous_map_type (numDoFs, rowIndexList.data(), nRows, indexBase, rowComm));
      Teuchos::RCP<const contiguous_map_type> newColMap = rcp (new contiguous_map_type (numDoFs, colIndexList.data(), nCols, indexBase, colComm));

      // Build Matrix with new Maps, 
      auto lclMatrix = this->mat_->getLocalMatrixDevice();
      RCP<matrix_t> contiguous_t_mat = rcp( new matrix_t(newRowMap, newColMap, lclMatrix));

      return rcp (new ConcreteMatrixAdapter<matrix_t> (contiguous_t_mat));
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
      Tpetra::MatrixMarket::Writer<matrix_t>::writeSparseFile ("matA.dat", this->mat_);
    }
} // end namespace Amesos2

#endif  // AMESOS2_TPETRACRSMATRIX_MATRIXADAPTER_DEF_HPP
