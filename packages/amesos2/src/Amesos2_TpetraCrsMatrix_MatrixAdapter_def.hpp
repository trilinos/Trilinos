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

} // end namespace Amesos2

#endif  // AMESOS2_TPETRACRSMATRIX_MATRIXADAPTER_DEF_HPP
