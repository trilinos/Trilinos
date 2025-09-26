// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifndef AMESOS2_EPETRACRSMATRIX_MATRIXADAPTER_DEF_HPP
#define AMESOS2_EPETRACRSMATRIX_MATRIXADAPTER_DEF_HPP

#include <Epetra_LocalMap.h>
#include <Epetra_Import.h>

#include "Amesos2_EpetraCrsMatrix_MatrixAdapter_decl.hpp"
#include "Amesos2_EpetraRowMatrix_AbstractMatrixAdapter_def.hpp"
#include "Amesos2_MatrixAdapter_def.hpp"


namespace Amesos2 {

  ConcreteMatrixAdapter<Epetra_CrsMatrix>::ConcreteMatrixAdapter(RCP<Epetra_CrsMatrix> m)
    : AbstractConcreteMatrixAdapter<Epetra_RowMatrix,Epetra_CrsMatrix>(m) // CrsMatrix inherits from RowMatrix virtually, so a dynamic cast is necessary
    {}

  Teuchos::RCP<const MatrixAdapter<Epetra_CrsMatrix> >
  ConcreteMatrixAdapter<Epetra_CrsMatrix>::get_impl(const Teuchos::Ptr<const map_t> map, EDistribution distribution) const
    {
      using Teuchos::as;
      using Teuchos::rcp;
      using Teuchos::rcpFromRef;

      RCP<const Epetra_Map> o_map, t_map;
      o_map = rcpFromRef(this->mat_->RowMap());
      t_map = Util::tpetra_map_to_epetra_map<local_ordinal_t,global_ordinal_t,global_size_t,node_t>(*map);

      const int maxRowNNZ = 0;
      RCP<Epetra_CrsMatrix> t_mat = rcp(new Epetra_CrsMatrix(Copy, *t_map, maxRowNNZ));

      Epetra_Import importer(*t_map, *o_map);
      t_mat->Import(*(this->mat_), importer, Insert);
      t_mat->FillComplete();

      // Case for non-contiguous GIDs
      if ( distribution == CONTIGUOUS_AND_ROOTED ) {

        auto myRank = map->getComm()->getRank();

        const int global_num_contiguous_entries = t_mat->NumGlobalRows();
        const int local_num_contiguous_entries = (myRank == 0) ? t_mat->NumGlobalRows() : 0;

        RCP<const Epetra_Map> contiguousRowMap = rcp( new Epetra_Map(global_num_contiguous_entries, local_num_contiguous_entries, 0, (t_mat->Comm() ) ) );
        RCP<const Epetra_Map> contiguousColMap = rcp( new Epetra_Map(global_num_contiguous_entries, local_num_contiguous_entries, 0, (t_mat->Comm() ) ) );
        RCP<const Epetra_Map> contiguousDomainMap = rcp( new Epetra_Map(global_num_contiguous_entries, local_num_contiguous_entries, 0, (t_mat->Comm() ) ) );
        RCP<const Epetra_Map> contiguousRangeMap  = rcp( new Epetra_Map(global_num_contiguous_entries, local_num_contiguous_entries, 0, (t_mat->Comm() ) ) );

        RCP<Epetra_CrsMatrix> contiguous_t_mat = rcp( new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *contiguousRowMap, *contiguousColMap, t_mat->MaxNumEntries()) );

        // fill local sparse matrix on rank zero
        if(myRank == 0) {
          int num_entries;
          int *indices;
          double *values;
          for (int row = 0; row < t_mat->NumMyRows(); row++) {
            t_mat->ExtractMyRowView(row, num_entries, values, indices);
            contiguous_t_mat->InsertMyValues(row, num_entries, values, indices);
          }
        }

        contiguous_t_mat->FillComplete(*contiguousDomainMap, *contiguousRangeMap);

        return rcp (new ConcreteMatrixAdapter<Epetra_CrsMatrix> (contiguous_t_mat));
      }

      return( rcp(new ConcreteMatrixAdapter<Epetra_CrsMatrix>(t_mat)) );
    }

  Teuchos::RCP<const MatrixAdapter<Epetra_CrsMatrix> >
  ConcreteMatrixAdapter<Epetra_CrsMatrix>::reindex_impl(Teuchos::RCP<const map_t> &contigRowMap,
                                                        Teuchos::RCP<const map_t> &contigColMap,
                                                        const EPhase /*current_phase*/) const
    {
      #if defined(HAVE_AMESOS2_EPETRAEXT)
      using Teuchos::RCP;
      using Teuchos::rcp;
      using Teuchos::rcpFromRef;
      auto CrsMatrix = const_cast<Epetra_CrsMatrix *>(this->mat_.getRawPtr());
      if(!CrsMatrix) {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Amesos2_EpetraCrsMatrix_MatrixAdapter requires CsrMatrix to reindex matrices.");
      }

      // Map
      RCP<const Epetra_Map> OriginalMap = rcpFromRef(CrsMatrix->RowMap());
      int NumGlobalElements = OriginalMap->NumGlobalElements();
      int NumMyElements = OriginalMap->NumMyElements();
      auto ReindexMap = rcp( new Epetra_Map( NumGlobalElements, NumMyElements, 0, OriginalMap->Comm() ) );

      // Matrix
      StdIndex_ = rcp( new EpetraExt::CrsMatrix_Reindex( *ReindexMap ) );
      ContigMat_ = rcpFromRef((*StdIndex_)( *CrsMatrix ));
      if(!ContigMat_) {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Amesos2_EpetraCrsMatrix_MatrixAdapter reindexing failed.");
      }

      auto reindexMat = rcp( new ConcreteMatrixAdapter<Epetra_CrsMatrix>(ContigMat_));
      contigRowMap = reindexMat->getRowMap();
      contigColMap = reindexMat->getColMap();

      return reindexMat;
      #else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "ConcreteMatrixAdapter<Epetra_CrsMatrix> requires EpetraExt to reindex matrices.");
      #endif
    }


  void
  ConcreteMatrixAdapter<Epetra_CrsMatrix>::describe (Teuchos::FancyOStream& os,
                 const Teuchos::EVerbosityLevel verbLevel) const
    {
      this->mat_->Print(*(os.getOStream()));
    }
} // end namespace Amesos2

#endif  // AMESOS2_EPETRACRSMATRIX_MATRIXADAPTER_DEF_HPP
