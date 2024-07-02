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
  ConcreteMatrixAdapter<Epetra_CrsMatrix>::get_impl(const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > map, EDistribution /* distribution */) const
    {
      using Teuchos::as;
      using Teuchos::rcp;
      using Teuchos::rcpFromRef;

      RCP<const Epetra_Map> o_map, t_map;
      o_map = rcpFromRef(this->mat_->RowMap());
      t_map = Util::tpetra_map_to_epetra_map<local_ordinal_t,global_ordinal_t,global_size_t,node_t>(*map);

      RCP<Epetra_CrsMatrix> t_mat = rcp(new Epetra_CrsMatrix(Copy, *t_map, this->getMaxRowNNZ()));

      Epetra_Import importer(*t_map, *o_map);
      t_mat->Import(*(this->mat_), importer, Insert);

      t_mat->FillComplete();    // Must be in local form for later extraction of rows

      return( rcp(new ConcreteMatrixAdapter<Epetra_CrsMatrix>(t_mat)) );
    }

} // end namespace Amesos2

#endif  // AMESOS2_EPETRACRSMATRIX_MATRIXADAPTER_DEF_HPP
