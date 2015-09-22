// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package 
//                  Copyright 2011 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//
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
  ConcreteMatrixAdapter<Epetra_CrsMatrix>::get_impl(const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > map) const
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
