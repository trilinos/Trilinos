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


#ifndef AMESOS2_TPETRACRSMATRIX_MATRIXADAPTER_DEF_HPP
#define AMESOS2_TPETRACRSMATRIX_MATRIXADAPTER_DEF_HPP

#include "Amesos2_TpetraCrsMatrix_MatrixAdapter_decl.hpp"
#include "Amesos2_TpetraRowMatrix_AbstractMatrixAdapter_def.hpp"
#include "Amesos2_MatrixAdapter_def.hpp"

namespace Amesos2 {

  template <typename Scalar,
	    typename LocalOrdinal,
	    typename GlobalOrdinal,
	    typename Node,
	    typename LocalMatOps >
  ConcreteMatrixAdapter<
    Tpetra::CrsMatrix<Scalar,
		      LocalOrdinal,
		      GlobalOrdinal,
		      Node,
		      LocalMatOps>
    >::ConcreteMatrixAdapter(Teuchos::RCP<matrix_t> m) 
      : AbstractConcreteMatrixAdapter<Tpetra::RowMatrix<Scalar,
							LocalOrdinal,
							GlobalOrdinal,
							Node>,
				      Tpetra::CrsMatrix<Scalar,
							LocalOrdinal,
							GlobalOrdinal,
							Node,
							LocalMatOps> >(m) // with implicit cast
    {}

  template <typename Scalar,
	    typename LocalOrdinal,
	    typename GlobalOrdinal,
	    typename Node,
	    typename LocalMatOps >
  Teuchos::RCP<const MatrixAdapter<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > >
  ConcreteMatrixAdapter<
    Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>
    >::get_impl(const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > map) const
    {
      using Teuchos::RCP;
      using Teuchos::rcp;
      using Teuchos::rcpFromPtr;
      typedef Tpetra::Import<local_ordinal_t, global_ordinal_t, node_t> import_t;

      RCP<matrix_t> t_mat;
      t_mat = rcp (new matrix_t (rcpFromPtr (map), this->getMaxRowNNZ()));
      
      RCP<import_t> importer = 
	rcp (new import_t (this->getRowMap(), rcpFromPtr (map)));
      
      // mfh 27 Mar 2012: INSERT is correct in this case, because
      // we're Importing into an empty matrix with a dynamic graph.
      t_mat->doImport (*(this->mat_), *importer, Tpetra::INSERT);
      
      return rcp (new ConcreteMatrixAdapter<matrix_t> (t_mat));
    }

} // end namespace Amesos2

#endif	// AMESOS2_TPETRACRSMATRIX_MATRIXADAPTER_DEF_HPP
