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


/**
 * \file   Amesos2_EpetraRowMatrix_AbstractMatrixAdapter_def.hpp
 * \author Eric Bavier <etbavie@sandia.gov>
 * \date   Tue Jun 14 17:21:32 2011
 *
 * \brief  Definitions for the Epetra_RowMatrix abstract adapter.
 */

#ifndef AMESOS2_EPETRAROWMATRIX_ABSTRACTMATRIXADAPTER_DEF_HPP
#define AMESOS2_EPETRAROWMATRIX_ABSTRACTMATRIXADAPTER_DEF_HPP

#include <Epetra_RowMatrix.h>
#include <Epetra_Map.h>
#include <Epetra_Comm.h>

#include "Amesos2_EpetraRowMatrix_AbstractMatrixAdapter_decl.hpp"


namespace Amesos2 {

  using Teuchos::RCP;
  using Teuchos::ArrayView;

  template <class DerivedMat>
  AbstractConcreteMatrixAdapter<Epetra_RowMatrix,DerivedMat>::AbstractConcreteMatrixAdapter(RCP<DerivedMat> m)
      : MatrixAdapter<DerivedMat>(m)
  {
    // anything else? probs not
  }

  // implementation functions
  template <class DerivedMat>
  void
  AbstractConcreteMatrixAdapter<
    Epetra_RowMatrix,
    DerivedMat>::getGlobalRowCopy_impl(global_ordinal_t row,
				       const ArrayView<global_ordinal_t>& indices,
				       const ArrayView<scalar_t>& vals,
				       size_t& nnz) const
  {
    using Teuchos::as;
    
    local_ordinal_t local_row = this->row_map_->getLocalElement(row);
    int nnz_ret = 0;
    int rowmatrix_return_val
      = this->mat_->ExtractMyRowCopy(as<int>(local_row),
				     as<int>(std::min(indices.size(), vals.size())),
				     nnz_ret,
				     vals.getRawPtr(),
				     indices.getRawPtr());
    TEST_FOR_EXCEPTION( rowmatrix_return_val != 0,
			std::runtime_error,
			"Epetra_RowMatrix object returned error code "
			<< rowmatrix_return_val << " from ExtractMyRowCopy." );
    nnz = as<size_t>(nnz_ret);

    // Epetra_CrsMatrix::ExtractMyRowCopy returns local column
    // indices, so transform these into global indices
    for( size_t i = 0; i < nnz; ++i ){
      indices[i] = this->col_map_->getGlobalElement(indices[i]);
    }
  }

  template <class DerivedMat>
  void
  AbstractConcreteMatrixAdapter<
    Epetra_RowMatrix,
    DerivedMat>::getGlobalColCopy_impl(global_ordinal_t col,
				       const ArrayView<global_ordinal_t>& indices,
				       const ArrayView<scalar_t>& vals,
				       size_t& nnz) const
  {
    TEST_FOR_EXCEPTION( true,
			std::runtime_error,
			"Column access to row-based object not yet supported.  "
			"Please contact the Amesos2 developers." );
  }


  template <class DerivedMat>
  typename AbstractConcreteMatrixAdapter<
    Epetra_RowMatrix,
    DerivedMat>::global_size_t
  AbstractConcreteMatrixAdapter<
    Epetra_RowMatrix,
    DerivedMat>::getGlobalNNZ_impl() const
  {
    return Teuchos::as<global_size_t>(this->mat_->NumGlobalNonzeros());
  }

  template <class DerivedMat>
  size_t
  AbstractConcreteMatrixAdapter<
    Epetra_RowMatrix,
    DerivedMat>::getLocalNNZ_impl() const
  {
    return Teuchos::as<size_t>(this->mat_->NumMyNonzeros());
  }

  template <class DerivedMat>
  size_t
  AbstractConcreteMatrixAdapter<
    Epetra_RowMatrix,
    DerivedMat>::getMaxRowNNZ_impl() const
  {
    return Teuchos::as<size_t>(this->mat_->MaxNumEntries());
  }

  template <class DerivedMat>
  size_t
  AbstractConcreteMatrixAdapter<
    Epetra_RowMatrix,
    DerivedMat>::getMaxColNNZ_impl() const
  {
    TEST_FOR_EXCEPTION( true,
			std::runtime_error,
			"Column access to row-based object not yet supported.  "
			"Please contact the Amesos2 developers." );
    return 0;
  }

  template <class DerivedMat>
  size_t
  AbstractConcreteMatrixAdapter<
    Epetra_RowMatrix,
    DerivedMat>::getGlobalRowNNZ_impl(global_ordinal_t row) const
  {
    // check whether row is local, then transform to local index
    Epetra_Map rowmap = this->mat_->RowMatrixRowMap();
    int gid = Teuchos::as<int>(row);
    TEST_FOR_EXCEPTION( !rowmap.MyGID(gid),
			std::invalid_argument,
			"The specified global row id does not belong to me" );
    int lid = rowmap.LID(gid);
    int nnz = 0;
    this->mat_->NumMyRowEntries(lid, nnz);
    return nnz;
  }

  template <class DerivedMat>
  size_t
  AbstractConcreteMatrixAdapter<
    Epetra_RowMatrix,
    DerivedMat>::getLocalRowNNZ_impl(local_ordinal_t row) const
  {
    Epetra_Map rowmap = this->mat_->RowMatrixRowMap();
    int lid = Teuchos::as<int>(row);
    TEST_FOR_EXCEPTION( !rowmap.MyLID(lid),
			std::invalid_argument,
			"The specified local row id does not beloing to me" );
    int num_entries = 0;
    this->mat_->NumMyRowEntries(row, num_entries);
    return num_entries;
  }

  template <class DerivedMat>
  size_t
  AbstractConcreteMatrixAdapter<
    Epetra_RowMatrix,
    DerivedMat>::getGlobalColNNZ_impl(global_ordinal_t col) const
  {
    TEST_FOR_EXCEPTION( true,
			std::runtime_error,
			"Column access to row-based object not yet supported.  "
			"Please contact the Amesos2 developers." );
    return 0;
  }

  template <class DerivedMat>
  size_t
  AbstractConcreteMatrixAdapter<
    Epetra_RowMatrix,
    DerivedMat>::getLocalColNNZ_impl(local_ordinal_t col) const
  {
    TEST_FOR_EXCEPTION( true,
			std::runtime_error,
			"Column access to row-based object not yet supported.  "
			"Please contact the Amesos2 developers." );
    return 0;
  }

  template <class DerivedMat>
  const RCP<const Tpetra::Map<MatrixTraits<Epetra_RowMatrix>::local_ordinal_t,
			      MatrixTraits<Epetra_RowMatrix>::global_ordinal_t,
			      MatrixTraits<Epetra_RowMatrix>::node_t> >
  AbstractConcreteMatrixAdapter<Epetra_RowMatrix, DerivedMat>::getRowMap_impl() const
  {
    // Must transform to a Tpetra::Map
    const Epetra_Map rowmap = this->mat_->RowMatrixRowMap();
    return( Util::epetra_map_to_tpetra_map<local_ordinal_t,global_ordinal_t,global_size_t,node_t>(rowmap) );
  }

  template <class DerivedMat>
  const RCP<const Tpetra::Map<MatrixTraits<Epetra_RowMatrix>::local_ordinal_t,
			      MatrixTraits<Epetra_RowMatrix>::global_ordinal_t,
			      MatrixTraits<Epetra_RowMatrix>::node_t> >
  AbstractConcreteMatrixAdapter<Epetra_RowMatrix, DerivedMat>::getColMap_impl() const
  {
    // Must transform this matrix' Epetra_Map to a Tpetra::Map
    const Epetra_Map colmap = this->mat_->RowMatrixColMap();
    return( Util::epetra_map_to_tpetra_map<local_ordinal_t,global_ordinal_t,global_size_t,node_t>(colmap) );
  }

  template <class DerivedMat>
  const RCP<const Teuchos::Comm<int> >
  AbstractConcreteMatrixAdapter<Epetra_RowMatrix, DerivedMat>::getComm_impl() const
  {
    return Util::to_teuchos_comm(Teuchos::rcpFromRef(this->mat_->Comm()));
  }

  template <class DerivedMat>
  bool
  AbstractConcreteMatrixAdapter<Epetra_RowMatrix, DerivedMat>::isLocallyIndexed_impl() const
  {
    return this->mat_->IndicesAreLocal();
  }

  template <class DerivedMat>
  bool
  AbstractConcreteMatrixAdapter<Epetra_RowMatrix, DerivedMat>::isGloballyIndexed_impl() const
  {
    return this->mat_->IndicesAreGlobal();
  }

  template <class DerivedMat>
  RCP<const MatrixAdapter<DerivedMat> >
  AbstractConcreteMatrixAdapter<Epetra_RowMatrix, DerivedMat>::get_impl(const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > map) const
  {
    // Delegate implementation to subclass
    return static_cast<ConcreteMatrixAdapter<DerivedMat>*>(this)->get_impl(map);
  }

} // end namespace Amesos2

#endif  // AMESOS2_EPETRAROWMATRIX_ABSTRACTMATRIXADAPTER_DEF_HPP
