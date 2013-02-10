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


#ifndef AMESOS2_TPETRAROWMATRIX_ABSTRACTMATRIXADAPTER_DEF_HPP
#define AMESOS2_TPETRAROWMATRIX_ABSTRACTMATRIXADAPTER_DEF_HPP

#include "Amesos2_ConcreteMatrixAdapter_decl.hpp"
#include "Amesos2_TpetraRowMatrix_AbstractMatrixAdapter_decl.hpp"

namespace Amesos2 {

  using Teuchos::RCP;
  using Teuchos::ArrayView;

  template <typename Scalar,
	    typename LocalOrdinal,
	    typename GlobalOrdinal,
	    typename Node,
	    class DerivedMat>
  AbstractConcreteMatrixAdapter<
    Tpetra::RowMatrix<Scalar,
		      LocalOrdinal,
		      GlobalOrdinal,
		      Node>,
    DerivedMat>::AbstractConcreteMatrixAdapter(RCP<Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > m)
      : MatrixAdapter<DerivedMat>(Teuchos::rcp_static_cast<DerivedMat>(m)) 
  {
    // anything else? probs not
  }

  // implementation functions
  template <typename Scalar,
	    typename LocalOrdinal,
	    typename GlobalOrdinal,
	    typename Node,
	    class DerivedMat>
  void
  AbstractConcreteMatrixAdapter<
    Tpetra::RowMatrix<Scalar,
		      LocalOrdinal,
		      GlobalOrdinal,
		      Node>,
    DerivedMat>::getGlobalRowCopy_impl(global_ordinal_t row,
				       const ArrayView<global_ordinal_t>& indices,
				       const ArrayView<scalar_t>& vals,
				       size_t& nnz) const
    {
      this->mat_->getGlobalRowCopy(row, indices, vals, nnz);
    }

  template <typename Scalar,
	    typename LocalOrdinal,
	    typename GlobalOrdinal,
	    typename Node,
	    class DerivedMat>
  void
  AbstractConcreteMatrixAdapter<
    Tpetra::RowMatrix<Scalar,
		      LocalOrdinal,
		      GlobalOrdinal,
		      Node>,
    DerivedMat>::getGlobalColCopy_impl(global_ordinal_t col,
			     const ArrayView<global_ordinal_t>& indices,
			     const ArrayView<scalar_t>& vals,
			     size_t& nnz) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION( true,
			std::runtime_error,
			"Column access to row-based object not yet supported.  "
			"Please contact the Amesos2 developers." );
  }


  template <typename Scalar,
	    typename LocalOrdinal,
	    typename GlobalOrdinal,
	    typename Node,
	    class DerivedMat>
  typename AbstractConcreteMatrixAdapter<
    Tpetra::RowMatrix<Scalar,
		      LocalOrdinal,
		      GlobalOrdinal,
		      Node>,
    DerivedMat>::global_size_t
  AbstractConcreteMatrixAdapter<
    Tpetra::RowMatrix<Scalar,
		      LocalOrdinal,
		      GlobalOrdinal,
		      Node>,
    DerivedMat>::getGlobalNNZ_impl() const
  {
    return this->mat_->getGlobalNumEntries();
  }

  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node, class DerivedMat>
  size_t
  AbstractConcreteMatrixAdapter<
    Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>,
    DerivedMat>::getLocalNNZ_impl() const
  {
    return this->mat_->getNodeNumEntries();
  }

  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node, class DerivedMat>
  typename AbstractConcreteMatrixAdapter<
    Tpetra::RowMatrix<Scalar,
		      LocalOrdinal,
		      GlobalOrdinal,
		      Node>,
    DerivedMat>::global_size_t
  AbstractConcreteMatrixAdapter<
    Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>,
    DerivedMat>::getGlobalNumRows_impl() const
  {
    return this->mat_->getGlobalNumRows();
  }

  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node, class DerivedMat>
  typename AbstractConcreteMatrixAdapter<
    Tpetra::RowMatrix<Scalar,
		      LocalOrdinal,
		      GlobalOrdinal,
		      Node>,
    DerivedMat>::global_size_t
  AbstractConcreteMatrixAdapter<
    Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>,
    DerivedMat>::getGlobalNumCols_impl() const
  {
    return this->mat_->getGlobalNumCols();
  }

  template <typename Scalar,
	    typename LocalOrdinal,
	    typename GlobalOrdinal,
	    typename Node,
	    class DerivedMat>
  size_t
  AbstractConcreteMatrixAdapter<
    Tpetra::RowMatrix<Scalar,
		      LocalOrdinal,
		      GlobalOrdinal,
		      Node>,
    DerivedMat>::getMaxRowNNZ_impl() const
  {
    return this->mat_->getGlobalMaxNumRowEntries();
  }

  template <typename Scalar,
	    typename LocalOrdinal,
	    typename GlobalOrdinal,
	    typename Node,
	    class DerivedMat>
  size_t
  AbstractConcreteMatrixAdapter<
    Tpetra::RowMatrix<Scalar,
		      LocalOrdinal,
		      GlobalOrdinal,
		      Node>,
    DerivedMat>::getMaxColNNZ_impl() const
  {
    TEUCHOS_TEST_FOR_EXCEPTION( true,
			std::runtime_error,
			"Column access to row-based object not yet supported.  "
			"Please contact the Amesos2 developers." );
    return 0;
  }

  template <typename Scalar,
	    typename LocalOrdinal,
	    typename GlobalOrdinal,
	    typename Node,
	    class DerivedMat>
  size_t
  AbstractConcreteMatrixAdapter<
    Tpetra::RowMatrix<Scalar,
		      LocalOrdinal,
		      GlobalOrdinal,
		      Node>,
    DerivedMat>::getGlobalRowNNZ_impl(global_ordinal_t row) const
  {
    return this->mat_->getNumEntriesInGlobalRow(row);
  }

  template <typename Scalar,
	    typename LocalOrdinal,
	    typename GlobalOrdinal,
	    typename Node,
	    class DerivedMat>
  size_t
  AbstractConcreteMatrixAdapter<
    Tpetra::RowMatrix<Scalar,
		      LocalOrdinal,
		      GlobalOrdinal,
		      Node>,
    DerivedMat>::getLocalRowNNZ_impl(local_ordinal_t row) const
  {
    return this->mat_->getNumEntriesInLocalRow(row);
  }

  template <typename Scalar,
	    typename LocalOrdinal,
	    typename GlobalOrdinal,
	    typename Node,
	    class DerivedMat>
  size_t
  AbstractConcreteMatrixAdapter<
    Tpetra::RowMatrix<Scalar,
		      LocalOrdinal,
		      GlobalOrdinal,
		      Node>,
    DerivedMat>::getGlobalColNNZ_impl(global_ordinal_t col) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION( true,
			std::runtime_error,
			"Column access to row-based object not yet supported.  "
			"Please contact the Amesos2 developers." );
    return 0;
  }

  template <typename Scalar,
	    typename LocalOrdinal,
	    typename GlobalOrdinal,
	    typename Node,
	    class DerivedMat>
  size_t
  AbstractConcreteMatrixAdapter<
    Tpetra::RowMatrix<Scalar,
		      LocalOrdinal,
		      GlobalOrdinal,
		      Node>,
    DerivedMat>::getLocalColNNZ_impl(local_ordinal_t col) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION( true,
			std::runtime_error,
			"Column access to row-based object not yet supported.  "
			"Please contact the Amesos2 developers." );
    return 0;
  }

  template <typename Scalar,
	    typename LocalOrdinal,
	    typename GlobalOrdinal,
	    typename Node,
	    class DerivedMat>
  const RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
  AbstractConcreteMatrixAdapter<
    Tpetra::RowMatrix<Scalar,
		      LocalOrdinal,
		      GlobalOrdinal,
		      Node>,
    DerivedMat>:: getRowMap_impl() const
  {
    return this->mat_->getRowMap();
  }
  
  template <typename Scalar,
	    typename LocalOrdinal,
	    typename GlobalOrdinal,
	    typename Node,
	    class DerivedMat>
  const RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
  AbstractConcreteMatrixAdapter<
    Tpetra::RowMatrix<Scalar,
		      LocalOrdinal,
		      GlobalOrdinal,
		      Node>,
    DerivedMat>::getColMap_impl() const
  {
    return this->mat_->getColMap();
  }

  template <typename Scalar,
	    typename LocalOrdinal,
	    typename GlobalOrdinal,
	    typename Node,
	    class DerivedMat>
  const RCP<const Teuchos::Comm<int> >
  AbstractConcreteMatrixAdapter<
    Tpetra::RowMatrix<Scalar,
		      LocalOrdinal,
		      GlobalOrdinal,
		      Node>,
    DerivedMat>::getComm_impl() const
  {
    return this->mat_->getComm();
  }

  template <typename Scalar,
	    typename LocalOrdinal,
	    typename GlobalOrdinal,
	    typename Node,
	    class DerivedMat>
  bool
  AbstractConcreteMatrixAdapter<
    Tpetra::RowMatrix<Scalar,
		      LocalOrdinal,
		      GlobalOrdinal,
		      Node>,
    DerivedMat>::isLocallyIndexed_impl() const
  {
    return this->mat_->isLocallyIndexed();
  }
  
  template <typename Scalar,
	    typename LocalOrdinal,
	    typename GlobalOrdinal,
	    typename Node,
	    class DerivedMat>
  bool
  AbstractConcreteMatrixAdapter<
    Tpetra::RowMatrix<Scalar,
		      LocalOrdinal,
		      GlobalOrdinal,
		      Node>,
    DerivedMat>::isGloballyIndexed_impl() const
  {
    return this->mat_->isGloballyIndexed();
  }

  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node, class DerivedMat>
  RCP<const MatrixAdapter<DerivedMat> >
  AbstractConcreteMatrixAdapter<
    Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>, DerivedMat
    >::get_impl(const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > map) const
  {
    return static_cast<ConcreteMatrixAdapter<DerivedMat>*>(this)->get_impl(map);
  }

} // end namespace Amesos2

#endif	// AMESOS2_TPETRAROWMATRIX_ABSTRACTMATRIXADAPTER_DEF_HPP
