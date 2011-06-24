#ifndef AMESOS2_TPETRAROWMATRIX_ABSTRACTMATRIXADAPTER_DEF_HPP
#define AMESOS2_TPETRAROWMATRIX_ABSTRACTMATRIXADAPTER_DEF_HPP

namespace Amesos {

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
    TEST_FOR_EXCEPTION( true,
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
    TEST_FOR_EXCEPTION( true,
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
    TEST_FOR_EXCEPTION( true,
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
    TEST_FOR_EXCEPTION( true,
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

} // end namespace Amesos

#endif	// AMESOS2_TPETRAROWMATRIX_ABSTRACTMATRIXADAPTER_DEF_HPP
