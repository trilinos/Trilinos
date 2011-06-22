#ifndef AMESOS2_MATRIXADAPTER_DEF_HPP
#define AMESOS2_MATRIXADAPTER_DEF_HPP

namespace Amesos {

  
  template < class Matrix >
  MatrixAdapter<Matrix>::MatrixAdapter(RCP<Matrix> m)
    : mat_(m)
  {
    comm_ = static_cast<const adapter_t*>(this)->getComm_impl();
    col_map_ = static_cast<const adapter_t*>(this)->getColMap_impl();
    row_map_ = static_cast<const adapter_t*>(this)->getRowMap_impl();
  }

  // implement virtual base class functions
  template < class Matrix >
  void
  MatrixAdapter<Matrix>::getCrs(const Teuchos::ArrayView<scalar_t> nzval,
				const Teuchos::ArrayView<global_ordinal_t> colind,
				const Teuchos::ArrayView<typename MatrixAdapter<Matrix>::global_size_t> rowptr,
				typename MatrixAdapter<Matrix>::global_size_t& nnz,
				Util::EDistribution distribution,
				Util::EStorage_Ordering ordering) const
  {
    help_getCrs(nzval, colind, rowptr,
		nnz, distribution, ordering,
		typename adapter_t::get_crs_spec());
  }

  template < class Matrix >
  void
  MatrixAdapter<Matrix>::getCcs(const Teuchos::ArrayView<scalar_t> nzval,
				const Teuchos::ArrayView<global_ordinal_t> colind,
				const Teuchos::ArrayView<typename MatrixAdapter<Matrix>::global_size_t> rowptr,
				typename MatrixAdapter<Matrix>::global_size_t& nnz,
				Util::EDistribution distribution,
				Util::EStorage_Ordering ordering) const
  {
    help_getCcs(nzval, colind, rowptr,
		nnz, distribution, ordering,
		typename adapter_t::get_ccs_spec());
  }

  template < class Matrix >
  typename MatrixAdapter<Matrix>::global_size_t
  MatrixAdapter<Matrix>::getGlobalNumRows() const
  {
    return row_map_->getMaxAllGlobalIndex() + 1;
  }

  template < class Matrix >
  typename MatrixAdapter<Matrix>::global_size_t
  MatrixAdapter<Matrix>::getGlobalNumCols() const
  {
    return col_map_->getMaxAllGlobalIndex() + 1;
  }
  
  template < class Matrix >
  typename MatrixAdapter<Matrix>::global_size_t
  MatrixAdapter<Matrix>::getGlobalNNZ() const
  {
    return static_cast<const adapter_t*>(this)->getGlobalNNZ_impl();
  }
  
  template < class Matrix >
  std::string
  MatrixAdapter<Matrix>::description() const
  {
    std::ostringstream oss;
    oss << "Amesos2::MatrixAdapter wrapping: ";
    oss << mat_->description();
    return oss.str();
  }
  
  template < class Matrix >
  void
  MatrixAdapter<Matrix>::describe(Teuchos::FancyOStream &out,
				  const Teuchos::EVerbosityLevel verbLevel) const
  {}
  


  /******************************
   * Private method definitions *
   ******************************/

  // Of course, all these helper functions have more parameters than this.
  template < class Matrix >
  void
  MatrixAdapter<Matrix>::help_getCrs(const Teuchos::ArrayView<scalar_t> nzval,
				     const Teuchos::ArrayView<global_ordinal_t> colind,
				     const Teuchos::ArrayView<typename MatrixAdapter<Matrix>::global_size_t> rowptr,
				     typename MatrixAdapter<Matrix>::global_size_t& nnz,
				     Util::EDistribution distribution,
				     Util::EStorage_Ordering ordering,
				     has_special_impl hsi) const
  {
    static_cast<const adapter_t*>(this)->getCrs_spec(nzval, colind, rowptr,
						     nnz, distribution, ordering);
  }

  template < class Matrix >
  void
  MatrixAdapter<Matrix>::help_getCrs(const Teuchos::ArrayView<scalar_t> nzval,
				     const Teuchos::ArrayView<global_ordinal_t> colind,
				     const Teuchos::ArrayView<typename MatrixAdapter<Matrix>::global_size_t> rowptr,
				     typename MatrixAdapter<Matrix>::global_size_t& nnz,
				     Util::EDistribution distribution,
				     Util::EStorage_Ordering ordering,
				     no_special_impl nsi) const
  {
    do_getCrs(nzval, colind, rowptr,
	      nnz, distribution, ordering,
	      typename adapter_t::major_access());
  }

  template < class Matrix >
  void
  MatrixAdapter<Matrix>::do_getCrs(const Teuchos::ArrayView<scalar_t> nzval,
				   const Teuchos::ArrayView<global_ordinal_t> colind,
				   const Teuchos::ArrayView<typename MatrixAdapter<Matrix>::global_size_t> rowptr,
				   typename MatrixAdapter<Matrix>::global_size_t& nnz,
				   Util::EDistribution distribution,
				   Util::EStorage_Ordering ordering,
				   row_access ra) const
  {
    using Teuchos::RCP;
    using Teuchos::ArrayView;
    using Teuchos::OrdinalTraits;
    
    RCP<const type> get_mat = get(distribution);

    RCP<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > rmap = get_mat->getRowMap();

    ArrayView<const global_ordinal_t> node_elements = rmap->getNodeElementList();
    typename ArrayView<const global_ordinal_t>::iterator row_it, row_end;
    row_end = node_elements.end();

    size_t rowptr_ind = OrdinalTraits<size_t>::zero();
    global_ordinal_t rowInd = OrdinalTraits<global_ordinal_t>::zero();
    for( row_it = node_elements.begin(); row_it != row_end; ++row_it ){
      rowptr[rowptr_ind++] = rowInd;
      size_t rowNNZ = get_mat->getGlobalRowNNZ(*row_it);
      size_t nnzRet = OrdinalTraits<size_t>::zero();
      ArrayView<global_ordinal_t> colind_view = colind.view(rowInd,rowNNZ);
      ArrayView<scalar_t> nzval_view = nzval.view(rowInd,rowNNZ);
      
      get_mat->getGlobalRowCopy(*row_it, colind_view, nzval_view, nnzRet);

      // It was suggested that instead of sorting each row's indices
      // individually, that we instead do a double-transpose at the
      // end, which would also lead to the indices being sorted.
      if( ordering == Util::Sorted_Indices ){
	Tpetra::sort2(colind_view.begin(), colind_view.end(), nzval_view.begin());
      }
      
      TEST_FOR_EXCEPTION( rowNNZ != nnzRet,
			  std::runtime_error,
			  "Number of values returned different from "
                          "number of values reported");
      rowInd += rowNNZ;
    }
    rowptr[rowptr_ind] = rowInd;
    Teuchos::reduceAll(*comm_,
		       Teuchos::REDUCE_SUM,
		       Teuchos::as<global_size_t>(rowInd),
		       Teuchos::ptrFromRef(nnz));
  }

  template < class Matrix >
  void
  MatrixAdapter<Matrix>::do_getCrs(const Teuchos::ArrayView<scalar_t> nzval,
				   const Teuchos::ArrayView<global_ordinal_t> colind,
				   const Teuchos::ArrayView<typename MatrixAdapter<Matrix>::global_size_t> rowptr,
				   typename MatrixAdapter<Matrix>::global_size_t& nnz,
				   Util::EDistribution distribution,
				   Util::EStorage_Ordering ordering,
				   col_access ca) const
  {
    using Teuchos::Array;
    // get the ccs and transpose

    Array<scalar_t> nzval_tmp(nzval.size(), 0);
    Array<global_ordinal_t> rowind(colind.size(), 0);
    Array<global_size_t> colptr(this->getGlobalNumCols() + 1);
    this->getCcs(nzval_tmp(), rowind(), colptr(), nnz, distribution, ordering);
    
    Util::transpose(nzval_tmp(), rowind(), colptr(), nzval, colind, rowptr);
  }

  template < class Matrix >
  void
  MatrixAdapter<Matrix>::help_getCcs(const Teuchos::ArrayView<scalar_t> nzval,
				     const Teuchos::ArrayView<global_ordinal_t> rowind,
				     const Teuchos::ArrayView<typename MatrixAdapter<Matrix>::global_size_t> colptr,
				     typename MatrixAdapter<Matrix>::global_size_t& nnz,
				     Util::EDistribution distribution,
				     Util::EStorage_Ordering ordering,
				     has_special_impl hsi) const
  {
    static_cast<const adapter_t*>(this)->getCcs_spec(nzval, rowind, colptr,
						     nnz, distribution, ordering);
  }

  template < class Matrix >
  void
  MatrixAdapter<Matrix>::help_getCcs(const Teuchos::ArrayView<scalar_t> nzval,
				     const Teuchos::ArrayView<global_ordinal_t> rowind,
				     const Teuchos::ArrayView<typename MatrixAdapter<Matrix>::global_size_t> colptr,
				     typename MatrixAdapter<Matrix>::global_size_t& nnz,
				     Util::EDistribution distribution,
				     Util::EStorage_Ordering ordering,
				     no_special_impl nsi) const
  {
    do_getCcs(nzval, rowind, colptr,
	      nnz, distribution, ordering,
	      typename adapter_t::major_access());
  }

  template < class Matrix >
  void
  MatrixAdapter<Matrix>::do_getCcs(const Teuchos::ArrayView<scalar_t> nzval,
				   const Teuchos::ArrayView<global_ordinal_t> rowind,
				   const Teuchos::ArrayView<typename MatrixAdapter<Matrix>::global_size_t> colptr,
				   typename MatrixAdapter<Matrix>::global_size_t& nnz,
				   Util::EDistribution distribution,
				   Util::EStorage_Ordering ordering,
				   row_access ra) const
  {
    using Teuchos::Array;
    // get the crs and transpose

    Array<scalar_t> nzval_tmp(nzval.size(), 0);
    Array<global_ordinal_t> colind(rowind.size(), 0);
    Array<global_size_t> rowptr(this->getGlobalNumRows() + 1);
    this->getCrs(nzval_tmp(), colind(), rowptr(), nnz, distribution, ordering);
    
    Util::transpose(nzval_tmp(), colind(), rowptr(), nzval, rowind, colptr);
  }

  template < class Matrix >
  void
  MatrixAdapter<Matrix>::do_getCcs(const Teuchos::ArrayView<scalar_t> nzval,
				   const Teuchos::ArrayView<global_ordinal_t> rowind,
				   const Teuchos::ArrayView<typename MatrixAdapter<Matrix>::global_size_t> colptr,
				   typename MatrixAdapter<Matrix>::global_size_t& nnz,
				   Util::EDistribution distribution,
				   Util::EStorage_Ordering ordering,
				   col_access ca) const
  {
    using Teuchos::RCP;
    using Teuchos::ArrayView;
    using Teuchos::OrdinalTraits;
    
    RCP<const type> get_mat = get(distribution);

    Tpetra::Map<scalar_t,local_ordinal_t,global_ordinal_t> rmap = get_mat->getColMap();

    ArrayView<global_ordinal_t> node_elements = rmap->getNodeElementList();
    typename ArrayView<global_ordinal_t>::iterator col_it, col_end;
    col_end = node_elements.end();

    size_t colptr_ind = OrdinalTraits<size_t>::zero();
    global_ordinal_t colInd = OrdinalTraits<global_ordinal_t>::zero();
    for( col_it = node_elements.begin(); col_it != col_end; ++col_it ){
      colptr[colptr_ind++] = colInd;
      size_t colNNZ = getGlobalColNNZ(*col_it);
      size_t nnzRet = 0;
      ArrayView<global_ordinal_t> rowind_view = rowind.view(colInd,colNNZ);
      ArrayView<scalar_t> nzval_view = nzval.view(colInd,colNNZ);
      getGlobalColCopy(*col_it, rowind_view, nzval_view, nnzRet);
      
      // It was suggested that instead of sorting each row's indices
      // individually, that we instead do a double-transpose at the
      // end, which would also lead to the indices being sorted.
      if( ordering == Util::Sorted_Indices ){
	Tpetra::sort2(rowind_view.begin(), rowind_view.end(), nzval_view.begin());
      }
      
      TEST_FOR_EXCEPTION( colNNZ != nnzRet,
			  std::runtime_error,
			  "Number of values returned different from "
                          "number of values reported");
      colInd += colNNZ;
    }
    colptr[colptr_ind] = colInd;
  }

  
  // These will link to concrete implementations
  template < class Matrix >
  void
  MatrixAdapter<Matrix>::getGlobalRowCopy(global_ordinal_t row,
					  const ArrayView<global_ordinal_t>& indices,
					  const ArrayView<scalar_t>& vals,
					  size_t& nnz) const
  {
    static_cast<const adapter_t*>(this)->getGlobalRowCopy_impl(row, indices, vals, nnz);
  }
  
  template < class Matrix >
  void
  MatrixAdapter<Matrix>::getGlobalColCopy(global_ordinal_t col,
					  const ArrayView<global_ordinal_t>& indices,
					  const ArrayView<scalar_t>& vals,
					  size_t& nnz) const
  {
    static_cast<const adapter_t*>(this)->getGlobalColCopy_impl(col, indices, vals, nnz);
  }

  template < class Matrix >
  size_t
  MatrixAdapter<Matrix>::getMaxRowNNZ() const
  {
    return static_cast<const adapter_t*>(this)->getMaxRowNNZ_impl();
  }

  template < class Matrix >
  size_t
  MatrixAdapter<Matrix>::getMaxColNNZ() const
  {
    return static_cast<const adapter_t*>(this)->getMaxColNNZ_impl();
  }
    
  template < class Matrix >
  size_t
  MatrixAdapter<Matrix>::getGlobalRowNNZ(global_ordinal_t row) const
  {
    return static_cast<const adapter_t*>(this)->getGlobalRowNNZ_impl(row);
  }

  template < class Matrix >
  size_t
  MatrixAdapter<Matrix>::getLocalRowNNZ(local_ordinal_t row) const
  {
    return static_cast<const adapter_t*>(this)->getLocalRowNNZ_impl(row);
  }

  template < class Matrix >
  size_t
  MatrixAdapter<Matrix>::getGlobalColNNZ(global_ordinal_t col) const
  {
    return static_cast<const adapter_t*>(this)->getGlobalColNNZ_impl(col);
  }

  template < class Matrix >
  size_t
  MatrixAdapter<Matrix>::getLocalColNNZ(local_ordinal_t col) const
  {
    return static_cast<const adapter_t*>(this)->getLocalColNNZ_impl(col);
  }

  template < class Matrix >
  bool
  MatrixAdapter<Matrix>::isLocallyIndexed() const
  {
    return static_cast<const adapter_t*>(this)->isLocallyIndexed_impl();
  }
  
  template < class Matrix >
  bool
  MatrixAdapter<Matrix>::isGloballyIndexed() const
  {
    return static_cast<const adapter_t*>(this)->isGloballyIndexed_impl();
  }
  
  template < class Matrix >
  RCP<const MatrixAdapter<Matrix> >
  MatrixAdapter<Matrix>::get(Util::EDistribution d) const
  {
    return static_cast<const adapter_t*>(this)->get_impl(d);
  }
  

} // end namespace Amesos

#endif	// AMESOS2_MATRIXADAPTER_DEF_HPP
