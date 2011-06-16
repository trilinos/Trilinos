/**
 * \file   Amesos2_EpetraRowMatrix_AbstractMatrixAdapter_def.hpp
 * \author Eric Bavier <etbavie@sandia.gov>
 * \date   Tue Jun 14 17:21:32 2011
 *
 * \brief  Definitions for the Epetra_RowMatrix abstract adapter.
 */

#ifndef AMESOS2_EPETRAROWMATRIX_ABSTRACTMATRIXADAPTER_DEF_HPP
#define AMESOS2_EPETRAROWMATRIX_ABSTRACTMATRIXADAPTER_DEF_HPP

#include <Epetra_Map.h>

namespace Amesos {

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
    
    int nnz_ret = as<int>(nnz);
    int rowmatrix_return_val
      = this->mat_->ExtractMyRowCopy(as<int>(row),
				     as<int>(std::min(indices.size(), vals.size())),
				     nnz_ret,
				     vals.getRawPtr(),
				     indices.getRawPtr());
    TEST_FOR_EXCEPTION( rowmatrix_return_val != 0,
			std::runtime_error,
			"Epetra_RowMatrix object returned error code "
			<< rowmatrix_return_val << " from ExtractMyRowCopy." );
    nnz = as<size_t>(nnz_ret);
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
    // not used for right now
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
    // not used
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
    // not used
  }

  template <class DerivedMat>
  size_t
  AbstractConcreteMatrixAdapter<
    Epetra_RowMatrix,
    DerivedMat>::getLocalColNNZ_impl(local_ordinal_t col) const
  {
    // not used
  }

  template <class DerivedMat>
  const RCP<const Tpetra::Map<MatrixTraits<Epetra_RowMatrix>::local_ordinal_t,
			      MatrixTraits<Epetra_RowMatrix>::global_ordinal_t,
			      MatrixTraits<Epetra_RowMatrix>::node_t> >
  AbstractConcreteMatrixAdapter<Epetra_RowMatrix, DerivedMat>::getRowMap_impl() const
  {
    // Must transform to a Tpetra::Map
    Epetra_Map rowmap = this->mat_->RowMatrixRowMap();

    int num_my_elements = rowmap.NumMyElements();
    Teuchos::Array<int> my_global_elements(num_my_elements);
    rowmap.MyGlobalElements(my_global_elements.getRawPtr());

    using Teuchos::as;
    typedef Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> map_t;
    RCP<map_t> tmap = rcp(new map_t(Teuchos::OrdinalTraits<global_size_t>::invalid(),
				    my_global_elements(),
				    as<global_ordinal_t>(rowmap.IndexBase()),
				    this->getComm()));
    return tmap;
  }

  template <class DerivedMat>
  const RCP<const Tpetra::Map<MatrixTraits<Epetra_RowMatrix>::local_ordinal_t,
			      MatrixTraits<Epetra_RowMatrix>::global_ordinal_t,
			      MatrixTraits<Epetra_RowMatrix>::node_t> >
  AbstractConcreteMatrixAdapter<Epetra_RowMatrix, DerivedMat>::getColMap_impl() const
  {
    // Must transform this matrix' Epetra_Map to a Tpetra::Map
    Epetra_Map colmap = this->mat_->RowMatrixColMap();

    int num_my_elements = colmap.NumMyElements();
    Teuchos::Array<int> my_global_elements(num_my_elements);
    colmap.MyGlobalElements(my_global_elements.getRawPtr());

    using Teuchos::as;
    typedef Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> map_t;
    RCP<map_t> tmap = rcp(new map_t(Teuchos::OrdinalTraits<global_size_t>::invalid(),
				    my_global_elements(),
				    as<global_ordinal_t>(colmap.IndexBase()),
				    this->getComm()));
    return tmap;
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
  AbstractConcreteMatrixAdapter<Epetra_RowMatrix, DerivedMat>::get_impl(EDistribution d) const
  {
    // Delegate implementation to subclass
    return static_cast<ConcreteMatrixAdapter<DerivedMat>*>(this)->get_impl(d);
  }

} // end namespace Amesos

#endif  // AMESOS2_EPETRAROWMATRIX_ABSTRACTMATRIXADAPTER_DEF_HPP
