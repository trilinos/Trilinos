// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
    const int local_row = this->row_map_->getLocalElement(row);
    bool threw = false;

    Teuchos::Array<local_ordinal_t> epetra_lcl_inds_buf;
    Teuchos::ArrayView<local_ordinal_t> epetra_lcl_inds;
    if (! std::is_same<global_ordinal_t, local_ordinal_t>::value) {
      int num_ent = 0;
      int err = 0;
      try {
        err = this->mat_->NumMyRowEntries (local_row, num_ent);
      }
      catch (int integer_exception) {
        threw = true;
        err = integer_exception;
      }
      TEUCHOS_TEST_FOR_EXCEPTION
        (threw && err != 0, std::runtime_error, "Epetra_RowMatrix::"
         "NumMyRowEntries, called on local row " << local_row << ", threw "
         "an integer exception " << err << ".");
      TEUCHOS_TEST_FOR_EXCEPTION
        (! threw && err != 0, std::runtime_error, "Epetra_RowMatrix returned "
         "error code " << err << " from NumMyRowEntries for local row "
         << local_row << ".");
      epetra_lcl_inds_buf.resize (num_ent);
      epetra_lcl_inds = epetra_lcl_inds_buf ();
    }
    else { // local_ordinal_t == global_ordinal_t
      using Teuchos::av_reinterpret_cast;
      epetra_lcl_inds = av_reinterpret_cast<int> (indices);
    }

    int nnz_ret = 0;
    int rowmatrix_return_val = 0;
    try {
      rowmatrix_return_val =
        this->mat_->ExtractMyRowCopy(local_row,
                                     as<int>(std::min(epetra_lcl_inds.size(), vals.size())),
                                     nnz_ret,
                                     vals.getRawPtr(),
                                     epetra_lcl_inds.getRawPtr());
    }
    catch (int integer_exception) {
      threw = true;
      rowmatrix_return_val = integer_exception;
    }
    TEUCHOS_TEST_FOR_EXCEPTION
      (threw && rowmatrix_return_val != 0, std::runtime_error,
       "Epetra_RowMatrix::ExtractMyRowCopy, called on local row " << local_row
       << ", threw an integer exception " << rowmatrix_return_val << ".");
    TEUCHOS_TEST_FOR_EXCEPTION
      (! threw && rowmatrix_return_val != 0, std::runtime_error,
       "Epetra_RowMatrix object returned error code "
       << rowmatrix_return_val << " from ExtractMyRowCopy." );
    nnz = as<size_t>(nnz_ret);

    // Epetra_CrsMatrix::ExtractMyRowCopy returns local column
    // indices, so transform these into global indices
    for( size_t i = 0; i < nnz; ++i ){
      indices[i] = this->col_map_->getGlobalElement(epetra_lcl_inds[i]);
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
    TEUCHOS_TEST_FOR_EXCEPTION( true,
                        std::runtime_error,
                        "Column access to row-based object not yet supported.  "
                        "Please contact the Amesos2 developers." );
  }


  template <class DerivedMat>
  template<typename KV_GO, typename KV_S>
  void
  AbstractConcreteMatrixAdapter<
    Epetra_RowMatrix,
    DerivedMat>::getGlobalRowCopy_kokkos_view_impl(global_ordinal_t row,
                                                   KV_GO & indices,
                                                   KV_S & vals,
                                                   size_t& nnz) const
  {
    using index_t = typename KV_GO::value_type;
    using value_t = typename KV_S::value_type;
    ArrayView<value_t>  vals_array    (vals.data(),    vals.extent(0));
    ArrayView<index_t>  indices_array (indices.data(), indices.extent(0));

    this->getGlobalRowCopy_impl(row, indices_array, vals_array, nnz);
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
  typename AbstractConcreteMatrixAdapter<
    Epetra_RowMatrix,
    DerivedMat>::global_size_t
  AbstractConcreteMatrixAdapter<
    Epetra_RowMatrix,
    DerivedMat>::getGlobalNumRows_impl() const
  {
    return Teuchos::as<global_size_t>(this->mat_->NumGlobalRows());
  }

  template <class DerivedMat>
  typename AbstractConcreteMatrixAdapter<
    Epetra_RowMatrix,
    DerivedMat>::global_size_t
  AbstractConcreteMatrixAdapter<
    Epetra_RowMatrix,
    DerivedMat>::getGlobalNumCols_impl() const
  {
    return Teuchos::as<global_size_t>(this->mat_->NumGlobalCols());
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
    TEUCHOS_TEST_FOR_EXCEPTION( true,
                        std::runtime_error,
                        "Column access to row-based object not yet supported.  "
                        "Please contact the Amesos2 developers." );
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
    TEUCHOS_TEST_FOR_EXCEPTION( !rowmap.MyGID(gid),
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
    TEUCHOS_TEST_FOR_EXCEPTION( !rowmap.MyLID(lid),
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
    TEUCHOS_TEST_FOR_EXCEPTION( true,
                        std::runtime_error,
                        "Column access to row-based object not yet supported.  "
                        "Please contact the Amesos2 developers." );
  }

  template <class DerivedMat>
  size_t
  AbstractConcreteMatrixAdapter<
    Epetra_RowMatrix,
    DerivedMat>::getLocalColNNZ_impl(local_ordinal_t col) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION( true,
                        std::runtime_error,
                        "Column access to row-based object not yet supported.  "
                        "Please contact the Amesos2 developers." );
  }

  template <class DerivedMat>
  const RCP<const Tpetra::Map<MatrixTraits<Epetra_RowMatrix>::local_ordinal_t,
                              MatrixTraits<Epetra_RowMatrix>::global_ordinal_t,
                              MatrixTraits<Epetra_RowMatrix>::node_t> >
  AbstractConcreteMatrixAdapter<Epetra_RowMatrix, DerivedMat>::getMap_impl() const
  {
    // Should Map() be used (returns Epetra_BlockMap)
    /*
    printf("Amesos2_EpetraRowMatrix_AbstractMatrixAdapter: Epetra does not support a 'getMap()' method. Returning rcp(Teuchos::null). \
        If your map contains non-contiguous GIDs please use Tpetra instead of Epetra. \n");
    */
    return( Teuchos::null );
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
  AbstractConcreteMatrixAdapter<Epetra_RowMatrix, DerivedMat>::get_impl(const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > map, EDistribution distribution) const
  {
    // Delegate implementation to subclass
#ifdef __CUDACC__
    // NVCC doesn't seem to like the static_cast, even though it is valid
    return dynamic_cast<ConcreteMatrixAdapter<DerivedMat>*>(this)->get_impl(map, distribution);
#else
    return static_cast<ConcreteMatrixAdapter<DerivedMat>*>(this)->get_impl(map, distribution);
#endif
  }

  template <class DerivedMat>
  typename AbstractConcreteMatrixAdapter<Epetra_RowMatrix,DerivedMat>
  ::spmtx_ptr_t
  AbstractConcreteMatrixAdapter<Epetra_RowMatrix, DerivedMat>::getSparseRowPtr() const
  {
    typename AbstractConcreteMatrixAdapter<Epetra_RowMatrix,DerivedMat>::spmtx_ptr_t  sp_rowptr = nullptr;
    typename AbstractConcreteMatrixAdapter<Epetra_RowMatrix,DerivedMat>::spmtx_idx_t  sp_colind = nullptr;
    typename AbstractConcreteMatrixAdapter<Epetra_RowMatrix,DerivedMat>::spmtx_val_t  sp_values = nullptr;

    this->mat_->ExtractCrsDataPointers(sp_rowptr, sp_colind, sp_values);

    return sp_rowptr;
  }

  template <class DerivedMat>
  typename AbstractConcreteMatrixAdapter<Epetra_RowMatrix,DerivedMat>
  ::spmtx_idx_t
  AbstractConcreteMatrixAdapter<Epetra_RowMatrix, DerivedMat>::getSparseColInd() const
  {
    typename AbstractConcreteMatrixAdapter<Epetra_RowMatrix,DerivedMat>::spmtx_ptr_t  sp_rowptr = nullptr;
    typename AbstractConcreteMatrixAdapter<Epetra_RowMatrix,DerivedMat>::spmtx_idx_t  sp_colind = nullptr;
    typename AbstractConcreteMatrixAdapter<Epetra_RowMatrix,DerivedMat>::spmtx_val_t  sp_values = nullptr;

    this->mat_->ExtractCrsDataPointers(sp_rowptr, sp_colind, sp_values);

    return sp_colind;
  }

  template <class DerivedMat>
  typename AbstractConcreteMatrixAdapter<Epetra_RowMatrix,DerivedMat>
  ::spmtx_val_t
  AbstractConcreteMatrixAdapter<Epetra_RowMatrix, DerivedMat>::getSparseValues() const
  {
    typename AbstractConcreteMatrixAdapter<Epetra_RowMatrix,DerivedMat>::spmtx_ptr_t  sp_rowptr = nullptr;
    typename AbstractConcreteMatrixAdapter<Epetra_RowMatrix,DerivedMat>::spmtx_idx_t  sp_colind = nullptr;
    typename AbstractConcreteMatrixAdapter<Epetra_RowMatrix,DerivedMat>::spmtx_val_t  sp_values = nullptr;

    this->mat_->ExtractCrsDataPointers(sp_rowptr, sp_colind, sp_values);

    return sp_values;
  }

} // end namespace Amesos2

#endif  // AMESOS2_EPETRAROWMATRIX_ABSTRACTMATRIXADAPTER_DEF_HPP
