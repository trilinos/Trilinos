// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifndef AMESOS2_MATRIXADAPTER_DEF_HPP
#define AMESOS2_MATRIXADAPTER_DEF_HPP
#include <Tpetra_Util.hpp>
#include "Amesos2_MatrixAdapter_decl.hpp"
#include "Amesos2_ConcreteMatrixAdapter_def.hpp"
//#include "Amesos2_ConcreteMatrixAdapter.hpp"

#define TESTING_AMESOS2_WITH_TPETRA_REMOVE_UVM
#if defined(TESTING_AMESOS2_WITH_TPETRA_REMOVE_UVM)
#include "KokkosSparse_Utils.hpp"
#include "KokkosSparse_SortCrs.hpp"
#include "KokkosKernels_Sorting.hpp"
#endif

namespace Amesos2 {

  
  template < class Matrix >
  MatrixAdapter<Matrix>::MatrixAdapter(const Teuchos::RCP<Matrix> m)
    : mat_(m)
  {
    comm_ = static_cast<const adapter_t*>(this)->getComm_impl();
    col_map_ = static_cast<const adapter_t*>(this)->getColMap_impl();
    row_map_ = static_cast<const adapter_t*>(this)->getRowMap_impl();
  }

  template < class Matrix >
  template<typename KV_S, typename KV_GO, typename KV_GS>
  void
  MatrixAdapter<Matrix>::getCrs_kokkos_view(KV_S & nzval,
        KV_GO & colind,
        KV_GS & rowptr,
        typename MatrixAdapter<Matrix>::global_size_t& nnz,
        const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t, global_ordinal_t, node_t> > rowmap,
        EStorage_Ordering ordering,
        EDistribution distribution) const
  {
    help_getCrs_kokkos_view(nzval, colind, rowptr,
      nnz, rowmap, distribution, ordering,
    typename adapter_t::get_crs_spec());
  }

  template < class Matrix >
  template<typename KV_S, typename KV_GO, typename KV_GS>
  void
  MatrixAdapter<Matrix>::getCrs_kokkos_view(KV_S & nzval,
        KV_GO & colind,
        KV_GS & rowptr,
        typename MatrixAdapter<Matrix>::global_size_t& nnz,
        EDistribution distribution,
        EStorage_Ordering ordering) const
  {
    const Teuchos::RCP<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > rowmap
      = Util::getDistributionMap<local_ordinal_t,global_ordinal_t,global_size_t,node_t>(distribution,
                                                                                        this->getGlobalNumRows(),
                                                                                        this->getComm());
    this->getCrs_kokkos_view(nzval, colind, rowptr, nnz, Teuchos::ptrInArg(*rowmap), ordering, distribution);
  }

  template < class Matrix >
  template<typename KV_S, typename KV_GO, typename KV_GS>
  void
  MatrixAdapter<Matrix>::getCcs_kokkos_view(KV_S & nzval,
        KV_GO & rowind,
        KV_GS & colptr,
        typename MatrixAdapter<Matrix>::global_size_t& nnz,
        const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t, global_ordinal_t, node_t> > colmap,
        EStorage_Ordering ordering,
        EDistribution distribution) const
  {
    help_getCcs_kokkos_view(nzval, rowind, colptr,
      nnz, colmap, distribution, ordering,
    typename adapter_t::get_ccs_spec());
  }

  template < class Matrix >
  template<typename KV_S, typename KV_GO, typename KV_GS>
  void
  MatrixAdapter<Matrix>::getCcs_kokkos_view(KV_S & nzval,
        KV_GO & rowind,
        KV_GS & colptr,
        typename MatrixAdapter<Matrix>::global_size_t& nnz,
        EDistribution distribution,
        EStorage_Ordering ordering) const
  {
    const Teuchos::RCP<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > colmap
      = Util::getDistributionMap<local_ordinal_t,global_ordinal_t,global_size_t,node_t>(distribution,
                                                                                        this->getGlobalNumCols(),
                                                                                        this->getComm());
    this->getCcs_kokkos_view(nzval, rowind, colptr, nnz, Teuchos::ptrInArg(*colmap), ordering, distribution);
  }


  template < class Matrix >
  typename MatrixAdapter<Matrix>::global_size_t
  MatrixAdapter<Matrix>::getGlobalNumRows() const
  {
    return static_cast<const adapter_t*>(this)->getGlobalNumRows_impl();
  }

  template < class Matrix >
  typename MatrixAdapter<Matrix>::global_size_t
  MatrixAdapter<Matrix>::getGlobalNumCols() const
  {
    return static_cast<const adapter_t*>(this)->getGlobalNumCols_impl();
  }
  
  template < class Matrix >
  typename MatrixAdapter<Matrix>::global_size_t
  MatrixAdapter<Matrix>::getRowIndexBase() const
  {
    // Kokkos adapter is for serial only testing right now and will not
    // create row_map_
    return (row_map_ != Teuchos::null) ? row_map_->getIndexBase() : 0;
  }

  template < class Matrix >
  typename MatrixAdapter<Matrix>::global_size_t
  MatrixAdapter<Matrix>::getColumnIndexBase() const
  {
    // Kokkos adapter is for serial only testing right now and will not
    // create col_map_
    return (col_map_ != Teuchos::null) ? col_map_->getIndexBase() : 0;
  }

  template < class Matrix >
  typename MatrixAdapter<Matrix>::global_size_t
  MatrixAdapter<Matrix>::getGlobalNNZ() const
  {
    return static_cast<const adapter_t*>(this)->getGlobalNNZ_impl();
  }
  
  template < class Matrix >
  size_t
  MatrixAdapter<Matrix>::getLocalNumRows() const
  {
    return row_map_->getLocalNumElements(); // TODO: verify
  }

  template < class Matrix >
  size_t
  MatrixAdapter<Matrix>::getLocalNumCols() const
  {
    return col_map_->getLocalNumElements();
  }
  
  template < class Matrix >
  size_t
  MatrixAdapter<Matrix>::getLocalNNZ() const
  {
    return static_cast<const adapter_t*>(this)->getLocalNNZ_impl();
  }

  // NDE: This is broken for Epetra_CrsMatrix
  template < class Matrix >
  std::string
  MatrixAdapter<Matrix>::description() const
  {
    std::ostringstream oss;
    oss << "Amesos2::MatrixAdapter wrapping: ";
    oss << mat_->description(); // NDE: This is not defined in Epetra_CrsMatrix, only in Tpetra::CrsMatrix
    return oss.str();
  }
  
  template < class Matrix >
  void
  MatrixAdapter<Matrix>::describe(Teuchos::FancyOStream &out,
                                  const Teuchos::EVerbosityLevel verbLevel) const
  {}

  template < class Matrix >
  template < class KV >
  void MatrixAdapter<Matrix>::returnRowPtr_kokkos_view(KV & view) const
  {
    return static_cast<const adapter_t*>(this)->getSparseRowPtr_kokkos_view(view);
  }

  template < class Matrix >
  template < class KV >
  void MatrixAdapter<Matrix>::returnColInd_kokkos_view(KV & view) const
  {
    return static_cast<const adapter_t*>(this)->getSparseColInd_kokkos_view(view);
  }

  template < class Matrix >
  template < class KV >
  void MatrixAdapter<Matrix>::returnValues_kokkos_view(KV & view) const
  {
    return static_cast<const adapter_t*>(this)->getSparseValues_kokkos_view(view);
  }

  /******************************
   * Private method definitions *
   ******************************/
  template < class Matrix >
  template<typename KV_S, typename KV_GO, typename KV_GS>
  void
  MatrixAdapter<Matrix>::help_getCrs_kokkos_view(KV_S & nzval,
             KV_GO & colind,
             KV_GS & rowptr,
             typename MatrixAdapter<Matrix>::global_size_t& nnz,
             const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > rowmap,
             EDistribution distribution,
             EStorage_Ordering ordering,
             no_special_impl nsi) const
  {

    //Added void to remove parameter not used warning
    ((void)nsi);
    do_getCrs_kokkos_view(nzval, colind, rowptr,
      nnz, rowmap, distribution, ordering,
      typename adapter_t::major_access());
  }

  template < class Matrix >
  template<typename KV_S, typename KV_GO, typename KV_GS>
  void
  MatrixAdapter<Matrix>::do_getCrs_kokkos_view(KV_S & nzval,
           KV_GO & colind,
           KV_GS & rowptr,
           typename MatrixAdapter<Matrix>::global_size_t& nnz,
           const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > rowmap,
           EDistribution distribution,
           EStorage_Ordering ordering,
           row_access ra) const
  {
    // Kokkos adapter will be serial and won't have the rowmap.
    // Tacho for example wouldn't ever call this in serial but Cholmod will
    // call ccs and want to convert using this.
    // If the kokkos adapter is extended to multiple ranks then this will
    // need to change.
    if(this->row_map_ == Teuchos::null) {
      this->returnValues_kokkos_view(nzval);
      this->returnRowPtr_kokkos_view(rowptr);
      this->returnColInd_kokkos_view(colind);
      nnz = nzval.size();
      return;
    }

    using Teuchos::rcp;
    using Teuchos::RCP;
    using Teuchos::ArrayView;
    using Teuchos::OrdinalTraits;

    ((void) ra);

    RCP<const type> get_mat;
    if( *rowmap == *this->row_map_ && distribution != CONTIGUOUS_AND_ROOTED ){
      // No need to redistribute
      get_mat = rcp(this,false); // non-owning
    } else {
      get_mat = get(rowmap, distribution);
    }
    // RCP<const type> get_mat = get(rowmap);

    // rmap may not necessarily check same as rowmap because rmap may
    // have been constructued with Tpetra's "expert" constructor,
    // which assumes that the map points are non-contiguous.
    //
    // TODO: There may be some more checking between the row map
    // compatibility, but things are working fine now.

    RCP<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > rmap = get_mat->getRowMap();
    ArrayView<const global_ordinal_t> node_elements = rmap->getLocalElementList();
    //if( node_elements.size() == 0 ) return; // no more contribution
    typename ArrayView<const global_ordinal_t>::iterator row_it, row_end;
    row_end = node_elements.end();

    size_t rowptr_ind = OrdinalTraits<size_t>::zero();
    global_ordinal_t rowInd = OrdinalTraits<global_ordinal_t>::zero();

    // For rowptr we can just make a mirror and deep_copy at the end
    typename KV_GS::HostMirror host_rowptr = Kokkos::create_mirror_view(rowptr);

    #if !defined(TESTING_AMESOS2_WITH_TPETRA_REMOVE_UVM)
    // Note nzval, colind, and rowptr will not all be in the same memory space.
    // Currently only Cholmod exercises this code which has all the arrays on host,
    // so this will need extension and testing when we have a solver using device here.
    Kokkos::View<scalar_t*, Kokkos::HostSpace>
      mat_nzval(Kokkos::ViewAllocateWithoutInitializing("mat_nzval"), nzval.size());

    Kokkos::View<global_ordinal_t*, Kokkos::HostSpace>
      mat_colind(Kokkos::ViewAllocateWithoutInitializing("mat_colind"), colind.size());

    ArrayView<scalar_t> nzval_arrayview(mat_nzval.data(), nzval.size());
    ArrayView<global_ordinal_t> colind_arrayview(mat_colind.data(), colind.size());

    for( row_it = node_elements.begin(); row_it != row_end; ++row_it ){
      host_rowptr(rowptr_ind++) = rowInd;
      size_t rowNNZ = get_mat->getGlobalRowNNZ(*row_it);
      size_t nnzRet = OrdinalTraits<size_t>::zero();
      ArrayView<global_ordinal_t> colind_view = colind_arrayview.view(rowInd,rowNNZ);
      ArrayView<scalar_t> nzval_view = nzval_arrayview.view(rowInd,rowNNZ);

      get_mat->getGlobalRowCopy(*row_it, colind_view, nzval_view, nnzRet);

      for (size_t rr = 0; rr < nnzRet ; rr++) {
        colind_view[rr] -= rmap->getIndexBase();
      }

      // It was suggested that instead of sorting each row's indices
      // individually, that we instead do a double-transpose at the
      // end, which would also lead to the indices being sorted.
      if( ordering == SORTED_INDICES ) {
        Tpetra::sort2(colind_view.begin(), colind_view.end(), nzval_view.begin());
      }

      TEUCHOS_TEST_FOR_EXCEPTION( rowNNZ != nnzRet,
        std::runtime_error,
        "Number of values returned different from "
                          "number of values reported");
      rowInd += rowNNZ;
    }
    host_rowptr(rowptr_ind) = nnz = rowInd;

    deep_copy_or_assign_view(nzval, mat_nzval);
    deep_copy_or_assign_view(colind, mat_colind);
    deep_copy_or_assign_view(rowptr, host_rowptr);
    #else
    // create temporary views to hold colind and nvals (TODO: allocate as much as needed, also for rowptr)
    global_host_idx_t mat_colind(Kokkos::ViewAllocateWithoutInitializing("mat_colind"), nzval.size());
    global_host_val_t mat_nzvals(Kokkos::ViewAllocateWithoutInitializing("mat_nzvals"), colind.size());

    auto host_colind = Kokkos::create_mirror_view(colind);
    auto host_nzval = Kokkos::create_mirror_view(nzval);

    // load crs (on host)
    for( row_it = node_elements.begin(); row_it != row_end; ++row_it ){
      size_t rowNNZ = get_mat->getGlobalRowNNZ(*row_it);
      size_t nnzRet = OrdinalTraits<size_t>::zero();
      //using range_type = Kokkos::pair<int, int>;
      //auto colind_view = Kokkos::subview(mat_colind, range_type(rowInd, rowInd+rowNNZ));
      //auto nzval_view = Kokkos::subview(mat_nzvals, range_type(rowInd, rowInd+rowNNZ));
      global_host_idx_t colind_view (&(mat_colind(rowInd)), rowNNZ);
      global_host_val_t nzvals_view (&(mat_nzvals(rowInd)), rowNNZ);

      global_ordinal_t row_id = *row_it;
      get_mat->getGlobalRowCopy_kokkos_view(row_id, colind_view, nzvals_view, nnzRet);

      TEUCHOS_TEST_FOR_EXCEPTION( rowNNZ != nnzRet,
        std::runtime_error,
        "Number of values returned different from "
                          "number of values reported");
      host_rowptr(rowptr_ind++) = rowInd;
      rowInd += rowNNZ;
    }
    host_rowptr(rowptr_ind) = nnz = rowInd;

    // fix index-base
    if (rmap->getIndexBase() != 0) {
      for (size_t k = 0; k < mat_colind.extent(0); k++) {
        mat_colind(k) -= rmap->getIndexBase();
      }
    }

    // copy to device (note: everything in the vectors are copied, though they may not be used)
    deep_copy_or_assign_view(nzval,  mat_nzvals);
    deep_copy_or_assign_view(colind, mat_colind);
    deep_copy_or_assign_view(rowptr, host_rowptr);

    // sort
    if( ordering == SORTED_INDICES ) {
      using execution_space = typename KV_GS::execution_space;
      KokkosSparse::sort_crs_matrix <execution_space, KV_GS, KV_GO, KV_S>
        (rowptr, colind, nzval);
    }
    #endif
  }

  template < class Matrix >
  template<typename KV_S, typename KV_GO, typename KV_GS>
  void
  MatrixAdapter<Matrix>::help_getCcs_kokkos_view(KV_S & nzval,
             KV_GO & rowind,
             KV_GS & colptr,
             typename MatrixAdapter<Matrix>::global_size_t& nnz,
             const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > colmap,
             EDistribution distribution,
             EStorage_Ordering ordering,
             no_special_impl nsi) const
  {

    //Added void to remove parameter not used warning
    ((void)nsi);
    do_getCcs_kokkos_view(nzval, rowind, colptr,
      nnz, colmap, distribution, ordering,
      typename adapter_t::major_access());
  }

  template < class Matrix >
  template<typename KV_S, typename KV_GO, typename KV_GS>
  void
  MatrixAdapter<Matrix>::do_getCcs_kokkos_view(KV_S & nzval,
           KV_GO & rowind,
           KV_GS & colptr,
           typename MatrixAdapter<Matrix>::global_size_t& nnz,
           const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > colmap,
           EDistribution distribution,
           EStorage_Ordering ordering,
           row_access ra) const
  {
    using Teuchos::ArrayView;
    // get the crs and transpose

    ((void) ra);

    KV_S nzval_tmp(Kokkos::ViewAllocateWithoutInitializing("nzval_tmp"), nzval.size());
    KV_GO colind(Kokkos::ViewAllocateWithoutInitializing("colind"), rowind.size());
    KV_GS rowptr(Kokkos::ViewAllocateWithoutInitializing("rowptr"), this->getGlobalNumRows() + 1);

    this->getCrs_kokkos_view(nzval_tmp, colind, rowptr, nnz, colmap, ordering, distribution);

    if(nnz > 0) {
      // This is currently just used by Cholmod in which case the views will be
      // host, even if Cholmod is using GPU. Will need to upgrade this section
      // to properly handle device when we have a solver that needs it.
      ArrayView<typename KV_S::value_type> av_nzval_tmp(nzval_tmp.data(), nzval_tmp.size());
      ArrayView<typename KV_GO::value_type> av_colind(colind.data(), colind.size());
      ArrayView<typename KV_GS::value_type> av_rowptr(rowptr.data(), rowptr.size());
      ArrayView<typename KV_S::value_type> av_nzval(nzval.data(), nzval.size());
      ArrayView<typename KV_GO::value_type> av_rowind(rowind.data(), rowind.size());
      ArrayView<typename KV_GS::value_type> av_colptr(colptr.data(), colptr.size());
      Util::transpose(av_nzval_tmp, av_colind, av_rowptr, av_nzval, av_rowind, av_colptr);
    }
  }
  
  // These will link to concrete implementations
  template < class Matrix >
  template<typename KV_GO, typename KV_S>
  void
  MatrixAdapter<Matrix>::getGlobalRowCopy_kokkos_view(global_ordinal_t row,
                                                      KV_GO & indices,
                                                      KV_S & vals,
                                                      size_t& nnz) const
  {
    static_cast<const adapter_t*>(this)->getGlobalRowCopy_kokkos_view_impl(row, indices, vals, nnz);
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
  Teuchos::RCP<const MatrixAdapter<Matrix> >
  MatrixAdapter<Matrix>::get(const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > map, EDistribution distribution) const
  {
    return static_cast<const adapter_t*>(this)->get_impl(map, distribution);
  }


  template <class Matrix>
  Teuchos::RCP<MatrixAdapter<Matrix> >
  createMatrixAdapter(Teuchos::RCP<Matrix> m){
    using Teuchos::rcp;
    using Teuchos::rcp_const_cast;
    
    if(m.is_null()) return Teuchos::null;
    return( rcp(new ConcreteMatrixAdapter<Matrix>(m)) );
  }

  template <class Matrix>
  Teuchos::RCP<const MatrixAdapter<Matrix> >
  createConstMatrixAdapter(Teuchos::RCP<const Matrix> m){
    using Teuchos::rcp;
    using Teuchos::rcp_const_cast;
    
    if(m.is_null()) return Teuchos::null;
    return( rcp(new ConcreteMatrixAdapter<Matrix>(rcp_const_cast<Matrix,const Matrix>(m))).getConst() );
  }

} // end namespace Amesos2

#endif        // AMESOS2_MATRIXADAPTER_DEF_HPP
