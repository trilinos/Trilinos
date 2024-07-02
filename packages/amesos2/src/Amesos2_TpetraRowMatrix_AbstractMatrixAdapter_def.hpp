// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifndef AMESOS2_TPETRAROWMATRIX_ABSTRACTMATRIXADAPTER_DEF_HPP
#define AMESOS2_TPETRAROWMATRIX_ABSTRACTMATRIXADAPTER_DEF_HPP

#include "Amesos2_ConcreteMatrixAdapter_decl.hpp"
#include "Amesos2_TpetraRowMatrix_AbstractMatrixAdapter_decl.hpp"

namespace Amesos2 {

  using Teuchos::RCP;

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
  template <typename KV_GO, typename KV_S>
  void
  AbstractConcreteMatrixAdapter<
    Tpetra::RowMatrix<Scalar,
                      LocalOrdinal,
                      GlobalOrdinal,
                      Node>,
    DerivedMat>::getGlobalRowCopy_kokkos_view_impl(global_ordinal_t row,
                                                   KV_GO & indices,
                                                   KV_S & vals,
                                                   size_t& nnz) const
    {
      this->mat_->getGlobalRowCopy(row, indices, vals, nnz);
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
    return this->mat_->getLocalNumEntries();
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
    DerivedMat>:: getMap_impl() const
  {
    return this->mat_->getMap();
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
    >::get_impl(const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > map, EDistribution distribution) const
  {
#ifdef __CUDACC__
    // NVCC doesn't seem to like the static_cast, even though it is valid
    return dynamic_cast<ConcreteMatrixAdapter<DerivedMat>*>(this)->get_impl(map, distribution);
#else
    return static_cast<ConcreteMatrixAdapter<DerivedMat>*>(this)->get_impl(map, distribution);
#endif
  }

} // end namespace Amesos2

#endif  // AMESOS2_TPETRAROWMATRIX_ABSTRACTMATRIXADAPTER_DEF_HPP
