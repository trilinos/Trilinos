// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifndef AMESOS2_KOKKOS_CRSMATRIX_MATRIXADAPTER_DEF_HPP
#define AMESOS2_KOKKOS_CRSMATRIX_MATRIXADAPTER_DEF_HPP

#include "Amesos2_KokkosCrsMatrix_MatrixAdapter_decl.hpp"
#include "Amesos2_MatrixAdapter_def.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include <Tpetra_Core.hpp>

namespace Amesos2 {

  template <typename Scalar, typename LocalOrdinal, typename ExecutionSpace>
  ConcreteMatrixAdapter<KokkosSparse::CrsMatrix<Scalar,LocalOrdinal,ExecutionSpace>>::ConcreteMatrixAdapter(Teuchos::RCP<matrix_t> m)
    : MatrixAdapter<KokkosSparse::CrsMatrix<Scalar,LocalOrdinal,ExecutionSpace>>(m)
  {

  }

  template <typename Scalar, typename LocalOrdinal, typename ExecutionSpace>
  typename ConcreteMatrixAdapter<KokkosSparse::CrsMatrix<Scalar,LocalOrdinal,ExecutionSpace>>::global_size_t
  ConcreteMatrixAdapter<KokkosSparse::CrsMatrix<Scalar,LocalOrdinal,ExecutionSpace>>::getRowIndexBase() const
  {
    return 0;
  }

  template <typename Scalar, typename LocalOrdinal, typename ExecutionSpace>
  typename ConcreteMatrixAdapter<KokkosSparse::CrsMatrix<Scalar,LocalOrdinal,ExecutionSpace>>::global_size_t
  ConcreteMatrixAdapter<KokkosSparse::CrsMatrix<Scalar,LocalOrdinal,ExecutionSpace>>::getColumnIndexBase() const
  {
    return 0;
  }

  template <typename Scalar, typename LocalOrdinal, typename ExecutionSpace>
  const Teuchos::RCP<const Teuchos::Comm<int> >
  ConcreteMatrixAdapter<KokkosSparse::CrsMatrix<Scalar,LocalOrdinal,ExecutionSpace>>::getComm_impl() const
  {
    return Tpetra::getDefaultComm(); // Kokkos CrsMatrix currently is just serial
  }

  template <typename Scalar, typename LocalOrdinal, typename ExecutionSpace>
  const RCP<const Tpetra::Map<typename MatrixTraits<KokkosSparse::CrsMatrix<Scalar,LocalOrdinal,ExecutionSpace>>::local_ordinal_t,
                              typename MatrixTraits<KokkosSparse::CrsMatrix<Scalar,LocalOrdinal,ExecutionSpace>>::global_ordinal_t,
                              typename MatrixTraits<KokkosSparse::CrsMatrix<Scalar,LocalOrdinal,ExecutionSpace>>::node_t> >
  ConcreteMatrixAdapter<KokkosSparse::CrsMatrix<Scalar,LocalOrdinal,ExecutionSpace>>::getRowMap_impl() const
  {
    return Teuchos::null; // not going to use this right now - serial
  }

  template <typename Scalar, typename LocalOrdinal, typename ExecutionSpace>
  const RCP<const Tpetra::Map<typename MatrixTraits<KokkosSparse::CrsMatrix<Scalar,LocalOrdinal,ExecutionSpace>>::local_ordinal_t,
                              typename MatrixTraits<KokkosSparse::CrsMatrix<Scalar,LocalOrdinal,ExecutionSpace>>::global_ordinal_t,
                              typename MatrixTraits<KokkosSparse::CrsMatrix<Scalar,LocalOrdinal,ExecutionSpace>>::node_t> >
  ConcreteMatrixAdapter<KokkosSparse::CrsMatrix<Scalar,LocalOrdinal,ExecutionSpace>>::getColMap_impl() const
  {
    return Teuchos::null; // not going to use this right now - serial
  }

  template <typename Scalar, typename LocalOrdinal, typename ExecutionSpace>
  Teuchos::RCP<const MatrixAdapter<KokkosSparse::CrsMatrix<Scalar,LocalOrdinal,ExecutionSpace> > >
  ConcreteMatrixAdapter<
    KokkosSparse::CrsMatrix<Scalar,LocalOrdinal,ExecutionSpace>
  >::get_impl(
    [[maybe_unused]] const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > map,
    [[maybe_unused]] EDistribution distribution
  ) const {
    TEUCHOS_TEST_FOR_EXCEPTION( true,
                        std::runtime_error,
                        "get_impl() not implemented for the Kokkos CrsMatrix adapter yet.  "
                        "Please contact the Amesos2 developers." );
  }

  template <typename Scalar, typename LocalOrdinal, typename ExecutionSpace>
  typename ConcreteMatrixAdapter<KokkosSparse::CrsMatrix<Scalar,LocalOrdinal,ExecutionSpace>>::global_size_t
  ConcreteMatrixAdapter<
    KokkosSparse::CrsMatrix<Scalar,LocalOrdinal,ExecutionSpace>>::getGlobalNumRows_impl() const
  {
    return this->mat_->numRows();
  }

  template <typename Scalar, typename LocalOrdinal, typename ExecutionSpace>
  typename ConcreteMatrixAdapter<KokkosSparse::CrsMatrix<Scalar,LocalOrdinal,ExecutionSpace>>::global_size_t
  ConcreteMatrixAdapter<
    KokkosSparse::CrsMatrix<Scalar,LocalOrdinal,ExecutionSpace>>::getGlobalNumCols_impl() const
  {
    return this->mat_->numCols();
  }

  template <typename Scalar, typename LocalOrdinal, typename ExecutionSpace>
  typename ConcreteMatrixAdapter<KokkosSparse::CrsMatrix<Scalar,LocalOrdinal,ExecutionSpace>>::global_size_t
  ConcreteMatrixAdapter<
    KokkosSparse::CrsMatrix<Scalar,LocalOrdinal,ExecutionSpace>>::getGlobalNNZ_impl() const
  {
    return this->mat_->nnz();
  }

  template <typename Scalar, typename LocalOrdinal, typename ExecutionSpace>
  const Teuchos::RCP<const Tpetra::Map<typename MatrixTraits<KokkosSparse::CrsMatrix<Scalar,LocalOrdinal,ExecutionSpace>>::local_ordinal_t,
                              typename MatrixTraits<KokkosSparse::CrsMatrix<Scalar,LocalOrdinal,ExecutionSpace>>::global_ordinal_t,
                              typename MatrixTraits<KokkosSparse::CrsMatrix<Scalar,LocalOrdinal,ExecutionSpace>>::node_t> >
  ConcreteMatrixAdapter<
    KokkosSparse::CrsMatrix<Scalar,LocalOrdinal,ExecutionSpace>>::getMap_impl() const
  {
    return( Teuchos::null );
  }

  template <typename Scalar, typename LocalOrdinal, typename ExecutionSpace>
  size_t
  ConcreteMatrixAdapter<
    KokkosSparse::CrsMatrix<Scalar,LocalOrdinal,ExecutionSpace>
  >::getGlobalRowNNZ_impl(
    [[maybe_unused]] global_ordinal_t row
  ) const {
    TEUCHOS_TEST_FOR_EXCEPTION( true,
                        std::runtime_error,
                        "getGlobalRowNNZ_impl() not implemented for the Kokkos CrsMatrix adapter yet.  "
                        "Please contact the Amesos2 developers." );
  }

  template <typename Scalar, typename LocalOrdinal, typename ExecutionSpace>
  size_t
  ConcreteMatrixAdapter<
    KokkosSparse::CrsMatrix<Scalar,LocalOrdinal,ExecutionSpace>>::getLocalRowNNZ_impl(local_ordinal_t row) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION( true,
                        std::runtime_error,
                        "getLocalRowNNZ_impl() not implemented for the Kokkos CrsMatrix adapter yet.  "
                        "Please contact the Amesos2 developers." );
  }

  template <typename Scalar, typename LocalOrdinal, typename ExecutionSpace>
  size_t
  ConcreteMatrixAdapter<
    KokkosSparse::CrsMatrix<Scalar,LocalOrdinal,ExecutionSpace>>::getGlobalColNNZ_impl(global_ordinal_t col) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION( true,
                        std::runtime_error,
                        "Column access to row-based object not yet supported.  "
                        "Please contact the Amesos2 developers." );
  }

  template <typename Scalar, typename LocalOrdinal, typename ExecutionSpace>
  size_t
  ConcreteMatrixAdapter<
    KokkosSparse::CrsMatrix<Scalar,LocalOrdinal,ExecutionSpace>>::getLocalColNNZ_impl(local_ordinal_t col) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION( true,
                        std::runtime_error,
                        "Column access to row-based object not yet supported.  "
                        "Please contact the Amesos2 developers." );
  }

  // implementation functions
  template <typename Scalar, typename LocalOrdinal, typename ExecutionSpace>
  template <typename KV_GO, typename KV_S>
  void
  ConcreteMatrixAdapter<
    KokkosSparse::CrsMatrix<Scalar,LocalOrdinal,ExecutionSpace>
  >::getGlobalRowCopy_kokkos_view_impl(
    [[maybe_unused]] global_ordinal_t row,
    [[maybe_unused]] KV_GO & indices,
    [[maybe_unused]] KV_S & vals,
    [[maybe_unused]] size_t& nnz
  ) const {
    TEUCHOS_TEST_FOR_EXCEPTION( true,
                      std::runtime_error,
                      "getGlobalRowCopy_kokkos_view_impl not implemented for Kokkos CrsMatrix yet.  "
                      "Please contact the Amesos2 developers." );
  }

} // end namespace Amesos2

#endif  // AMESOS2_KOKKOS_CRSMATRIX_MATRIXADAPTER_DEF_HPP
