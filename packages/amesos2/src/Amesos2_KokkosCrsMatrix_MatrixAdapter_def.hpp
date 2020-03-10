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
    >::get_impl(const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > map, EDistribution distribution) const
  {
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
    KokkosSparse::CrsMatrix<Scalar,LocalOrdinal,ExecutionSpace>>::getGlobalRowNNZ_impl(global_ordinal_t row) const
  {
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
  void
  ConcreteMatrixAdapter<
    KokkosSparse::CrsMatrix<Scalar,LocalOrdinal,ExecutionSpace>>::getGlobalRowCopy_impl(global_ordinal_t row,
                                       const ArrayView<global_ordinal_t>& indices,
                                       const ArrayView<scalar_t>& vals,
                                       size_t& nnz) const
    {
      TEUCHOS_TEST_FOR_EXCEPTION( true,
                        std::runtime_error,
                        "getGlobalRowCopy_impl not implemented for Kokkos CrsMatrix yet.  "
                        "Please contact the Amesos2 developers." );
    }

  template <typename Scalar, typename LocalOrdinal, typename ExecutionSpace>
  void
  ConcreteMatrixAdapter<
    KokkosSparse::CrsMatrix<Scalar,LocalOrdinal,ExecutionSpace>>::getGlobalColCopy_impl(global_ordinal_t col,
                             const ArrayView<global_ordinal_t>& indices,
                             const ArrayView<scalar_t>& vals,
                             size_t& nnz) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION( true,
                        std::runtime_error,
                        "Column access to row-based object not yet supported.  "
                        "Please contact the Amesos2 developers." );
  }

  template <typename Scalar, typename LocalOrdinal, typename ExecutionSpace>
  typename ConcreteMatrixAdapter<KokkosSparse::CrsMatrix<Scalar,LocalOrdinal,ExecutionSpace>>::spmtx_ptr_t
  ConcreteMatrixAdapter<KokkosSparse::CrsMatrix<Scalar,LocalOrdinal,ExecutionSpace>>::getSparseRowPtr() const
  {
    return this->mat_->graph.row_map.data();
  }

  template <typename Scalar, typename LocalOrdinal, typename ExecutionSpace>
  typename ConcreteMatrixAdapter<KokkosSparse::CrsMatrix<Scalar,LocalOrdinal,ExecutionSpace>>::spmtx_idx_t
  ConcreteMatrixAdapter<KokkosSparse::CrsMatrix<Scalar,LocalOrdinal,ExecutionSpace>>::getSparseColInd() const
  {
    return this->mat_->graph.entries.data();
  }

  template <typename Scalar, typename LocalOrdinal, typename ExecutionSpace>
  typename ConcreteMatrixAdapter<KokkosSparse::CrsMatrix<Scalar,LocalOrdinal,ExecutionSpace>>::spmtx_vals_t
  ConcreteMatrixAdapter<KokkosSparse::CrsMatrix<Scalar,LocalOrdinal,ExecutionSpace>>::getSparseValues() const
  {
    return this->mat_->values.data();
  }

} // end namespace Amesos2

#endif  // AMESOS2_KOKKOS_CRSMATRIX_MATRIXADAPTER_DEF_HPP
