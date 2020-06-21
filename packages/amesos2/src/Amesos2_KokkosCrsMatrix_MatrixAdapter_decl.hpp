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
 * \file   Amesos2_KokkosCrsMatrix_MatrixAdapter_decl.hpp
 * \author
 * \date
 *
 * \brief Specialization of the ConcreteMatrixAdapter for
 * KokkosSparse::CrsMatrix.
 */

#ifndef AMESOS2_KOKKOSCRSMATRIX_MATRIXADAPTER_DECL_HPP
#define AMESOS2_KOKKOSCRSMATRIX_MATRIXADAPTER_DECL_HPP

#include "Amesos2_config.h"

#include "Amesos2_MatrixAdapter_decl.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include "Amesos2_Kokkos_View_Copy_Assign.hpp"

namespace Amesos2 {

  /**
   * \brief MatrixAdapter definitions for KokkosSparse::CrsMatrix objects.
   *
   * All other significant functionality is inherited from this
   * class's superclass.
   *
   * \ingroup amesos2_matrix_adapters
   */
  template <typename Scalar,
            typename LocalOrdinal,
            typename ExecutionSpace>
  class ConcreteMatrixAdapter<KokkosSparse::CrsMatrix<Scalar, LocalOrdinal, ExecutionSpace>> :
    public MatrixAdapter<KokkosSparse::CrsMatrix<Scalar, LocalOrdinal, ExecutionSpace>>
  {
  public:
    typedef KokkosSparse::CrsMatrix<Scalar, LocalOrdinal, ExecutionSpace> matrix_t;

  public:
    typedef typename MatrixTraits<matrix_t>::scalar_t                 scalar_t;
    typedef typename MatrixTraits<matrix_t>::local_ordinal_t   local_ordinal_t;
    typedef typename MatrixTraits<matrix_t>::global_ordinal_t global_ordinal_t;
    typedef typename MatrixTraits<matrix_t>::node_t                     node_t;
    typedef typename MatrixTraits<matrix_t>::global_size_t       global_size_t;
    typedef typename MatrixTraits<matrix_t>::sparse_ptr_type       spmtx_ptr_t;
    typedef typename MatrixTraits<matrix_t>::sparse_idx_type       spmtx_idx_t;
    typedef typename MatrixTraits<matrix_t>::sparse_values_type   spmtx_vals_t;
    typedef no_special_impl                                       get_crs_spec;
    typedef no_special_impl                                       get_ccs_spec;
    typedef ConcreteMatrixAdapter<matrix_t>                               type;
    typedef typename MatrixTraits<matrix_t>::major_access         major_access;

    ConcreteMatrixAdapter(Teuchos::RCP<matrix_t> m);

    Teuchos::RCP<const MatrixAdapter<matrix_t> > get_impl(const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > map, EDistribution distribution = ROOTED) const;

    const Teuchos::RCP<const Teuchos::Comm<int> > getComm_impl() const;
    global_size_t getGlobalNumRows_impl() const;
    global_size_t getGlobalNumCols_impl() const;
    global_size_t getGlobalNNZ_impl() const;
    spmtx_ptr_t  getSparseRowPtr() const;
    spmtx_idx_t  getSparseColInd() const;
    spmtx_vals_t getSparseValues() const;

    template<class KV>
    void getSparseRowPtr_kokkos_view(KV & view) const {
      deep_copy_or_assign_view(view, this->mat_->graph.row_map);
    }

    template<class KV>
    void getSparseColInd_kokkos_view(KV & view) const {
      deep_copy_or_assign_view(view, this->mat_->graph.entries);
    }

    template<class KV>
    void getSparseValues_kokkos_view(KV & view) const {
      deep_copy_or_assign_view(view, this->mat_->values);
    }

    size_t getGlobalRowNNZ_impl(global_ordinal_t row) const;
    size_t getLocalRowNNZ_impl(local_ordinal_t row) const;
    size_t getGlobalColNNZ_impl(global_ordinal_t col) const;
    size_t getLocalColNNZ_impl(local_ordinal_t col) const;

    global_size_t getRowIndexBase() const;
    global_size_t getColumnIndexBase() const;

    const Teuchos::RCP<const Tpetra::Map<local_ordinal_t,
        global_ordinal_t,
        node_t> >
    getMap_impl() const;

    const Teuchos::RCP<const Tpetra::Map<local_ordinal_t,
        global_ordinal_t,
        node_t> >
    getRowMap_impl() const;

    const Teuchos::RCP<const Tpetra::Map<local_ordinal_t,
        global_ordinal_t,
        node_t> >
    getColMap_impl() const;

    // implementation functions
    void getGlobalRowCopy_impl(global_ordinal_t row,
             const Teuchos::ArrayView<global_ordinal_t>& indices,
             const Teuchos::ArrayView<scalar_t>& vals,
             size_t& nnz) const;

    void getGlobalColCopy_impl(global_ordinal_t col,
             const Teuchos::ArrayView<global_ordinal_t>& indices,
             const Teuchos::ArrayView<scalar_t>& vals,
             size_t& nnz) const;

  };

} // end namespace Amesos2

#endif  // AMESOS2_KOKKOSCRSMATRIX_MATRIXADAPTER_DECL_HPP
