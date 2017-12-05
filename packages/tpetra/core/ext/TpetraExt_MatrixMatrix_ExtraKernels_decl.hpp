// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// ************************************************************************
// @HEADER

#ifndef TPETRA_MATRIXMATRIX_EXTRAKERNELS_DECL_HPP
#define TPETRA_MATRIXMATRIX_EXTRAKERNELS_DECL_HPP

#ifdef HAVE_KOKKOSKERNELS_EXPERIMENTAL
#include "KokkosSparse_spgemm.hpp"
#endif
namespace Tpetra {

namespace MatrixMatrix {


#ifdef HAVE_KOKKOSKERNELS_EXPERIMENTAL
  template <typename KernelHandle,
  typename alno_row_view_t_,
  typename alno_nnz_view_t_,
  typename blno_row_view_t_,
  typename blno_nnz_view_t_,
  typename clno_row_view_t_>
  void jacobi_mm_symbolic(
      KernelHandle *handle,
      typename KernelHandle::nnz_lno_t m,
      typename KernelHandle::nnz_lno_t n,
      typename KernelHandle::nnz_lno_t k,
      alno_row_view_t_ row_mapA,
      alno_nnz_view_t_ entriesA,
      bool transposeA,
      blno_row_view_t_ row_mapB,
      blno_nnz_view_t_ entriesB,
      bool transposeB,
      clno_row_view_t_ row_mapC);



 template <typename KernelHandle,
    typename alno_row_view_t_,
    typename alno_nnz_view_t_,
    typename ascalar_nnz_view_t_,
    typename adiag_nnz_view_t_,
    typename ascalar_t_;
    typename blno_row_view_t_,
    typename blno_nnz_view_t_,
    typename bscalar_nnz_view_t_,
    typename clno_row_view_t_,
    typename clno_nnz_view_t_,
    typename cscalar_nnz_view_t_>
  void jacobi_mm_numeric(
      KernelHandle *handle,
      typename KernelHandle::nnz_lno_t m,
      typename KernelHandle::nnz_lno_t n,
      typename KernelHandle::nnz_lno_t k,
      alno_row_view_t_ row_mapA,
      alno_nnz_view_t_ entriesA,
      ascalar_nnz_view_t_ valuesA,
      adiag_nnz_view_t_ diagA,
      ascalar_t_ omegaA,
      bool transposeA,
      blno_row_view_t_ row_mapB,
      blno_nnz_view_t_ entriesB,
      bscalar_nnz_view_t_ valuesB,
      bool transposeB,
      clno_row_view_t_ row_mapC,
      clno_nnz_view_t_ &entriesC,
      cscalar_nnz_view_t_ &valuesC);
#endif                         
  


}//end namespace MM


} // end of Tpetra namespace


#endif
