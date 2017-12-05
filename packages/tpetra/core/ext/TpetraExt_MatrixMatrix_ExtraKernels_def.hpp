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

#ifndef TPETRA_MATRIXMATRIX_EXTRAKERNELS_DEF_HPP
#define TPETRA_MATRIXMATRIX_EXTRAKERNELS_DEF_HPP
#include "TpetraExt_MatrixMatrix_ExtraKernels_decl.hpp"

namespace Tpetra {

namespace MatrixMatrix{

#ifdef HAVE_KOKKOSKERNELS_EXPERIMENTAL
template <typename KernelHandle,
          typename alno_row_view_t_,
          typename alno_nnz_view_t_,
          typename blno_row_view_t_,
          typename blno_nnz_view_t_,
          typename clno_row_view_t_>
void jacobi_mm_symbolic(KernelHandle *handle,                                            
                        typename KernelHandle::nnz_lno_t m,
                        typename KernelHandle::nnz_lno_t n,
                        typename KernelHandle::nnz_lno_t k,
                        alno_row_view_t_ row_mapA,
                        alno_nnz_view_t_ entriesA,
                        bool transposeA,
                        blno_row_view_t_ row_mapB,
                        blno_nnz_view_t_ entriesB,
                        bool transposeB,
                        clno_row_view_t_ row_mapC) {
  // Just use the SPGEMM symbolic for now
  spgemm_symbolic(handle,m,n,k,riow_mapA,entriesA,transposeA,row_mapB,entriesB,transposeB,row_mapC);

}

  

template <typename KernelHandle,
          typename alno_row_view_t_,
          typename alno_nnz_view_t_,
          typename ascalar_nnz_view_t_,
          typename ascalar_t_,
          typename ddiag_nnz_view_,
          typename blno_row_view_t_,
          typename blno_nnz_view_t_,
          typename bscalar_nnz_view_t_,
          typename clno_row_view_t_,
          typename clno_nnz_view_t_,
          typename cscalar_nnz_view_t_>
void jacobi_mm_numeric(KernelHandle *handle,                       
                       typename KernelHandle::nnz_lno_t m,
                       typename KernelHandle::nnz_lno_t n,
                       typename KernelHandle::nnz_lno_t k,
                       alno_row_view_t_ row_mapA,
                       alno_nnz_view_t_ entriesA,
                       ascalar_nnz_view_t_ valuesA,
                       ascalar_t_ omegaA,
                       ddiag_nnz_view_t_ diaginvD,
                       bool transposeA,
                       blno_row_view_t_ row_mapB,
                       blno_nnz_view_t_ entriesB,
                       bscalar_nnz_view_t_ valuesB,
                       bool transposeB,
                       clno_row_view_t_ row_mapC,
                       clno_nnz_view_t_ &entriesC,
                       cscalar_nnz_view_t_ &valuesC) {


  using namespace KokkosSparse::Impl;
  
  typedef typename KernelHandle::SPGEMMHandleType spgemmHandleType;
  spgemmHandleType *sh = handle->get_spgemm_handle();
  if (!sh->is_symbolic_called()){
    spgemm_symbolic<KernelHandle,
      alno_row_view_t_, alno_nnz_view_t_,
      blno_row_view_t_, blno_nnz_view_t_,
      clno_row_view_t_>(handle, m, n, k,                        
                        row_mapA, entriesA, transposeA,
                        row_mapB, entriesB, transposeB,
                        row_mapC
                        );
    typename clno_row_view_t_::value_type c_nnz_size = handle->get_spgemm_handle()->get_c_nnz();
    if (c_nnz_size){
      entriesC = clno_nnz_view_t_ (Kokkos::ViewAllocateWithoutInitializing("entriesC"), c_nnz_size);
      valuesC = cscalar_nnz_view_t_ (Kokkos::ViewAllocateWithoutInitializing("valuesC"), c_nnz_size);
    }
  }
  

  switch (sh->get_algorithm_type()){
  case SPGEMM_DEFAULT:
  case SPGEMM_KK_MEMSPEED:
  case SPGEMM_KK_SPEED:
  case SPGEMM_KK_MEMORY:
  case SPGEMM_KK_MEMORY2:
  case SPGEMM_KK_COLOR:
  case SPGEMM_KK_MULTICOLOR:
  case SPGEMM_KK_MULTICOLOR2:
  case SPGEMM_KK_MULTIMEM:
  case SPGEMM_KK_OUTERMULTIMEM:
    {
      KokkosJacobi
      <KernelHandle,
      alno_row_view_t_, alno_nnz_view_t_, ascalar_nnz_view_t_,
         ascalar_t_, ddiag_nnz_view_t_,
        blno_row_view_t_, blno_nnz_view_t_,  bscalar_nnz_view_t_>
        kjacobi (handle,m,n,k,row_mapA, entriesA, valuesA, omegaA, diaginvD, transposeA, row_mapB, entriesB, valuesB, transposeB);
      kspgemm.KokkosJacobi_numeric(row_mapC, entriesC, valuesC);
    }
  default:
    throw std::runtime_error("Invalid spgemm algorithm type for Jacobi");
  }

}
 
#endif //KokkosKernels_Experimental

}//MatrixMatirx
}//Tpetra
                        

#endif
