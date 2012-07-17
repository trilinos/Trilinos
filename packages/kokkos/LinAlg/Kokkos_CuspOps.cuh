//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
//@HEADER

#ifndef KOKKOS_CUSPOPS_CUH
#define KOKKOS_CUSPOPS_CUH

#include <cusp/csr_matrix.h>

namespace Kokkos {

  namespace Cuspdetails {

    template <class O, class S>
    struct types {
      typedef  thrust::device_ptr<O>                                              ind_ptr;
      typedef  thrust::device_ptr<S>                                              val_ptr;
      typedef  cusp::array1d_view< ind_ptr >                             ind_arr_view;
      typedef  cusp::array1d_view< val_ptr >                             val_arr_view;
      typedef  cusp::array2d_view< val_arr_view, cusp::column_major>              val_arr_view2d;
      typedef  cusp::csr_matrix_view<ind_arr_view, ind_arr_view, val_arr_view>    csr_view;
    };

    // forward declaration for cusp adaptors
    template <class Ordinal, class Scalar> 
    void cuspCrsMultiply(const typename types<const Ordinal,const Scalar>::csr_view &A, 
                         const typename types<const Ordinal,const Scalar>::val_arr_view2d &X, 
                         const typename types<const Ordinal,      Scalar>::val_arr_view2d &Y) {}
  } // end of namespace Cuspdetails

}

#endif // KOKKOS_CUSPOPS_CUH
