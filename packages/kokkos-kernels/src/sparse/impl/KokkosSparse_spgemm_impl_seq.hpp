/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
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
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

namespace KokkosSparse{

namespace Impl{

template <typename KernelHandle,
  typename alno_row_view_t_,
  typename alno_nnz_view_t_,
  typename ascalar_nnz_view_t_,
  typename blno_row_view_t_,
  typename blno_nnz_view_t_,
  typename bscalar_nnz_view_t_,
  typename clno_row_view_t_,
  typename clno_nnz_view_t_,
  typename cscalar_nnz_view_t_>
void spgemm_debug(
    KernelHandle *handle,
    typename KernelHandle::nnz_lno_t m,
    typename KernelHandle::nnz_lno_t n,
    typename KernelHandle::nnz_lno_t k,
    alno_row_view_t_ row_mapA,
    alno_nnz_view_t_ entriesA,
    ascalar_nnz_view_t_ valuesA,

    bool transposeA,
    blno_row_view_t_ row_mapB,
    blno_nnz_view_t_ entriesB,
    bscalar_nnz_view_t_ valuesB,
    bool transposeB,
    clno_row_view_t_ row_mapC,
    clno_nnz_view_t_ &entriesC,
    cscalar_nnz_view_t_ &valuesC
    ){
  typename alno_row_view_t_::HostMirror h_rma = Kokkos::create_mirror_view (row_mapA);
  Kokkos::deep_copy (h_rma, row_mapA);
  typename alno_nnz_view_t_::HostMirror h_enta = Kokkos::create_mirror_view (entriesA);
  Kokkos::deep_copy (h_enta, entriesA);
  typename ascalar_nnz_view_t_::HostMirror h_vala = Kokkos::create_mirror_view (valuesA);
  Kokkos::deep_copy (h_vala, valuesA);

  typename blno_row_view_t_::HostMirror h_rmb = Kokkos::create_mirror_view (row_mapB);
  Kokkos::deep_copy (h_rmb, row_mapB);
  typename blno_nnz_view_t_::HostMirror h_entb = Kokkos::create_mirror_view (entriesB);
  Kokkos::deep_copy (h_entb, entriesB);
  typename bscalar_nnz_view_t_::HostMirror h_valb = Kokkos::create_mirror_view (valuesB);
  Kokkos::deep_copy (h_valb, valuesB);
  typename clno_row_view_t_::HostMirror h_rmc = Kokkos::create_mirror_view (row_mapC);
  Kokkos::fence();

  typedef typename KernelHandle::nnz_lno_t lno_t;
  typedef typename KernelHandle::size_type size_type;
  typedef typename KernelHandle::nnz_scalar_t scalar_t;

  std::vector<scalar_t> accumulator(k, 0);
  std::vector<bool> acc_flag(k, false);

  std::vector<lno_t> result_c_col_indices(k);
  std::vector<scalar_t> result_c_col_values(k);

  const size_type alloc_step = 10000;
  size_type result_index = 0;
  size_type current_c_size = k;

  h_rmc(0) = 0;
  for (lno_t i = 0; i < m; ++i){
    const size_type a_row_begin = h_rma(i);
    const size_type a_row_end = h_rma(i + 1);
    lno_t a_row_size = a_row_end - a_row_begin;

    for (lno_t j = 0; j < a_row_size; ++j){
      size_type ind = a_row_begin + j;
      lno_t col = h_enta(ind);
      scalar_t val = h_vala(ind);
      //if (i == 0) std::cout << "a row:" <<  col << std::endl;
      const size_type b_row_begin = h_rmb(col);
      const size_type b_row_end = h_rmb(col + 1);
      lno_t b_row_size = b_row_end - b_row_begin;
      for (lno_t z = 0; z < b_row_size; ++z){
        size_type ind_ = b_row_begin + z;
        lno_t b_col = h_entb(ind_);
        scalar_t b_val = h_valb(ind_);
        //if (i == 0) std::cout << "\tb col:" <<  b_col << std::endl;
        if (acc_flag[b_col] == false){
          acc_flag[b_col] = true;
          result_c_col_indices[result_index++] = b_col;
          if (current_c_size == result_index){
            current_c_size += alloc_step;
            result_c_col_indices.resize(current_c_size);
            result_c_col_values.resize(current_c_size);
          }
        }
        accumulator[b_col] += b_val * val;
      }
    }
    h_rmc(i+1) = result_index;
    size_type c_row_begin = h_rmc(i);
    lno_t c_row_size =  result_index - c_row_begin;

    //if (i == 0) std::cout << "result_cols" << std::endl;

    for (lno_t j = 0; j < c_row_size; ++j){

      size_type ind = c_row_begin + j;
      lno_t result_col = result_c_col_indices[ind];
      //if (i == 0) std::cout << result_col << std::endl;
      result_c_col_values[ind] = accumulator[result_col];
      accumulator[result_col] = 0;
      acc_flag[result_col] = false;
    }

  }

  entriesC = clno_nnz_view_t_("entriesC", result_index);
  valuesC = cscalar_nnz_view_t_("entriesC", result_index);

  typename clno_nnz_view_t_::HostMirror h_entc = Kokkos::create_mirror_view (entriesC);
  typename cscalar_nnz_view_t_::HostMirror h_valc = Kokkos::create_mirror_view (valuesC);

  for (size_type i = 0; i < result_index; ++i){
    h_entc(i) = result_c_col_indices[i];
    h_valc(i) = result_c_col_values[i];
  }

  Kokkos::deep_copy (row_mapC, h_rmc);
  Kokkos::deep_copy (entriesC, h_entc);
  Kokkos::deep_copy (valuesC, h_valc);
  Kokkos::fence();


}

}
}
