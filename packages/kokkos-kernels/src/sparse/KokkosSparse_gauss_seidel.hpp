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
#ifndef _KOKKOS_GAUSSSEIDEL_HPP
#define _KOKKOS_GAUSSSEIDEL_HPP

#include "KokkosSparse_gauss_seidel_impl.hpp"
#include "KokkosKernels_Handle.hpp"

namespace KokkosSparse{

namespace Experimental{

  template <typename KernelHandle, typename lno_row_view_t_, typename lno_nnz_view_t_>
  void gauss_seidel_symbolic(
      KernelHandle *handle,
      typename KernelHandle::nnz_lno_t num_rows,
      typename KernelHandle::nnz_lno_t num_cols,
      lno_row_view_t_ row_map,
      lno_nnz_view_t_ entries,
      bool is_graph_symmetric = true){
    typedef typename Impl::GaussSeidel<KernelHandle, lno_row_view_t_,
          lno_nnz_view_t_, typename KernelHandle::in_scalar_nnz_view_t> SGS;
    SGS sgs(handle,num_rows, num_cols, row_map, entries, is_graph_symmetric);
    sgs.initialize_symbolic();
  }

  template <typename KernelHandle,
            typename lno_row_view_t_,
            typename lno_nnz_view_t_,
            typename scalar_nnz_view_t_>
  void gauss_seidel_numeric(KernelHandle *handle,
      typename KernelHandle::nnz_lno_t num_rows,
      typename KernelHandle::nnz_lno_t num_cols,
      lno_row_view_t_ row_map,
      lno_nnz_view_t_ entries,
      scalar_nnz_view_t_ values,
      bool is_graph_symmetric = true
      ){
    typedef typename Impl::GaussSeidel
        <KernelHandle,lno_row_view_t_,
        lno_nnz_view_t_,scalar_nnz_view_t_> SGS;
    SGS sgs(handle, num_rows, num_cols, row_map, entries, values, is_graph_symmetric);
    sgs.initialize_numeric();
  }

  template <typename KernelHandle,
    typename lno_row_view_t_,
    typename lno_nnz_view_t_,
    typename scalar_nnz_view_t_,
    typename x_scalar_view_t,
    typename y_scalar_view_t>
  void symmetric_gauss_seidel_apply(KernelHandle *handle,
      typename KernelHandle::nnz_lno_t num_rows,
      typename KernelHandle::nnz_lno_t num_cols,
      lno_row_view_t_ row_map,
      lno_nnz_view_t_ entries,
      scalar_nnz_view_t_ values,
      x_scalar_view_t x_lhs_output_vec,
      y_scalar_view_t y_rhs_input_vec,
      bool init_zero_x_vector = false,
      bool update_y_vector = true,
      int numIter = 1){
    typedef
        typename Impl::GaussSeidel
            <KernelHandle,
            lno_row_view_t_, lno_nnz_view_t_,scalar_nnz_view_t_ > SGS;

    SGS sgs(handle, num_rows, num_cols, row_map, entries, values);
    sgs.apply(
        x_lhs_output_vec,
        y_rhs_input_vec,
        init_zero_x_vector,
        numIter,
        true,
        true,update_y_vector);

  }

  template <class KernelHandle,
  typename lno_row_view_t_,
  typename lno_nnz_view_t_,
  typename scalar_nnz_view_t_,
  typename x_scalar_view_t, typename y_scalar_view_t>
  void forward_sweep_gauss_seidel_apply(
      KernelHandle *handle,
      typename KernelHandle::nnz_lno_t num_rows,
      typename KernelHandle::nnz_lno_t num_cols,
      lno_row_view_t_ row_map,
      lno_nnz_view_t_ entries,
      scalar_nnz_view_t_ values,
      x_scalar_view_t x_lhs_output_vec,
      y_scalar_view_t y_rhs_input_vec,
      bool init_zero_x_vector = false,
      bool update_y_vector = true,
      int numIter = 1){
    typedef typename Impl::GaussSeidel <KernelHandle,
            lno_row_view_t_, lno_nnz_view_t_,scalar_nnz_view_t_ > SGS;

    SGS sgs(handle, num_rows, num_cols, row_map, entries, values);
    sgs.apply(
        x_lhs_output_vec,
        y_rhs_input_vec,
        init_zero_x_vector,
        numIter,
        true,
        false,update_y_vector);

  }
  template <class KernelHandle,
  typename lno_row_view_t_,
  typename lno_nnz_view_t_,
  typename scalar_nnz_view_t_,
  typename x_scalar_view_t, typename y_scalar_view_t>
  void backward_sweep_gauss_seidel_apply(
      KernelHandle *handle,
      typename KernelHandle::nnz_lno_t num_rows,
      typename KernelHandle::nnz_lno_t num_cols,
      lno_row_view_t_ row_map,
      lno_nnz_view_t_ entries,
      scalar_nnz_view_t_ values,
      x_scalar_view_t x_lhs_output_vec,
      y_scalar_view_t y_rhs_input_vec,
      bool init_zero_x_vector = false,
      bool update_y_vector = true,
      int numIter = 1){
    typedef typename Impl::GaussSeidel <KernelHandle,
            lno_row_view_t_, lno_nnz_view_t_,scalar_nnz_view_t_ > SGS;
    SGS sgs(handle, num_rows, num_cols, row_map, entries, values);
    sgs.apply(
        x_lhs_output_vec,
        y_rhs_input_vec,
        init_zero_x_vector,
        numIter,
        false,
        true, update_y_vector);

  }
}
}
#endif
