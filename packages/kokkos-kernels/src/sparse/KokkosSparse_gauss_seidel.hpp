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

#include "KokkosSparse_gauss_seidel_spec.hpp"
#include "KokkosKernels_Handle.hpp"
#include "KokkosKernels_helpers.hpp"

namespace KokkosSparse{

namespace Experimental{

  template <typename KernelHandle, typename lno_row_view_t_, typename lno_nnz_view_t_>
  void gauss_seidel_symbolic(
      KernelHandle *handle,
      typename KernelHandle::const_nnz_lno_t num_rows,
      typename KernelHandle::const_nnz_lno_t num_cols,
      lno_row_view_t_ row_map,
      lno_nnz_view_t_ entries,
	  bool is_graph_symmetric = true){



	  static_assert (std::is_same<typename KernelHandle::const_size_type,
			  typename lno_row_view_t_::const_value_type>::value,
			  "KokkosSparse::gauss_seidel_symbolic: Size type of the matrix should be same as kernelHandle sizetype.");

	  static_assert (std::is_same<typename KernelHandle::const_nnz_lno_t,
			  typename lno_nnz_view_t_::const_value_type>::value,
			  "KokkosSparse::gauss_seidel_symbolic: lno type of the matrix should be same as kernelHandle lno_t.");


	  typedef typename KernelHandle::const_size_type c_size_t;
	  typedef typename KernelHandle::const_nnz_lno_t c_lno_t;
	  typedef typename KernelHandle::const_nnz_scalar_t c_scalar_t;

	  typedef typename KernelHandle::HandleExecSpace c_exec_t;
	  typedef typename KernelHandle::HandleTempMemorySpace c_temp_t;
	  typedef typename KernelHandle::HandlePersistentMemorySpace c_persist_t;

	  typedef typename  KokkosKernels::Experimental::KokkosKernelsHandle<c_size_t, c_lno_t, c_scalar_t, c_exec_t, c_temp_t, c_persist_t> const_handle_type;
	  //const_handle_type tmp_handle = *handle;
	  const_handle_type tmp_handle (*handle);

	  typedef Kokkos::View<
			  typename lno_row_view_t_::const_value_type*,
			  typename KokkosKernels::Impl::GetUnifiedLayout<lno_row_view_t_>::array_layout,
			  typename lno_row_view_t_::device_type,
			  Kokkos::MemoryTraits<Kokkos::Unmanaged> > Internal_alno_row_view_t_;

	  typedef Kokkos::View<
			  typename lno_nnz_view_t_::const_value_type*,
			  typename KokkosKernels::Impl::GetUnifiedLayout<lno_nnz_view_t_>::array_layout,
			  typename lno_nnz_view_t_::device_type,
			  Kokkos::MemoryTraits<Kokkos::Unmanaged> > Internal_alno_nnz_view_t_;


	  //Internal_alno_row_view_t_ const_a_r  = row_map;
	  //Internal_alno_nnz_view_t_ const_a_l  = entries;
          Internal_alno_row_view_t_ const_a_r (row_map.data(), row_map.extent(0));
          Internal_alno_nnz_view_t_ const_a_l (entries.data(), entries.extent(0));

	  using namespace KokkosSparse::Impl;

	  GAUSS_SEIDEL_SYMBOLIC<const_handle_type, Internal_alno_row_view_t_, Internal_alno_nnz_view_t_>::gauss_seidel_symbolic
	  (&tmp_handle, num_rows, num_cols, const_a_r, const_a_l, is_graph_symmetric);

  }

  template <typename KernelHandle, typename lno_row_view_t_, typename lno_nnz_view_t_>
  void block_gauss_seidel_symbolic(
      KernelHandle *handle,
      typename KernelHandle::const_nnz_lno_t num_rows,
      typename KernelHandle::const_nnz_lno_t num_cols,
	  typename KernelHandle::const_nnz_lno_t block_size,
      lno_row_view_t_ row_map,
      lno_nnz_view_t_ entries,
	  bool is_graph_symmetric = true){


	  handle->get_gs_handle()->set_block_size(block_size);

	  gauss_seidel_symbolic(handle,
			  num_rows, num_cols,
			  row_map, entries, is_graph_symmetric);

  }


  template <typename KernelHandle,
            typename lno_row_view_t_,
            typename lno_nnz_view_t_,
            typename scalar_nnz_view_t_>
  void gauss_seidel_numeric(KernelHandle *handle,
      typename KernelHandle::const_nnz_lno_t num_rows,
      typename KernelHandle::const_nnz_lno_t num_cols,
      lno_row_view_t_ row_map,
      lno_nnz_view_t_ entries,
      scalar_nnz_view_t_ values,
      bool is_graph_symmetric = true
      ){

	  static_assert (std::is_same<typename KernelHandle::const_size_type,
			  typename lno_row_view_t_::const_value_type>::value,
			  "KokkosSparse::gauss_seidel_symbolic: Size type of the matrix should be same as kernelHandle sizetype.");

	  static_assert (std::is_same<typename KernelHandle::const_nnz_lno_t,
			  typename lno_nnz_view_t_::const_value_type>::value,
			  "KokkosSparse::gauss_seidel_symbolic: lno type of the matrix should be same as kernelHandle lno_t.");

	  static_assert (std::is_same<typename KernelHandle::const_nnz_scalar_t,
	  			  typename scalar_nnz_view_t_::const_value_type>::value,
	  			  "KokkosSparse::gauss_seidel_symbolic: scalar type of the matrix should be same as kernelHandle scalar_t.");


	  typedef typename KernelHandle::const_size_type c_size_t;
	  typedef typename KernelHandle::const_nnz_lno_t c_lno_t;
	  typedef typename KernelHandle::const_nnz_scalar_t c_scalar_t;

	  typedef typename KernelHandle::HandleExecSpace c_exec_t;
	  typedef typename KernelHandle::HandleTempMemorySpace c_temp_t;
	  typedef typename KernelHandle::HandlePersistentMemorySpace c_persist_t;

	  typedef typename  KokkosKernels::Experimental::KokkosKernelsHandle<c_size_t, c_lno_t, c_scalar_t, c_exec_t, c_temp_t, c_persist_t> const_handle_type;
	  //const_handle_type tmp_handle = *handle;
	  const_handle_type tmp_handle (*handle);

	  typedef Kokkos::View<
			  typename lno_row_view_t_::const_value_type*,
			  typename KokkosKernels::Impl::GetUnifiedLayout<lno_row_view_t_>::array_layout,
			  typename lno_row_view_t_::device_type,
			  Kokkos::MemoryTraits<Kokkos::Unmanaged> > Internal_alno_row_view_t_;

	  typedef Kokkos::View<
			  typename lno_nnz_view_t_::const_value_type*,
			  typename KokkosKernels::Impl::GetUnifiedLayout<lno_nnz_view_t_>::array_layout,
			  typename lno_nnz_view_t_::device_type,
			  Kokkos::MemoryTraits<Kokkos::Unmanaged> > Internal_alno_nnz_view_t_;

	  typedef Kokkos::View<
			  typename scalar_nnz_view_t_::const_value_type*,
			  typename KokkosKernels::Impl::GetUnifiedLayout<scalar_nnz_view_t_>::array_layout,
			  typename scalar_nnz_view_t_::device_type,
			  Kokkos::MemoryTraits<Kokkos::Unmanaged> > Internal_ascalar_nnz_view_t_;

	  Internal_alno_row_view_t_ const_a_r (row_map.data(), row_map.extent(0));
	  Internal_alno_nnz_view_t_ const_a_l (entries.data(), entries.extent(0));
	  Internal_ascalar_nnz_view_t_ const_a_v (values.data(), values.extent(0));


	  using namespace KokkosSparse::Impl;

	  GAUSS_SEIDEL_NUMERIC<const_handle_type, Internal_alno_row_view_t_, Internal_alno_nnz_view_t_, Internal_ascalar_nnz_view_t_>::gauss_seidel_numeric
	  	  (&tmp_handle, num_rows, num_cols, const_a_r, const_a_l, const_a_v, is_graph_symmetric);

  }

  template <typename KernelHandle,
            typename lno_row_view_t_,
            typename lno_nnz_view_t_,
            typename scalar_nnz_view_t_>
  void block_gauss_seidel_numeric(KernelHandle *handle,
      typename KernelHandle::const_nnz_lno_t num_rows,
      typename KernelHandle::const_nnz_lno_t num_cols,
	  typename KernelHandle::const_nnz_lno_t block_size,
      lno_row_view_t_ row_map,
      lno_nnz_view_t_ entries,
      scalar_nnz_view_t_ values,
      bool is_graph_symmetric = true
      ){
	  handle->get_gs_handle()->set_block_size(block_size);

	  gauss_seidel_numeric(handle,
			  num_rows,num_cols,
			  row_map, entries, values,
			  is_graph_symmetric
	        );
  }

  template <typename KernelHandle,
    typename lno_row_view_t_,
    typename lno_nnz_view_t_,
    typename scalar_nnz_view_t_,
    typename x_scalar_view_t,
    typename y_scalar_view_t>
  void symmetric_gauss_seidel_apply(KernelHandle *handle,
      typename KernelHandle::const_nnz_lno_t num_rows,
      typename KernelHandle::const_nnz_lno_t num_cols,
      lno_row_view_t_ row_map,
      lno_nnz_view_t_ entries,
      scalar_nnz_view_t_ values,
      x_scalar_view_t x_lhs_output_vec,
      y_scalar_view_t y_rhs_input_vec,
      bool init_zero_x_vector = false,
      bool update_y_vector = true,
      int numIter = 1){

	  static_assert (std::is_same<typename KernelHandle::const_size_type,
			  typename lno_row_view_t_::const_value_type>::value,
			  "KokkosSparse::symmetric_gauss_seidel_apply: Size type of the matrix should be same as kernelHandle sizetype.");

	  static_assert (std::is_same<typename KernelHandle::const_nnz_lno_t,
			  typename lno_nnz_view_t_::const_value_type>::value,
			  "KokkosSparse::symmetric_gauss_seidel_apply: lno type of the matrix should be same as kernelHandle lno_t.");

	  static_assert (std::is_same<typename KernelHandle::const_nnz_scalar_t,
	  			  typename scalar_nnz_view_t_::const_value_type>::value,
	  			  "KokkosSparse::symmetric_gauss_seidel_apply: scalar type of the matrix should be same as kernelHandle scalar_t.");

	  static_assert (std::is_same<typename KernelHandle::const_nnz_scalar_t,
	  			  typename y_scalar_view_t::const_value_type>::value,
	  			  "KokkosSparse::symmetric_gauss_seidel_apply: scalar type of the y-vector should be same as kernelHandle scalar_t.");

	  static_assert (std::is_same<typename KernelHandle::nnz_scalar_t,
	  			  typename x_scalar_view_t::value_type>::value,
	  			  "KokkosSparse::symmetric_gauss_seidel_apply: scalar type of the x-vector should be same as kernelHandle non-const scalar_t.");




	  typedef typename KernelHandle::const_size_type c_size_t;
	  typedef typename KernelHandle::const_nnz_lno_t c_lno_t;
	  typedef typename KernelHandle::const_nnz_scalar_t c_scalar_t;

	  typedef typename KernelHandle::HandleExecSpace c_exec_t;
	  typedef typename KernelHandle::HandleTempMemorySpace c_temp_t;
	  typedef typename KernelHandle::HandlePersistentMemorySpace c_persist_t;

	  typedef typename  KokkosKernels::Experimental::KokkosKernelsHandle<c_size_t, c_lno_t, c_scalar_t, c_exec_t, c_temp_t, c_persist_t> const_handle_type;
	  //const_handle_type tmp_handle = *handle;
	  const_handle_type tmp_handle (*handle);

	  typedef Kokkos::View<
			  typename lno_row_view_t_::const_value_type*,
			  typename KokkosKernels::Impl::GetUnifiedLayout<lno_row_view_t_>::array_layout,
			  typename lno_row_view_t_::device_type,
			  Kokkos::MemoryTraits<Kokkos::Unmanaged> > Internal_alno_row_view_t_;

	  typedef Kokkos::View<
			  typename lno_nnz_view_t_::const_value_type*,
			  typename KokkosKernels::Impl::GetUnifiedLayout<lno_nnz_view_t_>::array_layout,
			  typename lno_nnz_view_t_::device_type,
			  Kokkos::MemoryTraits<Kokkos::Unmanaged> > Internal_alno_nnz_view_t_;

	  typedef Kokkos::View<
			  typename scalar_nnz_view_t_::const_value_type*,
			  typename KokkosKernels::Impl::GetUnifiedLayout<scalar_nnz_view_t_>::array_layout,
			  typename scalar_nnz_view_t_::device_type,
			  Kokkos::MemoryTraits<Kokkos::Unmanaged> > Internal_ascalar_nnz_view_t_;


	  typedef Kokkos::View<
			  typename y_scalar_view_t::const_value_type*,
			  typename KokkosKernels::Impl::GetUnifiedLayout<y_scalar_view_t>::array_layout,
			  typename y_scalar_view_t::device_type,
			  Kokkos::MemoryTraits<Kokkos::Unmanaged> > Internal_yscalar_nnz_view_t_;

	  typedef Kokkos::View<
			  typename x_scalar_view_t::non_const_value_type*,
			  typename KokkosKernels::Impl::GetUnifiedLayout<x_scalar_view_t>::array_layout,
			  typename x_scalar_view_t::device_type,
			  Kokkos::MemoryTraits<Kokkos::Unmanaged> > Internal_xscalar_nnz_view_t_;

	  /*

	  Internal_alno_row_view_t_ const_a_r  = row_map;
	  Internal_alno_nnz_view_t_ const_a_l  = entries;
	  Internal_ascalar_nnz_view_t_ const_a_v  = values;
	  Internal_xscalar_nnz_view_t_ nonconst_x_v = x_lhs_output_vec;
	  Internal_yscalar_nnz_view_t_ const_y_v = y_rhs_input_vec;
	  */

          Internal_alno_row_view_t_ const_a_r (row_map.data(), row_map.extent(0));
          Internal_alno_nnz_view_t_ const_a_l (entries.data(), entries.extent(0));
          Internal_ascalar_nnz_view_t_ const_a_v (values.data(), values.extent(0));

          Internal_xscalar_nnz_view_t_ nonconst_x_v (x_lhs_output_vec.data(), x_lhs_output_vec.extent(0));
          Internal_yscalar_nnz_view_t_ const_y_v (y_rhs_input_vec.data(), y_rhs_input_vec.extent(0));

	  using namespace KokkosSparse::Impl;

	  GAUSS_SEIDEL_APPLY<const_handle_type,
	  Internal_alno_row_view_t_, Internal_alno_nnz_view_t_, Internal_ascalar_nnz_view_t_,
	  Internal_xscalar_nnz_view_t_, Internal_yscalar_nnz_view_t_>::gauss_seidel_apply (
	  		  &tmp_handle, num_rows, num_cols,
			  const_a_r, const_a_l, const_a_v,
			  nonconst_x_v, const_y_v,
			  init_zero_x_vector,
			  update_y_vector,
			  numIter, true, true);


  }

  template <typename KernelHandle,
    typename lno_row_view_t_,
    typename lno_nnz_view_t_,
    typename scalar_nnz_view_t_,
    typename x_scalar_view_t,
    typename y_scalar_view_t>
  void symmetric_block_gauss_seidel_apply(KernelHandle *handle,
      typename KernelHandle::const_nnz_lno_t num_rows,
      typename KernelHandle::const_nnz_lno_t num_cols,
	  typename KernelHandle::const_nnz_lno_t block_size,

      lno_row_view_t_ row_map,
      lno_nnz_view_t_ entries,
      scalar_nnz_view_t_ values,
      x_scalar_view_t x_lhs_output_vec,
      y_scalar_view_t y_rhs_input_vec,
      bool init_zero_x_vector = false,
      bool update_y_vector = true,
      int numIter = 1){
	  handle->get_gs_handle()->set_block_size(block_size);
	  symmetric_gauss_seidel_apply(handle,num_rows,num_cols,
			  row_map,
		      entries,
		      values,
		      x_lhs_output_vec,
		      y_rhs_input_vec,
		      init_zero_x_vector,
		      update_y_vector,
		      numIter);
  }
  template <class KernelHandle,
  typename lno_row_view_t_,
  typename lno_nnz_view_t_,
  typename scalar_nnz_view_t_,
  typename x_scalar_view_t, typename y_scalar_view_t>
  void forward_sweep_gauss_seidel_apply(
      KernelHandle *handle,
      typename KernelHandle::const_nnz_lno_t num_rows,
      typename KernelHandle::const_nnz_lno_t num_cols,
      lno_row_view_t_ row_map,
      lno_nnz_view_t_ entries,
      scalar_nnz_view_t_ values,
      x_scalar_view_t x_lhs_output_vec,
      y_scalar_view_t y_rhs_input_vec,
      bool init_zero_x_vector = false,
      bool update_y_vector = true,
      int numIter = 1){

	  static_assert (std::is_same<typename KernelHandle::const_size_type,
			  typename lno_row_view_t_::const_value_type>::value,
			  "KokkosSparse::symmetric_gauss_seidel_apply: Size type of the matrix should be same as kernelHandle sizetype.");

	  static_assert (std::is_same<typename KernelHandle::const_nnz_lno_t,
			  typename lno_nnz_view_t_::const_value_type>::value,
			  "KokkosSparse::symmetric_gauss_seidel_apply: lno type of the matrix should be same as kernelHandle lno_t.");

	  static_assert (std::is_same<typename KernelHandle::const_nnz_scalar_t,
	  			  typename scalar_nnz_view_t_::const_value_type>::value,
	  			  "KokkosSparse::symmetric_gauss_seidel_apply: scalar type of the matrix should be same as kernelHandle scalar_t.");

	  static_assert (std::is_same<typename KernelHandle::const_nnz_scalar_t,
	  			  typename y_scalar_view_t::const_value_type>::value,
	  			  "KokkosSparse::symmetric_gauss_seidel_apply: scalar type of the y-vector should be same as kernelHandle scalar_t.");

	  static_assert (std::is_same<typename KernelHandle::nnz_scalar_t,
	  			  typename x_scalar_view_t::value_type>::value,
	  			  "KokkosSparse::symmetric_gauss_seidel_apply: scalar type of the x-vector should be same as kernelHandle non-const scalar_t.");




	  typedef typename KernelHandle::const_size_type c_size_t;
	  typedef typename KernelHandle::const_nnz_lno_t c_lno_t;
	  typedef typename KernelHandle::const_nnz_scalar_t c_scalar_t;

	  typedef typename KernelHandle::HandleExecSpace c_exec_t;
	  typedef typename KernelHandle::HandleTempMemorySpace c_temp_t;
	  typedef typename KernelHandle::HandlePersistentMemorySpace c_persist_t;

	  typedef typename  KokkosKernels::Experimental::KokkosKernelsHandle<c_size_t, c_lno_t, c_scalar_t, c_exec_t, c_temp_t, c_persist_t> const_handle_type;
	  //const_handle_type tmp_handle = *handle;
	  const_handle_type tmp_handle (*handle);

	  typedef Kokkos::View<
			  typename lno_row_view_t_::const_value_type*,
			  typename KokkosKernels::Impl::GetUnifiedLayout<lno_row_view_t_>::array_layout,
			  typename lno_row_view_t_::device_type,
			  Kokkos::MemoryTraits<Kokkos::Unmanaged> > Internal_alno_row_view_t_;

	  typedef Kokkos::View<
			  typename lno_nnz_view_t_::const_value_type*,
			  typename KokkosKernels::Impl::GetUnifiedLayout<lno_nnz_view_t_>::array_layout,
			  typename lno_nnz_view_t_::device_type,
			  Kokkos::MemoryTraits<Kokkos::Unmanaged> > Internal_alno_nnz_view_t_;

	  typedef Kokkos::View<
			  typename scalar_nnz_view_t_::const_value_type*,
			  typename KokkosKernels::Impl::GetUnifiedLayout<scalar_nnz_view_t_>::array_layout,
			  typename scalar_nnz_view_t_::device_type,
			  Kokkos::MemoryTraits<Kokkos::Unmanaged> > Internal_ascalar_nnz_view_t_;


	  typedef Kokkos::View<
			  typename y_scalar_view_t::const_value_type*,
			  typename KokkosKernels::Impl::GetUnifiedLayout<y_scalar_view_t>::array_layout,
			  typename y_scalar_view_t::device_type,
			  Kokkos::MemoryTraits<Kokkos::Unmanaged> > Internal_yscalar_nnz_view_t_;

	  typedef Kokkos::View<
			  typename x_scalar_view_t::non_const_value_type*,
			  typename KokkosKernels::Impl::GetUnifiedLayout<x_scalar_view_t>::array_layout,
			  typename x_scalar_view_t::device_type,
			  Kokkos::MemoryTraits<Kokkos::Unmanaged> > Internal_xscalar_nnz_view_t_;

	  /*
	  Internal_alno_row_view_t_ const_a_r  = row_map;
	  Internal_alno_nnz_view_t_ const_a_l  = entries;
	  Internal_ascalar_nnz_view_t_ const_a_v  = values;
	  Internal_xscalar_nnz_view_t_ nonconst_x_v = x_lhs_output_vec;
	  Internal_yscalar_nnz_view_t_ const_y_v = y_rhs_input_vec;
	  */

          Internal_alno_row_view_t_ const_a_r (row_map.data(), row_map.extent(0));
          Internal_alno_nnz_view_t_ const_a_l (entries.data(), entries.extent(0));
          Internal_ascalar_nnz_view_t_ const_a_v (values.data(), values.extent(0));

          Internal_xscalar_nnz_view_t_ nonconst_x_v (x_lhs_output_vec.data(), x_lhs_output_vec.extent(0));
          Internal_yscalar_nnz_view_t_ const_y_v (y_rhs_input_vec.data(), y_rhs_input_vec.extent(0));



	  using namespace KokkosSparse::Impl;

	  GAUSS_SEIDEL_APPLY<const_handle_type,
	  Internal_alno_row_view_t_, Internal_alno_nnz_view_t_, Internal_ascalar_nnz_view_t_,
	  Internal_xscalar_nnz_view_t_, Internal_yscalar_nnz_view_t_>::gauss_seidel_apply(
	  		  &tmp_handle, num_rows, num_cols,
			  const_a_r, const_a_l, const_a_v,
			  nonconst_x_v, const_y_v,
			  init_zero_x_vector,
			  update_y_vector,
			  numIter, true, false);


  }


  template <typename KernelHandle,
    typename lno_row_view_t_,
    typename lno_nnz_view_t_,
    typename scalar_nnz_view_t_,
    typename x_scalar_view_t,
    typename y_scalar_view_t>
  void forward_sweep_block_gauss_seidel_apply(KernelHandle *handle,
      typename KernelHandle::const_nnz_lno_t num_rows,
      typename KernelHandle::const_nnz_lno_t num_cols,
	  typename KernelHandle::const_nnz_lno_t block_size,

      lno_row_view_t_ row_map,
      lno_nnz_view_t_ entries,
      scalar_nnz_view_t_ values,
      x_scalar_view_t x_lhs_output_vec,
      y_scalar_view_t y_rhs_input_vec,
      bool init_zero_x_vector = false,
      bool update_y_vector = true,
      int numIter = 1){
	  handle->get_gs_handle()->set_block_size(block_size);
	  forward_sweep_gauss_seidel_apply(handle,num_rows,num_cols,
			  row_map,
		      entries,
		      values,
		      x_lhs_output_vec,
		      y_rhs_input_vec,
		      init_zero_x_vector,
		      update_y_vector,
		      numIter);
  }
  template <class KernelHandle,
  typename lno_row_view_t_,
  typename lno_nnz_view_t_,
  typename scalar_nnz_view_t_,
  typename x_scalar_view_t, typename y_scalar_view_t>
  void backward_sweep_gauss_seidel_apply(
      KernelHandle *handle,
      typename KernelHandle::const_nnz_lno_t num_rows,
      typename KernelHandle::const_nnz_lno_t num_cols,
      lno_row_view_t_ row_map,
      lno_nnz_view_t_ entries,
      scalar_nnz_view_t_ values,
      x_scalar_view_t x_lhs_output_vec,
      y_scalar_view_t y_rhs_input_vec,
      bool init_zero_x_vector = false,
      bool update_y_vector = true,
      int numIter = 1){

	  static_assert (std::is_same<typename KernelHandle::const_size_type,
			  typename lno_row_view_t_::const_value_type>::value,
			  "KokkosSparse::symmetric_gauss_seidel_apply: Size type of the matrix should be same as kernelHandle sizetype.");

	  static_assert (std::is_same<typename KernelHandle::const_nnz_lno_t,
			  typename lno_nnz_view_t_::const_value_type>::value,
			  "KokkosSparse::symmetric_gauss_seidel_apply: lno type of the matrix should be same as kernelHandle lno_t.");

	  static_assert (std::is_same<typename KernelHandle::const_nnz_scalar_t,
	  			  typename scalar_nnz_view_t_::const_value_type>::value,
	  			  "KokkosSparse::symmetric_gauss_seidel_apply: scalar type of the matrix should be same as kernelHandle scalar_t.");

	  static_assert (std::is_same<typename KernelHandle::const_nnz_scalar_t,
	  			  typename y_scalar_view_t::const_value_type>::value,
	  			  "KokkosSparse::symmetric_gauss_seidel_apply: scalar type of the y-vector should be same as kernelHandle scalar_t.");

	  static_assert (std::is_same<typename KernelHandle::nnz_scalar_t,
	  			  typename x_scalar_view_t::value_type>::value,
	  			  "KokkosSparse::symmetric_gauss_seidel_apply: scalar type of the x-vector should be same as kernelHandle non-const scalar_t.");




	  typedef typename KernelHandle::const_size_type c_size_t;
	  typedef typename KernelHandle::const_nnz_lno_t c_lno_t;
	  typedef typename KernelHandle::const_nnz_scalar_t c_scalar_t;

	  typedef typename KernelHandle::HandleExecSpace c_exec_t;
	  typedef typename KernelHandle::HandleTempMemorySpace c_temp_t;
	  typedef typename KernelHandle::HandlePersistentMemorySpace c_persist_t;

	  typedef typename  KokkosKernels::Experimental::KokkosKernelsHandle<c_size_t, c_lno_t, c_scalar_t, c_exec_t, c_temp_t, c_persist_t> const_handle_type;
	  //const_handle_type tmp_handle = *handle;
	  const_handle_type tmp_handle (*handle);

	  typedef Kokkos::View<
			  typename lno_row_view_t_::const_value_type*,
			  typename KokkosKernels::Impl::GetUnifiedLayout<lno_row_view_t_>::array_layout,
			  typename lno_row_view_t_::device_type,
			  Kokkos::MemoryTraits<Kokkos::Unmanaged> > Internal_alno_row_view_t_;

	  typedef Kokkos::View<
			  typename lno_nnz_view_t_::const_value_type*,
			  typename KokkosKernels::Impl::GetUnifiedLayout<lno_nnz_view_t_>::array_layout,
			  typename lno_nnz_view_t_::device_type,
			  Kokkos::MemoryTraits<Kokkos::Unmanaged> > Internal_alno_nnz_view_t_;

	  typedef Kokkos::View<
			  typename scalar_nnz_view_t_::const_value_type*,
			  typename KokkosKernels::Impl::GetUnifiedLayout<scalar_nnz_view_t_>::array_layout,
			  typename scalar_nnz_view_t_::device_type,
			  Kokkos::MemoryTraits<Kokkos::Unmanaged> > Internal_ascalar_nnz_view_t_;


	  typedef Kokkos::View<
			  typename y_scalar_view_t::const_value_type*,
			  typename KokkosKernels::Impl::GetUnifiedLayout<y_scalar_view_t>::array_layout,
			  typename y_scalar_view_t::device_type,
			  Kokkos::MemoryTraits<Kokkos::Unmanaged> > Internal_yscalar_nnz_view_t_;

	  typedef Kokkos::View<
			  typename x_scalar_view_t::non_const_value_type*,
			  typename KokkosKernels::Impl::GetUnifiedLayout<x_scalar_view_t>::array_layout,
			  typename x_scalar_view_t::device_type,
			  Kokkos::MemoryTraits<Kokkos::Unmanaged> > Internal_xscalar_nnz_view_t_;

	  /*
	  Internal_alno_row_view_t_ const_a_r  = row_map;
	  Internal_alno_nnz_view_t_ const_a_l  = entries;
	  Internal_ascalar_nnz_view_t_ const_a_v  = values;
	  Internal_xscalar_nnz_view_t_ nonconst_x_v = x_lhs_output_vec;
	  Internal_yscalar_nnz_view_t_ const_y_v = y_rhs_input_vec;
	  */
          Internal_alno_row_view_t_ const_a_r (row_map.data(), row_map.extent(0));
          Internal_alno_nnz_view_t_ const_a_l (entries.data(), entries.extent(0));
          Internal_ascalar_nnz_view_t_ const_a_v (values.data(), values.extent(0));

          Internal_xscalar_nnz_view_t_ nonconst_x_v (x_lhs_output_vec.data(), x_lhs_output_vec.extent(0));
          Internal_yscalar_nnz_view_t_ const_y_v (y_rhs_input_vec.data(), y_rhs_input_vec.extent(0));


	  using namespace KokkosSparse::Impl;

	  GAUSS_SEIDEL_APPLY<const_handle_type,
	  Internal_alno_row_view_t_, Internal_alno_nnz_view_t_, Internal_ascalar_nnz_view_t_,
	  Internal_xscalar_nnz_view_t_, Internal_yscalar_nnz_view_t_>::gauss_seidel_apply (
	  		  &tmp_handle, num_rows, num_cols,
			  const_a_r, const_a_l, const_a_v,
			  nonconst_x_v, const_y_v,
			  init_zero_x_vector,
			  update_y_vector,
			  numIter, false, true);


  }

  template <typename KernelHandle,
    typename lno_row_view_t_,
    typename lno_nnz_view_t_,
    typename scalar_nnz_view_t_,
    typename x_scalar_view_t,
    typename y_scalar_view_t>
  void backward_sweep_block_gauss_seidel_apply(KernelHandle *handle,
      typename KernelHandle::const_nnz_lno_t num_rows,
      typename KernelHandle::const_nnz_lno_t num_cols,
	  typename KernelHandle::const_nnz_lno_t block_size,

      lno_row_view_t_ row_map,
      lno_nnz_view_t_ entries,
      scalar_nnz_view_t_ values,
      x_scalar_view_t x_lhs_output_vec,
      y_scalar_view_t y_rhs_input_vec,
      bool init_zero_x_vector = false,
      bool update_y_vector = true,
      int numIter = 1){
	  handle->get_gs_handle()->set_block_size(block_size);
	  backward_sweep_gauss_seidel_apply(handle,num_rows,num_cols,
			  row_map,
		      entries,
		      values,
		      x_lhs_output_vec,
		      y_rhs_input_vec,
		      init_zero_x_vector,
		      update_y_vector,
		      numIter);
  }
}
}
#endif
