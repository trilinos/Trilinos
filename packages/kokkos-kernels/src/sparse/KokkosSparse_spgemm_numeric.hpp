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
#ifndef _KOKKOS_SPGEMM_NUMERIC_HPP
#define _KOKKOS_SPGEMM_NUMERIC_HPP

/*
#include "KokkosSparse_spgemm_cuSPARSE_impl.hpp"
#include "KokkosSparse_spgemm_CUSP_impl.hpp"
#include "KokkosSparse_spgemm_impl.hpp"
#include "KokkosSparse_spgemm_impl_seq.hpp"
#include "KokkosSparse_spgemm_mkl_impl.hpp"
#include "KokkosSparse_spgemm_mkl2phase_impl.hpp"
#include "KokkosSparse_spgemm_viennaCL_impl.hpp"
#include "KokkosKernels_Handle.hpp"
*/
#include "KokkosKernels_helpers.hpp"
#include "KokkosSparse_spgemm_numeric_spec.hpp"


namespace KokkosSparse{

namespace Experimental{


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
void spgemm_numeric(
    KernelHandle *handle,
    typename KernelHandle::const_nnz_lno_t m,
    typename KernelHandle::const_nnz_lno_t n,
    typename KernelHandle::const_nnz_lno_t k,
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



  static_assert (std::is_same<typename clno_row_view_t_::value_type,
      typename clno_row_view_t_::non_const_value_type>::value,
      "KokkosSparse::spgemm_symbolic: Output matrix rowmap must be non-const.");

  static_assert (std::is_same<typename clno_nnz_view_t_::value_type,
      typename clno_nnz_view_t_::non_const_value_type>::value,
      "KokkosSparse::spgemm_symbolic: Output matrix entriesView must be non-const.");

  static_assert (std::is_same<typename cscalar_nnz_view_t_::value_type,
      typename cscalar_nnz_view_t_::non_const_value_type>::value,
      "KokkosSparse::spgemm_symbolic: Output matrix scalar view must be non-const.");

  static_assert (std::is_same<typename KernelHandle::const_size_type,
      typename alno_row_view_t_::const_value_type>::value,
      "KokkosSparse::spgemm_symbolic: Size type of left handside matrix should be same as kernelHandle sizetype.");

  static_assert (std::is_same<typename KernelHandle::const_size_type,
      typename blno_row_view_t_::const_value_type>::value,
      "KokkosSparse::spgemm_symbolic: Size type of right handside matrix should be same as kernelHandle sizetype.");

  static_assert (std::is_same<typename KernelHandle::const_size_type,
      typename clno_row_view_t_::const_value_type>::value,
      "KokkosSparse::spgemm_symbolic: Size type of output matrix should be same as kernelHandle sizetype.");


  static_assert (std::is_same<typename KernelHandle::const_nnz_lno_t,
      typename alno_nnz_view_t_::const_value_type>::value,
      "KokkosSparse::spgemm_symbolic: lno type of left handside matrix should be same as kernelHandle lno_t.");

  static_assert (std::is_same<typename KernelHandle::const_nnz_lno_t,
      typename blno_nnz_view_t_::const_value_type>::value,
      "KokkosSparse::spgemm_symbolic: lno type of right handside matrix should be same as kernelHandle lno_t.");

  static_assert (std::is_same<typename KernelHandle::const_nnz_lno_t,
      typename clno_nnz_view_t_::const_value_type>::value,
      "KokkosSparse::spgemm_symbolic: lno type of output matrix should be same as kernelHandle lno_t.");

  static_assert (std::is_same<typename KernelHandle::const_nnz_scalar_t,
      typename ascalar_nnz_view_t_::const_value_type>::value,
      "KokkosSparse::spgemm_symbolic: scalar type of left handside matrix should be same as kernelHandle scalar.");

  static_assert (std::is_same<typename KernelHandle::const_nnz_scalar_t,
      typename bscalar_nnz_view_t_::const_value_type>::value,
      "KokkosSparse::spgemm_symbolic: scalar type of right handside matrix should be same as kernelHandle scalar.");

  static_assert (std::is_same<typename KernelHandle::const_nnz_scalar_t,
      typename cscalar_nnz_view_t_::const_value_type>::value,
      "KokkosSparse::spgemm_symbolic: scalar type of output matrix should be same as kernelHandle scalar.");


  if (transposeA || transposeB){
    throw std::runtime_error ("SpGEMM is not implemented for Transposes yet. "
        "If you need this case please let kokkos-kernels developers know.\n");
  }


    typedef typename KernelHandle::const_size_type c_size_t;
    typedef typename KernelHandle::const_nnz_lno_t c_lno_t;
    typedef typename KernelHandle::const_nnz_scalar_t c_scalar_t;

    typedef typename KernelHandle::HandleExecSpace c_exec_t;
    typedef typename KernelHandle::HandleTempMemorySpace c_temp_t;
    typedef typename KernelHandle::HandlePersistentMemorySpace c_persist_t;
    typedef typename Kokkos::Device<c_exec_t, c_temp_t> UniformDevice_t;

    typedef typename  KokkosKernels::Experimental::KokkosKernelsHandle<c_size_t, c_lno_t, c_scalar_t, c_exec_t, c_temp_t, c_persist_t> const_handle_type;


//  /const_handle_type tmp_handle = *handle;
  const_handle_type tmp_handle (*handle);

  typedef Kokkos::View<
      typename alno_row_view_t_::const_value_type*,
      typename KokkosKernels::Impl::GetUnifiedLayout<alno_row_view_t_>::array_layout,
      UniformDevice_t, //typename alno_row_view_t_::device_type,
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > Internal_alno_row_view_t_;

  typedef Kokkos::View<
      typename alno_nnz_view_t_::const_value_type*,
      typename KokkosKernels::Impl::GetUnifiedLayout<alno_nnz_view_t_>::array_layout,
      UniformDevice_t, // typename alno_nnz_view_t_::device_type,
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > Internal_alno_nnz_view_t_;


  typedef Kokkos::View<
      typename ascalar_nnz_view_t_::const_value_type*,
      typename KokkosKernels::Impl::GetUnifiedLayout<ascalar_nnz_view_t_>::array_layout,
      UniformDevice_t, //       typename ascalar_nnz_view_t_::device_type,
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > Internal_ascalar_nnz_view_t_;


  typedef Kokkos::View<
      typename blno_row_view_t_::const_value_type*,
      typename KokkosKernels::Impl::GetUnifiedLayout<blno_row_view_t_>::array_layout,
      UniformDevice_t, //       typename blno_row_view_t_::device_type,
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > Internal_blno_row_view_t_;


  typedef Kokkos::View<
      typename blno_nnz_view_t_::const_value_type*,
      typename KokkosKernels::Impl::GetUnifiedLayout<blno_nnz_view_t_>::array_layout,
      UniformDevice_t, //       typename blno_nnz_view_t_::device_type,
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > Internal_blno_nnz_view_t_;

  typedef Kokkos::View<
      typename bscalar_nnz_view_t_::const_value_type*,
      typename KokkosKernels::Impl::GetUnifiedLayout<bscalar_nnz_view_t_>::array_layout,
      UniformDevice_t, //       typename bscalar_nnz_view_t_::device_type,
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > Internal_bscalar_nnz_view_t_;

  //static assert clno_row_view_t_ cannot be const type.
  typedef Kokkos::View<
      typename clno_row_view_t_::non_const_value_type*,
      typename KokkosKernels::Impl::GetUnifiedLayout<clno_row_view_t_>::array_layout,
      UniformDevice_t, //       typename clno_row_view_t_::device_type,
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > Internal_clno_row_view_t_;

  typedef Kokkos::View<
      typename clno_nnz_view_t_::non_const_value_type*,
      typename KokkosKernels::Impl::GetUnifiedLayout<clno_nnz_view_t_>::array_layout,
      UniformDevice_t, //       typename clno_nnz_view_t_::device_type,
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > Internal_clno_nnz_view_t_;

  typedef Kokkos::View<
      typename cscalar_nnz_view_t_::non_const_value_type*,
      typename KokkosKernels::Impl::GetUnifiedLayout<cscalar_nnz_view_t_>::array_layout,
      UniformDevice_t, //       typename cscalar_nnz_view_t_::device_type,
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > Internal_cscalar_nnz_view_t_;

  //const_handle_type *const_handle = (const_handle_type *)handle;
  /*
  Internal_alno_row_view_t_ const_a_r  = row_mapA;
  Internal_alno_nnz_view_t_ const_a_l  = entriesA;
  Internal_ascalar_nnz_view_t_ const_a_s = valuesA;
  Internal_blno_row_view_t_ const_b_r  = row_mapB;
  Internal_blno_nnz_view_t_ const_b_l  = entriesB;
  Internal_bscalar_nnz_view_t_ const_b_s = valuesB;
  Internal_clno_row_view_t_ nonconst_c_r  = row_mapC;
  Internal_clno_nnz_view_t_ nonconst_c_l  = entriesC;
  Internal_cscalar_nnz_view_t_ nonconst_c_s = valuesC;
  */

  Internal_alno_row_view_t_ const_a_r (row_mapA.data(), row_mapA.extent(0));
  Internal_alno_nnz_view_t_ const_a_l (entriesA.data(), entriesA.extent(0));
  Internal_ascalar_nnz_view_t_ const_a_s (valuesA.data(), valuesA.extent(0));
  Internal_blno_row_view_t_ const_b_r (row_mapB.data(), row_mapB.extent(0));
  Internal_blno_nnz_view_t_ const_b_l  (entriesB.data(), entriesB.extent(0));
  Internal_bscalar_nnz_view_t_ const_b_s ( valuesB.data(), valuesB.extent(0));
  Internal_clno_row_view_t_ nonconst_c_r  ( row_mapC.data(), row_mapC.extent(0));
  Internal_clno_nnz_view_t_ nonconst_c_l  ( entriesC.data(), entriesC.extent(0));
  Internal_cscalar_nnz_view_t_ nonconst_c_s ( valuesC.data(), valuesC.extent(0));


  KokkosSparse::Impl::SPGEMM_NUMERIC<
  const_handle_type, //KernelHandle,
  Internal_alno_row_view_t_, Internal_alno_nnz_view_t_, Internal_ascalar_nnz_view_t_,
  Internal_blno_row_view_t_, Internal_blno_nnz_view_t_, Internal_bscalar_nnz_view_t_,
  Internal_clno_row_view_t_, Internal_clno_nnz_view_t_, Internal_cscalar_nnz_view_t_>::
  spgemm_numeric(
      &tmp_handle, //handle,
      m,
      n,
      k,
      const_a_r,
      const_a_l,
      const_a_s,
      transposeA,
      const_b_r,
      const_b_l,
      const_b_s,
      transposeB,
      nonconst_c_r,
      nonconst_c_l,
      nonconst_c_s);
}


}
}
#endif
