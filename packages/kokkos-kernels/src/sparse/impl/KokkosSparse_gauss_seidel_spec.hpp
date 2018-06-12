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
#ifndef KOKKOSSPARSE_IMPL_GAUSS_SEIDEL_SPEC_HPP_
#define KOKKOSSPARSE_IMPL_GAUSS_SEIDEL_SPEC_HPP_

#include <KokkosKernels_config.h>

#include <Kokkos_Core.hpp>
//#include <Kokkos_ArithTraits.hpp>
#include "KokkosKernels_Handle.hpp"
// Include the actual functors
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include "KokkosSparse_gauss_seidel_impl.hpp"
#endif

namespace KokkosSparse {
  namespace Impl {
    // Specialization struct which defines whether a specialization exists
    template<class KernelHandle, class a_size_view_t_, class a_lno_view_t>
    struct gauss_seidel_symbolic_eti_spec_avail {
      enum : bool { value = false };
    };
    template<class KernelHandle, class a_size_view_t_, class a_lno_view_t, class a_scalar_view_t>
    struct gauss_seidel_numeric_eti_spec_avail {
      enum : bool { value = false };
    };
    template<class KernelHandle, class a_size_view_t_, class a_lno_view_t, class a_scalar_view_t, class x_scalar_view_t, class y_scalar_view_t>
    struct gauss_seidel_apply_eti_spec_avail {
      enum : bool { value = false };
    };
  }
}


#define KOKKOSSPARSE_GAUSS_SEIDEL_SYMBOLIC_ETI_SPEC_AVAIL( SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE, SLOW_MEM_SPACE ) \
  template<>                                                            \
  struct gauss_seidel_symbolic_eti_spec_avail<                          \
                                               KokkosKernels::Experimental::KokkosKernelsHandle< \
                                                                                                const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE, \
                                                                                                EXEC_SPACE_TYPE, MEM_SPACE_TYPE, SLOW_MEM_SPACE> , \
                                               Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, \
                                                            Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                                            Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                                               Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, \
                                                            Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                                            Kokkos::MemoryTraits<Kokkos::Unmanaged> > > \
  { enum : bool { value = true }; };

#define KOKKOSSPARSE_GAUSS_SEIDEL_NUMERIC_ETI_SPEC_AVAIL( SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE, SLOW_MEM_SPACE ) \
  template<>                                                            \
  struct gauss_seidel_numeric_eti_spec_avail<                           \
                                              KokkosKernels::Experimental::KokkosKernelsHandle< \
                                                                                               const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE, \
                                                                                               EXEC_SPACE_TYPE, MEM_SPACE_TYPE, SLOW_MEM_SPACE> , \
                                              Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, \
                                                           Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                                              Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, \
                                                           Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                                              Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE, \
                                                           Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                                           Kokkos::MemoryTraits<Kokkos::Unmanaged> > > \
  { enum : bool { value = true }; };

#define KOKKOSSPARSE_GAUSS_SEIDEL_APPLY_ETI_SPEC_AVAIL( SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE, SLOW_MEM_SPACE ) \
  template<>                                                            \
  struct gauss_seidel_apply_eti_spec_avail<                             \
                                            KokkosKernels::Experimental::KokkosKernelsHandle< \
                                                                                             const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE, \
                                                                                             EXEC_SPACE_TYPE, MEM_SPACE_TYPE, SLOW_MEM_SPACE> , \
                                            Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, \
                                                         Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                                            Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, \
                                                         Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                                            Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE, \
                                                         Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                                            Kokkos::View< SCALAR_TYPE *, LAYOUT_TYPE, \
                                                          Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                                            Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE, \
                                                         Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                                         Kokkos::MemoryTraits<Kokkos::Unmanaged> > > \
  { enum : bool { value = true }; };

// Include the actual specialization declarations
#include<KokkosSparse_gauss_seidel_tpl_spec_avail.hpp>
#include<generated_specializations_hpp/KokkosSparse_gauss_seidel_symbolic_eti_spec_avail.hpp>
#include<generated_specializations_hpp/KokkosSparse_gauss_seidel_numeric_eti_spec_avail.hpp>
#include<generated_specializations_hpp/KokkosSparse_gauss_seidel_apply_eti_spec_avail.hpp>

namespace KokkosSparse {
  namespace Impl {



    template<
      class KernelHandle,
      class a_size_view_t_, class a_lno_view_t,
      bool tpl_spec_avail =
      gauss_seidel_symbolic_tpl_spec_avail<
        KernelHandle, a_size_view_t_,  a_lno_view_t>::value,
      bool eti_spec_avail =
      gauss_seidel_symbolic_eti_spec_avail<
        KernelHandle,
        a_size_view_t_,  a_lno_view_t>::value >
    struct GAUSS_SEIDEL_SYMBOLIC{
      static void
      gauss_seidel_symbolic (
                             KernelHandle *handle,
                             typename KernelHandle::const_nnz_lno_t num_rows,
                             typename KernelHandle::const_nnz_lno_t num_cols,
                             a_size_view_t_ row_map,
                             a_lno_view_t entries,
                             bool is_graph_symmetric);
    };

    template<
      class KernelHandle,
      class a_size_view_t_, class a_lno_view_t, class  a_scalar_view_t,
      bool tpl_spec_avail =
      gauss_seidel_numeric_tpl_spec_avail<
        KernelHandle, a_size_view_t_,  a_lno_view_t, a_scalar_view_t>::value,
      bool eti_spec_avail =
      gauss_seidel_numeric_eti_spec_avail<
        KernelHandle,
        a_size_view_t_,  a_lno_view_t, a_scalar_view_t>::value >
    struct GAUSS_SEIDEL_NUMERIC{
      static void
      gauss_seidel_numeric (KernelHandle *handle,
                            typename KernelHandle::const_nnz_lno_t num_rows,
                            typename KernelHandle::const_nnz_lno_t num_cols,
                            a_size_view_t_ row_map,
                            a_lno_view_t entries,
                            a_scalar_view_t values,
                            bool is_graph_symmetric
                            );

      static void
      gauss_seidel_numeric (KernelHandle *handle,
                            typename KernelHandle::const_nnz_lno_t num_rows,
                            typename KernelHandle::const_nnz_lno_t num_cols,
                            a_size_view_t_ row_map,
                            a_lno_view_t entries,
                            a_scalar_view_t values,
                            a_scalar_view_t given_inverse_diagonal,
                            bool is_graph_symmetric
                            );
    };


    template<
      class KernelHandle,
      class a_size_view_t_, class a_lno_view_t, class  a_scalar_view_t, class x_scalar_view_t, class y_scalar_view_t,
      bool tpl_spec_avail =
      gauss_seidel_apply_tpl_spec_avail<
        KernelHandle, a_size_view_t_,  a_lno_view_t, a_scalar_view_t,x_scalar_view_t, y_scalar_view_t>::value,
      bool eti_spec_avail =
      gauss_seidel_apply_eti_spec_avail<
        KernelHandle,
        a_size_view_t_,  a_lno_view_t, a_scalar_view_t,x_scalar_view_t, y_scalar_view_t>::value >
    struct GAUSS_SEIDEL_APPLY{
      static void
      gauss_seidel_apply (
                          KernelHandle *handle,
                          typename KernelHandle::const_nnz_lno_t num_rows,
                          typename KernelHandle::const_nnz_lno_t num_cols,
                          a_size_view_t_ row_map,
                          a_lno_view_t entries,
                          a_scalar_view_t values,
                          x_scalar_view_t x_lhs_output_vec,
                          y_scalar_view_t y_rhs_input_vec,
                          bool init_zero_x_vector,
                          bool update_y_vector,
                          typename KernelHandle::nnz_scalar_t omega, int numIter, bool apply_forward, bool apply_backward);
    };


#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY


    template<class KernelHandle, class a_size_view_t_, class a_lno_view_t>
    struct GAUSS_SEIDEL_SYMBOLIC<KernelHandle,
                                 a_size_view_t_,  a_lno_view_t,
                                 false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY>{

      static void
      gauss_seidel_symbolic (
                             KernelHandle *handle,
                             typename KernelHandle::const_nnz_lno_t num_rows,
                             typename KernelHandle::const_nnz_lno_t num_cols,
                             a_size_view_t_ row_map,
                             a_lno_view_t entries,
                             bool is_graph_symmetric){

        typedef typename Impl::GaussSeidel<KernelHandle, a_size_view_t_,
                                           a_lno_view_t, typename KernelHandle::in_scalar_nnz_view_t> SGS;
        SGS sgs(handle,num_rows, num_cols, row_map, entries, is_graph_symmetric);
        sgs.initialize_symbolic();
      }
    };

    template<class KernelHandle, class a_size_view_t_, class a_lno_view_t, class a_scalar_view_t>
    struct GAUSS_SEIDEL_NUMERIC<KernelHandle,
                                a_size_view_t_,  a_lno_view_t, a_scalar_view_t,
                                false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY>{

      static void
      gauss_seidel_numeric(KernelHandle *handle,
                           typename KernelHandle::const_nnz_lno_t num_rows,
                           typename KernelHandle::const_nnz_lno_t num_cols,
                           a_size_view_t_ row_map,
                           a_lno_view_t entries,
                           a_scalar_view_t values,
                           bool is_graph_symmetric
                           ){
        typedef typename Impl::GaussSeidel
          <KernelHandle,a_size_view_t_,
           a_lno_view_t,a_scalar_view_t> SGS;
        SGS sgs(handle, num_rows, num_cols, row_map, entries, values, is_graph_symmetric);
        sgs.initialize_numeric();
      }

      static void
      gauss_seidel_numeric(KernelHandle *handle,
                           typename KernelHandle::const_nnz_lno_t num_rows,
                           typename KernelHandle::const_nnz_lno_t num_cols,
                           a_size_view_t_ row_map,
                           a_lno_view_t entries,
                           a_scalar_view_t values,
                           a_scalar_view_t given_inverse_diagonal,
                           bool is_graph_symmetric
                           ){
        typedef typename Impl::GaussSeidel
          <KernelHandle,a_size_view_t_,
           a_lno_view_t,a_scalar_view_t> SGS;
        SGS sgs(handle, num_rows, num_cols, row_map, entries, values, given_inverse_diagonal, is_graph_symmetric);
        sgs.initialize_numeric();
      }
    };

    template<class KernelHandle, class a_size_view_t_, class a_lno_view_t, class a_scalar_view_t, class x_scalar_view_t, class y_scalar_view_t>
    struct GAUSS_SEIDEL_APPLY<KernelHandle,
                              a_size_view_t_,  a_lno_view_t, a_scalar_view_t,x_scalar_view_t, y_scalar_view_t,
                              false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY>{

      static void
      gauss_seidel_apply(
                         KernelHandle *handle,
                         typename KernelHandle::const_nnz_lno_t num_rows,
                         typename KernelHandle::const_nnz_lno_t num_cols,
                         a_size_view_t_ row_map,
                         a_lno_view_t entries,
                         a_scalar_view_t values,
                         x_scalar_view_t x_lhs_output_vec,
                         y_scalar_view_t y_rhs_input_vec,
                         bool init_zero_x_vector,
                         bool update_y_vector,
                         typename KernelHandle::nnz_scalar_t omega, int numIter, bool apply_forward, bool apply_backward){

        typedef typename Impl::GaussSeidel <KernelHandle,
                                            a_size_view_t_, a_lno_view_t,a_scalar_view_t > SGS;
        SGS sgs(handle, num_rows, num_cols, row_map, entries, values);
        sgs.apply(
                  x_lhs_output_vec,
                  y_rhs_input_vec,
                  init_zero_x_vector,
                  numIter,
                  omega,
                  apply_forward,
                  apply_backward, update_y_vector);
      }
    };
#endif



  }
}

#define KOKKOSSPARSE_GAUSS_SEIDEL_SYMBOLIC_ETI_SPEC_DECL( SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE, SLOW_MEM_SPACE ) \
  extern template struct                                                \
  GAUSS_SEIDEL_SYMBOLIC<                                                \
                         KokkosKernels::Experimental::KokkosKernelsHandle< \
                                                                          const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE, \
                                                                          EXEC_SPACE_TYPE, MEM_SPACE_TYPE, SLOW_MEM_SPACE> , \
                         Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, \
                                      Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                      Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                         Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, \
                                      Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                      Kokkos::MemoryTraits<Kokkos::Unmanaged> > , \
                         false, true >;


#define KOKKOSSPARSE_GAUSS_SEIDEL_SYMBOLIC_ETI_SPEC_INST( SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE, SLOW_MEM_SPACE) \
  template struct                                                       \
  GAUSS_SEIDEL_SYMBOLIC<                                                \
                         KokkosKernels::Experimental::KokkosKernelsHandle< \
                                                                          const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE, \
                                                                          EXEC_SPACE_TYPE, MEM_SPACE_TYPE, SLOW_MEM_SPACE> , \
                         Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, \
                                      Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                      Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                         Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, \
                                      Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                      Kokkos::MemoryTraits<Kokkos::Unmanaged> > , \
                         false, true > ;

#define KOKKOSSPARSE_GAUSS_SEIDEL_NUMERIC_ETI_SPEC_DECL( SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE,SLOW_MEM_SPACE ) \
  extern template struct                                                \
  GAUSS_SEIDEL_NUMERIC<                                                 \
                        KokkosKernels::Experimental::KokkosKernelsHandle< \
                                                                         const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE, \
                                                                         EXEC_SPACE_TYPE, MEM_SPACE_TYPE, SLOW_MEM_SPACE> , \
                        Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE,  \
                                     Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                        Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, \
                                     Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                        Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE,  \
                                     Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged> > , \
                        false, true >;


#define KOKKOSSPARSE_GAUSS_SEIDEL_NUMERIC_ETI_SPEC_INST( SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE, SLOW_MEM_SPACE) \
  template struct                                                       \
  GAUSS_SEIDEL_NUMERIC<                                                 \
                        KokkosKernels::Experimental::KokkosKernelsHandle< \
                                                                         const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE, \
                                                                         EXEC_SPACE_TYPE, MEM_SPACE_TYPE, SLOW_MEM_SPACE> , \
                        Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE,  \
                                     Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                        Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, \
                                     Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                        Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE,  \
                                     Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged> > , \
                        false, true > ;

#define KOKKOSSPARSE_GAUSS_SEIDEL_APPLY_ETI_SPEC_DECL( SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE, SLOW_MEM_SPACE) \
  extern template struct                                                \
  GAUSS_SEIDEL_APPLY<                                                   \
                      KokkosKernels::Experimental::KokkosKernelsHandle< \
                                                                       const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE, \
                                                                       EXEC_SPACE_TYPE, MEM_SPACE_TYPE, SLOW_MEM_SPACE> , \
                      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE,    \
                                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE,   \
                                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE,    \
                                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                      Kokkos::View<SCALAR_TYPE *, LAYOUT_TYPE,          \
                                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE,    \
                                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged> > , \
                      false, true >;


#define KOKKOSSPARSE_GAUSS_SEIDEL_APPLY_ETI_SPEC_INST( SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE, SLOW_MEM_SPACE) \
  template struct                                                       \
  GAUSS_SEIDEL_APPLY<                                                   \
                      KokkosKernels::Experimental::KokkosKernelsHandle< \
                                                                       const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE, \
                                                                       EXEC_SPACE_TYPE, MEM_SPACE_TYPE, SLOW_MEM_SPACE> , \
                      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE,    \
                                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE,   \
                                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE,    \
                                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                      Kokkos::View<SCALAR_TYPE *, LAYOUT_TYPE,          \
                                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE,    \
                                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged> > , \
                      false, true > ;

#include<KokkosSparse_gauss_seidel_tpl_spec_decl.hpp>
#include<generated_specializations_hpp/KokkosSparse_gauss_seidel_symbolic_eti_spec_decl.hpp>
#include<generated_specializations_hpp/KokkosSparse_gauss_seidel_numeric_eti_spec_decl.hpp>
#include<generated_specializations_hpp/KokkosSparse_gauss_seidel_apply_eti_spec_decl.hpp>



#endif // KOKKOS_BLAS1_MV_IMPL_DOT_HPP_
