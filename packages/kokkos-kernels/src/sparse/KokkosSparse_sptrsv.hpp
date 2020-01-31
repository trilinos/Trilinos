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

/// \file KokkosSparse_sptrsv.hpp
/// \brief Parallel sparse triangular solve
///
/// This file provides KokkosSparse::sptrsv.  This function performs a
/// local (no MPI) sparse triangular solve on matrices stored in
/// compressed row sparse ("Crs") format.

#ifndef KOKKOSSPARSE_SPTRSV_HPP_
#define KOKKOSSPARSE_SPTRSV_HPP_

#include <type_traits>

//#include "KokkosSparse_sptrsv_handle.hpp"
#include "KokkosKernels_helpers.hpp"
#include "KokkosSparse_sptrsv_symbolic_spec.hpp"
#include "KokkosSparse_sptrsv_solve_spec.hpp"

#ifdef KOKKOSKERNELS_ENABLE_TPL_SUPERLU
 #include "slu_ddefs.h"
 #include "KokkosSparse_sptrsv_superlu.hpp"
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_CHOLMOD
 #include "cholmod.h"
 #include "KokkosSparse_sptrsv_cholmod.hpp"
#endif

#include "KokkosSparse_sptrsv_cuSPARSE_impl.hpp"

namespace KokkosSparse {
namespace Experimental {

#define KOKKOSKERNELS_SPTRSV_SAME_TYPE(A, B) std::is_same<typename std::remove_const<A>::type, typename std::remove_const<B>::type>::value

  template <typename KernelHandle,
            typename lno_row_view_t_,
            typename lno_nnz_view_t_>
  void sptrsv_symbolic(
      KernelHandle *handle, 
      lno_row_view_t_ rowmap,
      lno_nnz_view_t_ entries)
  {
    typedef typename KernelHandle::size_type size_type;
    typedef typename KernelHandle::nnz_lno_t ordinal_type;

    static_assert(KOKKOSKERNELS_SPTRSV_SAME_TYPE(typename lno_row_view_t_::non_const_value_type, size_type),
        "sptrsv_symbolic: A size_type must match KernelHandle size_type (const doesn't matter)");

    static_assert(KOKKOSKERNELS_SPTRSV_SAME_TYPE(typename lno_nnz_view_t_::non_const_value_type, ordinal_type),
        "sptrsv_symbolic: A entry type must match KernelHandle entry type (aka nnz_lno_t, and const doesn't matter)");


    typedef typename KernelHandle::const_size_type c_size_t;
    typedef typename KernelHandle::const_nnz_lno_t c_lno_t;
    typedef typename KernelHandle::const_nnz_scalar_t c_scalar_t;

    typedef typename KernelHandle::HandleExecSpace c_exec_t;
    typedef typename KernelHandle::HandleTempMemorySpace c_temp_t;
    typedef typename KernelHandle::HandlePersistentMemorySpace c_persist_t;

    typedef typename  KokkosKernels::Experimental::KokkosKernelsHandle<c_size_t, c_lno_t, c_scalar_t, c_exec_t, c_temp_t, c_persist_t> const_handle_type;
    const_handle_type tmp_handle (*handle);

    typedef Kokkos::View<
          typename lno_row_view_t_::const_value_type*,
          typename KokkosKernels::Impl::GetUnifiedLayout<lno_row_view_t_>::array_layout,
          typename lno_row_view_t_::device_type,
          Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> > RowMap_Internal;

    typedef Kokkos::View<
          typename lno_nnz_view_t_::const_value_type*,
          typename KokkosKernels::Impl::GetUnifiedLayout<lno_nnz_view_t_>::array_layout,
          typename lno_nnz_view_t_::device_type,
          Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> > Entries_Internal;


    RowMap_Internal rowmap_i = rowmap;
    Entries_Internal entries_i = entries;

    KokkosSparse::Impl::SPTRSV_SYMBOLIC<const_handle_type, RowMap_Internal, Entries_Internal>::sptrsv_symbolic (&tmp_handle, rowmap_i, entries_i);

  } // sptrsv_symbolic

  template <typename KernelHandle,
            typename lno_row_view_t_,
            typename lno_nnz_view_t_,
            typename scalar_nnz_view_t_>
  void sptrsv_symbolic(
      KernelHandle *handle, 
      lno_row_view_t_ rowmap,
      lno_nnz_view_t_ entries,
      scalar_nnz_view_t_ values)
  {
    typedef typename KernelHandle::size_type size_type;
    typedef typename KernelHandle::nnz_lno_t ordinal_type;
    typedef typename KernelHandle::nnz_scalar_t scalar_type;

    static_assert(KOKKOSKERNELS_SPTRSV_SAME_TYPE(typename lno_row_view_t_::non_const_value_type, size_type),
        "sptrsv_symbolic: A size_type must match KernelHandle size_type (const doesn't matter)");

    static_assert(KOKKOSKERNELS_SPTRSV_SAME_TYPE(typename lno_nnz_view_t_::non_const_value_type, ordinal_type),
        "sptrsv_symbolic: A entry type must match KernelHandle entry type (aka nnz_lno_t, and const doesn't matter)");

    static_assert(KOKKOSKERNELS_SPTRSV_SAME_TYPE(typename scalar_nnz_view_t_::value_type, scalar_type),
        "sptrsv_symbolic: A scalar type must match KernelHandle entry type (aka nnz_lno_t, and const doesn't matter)");

    typedef typename KernelHandle::const_size_type c_size_t;
    typedef typename KernelHandle::const_nnz_lno_t c_lno_t;
    typedef typename KernelHandle::const_nnz_scalar_t c_scalar_t;

    typedef typename KernelHandle::HandleExecSpace c_exec_t;
    typedef typename KernelHandle::HandleTempMemorySpace c_temp_t;
    typedef typename KernelHandle::HandlePersistentMemorySpace c_persist_t;

    typedef typename  KokkosKernels::Experimental::KokkosKernelsHandle<c_size_t, c_lno_t, c_scalar_t, c_exec_t, c_temp_t, c_persist_t> const_handle_type;
    const_handle_type tmp_handle (*handle);

    typedef Kokkos::View<
          typename lno_row_view_t_::const_value_type*,
          typename KokkosKernels::Impl::GetUnifiedLayout<lno_row_view_t_>::array_layout,
          typename lno_row_view_t_::device_type,
          Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> > RowMap_Internal;

    typedef Kokkos::View<
          typename lno_nnz_view_t_::const_value_type*,
          typename KokkosKernels::Impl::GetUnifiedLayout<lno_nnz_view_t_>::array_layout,
          typename lno_nnz_view_t_::device_type,
          Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> > Entries_Internal;

    typedef Kokkos::View<
          typename scalar_nnz_view_t_::const_value_type*,
          typename KokkosKernels::Impl::GetUnifiedLayout<scalar_nnz_view_t_>::array_layout,
          typename scalar_nnz_view_t_::device_type,
          Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> > Values_Internal;

    auto sptrsv_handle = handle->get_sptrsv_handle();
    if (sptrsv_handle->get_algorithm() == KokkosSparse::Experimental::SPTRSVAlgorithm::SPTRSV_CUSPARSE) {
      RowMap_Internal rowmap_i = rowmap;
      Entries_Internal entries_i = entries;
      Values_Internal values_i = values;


      typedef typename KernelHandle::SPTRSVHandleType sptrsvHandleType;
      sptrsvHandleType *sh = handle->get_sptrsv_handle();
      auto nrows = sh->get_nrows();

      KokkosSparse::Impl::sptrsvcuSPARSE_symbolic
      < sptrsvHandleType,
        RowMap_Internal,
        Entries_Internal,
        Values_Internal >
        (sh, nrows, rowmap_i, entries_i, values_i, false);

    }
    else {
      KokkosSparse::Experimental::sptrsv_symbolic (handle, rowmap, entries);
    }

  } // sptrsv_symbolic

  template <typename KernelHandle,
            typename lno_row_view_t_,
            typename lno_nnz_view_t_,
            typename scalar_nnz_view_t_,
            class BType,
            class XType>
  void sptrsv_solve(
      KernelHandle *handle, 
      lno_row_view_t_ rowmap,
      lno_nnz_view_t_ entries,
      scalar_nnz_view_t_ values,
      BType b,
      XType x)
  {
    typedef typename KernelHandle::size_type size_type;
    typedef typename KernelHandle::nnz_lno_t ordinal_type;
    typedef typename KernelHandle::nnz_scalar_t scalar_type;
    
    static_assert(KOKKOSKERNELS_SPTRSV_SAME_TYPE(typename lno_row_view_t_::non_const_value_type, size_type),
        "sptrsv_solve: A size_type must match KernelHandle size_type (const doesn't matter)");
    static_assert(KOKKOSKERNELS_SPTRSV_SAME_TYPE(typename lno_nnz_view_t_::non_const_value_type, ordinal_type),
        "sptrsv_solve: A entry type must match KernelHandle entry type (aka nnz_lno_t, and const doesn't matter)");
    static_assert(KOKKOSKERNELS_SPTRSV_SAME_TYPE(typename scalar_nnz_view_t_::value_type, scalar_type),
        "sptrsv_solve: A scalar type must match KernelHandle entry type (aka nnz_lno_t, and const doesn't matter)");

    static_assert (Kokkos::Impl::is_view<BType>::value,
        "sptrsv: b is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<XType>::value,
        "sptrsv: x is not a Kokkos::View.");
    static_assert ((int) BType::rank == (int) XType::rank,
        "sptrsv: The ranks of b and x do not match.");
    static_assert (BType::rank == 1,
        "sptrsv: b and x must both either have rank 1.");
    static_assert (std::is_same<typename XType::value_type,
        typename XType::non_const_value_type>::value,
        "sptrsv: The output x must be nonconst.");
    static_assert (std::is_same<typename BType::device_type, typename XType::device_type>::value,
        "sptrsv: Views BType and XType have different device_types.");
    static_assert (std::is_same<typename BType::device_type::execution_space, typename KernelHandle::SPTRSVHandleType::execution_space>::value,
        "sptrsv: KernelHandle and Views have different execution spaces.");
    static_assert (std::is_same<typename lno_row_view_t_::device_type, typename lno_nnz_view_t_::device_type>::value,
        "sptrsv: rowmap and entries have different device types.");
    static_assert (std::is_same<typename lno_row_view_t_::device_type, typename scalar_nnz_view_t_::device_type>::value,
        "sptrsv: rowmap and values have different device types.");


    typedef typename KernelHandle::const_size_type c_size_t;
    typedef typename KernelHandle::const_nnz_lno_t c_lno_t;
    typedef typename KernelHandle::const_nnz_scalar_t c_scalar_t;

    typedef typename KernelHandle::HandleExecSpace c_exec_t;
    typedef typename KernelHandle::HandleTempMemorySpace c_temp_t;
    typedef typename KernelHandle::HandlePersistentMemorySpace c_persist_t;

    typedef typename  KokkosKernels::Experimental::KokkosKernelsHandle<c_size_t, c_lno_t, c_scalar_t, c_exec_t, c_temp_t, c_persist_t> const_handle_type;
    const_handle_type tmp_handle (*handle);

    typedef Kokkos::View<
          typename lno_row_view_t_::const_value_type*,
          typename KokkosKernels::Impl::GetUnifiedLayout<lno_row_view_t_>::array_layout,
          typename lno_row_view_t_::device_type,
          Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> > RowMap_Internal;

    typedef Kokkos::View<
          typename lno_nnz_view_t_::const_value_type*,
          typename KokkosKernels::Impl::GetUnifiedLayout<lno_nnz_view_t_>::array_layout,
          typename lno_nnz_view_t_::device_type,
          Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> > Entries_Internal;

    typedef Kokkos::View<
          typename scalar_nnz_view_t_::const_value_type*,
          typename KokkosKernels::Impl::GetUnifiedLayout<scalar_nnz_view_t_>::array_layout,
          typename scalar_nnz_view_t_::device_type,
          Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> > Values_Internal;


    typedef Kokkos::View<
          typename BType::const_value_type*,
          typename KokkosKernels::Impl::GetUnifiedLayout<BType>::array_layout,
          typename BType::device_type,
          Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> > BType_Internal;

    typedef Kokkos::View<
          typename XType::non_const_value_type*,
          typename KokkosKernels::Impl::GetUnifiedLayout<XType>::array_layout,
          typename XType::device_type,
          Kokkos::MemoryTraits<Kokkos::Unmanaged> > XType_Internal;


    RowMap_Internal rowmap_i = rowmap;
    Entries_Internal entries_i = entries;
    Values_Internal values_i = values;

    BType_Internal b_i = b;
    XType_Internal x_i = x;

    auto sptrsv_handle = handle->get_sptrsv_handle();
    if (sptrsv_handle->get_algorithm() == KokkosSparse::Experimental::SPTRSVAlgorithm::SPTRSV_CUSPARSE) {
      typedef typename KernelHandle::SPTRSVHandleType sptrsvHandleType;
      sptrsvHandleType *sh = handle->get_sptrsv_handle();
      auto nrows = sh->get_nrows();

      KokkosSparse::Impl::sptrsvcuSPARSE_solve
      < sptrsvHandleType,
        RowMap_Internal,
        Entries_Internal,
        Values_Internal,
        BType_Internal,
        XType_Internal >
        (sh, nrows, rowmap_i, entries_i, values_i, b_i, x_i, false);

    }
    else {
      KokkosSparse::Impl::SPTRSV_SOLVE<const_handle_type, RowMap_Internal, Entries_Internal, Values_Internal, BType_Internal, XType_Internal>::sptrsv_solve (&tmp_handle, rowmap_i, entries_i, values_i, b_i, x_i);
    }

  } // sptrsv_solve


#ifdef KOKKOSKERNELS_ENABLE_TPL_SUPERLU
  // ---------------------------------------------------------------------
  template <typename scalar_type,
            typename ordinal_type,
            typename size_type,
            typename KernelHandle,
            typename execution_space      = Kokkos::DefaultExecutionSpace,
            typename host_execution_space = Kokkos::DefaultHostExecutionSpace>
  void sptrsv_symbolic(
      KernelHandle *kernelHandleL,
      KernelHandle *kernelHandleU,
      SuperMatrix &L,
      SuperMatrix &U)
  {
    typedef KokkosSparse::CrsMatrix<scalar_type, ordinal_type, host_execution_space, void, size_type> host_crsmat_t;
    typedef KokkosSparse::CrsMatrix<scalar_type, ordinal_type,      execution_space, void, size_type> crsmat_t;

    typedef typename host_crsmat_t::StaticCrsGraphType host_graph_t;
    typedef typename      crsmat_t::StaticCrsGraphType      graph_t;

    Kokkos::Timer timer;
    Kokkos::Timer tic;

    // ===================================================================
    // load sptrsv-handles
    auto *handleL = kernelHandleL->get_sptrsv_handle ();
    auto *handleU = kernelHandleU->get_sptrsv_handle ();
    int *etree = handleL->get_etree ();

    // ===================================================================
    // load options
    bool col_majorL = handleL->is_column_major ();
    bool col_majorU = handleU->is_column_major ();
    bool merge = handleL->get_merge_supernodes ();
    bool UinCSC = handleU->is_column_major ();
    bool needEtree = (handleL->get_algorithm () == SPTRSVAlgorithm::SUPERNODAL_SPMV ||
                      handleL->get_algorithm () == SPTRSVAlgorithm::SUPERNODAL_ETREE);
    if (needEtree && etree == nullptr) {
      std::cout << std::endl
                << " ** etree needs to be set before calling sptrsv_symbolic with SuperLU **"
                << std::endl << std::endl;
      return;
    }

    // ===================================================================
    // read CrsGraph from SuperLU factor
    #ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
    std::cout << " > Read SuperLU factor into KokkosSparse::CrsMatrix (invert diagonal and copy to device)" << std::endl;
    if (merge) {
      std::cout << " > Merge supernodes" << std::endl;
    }
    #endif
    bool cusparse = false; // pad diagonal blocks with zeros
    host_graph_t graphL_host;
    host_graph_t graphU_host;
    graph_t graphL;
    graph_t graphU;

    tic.reset();
    graphL_host = read_superlu_graphL<host_graph_t> (cusparse, merge, &L);
    if (UinCSC) {
      graphU_host = read_superlu_graphU_CSC<host_graph_t> (&L, &U); 
    } else {
      graphU_host = read_superlu_graphU<host_graph_t> (&L, &U); 
    }
    #ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
    double time_seconds = tic.seconds ();
    std::cout << "   Conversion Time (from SuperLU to CSR): " << time_seconds << std::endl;
    #endif

    // ===================================================================
    // setup supnodal info
    SCformat *Lstore = (SCformat*)(L.Store);
    int nsuper = 1 + Lstore->nsuper;
    int *supercols = Lstore->sup_to_col;

    // > make a copy of supercols (merge needs both original and merged supercols)
    typename KernelHandle::SPTRSVHandleType::integer_view_host_t supercols_view ("supercols view", 1+nsuper);
    int *supercols_merged = supercols_view.data ();
    for (int i = 0; i <= nsuper; i++) {
      supercols_merged[i] = supercols[i];
    }
    if (merge) {
      // =================================================================
      // merge supernodes
      tic.reset ();
      int nsuper_merged = nsuper;
      #ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
      int nrows = L.nrow;
      check_supernode_sizes("Original L-structure", nrows, nsuper, supercols_merged, graphL_host);
      check_supernode_sizes("Original U-structure", nrows, nsuper, supercols_merged, graphU_host);
      #endif
      // etree will be updated
      merge_supernodal_graph (&nsuper_merged, supercols_merged,
                              col_majorL, graphL_host, col_majorU, graphU_host,
                              etree);

      // =================================================================
      // generate merged graph for L-solve
      #ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
      int nnzL = graphL_host.row_map (nrows);
      #endif
      int nnzL_merged;
      bool lower = true;
      handleL->set_original_graph_host (graphL_host); // save graph before merge
      graphL_host = generate_merged_supernodal_graph<host_graph_t> (lower, nsuper, supercols,
                                                                    nsuper_merged, supercols_merged,
                                                                    graphL_host, &nnzL_merged);
      #ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
      time_seconds = tic.seconds ();
      check_supernode_sizes ("After Merge", nrows, nsuper_merged, supercols_merged, graphL_host);
      std::cout << " for L factor:" << std::endl;
      std::cout << "   Merge Supernodes Time: " << time_seconds << std::endl;
      std::cout << "   Number of nonzeros   : " << nnzL << " -> " << nnzL_merged
                << " : " << double(nnzL_merged) / double(nnzL) << "x" << std::endl;
      #endif

      // =================================================================
      // generate merged graph for U-solve
      #ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
      int nnzU = graphU_host.row_map (nrows);
      #endif
      int nnzU_merged;
      lower = (UinCSC ? false : true);
      tic.reset ();
      handleU->set_original_graph_host (graphU_host); // save graph before merge
      graphU_host = generate_merged_supernodal_graph<host_graph_t> (lower, nsuper, supercols,
                                                                    nsuper_merged, supercols_merged,
                                                                    graphU_host, &nnzU_merged);
      #ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
      time_seconds = tic.seconds ();
      check_supernode_sizes("After Merge", nrows, nsuper_merged, supercols_merged, graphU_host);
      std::cout << " for U factor:" << std::endl;
      std::cout << "   Merge Supernodes Time: " << time_seconds << std::endl;
      std::cout << "   Number of nonzeros   : " << nnzU << " -> " << nnzU_merged
                << " : " << double(nnzU_merged) / double(nnzU) << "x" << std::endl;
      #endif

      // update the number of supernodes
      nsuper = nsuper_merged;
    }
    // replace the supernodal info with the merged ones
    supercols = supercols_merged;

    // ===================================================================
    // copy graph to device
    graphL = deep_copy_graph<host_graph_t, graph_t> (graphL_host);
    graphU = deep_copy_graph<host_graph_t, graph_t> (graphU_host);

    // ===================================================================
    // save the supernodal info in the handles for L/U solves
    handleL->set_supernodes (nsuper, supercols_view, etree);
    handleU->set_supernodes (nsuper, supercols_view, etree);

    if (handleL->get_algorithm () == SPTRSVAlgorithm::SUPERNODAL_DAG ||
        handleL->get_algorithm () == SPTRSVAlgorithm::SUPERNODAL_SPMV_DAG) {
      // generate supernodal graphs for DAG scheduling
      auto supL = generate_supernodal_graph<host_graph_t> (!col_majorL, graphL_host, nsuper, supercols);
      auto supU = generate_supernodal_graph<host_graph_t> ( col_majorU, graphU_host, nsuper, supercols);

      auto dagL = generate_supernodal_dag<host_graph_t> (nsuper, supL, supU);
      auto dagU = generate_supernodal_dag<host_graph_t> (nsuper, supU, supL);
      handleL->set_supernodal_dag (dagL);
      handleU->set_supernodal_dag (dagU);
    }

    // ===================================================================
    // do symbolic for L solve on the host
    auto row_mapL = graphL.row_map;
    auto entriesL = graphL.entries;
    tic.reset();
    #ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
    std::cout << std::endl;
    #endif
    sptrsv_symbolic (kernelHandleL, row_mapL, entriesL);
    #ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
    time_seconds = tic.seconds ();
    std::cout << " > Lower-TRI: " << std::endl;
    std::cout << "   Symbolic Time: " << time_seconds << std::endl;
    #endif

    // ===================================================================
    // do symbolic for U solve on the host
    auto row_mapU = graphU.row_map;
    auto entriesU = graphU.entries;
    tic.reset ();
    sptrsv_symbolic (kernelHandleU, row_mapU, entriesU);
    #ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
    time_seconds = tic.seconds ();
    std::cout << " > Upper-TRI: " << std::endl;
    std::cout << "   Symbolic Time: " << time_seconds << std::endl;
    #endif

    // ===================================================================
    // save options
    handleL->set_merge_supernodes (merge);
    handleU->set_merge_supernodes (merge);

    // ===================================================================
    // save graphs
    handleL->set_graph (graphL);
    handleU->set_graph (graphU);
    // graph on host (merged)
    handleL->set_graph_host (graphL_host);
    handleU->set_graph_host (graphU_host);

    // ===================================================================
    handleL->set_symbolic_complete ();
    handleU->set_symbolic_complete ();
    handleL->set_etree (etree);
    handleU->set_etree (etree);
    #ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
    time_seconds = tic.seconds ();
    std::cout << "   Total Symbolic Time: " << time_seconds << std::endl << std::endl;
    #endif
  }


  // ---------------------------------------------------------------------
  template <typename scalar_type,
            typename ordinal_type,
            typename size_type,
            typename KernelHandle,
            typename execution_space      = Kokkos::DefaultExecutionSpace,
            typename host_execution_space = Kokkos::DefaultHostExecutionSpace>
  void sptrsv_compute(
      KernelHandle *kernelHandleL,
      KernelHandle *kernelHandleU,
      SuperMatrix &L,
      SuperMatrix &U)
  {
    typedef KokkosSparse::CrsMatrix<scalar_type, ordinal_type, host_execution_space, void, size_type> host_crsmat_t;
    typedef KokkosSparse::CrsMatrix<scalar_type, ordinal_type,      execution_space, void, size_type> crsmat_t;

    Kokkos::Timer tic;
    Kokkos::Timer timer;
    // ===================================================================
    // load sptrsv-handles
    auto *handleL = kernelHandleL->get_sptrsv_handle ();
    auto *handleU = kernelHandleU->get_sptrsv_handle ();

    if (!(handleL->is_symbolic_complete()) ||
        !(handleU->is_symbolic_complete())) {
      std::cout << std::endl
                << " ** needs to call sptrsv_symbolic before calling sptrsv_numeric **"
                << std::endl << std::endl;
      return;
    }

    // ===================================================================
    // load options
    bool merge = handleL->get_merge_supernodes ();
    bool invert_offdiag = handleL->get_invert_offdiagonal ();
    bool UinCSC = handleU->is_column_major ();
    bool useSpMV = (handleL->get_algorithm () == SPTRSVAlgorithm::SUPERNODAL_SPMV ||
                    handleL->get_algorithm () == SPTRSVAlgorithm::SUPERNODAL_SPMV_DAG);
    #ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
    double time_seconds = 0.0;
    if (merge)          std::cout << " >> merge\n" << std::endl;
    if (invert_offdiag) std::cout << " >> invert offdiag\n" << std::endl;
    #endif

    // ===================================================================
    // load graphs
    auto graphL = handleL->get_graph ();
    auto graphL_host = handleL->get_graph_host ();

    auto graphU = handleU->get_graph ();
    auto graphU_host = handleU->get_graph_host ();

    int nsuper = handleL->get_num_supernodes ();
    const int* supercols = handleL->get_supercols_host ();
    crsmat_t superluL, superluU;
    host_crsmat_t superluL_host, superluU_host;
    if (merge) {
      tic.reset ();
      // ========================================================
      // read in the numerical L-values into merged csc
      // NOTE: we first load into CRS, and then merge (should be combined)
      bool cusparse = false;
      bool invert_diag = false; // invert after merge
      // 1) load L into crs (offdiagonal not inverted, unless invert diag)
      auto original_graphL_host = handleL->get_original_graph_host ();
      superluL_host = read_superlu_valuesL<host_crsmat_t> (cusparse, merge, invert_diag, invert_offdiag, &L, original_graphL_host);
      // 2) re-load L into merged crs
      invert_diag = true;       // TODO: diagonals are always inverted
      bool lower = true;
      bool unit_diag = true;
      if (useSpMV) {
        superluL_host = read_merged_supernodes<host_crsmat_t> (nsuper, supercols,
                                                               lower, unit_diag, invert_diag, invert_offdiag,
                                                               superluL_host, graphL_host);
      } else {
        superluL = read_merged_supernodes<crsmat_t> (nsuper, supercols,
                                                     lower, unit_diag, invert_diag, invert_offdiag,
                                                     superluL_host, graphL);
      }

      // ========================================================
      // read in the numerical U-values into merged csr
      // 1) load U into crs
      invert_offdiag = handleU->get_invert_offdiagonal ();
      invert_diag = false;     // invert after merge
      auto original_graphU_host = handleU->get_original_graph_host ();
      if (UinCSC) {
        superluU_host = read_superlu_valuesU_CSC<host_crsmat_t> (invert_diag, invert_offdiag, &L, &U, original_graphU_host);
      } else {
        // NOTE: invert-offdiag not supported in CSR format
        superluU_host = read_superlu_valuesU<host_crsmat_t> (invert_diag, &L, &U, original_graphU_host);
      }
      invert_diag = true;      // TODO: diagonals are always inverted
      // 2) re-load U into merged crs
      lower = (UinCSC ? false : true);
      unit_diag = false;
      if (useSpMV) {
        superluU_host = read_merged_supernodes<host_crsmat_t> (nsuper, supercols,
                                                               lower, unit_diag, invert_diag, invert_offdiag,
                                                               superluU_host, graphU_host);
      } else {
        superluU = read_merged_supernodes<crsmat_t> (nsuper, supercols,
                                                     lower, unit_diag, invert_diag, invert_offdiag,
                                                     superluU_host, graphU);
      }
      #ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
      time_seconds = tic.seconds ();
      std::cout << "   Time to Merge and Copy to device: " << time_seconds << std::endl;
      #endif
    } else {
      // ========================================================
      // read in the numerical values into merged csc for L
      bool cusparse = false;
      bool invert_diag = true; // only, invert diag is supported for now
      tic.reset ();
      if (useSpMV) {
        superluL_host = read_superlu_valuesL<host_crsmat_t> (cusparse, merge, invert_diag, invert_offdiag, &L, graphL_host);
      } else {
        superluL = read_superlu_valuesL<crsmat_t> (cusparse, merge, invert_diag, invert_offdiag, &L, graphL);
      }
      #ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
      double timeL = tic.seconds ();
      #endif

      // ========================================================
      // read in the numerical values into merged csc/csr for U
      tic.reset ();
      if (useSpMV) {
        if (UinCSC) {
          superluU_host = read_superlu_valuesU_CSC<host_crsmat_t> (invert_diag, invert_offdiag, &L, &U, graphU_host);
        } else {
          superluU_host = read_superlu_valuesU<host_crsmat_t> (invert_diag, &L, &U, graphU_host);
        }
      } else {
        if (UinCSC) {
          superluU = read_superlu_valuesU_CSC<crsmat_t> (invert_diag, invert_offdiag, &L, &U, graphU);
        } else {
          superluU = read_superlu_valuesU<crsmat_t> (invert_diag, &L, &U, graphU);
        }
      }
      #ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
      double timeU = tic.seconds ();
      std::cout << "   Time to copy to device: " << std::endl;
      std::cout << "   > copy L to device: " << timeL << std::endl;
      std::cout << "   > copy U to device: " << timeU << std::endl;
      #endif
    }

    // ===================================================================
    if (useSpMV) {
      // ----------------------------------------------------
      // split the matrix into submatrices for spmv at each level
      tic.reset ();
      split_crsmat<crsmat_t> (kernelHandleL, superluL_host);
      split_crsmat<crsmat_t> (kernelHandleU, superluU_host);
      #ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
      time_seconds = tic.seconds ();
      std::cout << "   Time to Split to submatrix: " << time_seconds << std::endl;
      #endif
    }


    // ==============================================
    // save crsmat
    handleL->set_crsmat (superluL);
    handleU->set_crsmat (superluU);

    // ===================================================================
    handleL->set_numeric_complete ();
    handleU->set_numeric_complete ();
    #ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
    time_seconds = timer.seconds ();
    std::cout << "   Total Compute Time: " << time_seconds << std::endl << std::endl;
    #endif
  } // sptrsv_compute
#endif //KOKKOSKERNELS_ENABLE_TPL_SUPERLU


#ifdef KOKKOSKERNELS_ENABLE_TPL_CHOLMOD
  // ---------------------------------------------------------------------
  template <typename scalar_type,
            typename ordinal_type,
            typename size_type,
            typename KernelHandle,
            typename execution_space      = Kokkos::DefaultExecutionSpace,
            typename host_execution_space = Kokkos::DefaultHostExecutionSpace>
  void sptrsv_symbolic(
      KernelHandle *kernelHandleL,
      KernelHandle *kernelHandleU,
      cholmod_factor *L,
      cholmod_common *cm)
  {
    typedef KokkosSparse::CrsMatrix<scalar_type, ordinal_type, execution_space, void, size_type>  crsmat_t;
    typedef typename  crsmat_t::StaticCrsGraphType  graph_t;

    // ===================================================================
    // load sptrsv-handles
    auto *handleL = kernelHandleL->get_sptrsv_handle ();
    auto *handleU = kernelHandleU->get_sptrsv_handle ();

    // load options
    int *etree = handleL->get_etree ();
    bool needEtree = (handleL->get_algorithm () == SPTRSVAlgorithm::SUPERNODAL_SPMV ||
                      handleL->get_algorithm () == SPTRSVAlgorithm::SUPERNODAL_ETREE);
    if (needEtree && etree == nullptr) {
      std::cout << std::endl
                << " ** etree needs to be set before calling sptrsv_symbolic with SuperLU **"
                << std::endl << std::endl;
      return;
    }

    Kokkos::Timer timer;
    // ==============================================
    // extract CrsGraph from Cholmod
    bool cusparse = false; // pad diagonal blocks with zeros
    auto graph = read_cholmod_graphL<graph_t>(cusparse, L, cm);
    auto row_map = graph.row_map;
    auto entries = graph.entries;

    // ==============================================
    // setup supnodal info 
    int nsuper = (int)(L->nsuper);
    int *supercols = (int*)(L->super);
    handleL->set_supernodes (nsuper, supercols, etree);
    handleU->set_supernodes (nsuper, supercols, etree);

    // ==============================================
    // symbolic for L-solve on the host
    timer.reset();
    sptrsv_symbolic (kernelHandleL, row_map, entries);
    #ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
    double timeL = timer.seconds ();
    std::cout << " > Lower-TRI: " << std::endl;
    std::cout << "   Symbolic Time: " << timeL << std::endl;
    #endif

    // ==============================================
    // symbolic for L^T-solve on the host
    timer.reset ();
    sptrsv_symbolic (kernelHandleU, row_map, entries);
    #ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
    double timeU = timer.seconds ();
    std::cout << " > Upper-TRI: " << std::endl;
    std::cout << "   Symbolic Time: " << timeU << std::endl;
    #endif

    // ==============================================
    // save graphs
    handleL->set_graph (graph);
    handleU->set_graph (graph);

    // ===================================================================
    handleU->set_symbolic_complete ();
    handleU->set_symbolic_complete ();
  }


  // ---------------------------------------------------------------------
  template <typename scalar_type,
            typename ordinal_type,
            typename size_type,
            typename KernelHandle,
            typename execution_space      = Kokkos::DefaultExecutionSpace,
            typename host_execution_space = Kokkos::DefaultHostExecutionSpace>
  void sptrsv_compute(
      KernelHandle *kernelHandleL,
      KernelHandle *kernelHandleU,
      cholmod_factor *L,
      cholmod_common *cm)
  {
    typedef KokkosSparse::CrsMatrix<scalar_type, ordinal_type, execution_space, void, size_type> crsmat_t;

    // ===================================================================
    // load sptrsv-handles
    auto *handleL = kernelHandleL->get_sptrsv_handle ();
    auto *handleU = kernelHandleU->get_sptrsv_handle ();

    if (!(handleL->is_symbolic_complete()) ||
        !(handleU->is_symbolic_complete())) {
      std::cout << std::endl
                << " ** needs to call sptrsv_symbolic before calling sptrsv_numeric **"
                << std::endl << std::endl;
      return;
    }

    // ==============================================
    // load crsGraph
    auto graph = handleL->get_graph ();

    // ==============================================
    // read numerical values of L from Cholmod
    bool invert_diag = true;
    bool cusparse = false; // pad diagonal blocks with zeros
    auto cholmodL = read_cholmod_factor<crsmat_t> (cusparse, invert_diag, L, cm, graph);

    // ==============================================
    // save crsmat
    bool invert_offdiag = false;
    handleL->set_invert_offdiagonal(invert_offdiag);
    handleU->set_invert_offdiagonal(invert_offdiag);
    handleL->set_crsmat (cholmodL);
    handleU->set_crsmat (cholmodL);

    // ===================================================================
    handleL->set_numeric_complete ();
    handleU->set_numeric_complete ();
  }
#endif // KOKKOSKERNELS_ENABLE_TPL_CHOLMOD

#if defined(KOKKOSKERNELS_ENABLE_TPL_CHOLMOD) | defined(KOKKOSKERNELS_ENABLE_TPL_SUPERLU)
  // ---------------------------------------------------------------------
  template <typename KernelHandle,
            class XType>
  void sptrsv_solve(
      KernelHandle *handle, 
      XType x,
      XType b)
  {
    auto crsmat = handle->get_sptrsv_handle ()->get_crsmat ();
    auto values  = crsmat.values;
    auto graph   = crsmat.graph;
    auto row_map = graph.row_map;
    auto entries = graph.entries;

    if (!(handle->get_sptrsv_handle ()->is_numeric_complete())) {
      std::cout << std::endl
                << " ** needs to call sptrsv_compute before calling sptrsv_solve **"
                << std::endl << std::endl;
      return;
    }

    if (handle->is_sptrsv_lower_tri ()) {
      // apply forward pivoting
      Kokkos::deep_copy (x, b);

      // the fifth argument (i.e., first x) is not used
      sptrsv_solve (handle, row_map, entries, values, x, x);
    } else {
      // the fifth argument (i.e., first x) is not used
      sptrsv_solve (handle, row_map, entries, values, b, b);

      // apply backward pivoting
      Kokkos::deep_copy (x, b);
   }
  }

  // ---------------------------------------------------------------------
  template <typename KernelHandle,
            class XType>
  void sptrsv_solve(
      KernelHandle *handleL,
      KernelHandle *handleU,
      XType x,
      XType b)
  {
    // Lower-triangular solve
    sptrsv_solve(handleL, x, b);

    // copy the solution to rhs
    Kokkos::deep_copy (b, x);

    // uper-triangular solve
    sptrsv_solve(handleU, x, b);
  }
#endif

} // namespace Experimental
} // namespace KokkosSparse

#undef KOKKOSKERNELS_SPTRSV_SAME_TYPE

#endif // KOKKOSSPARSE_SPTRSV_HPP_

