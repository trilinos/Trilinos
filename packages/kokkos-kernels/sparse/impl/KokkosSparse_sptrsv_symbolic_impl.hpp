//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOSSPARSE_IMPL_SPTRSV_SYMBOLIC_HPP_
#define KOKKOSSPARSE_IMPL_SPTRSV_SYMBOLIC_HPP_
/// \file KokkosSparse_impl_sptrsv_symbolic.hpp
/// \brief Implementation(s) of sparse triangular solve.

#include <KokkosKernels_config.h>
#include <Kokkos_ArithTraits.hpp>
#include <KokkosSparse_sptrsv_handle.hpp>

// #define TRISOLVE_SYMB_TIMERS
// #define LVL_OUTPUT_INFO
// #define CHAIN_LVL_OUTPUT_INFO

// TODO Pass values array and store diagonal entries - should this always be
// done or optional?

namespace KokkosSparse {
namespace Impl {
namespace Experimental {

template <class ViewType>
void print_view1d_symbolic(const ViewType dv, size_t range = 0) {
  auto v = Kokkos::create_mirror_view(dv);
  Kokkos::deep_copy(v, dv);
  std::cout << "Output for view " << v.label() << std::endl;
  range = range == 0 ? dv.extent(0) : range;
  for (size_t i = 0; i < range; ++i) {
    std::cout << "v(" << i << ") = " << v(i) << " , ";
  }
  std::cout << std::endl;
}

// Usage:
// for c in [0, num_chain_entries)
//   s = h_chain_ptr(c); e = h_chain_ptr(c+1);
//   num_levels_in_current_chain = e - s;
//   if nlicc > 256
//     call current_alg
//   else
//     call single_block(s,e)

template <class TriSolveHandle, class NPLViewType>
void symbolic_chain_phase(TriSolveHandle& thandle, const NPLViewType& nodes_per_level) {
#ifdef TRISOLVE_SYMB_TIMERS
  Kokkos::Timer timer_sym_chain_total;
#endif
  typedef typename TriSolveHandle::size_type size_type;

  size_type nlevels = thandle.get_num_levels();

  // Create the chain now
  // FIXME Implementations will need to be templated on exec space it seems...
  auto cutoff_threshold = thandle.get_chain_threshold();
  if (thandle.algm_requires_symb_chain()) {
    auto h_chain_ptr            = thandle.get_host_chain_ptr();
    h_chain_ptr(0)              = 0;
    size_type chainlinks_length = 0;
    size_type num_chain_entries = 0;
    int chain_state             = 0;
    const int cutoff            = cutoff_threshold;
    for (size_type i = 0; i < nlevels; ++i) {
      auto cnpl = nodes_per_level(i);
      if (cnpl <= cutoff) {
        // this nlevels may be part of a chain passed to the "single_block"
        // solver to reduce kernel launches
        chainlinks_length += 1;
      } else {
        // Too many levels to run on single block...
        // If first lvl <= cutoff but next nlevels isn't, the two aren't
        // separately updated and info is lost... if chainlinks_length > 0, take
        // path so that chain-links updated, then current too large chain
        // updated (i.e. 2 updates); if chainlinks_length == 0, then no previous
        // chains and only one update required (npl too large for single-block
        chain_state = chainlinks_length > 0 ? 2 : 1;
      }

      // if we hit final nlevels before a trigger to update the chain, than
      // override it in this case, there was not a larger value to miss cutoff
      // and reset the update
      if (chain_state == 0 && i == nlevels - 1) {
        chain_state = 1;
      }

      if (chain_state == 1) {
        num_chain_entries += 1;
        if (chainlinks_length == 0) {
          h_chain_ptr(num_chain_entries) = h_chain_ptr(num_chain_entries - 1) + 1;
        } else {
          h_chain_ptr(num_chain_entries) = h_chain_ptr(num_chain_entries - 1) + chainlinks_length;
        }
        chainlinks_length = 0;  // reset
        chain_state       = 0;  // reset
      }
      // Two updates required - should only occur if chainlinks_length > 0
      // We have found two things: a non-one length chain, and a subsequent one
      // length chain
      if (chain_state == 2) {
        if (chainlinks_length == 0) {
          throw(std::runtime_error("MAJOR LOGIC ERROR! TERMINATE!"));
        }

        num_chain_entries += 1;
        h_chain_ptr(num_chain_entries) = h_chain_ptr(num_chain_entries - 1) + chainlinks_length;

        num_chain_entries += 1;
        h_chain_ptr(num_chain_entries) = h_chain_ptr(num_chain_entries - 1) + 1;

        chainlinks_length = 0;  // reset
        chain_state       = 0;  // reset
      }
    }
    thandle.set_num_chain_entries(num_chain_entries);

#ifdef CHAIN_LVL_OUTPUT_INFO
    std::cout << "  num_chain_entries = " << thandle.get_num_chain_entries() << std::endl;
    for (size_type i = 0; i < num_chain_entries + 1; ++i) {
      std::cout << "chain_ptr(" << i << "): " << h_chain_ptr(i) << std::endl;
    }
#endif
  }

#ifdef TRISOLVE_SYMB_TIMERS
  std::cout << "  Symbolic Chain Phase Total Time: " << timer_sym_chain_total.seconds() << std::endl;
  ;
#endif
}  // end symbolic_chain_phase

template <class ExecSpaceIn, class TriSolveHandle, class RowMapType, class EntriesType>
void lower_tri_symbolic(ExecSpaceIn& space, TriSolveHandle& thandle, const RowMapType drow_map,
                        const EntriesType dentries) {
#ifdef TRISOLVE_SYMB_TIMERS
  Kokkos::Timer timer_sym_lowertri_total;
  Kokkos::Timer timer;
#endif

  using namespace KokkosSparse::Experimental;
  if (thandle.get_algorithm() == SPTRSVAlgorithm::SEQLVLSCHD_RP ||
      thandle.get_algorithm() == SPTRSVAlgorithm::SEQLVLSCHD_TP1 ||
      /*thandle.get_algorithm () == SPTRSVAlgorithm::SEQLVLSCHED_TP2*/
      thandle.get_algorithm() == SPTRSVAlgorithm::SEQLVLSCHD_TP1CHAIN) {
    // Scheduling currently computes on host - need host copy of all views

    typedef typename TriSolveHandle::size_type size_type;

    typedef typename TriSolveHandle::nnz_lno_view_t DeviceEntriesType;

    typedef typename TriSolveHandle::signed_nnz_lno_view_t DeviceSignedEntriesType;
    typedef typename TriSolveHandle::signed_nnz_lno_view_t::HostMirror HostSignedEntriesType;

    typedef typename TriSolveHandle::signed_integral_t signed_integral_t;

    // Necessary for partitioned persisting sparse matrix
    size_type nrows = drow_map.extent(0) - 1;

    auto row_map = Kokkos::create_mirror_view(drow_map);
    Kokkos::deep_copy(space, row_map, drow_map);

    auto entries = Kokkos::create_mirror_view(dentries);
    Kokkos::deep_copy(space, entries, dentries);

    // get device view - will deep_copy to it at end of this host routine
    DeviceEntriesType dnodes_per_level = thandle.get_nodes_per_level();
    auto nodes_per_level               = thandle.get_host_nodes_per_level();

    // get device view - will deep_copy to it at end of this host routine
    DeviceEntriesType dnodes_grouped_by_level = thandle.get_nodes_grouped_by_level();
    auto nodes_grouped_by_level               = thandle.get_host_nodes_grouped_by_level();

    DeviceSignedEntriesType dlevel_list = thandle.get_level_list();
    HostSignedEntriesType level_list    = Kokkos::create_mirror_view(dlevel_list);
    Kokkos::deep_copy(space, level_list, dlevel_list);

    signed_integral_t level = 0;
    size_type node_count    = 0;

    space.fence();                                                      // wait for deep copy write to land
    typename DeviceEntriesType::HostMirror level_ptr("lp", nrows + 1);  // temp View used for index bookkeeping
    level_ptr(0) = 0;
    for (size_type i = 0; i < nrows; ++i) {
      signed_integral_t l = 0;
      size_type rowstart  = row_map(i);
      size_type rowend    = row_map(i + 1);
      for (size_type j = rowstart; j < rowend; j++) {
        size_type col = entries(j);
        l             = std::max(l, level_list(col));
      }
      level_list(i) = l + 1;
      nodes_per_level(l) += 1;  // 0-based indexing
      level_ptr(l + 1) += 1;
      level = std::max(level, l + 1);
      node_count++;
    }
    for (signed_integral_t i = 1; i <= level; ++i) {
      level_ptr(i) += level_ptr(i - 1);
    }
    for (size_type i = 0; i < nrows; i++) {
      nodes_grouped_by_level(level_ptr(level_list(i) - 1)) = i;
      level_ptr(level_list(i) - 1) += 1;
    }

    thandle.set_num_levels(level);

    // Create the chain now
    if (thandle.algm_requires_symb_chain()) {
      // No need to pass in space, chain phase runs on the host
      symbolic_chain_phase(thandle, nodes_per_level);
    }
    thandle.set_symbolic_complete();

    // Output check
#ifdef LVL_OUTPUT_INFO
    std::cout << "  set symbolic complete: " << thandle.is_symbolic_complete() << std::endl;
    std::cout << "  set num levels: " << thandle.get_num_levels() << std::endl;

    std::cout << "  lower_tri_symbolic result: " << std::endl;
    for (size_type i = 0; i < node_count; ++i) {
      std::cout << "node: " << i << "  level_list = " << level_list(i) << std::endl;
    }

    for (size_type i = 0; i < level; ++i) {
      std::cout << "level: " << i << "  nodes_per_level = " << nodes_per_level(i) << std::endl;
    }

    for (size_type i = 0; i < node_count; ++i) {
      std::cout << "i: " << i << "  nodes_grouped_by_level = " << nodes_grouped_by_level(i) << std::endl;
    }
#endif

    // Deep copy to device views
    Kokkos::deep_copy(space, dnodes_grouped_by_level, nodes_grouped_by_level);
    Kokkos::deep_copy(space, dnodes_per_level, nodes_per_level);
    Kokkos::deep_copy(space, dlevel_list, level_list);

    // Extra check:
#ifdef LVL_OUTPUT_INFO
    {
      std::cout << "  End symb - extra checks" << std::endl;
      std::cout << "  node_count = " << node_count << std::endl;
      std::cout << "  nlevel = " << level << std::endl;
      std::cout << "  npl.extent = " << nodes_per_level.extent(0) << std::endl;
      long check_count = 0;
      Kokkos::parallel_reduce(
          "check_count host", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, nodes_per_level.extent(0)),
          KOKKOS_LAMBDA(const long i, long& update) { update += nodes_per_level(i); }, check_count);
      std::cout << "  host check_count= " << check_count << std::endl;

      space.fence();    // wait for deep copy writes to land
      check_count = 0;  // reset
      Kokkos::parallel_reduce(
          "check_count device", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, dnodes_per_level.extent(0)),
          KOKKOS_LAMBDA(const long i, long& update) { update += dnodes_per_level(i); }, check_count);
      std::cout << "  devicecheck_count= " << check_count << std::endl;
    }
#endif
  }
#ifdef KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV
  else if (thandle.get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_NAIVE ||
           thandle.get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_ETREE ||
           thandle.get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_DAG ||
           thandle.get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_SPMV ||
           thandle.get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_SPMV_DAG) {

    using size_type           = typename TriSolveHandle::size_type;
    using signed_integral_t   = typename TriSolveHandle::signed_integral_t;
    using integer_view_t      = typename TriSolveHandle::integer_view_t;
    using integer_view_host_t = typename integer_view_t::HostMirror;

    // rowptr: pointer to begining of each row (CRS)
    auto row_map = Kokkos::create_mirror_view(drow_map);
    Kokkos::deep_copy(row_map, drow_map);

    // # of nodes per level
    auto dnodes_per_level = thandle.get_nodes_per_level();
    auto nodes_per_level  = thandle.get_host_nodes_per_level();

    // node ids in each level
    auto dnodes_grouped_by_level = thandle.get_nodes_grouped_by_level();
    auto nodes_grouped_by_level  = thandle.get_host_nodes_grouped_by_level();

    // map node id to level that this node belongs to
    auto dlevel_list = thandle.get_level_list();
    auto level_list  = Kokkos::create_mirror_view(dlevel_list);

    // type of kernels used at each level
    int size_unblocked = thandle.get_supernode_size_unblocked();
    // int size_blocked = thandle.get_supernode_size_blocked();
    auto kernel_type_by_level      = thandle.get_kernel_type_host();
    auto diag_kernel_type_by_level = thandle.get_diag_kernel_type_host();

    // # of supernodal columns
    size_type nsuper     = thandle.get_num_supernodes();
    const int* supercols = thandle.get_supercols_host();

    // workspace
    signed_integral_t max_lwork = 0;
    auto work_offset_host       = thandle.get_work_offset_host();
    if (thandle.get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_NAIVE) {
      // >> Naive (sequential) version: going through supernodal column one at a
      // time from 1 to nsuper Set number of level equal to be the number of
      // supernodal columns
      thandle.set_num_levels(nsuper);

      // Set up level sets: going through supernodal column one at a time from 1
      // to nsuper
      for (size_type s = 0; s < nsuper; s++) {
        nodes_per_level(s)        = 1;  // # of nodes per level
        nodes_grouped_by_level(s) = s;  // only one task per level (task id)
        level_list(s)             = s;  // map task id to level

        // local/max workspace size
        size_type row           = supercols[s];
        signed_integral_t lwork = row_map(row + 1) - row_map(row);
        if (max_lwork < lwork) {
          max_lwork = lwork;
        }

        // kernel type
        if (lwork < size_unblocked) {
          // batched unblocked
          kernel_type_by_level(s)      = 0;
          diag_kernel_type_by_level(s) = 0;
          //} else if (lwork < size_blocked) {
          //  // batched blocked
          //  kernel_type_by_level (s) = 1;
          //  diag_kernel_type_by_level (s) = 1;
        } else {
          // device
          kernel_type_by_level(s)      = 3;
          diag_kernel_type_by_level(s) = 3;
        }
        work_offset_host(s) = 0;
      }
    } else {
// #define profile_supernodal_etree
#ifdef profile_supernodal_etree
      // min, max, tot size of supernodes
      signed_integral_t max_nsrow = 0;
      signed_integral_t min_nsrow = 0;
      signed_integral_t tot_nsrow = 0;

      signed_integral_t max_nscol = 0;
      signed_integral_t min_nscol = 0;
      signed_integral_t tot_nscol = 0;

      // min, max, tot num of leaves
      signed_integral_t max_nleave = 0;
      signed_integral_t min_nleave = 0;
      signed_integral_t tot_nleave = 0;
#endif

      /* initialize the ready tasks with leaves */
      const int* parents = thandle.get_etree_parents();
      integer_view_host_t check("check", nsuper);
      Kokkos::deep_copy(check, 0);

      auto dag         = thandle.get_supernodal_dag();
      auto dag_row_map = dag.row_map;
      auto dag_entries = dag.entries;
      bool use_dag     = (thandle.get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_DAG ||
                      thandle.get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_SPMV_DAG);
      for (size_type s = 0; s < nsuper; s++) {
        if (use_dag) {
          for (size_type e = dag_row_map(s); e < dag_row_map(s + 1); e++) {
            check(dag_entries(e))++;
          }
        } else {
          if (parents[s] >= 0) {
            check(parents[s])++;
          }
        }
      }

      size_type num_done = 0;
      size_type level    = 0;
      while (num_done < nsuper) {
        nodes_per_level(level) = 0;
        // look for ready-tasks
        signed_integral_t lwork     = 0;
        signed_integral_t num_leave = 0;
        signed_integral_t avg_nscol = 0;
        signed_integral_t avg_nsrow = 0;
        for (size_type s = 0; s < nsuper; s++) {
          if (check(s) == 0) {
            nodes_per_level(level)++;
            nodes_grouped_by_level(num_done + num_leave) = s;
            level_list(s)                                = level;

            // work offset
            work_offset_host(s) = lwork;

            // update workspace size
            size_type row           = supercols[s];
            signed_integral_t nsrow = row_map(row + 1) - row_map(row);
            lwork += nsrow;
            // printf( " %d: nuum_leave=%d, level=%d, nsrow=%d, nscol=%d,
            // s=%d\n",(int)(num_done+num_leave), (int)num_leave, (int)level,
            // (int)nsrow, (int)(supercols[s+1]-supercols[s]),(int)s ); for (int
            // i = supercols[s]; i < supercols[s+1]; i++) printf("%d %d
            // %d\n",i,s,level );  // permute matrix based on scheduling

            // total supernode size
            avg_nsrow += row_map(row + 1) - row_map(row);
            avg_nscol += supercols[s + 1] - supercols[s];

#ifdef profile_supernodal_etree
            // gather static if requested
            signed_integral_t nscol = supercols[s + 1] - supercols[s];
            if (num_done + num_leave == 0) {
              // initialization at the initial step
              max_nscol = nscol;
              min_nscol = nscol;

              max_nsrow = nsrow;
              min_nsrow = nsrow;
            } else {
              if (max_nscol < nscol) {
                max_nscol = nscol;
              }
              if (min_nscol > nscol) {
                min_nscol = nscol;
              }

              if (max_nsrow < nsrow) {
                max_nsrow = nsrow;
              }
              if (min_nsrow > nsrow) {
                min_nsrow = nsrow;
              }
            }
            tot_nsrow += nsrow;
            tot_nscol += nscol;
#endif

            num_leave++;
          }
        }
        if (lwork > max_lwork) {
          max_lwork = lwork;
        }

        // average supernode size at this level
        avg_nsrow /= num_leave;
        avg_nscol /= num_leave;
        // kernel type
        if (avg_nscol < size_unblocked) {
          // batched unblocked
          kernel_type_by_level(level)      = 0;
          diag_kernel_type_by_level(level) = 0;
          //} else if (avg_nscol < size_blocked) {
          //  // batched blocked
          //  kernel_type_by_level (level) = 1;
          //  diag_kernel_type_by_level (level) = 1;
        } else {
          // device
          kernel_type_by_level(level)      = 3;
          diag_kernel_type_by_level(level) = 3;
        }
#ifdef profile_supernodal_etree
        std::cout << level << " : num_leave=" << num_leave << ", nsrow=" << min_nsrow << ", " << avg_nsrow << ", "
                  << max_nsrow << ", nscol=" << min_nscol << ", " << avg_nscol << ", " << max_nscol
                  << ", kernel_type=" << diag_kernel_type_by_level(level) << "(" << size_unblocked << ","
                  << thandle.get_supernode_size_blocked() << ")" << std::endl;
        if (level == 0) {
          max_nleave = num_leave;
          min_nleave = num_leave;
        } else {
          if (max_nleave < num_leave) {
            max_nleave = num_leave;
          }
          if (min_nleave > num_leave) {
            min_nleave = num_leave;
          }
        }
        tot_nleave += num_leave;
#endif

        // free the dependency
        for (signed_integral_t task = 0; task < num_leave; task++) {
          size_type s = nodes_grouped_by_level(num_done + task);
          check(s)    = -1;
          // printf( " %d: check[%d]=%d ",level,s,check[s]);
          if (use_dag) {
            for (size_type e = dag_row_map(s); e < dag_row_map(s + 1); e++) {
              check(dag_entries(e))--;
            }
          } else {
            if (parents[s] >= 0) {
              check(parents[s])--;
              // printf( " -> check[%d]=%d",parents[s],check (parents[s]));
            }
          }
          // printf( "\n" );
        }
        num_done += num_leave;
        // printf( " level=%d: num_done=%d / %d\n",level,num_done,nsuper );
        level++;
      }
#ifdef profile_supernodal_etree
      std::cout << "   * number of supernodes = " << nsuper << std::endl;
      std::cout << "   * supernodal rows: min = " << min_nsrow << "\t max = " << max_nsrow
                << "\t avg = " << tot_nsrow / nsuper << std::endl;
      std::cout << "   * supernodal cols: min = " << min_nscol << "\t max = " << max_nscol
                << "\t avg = " << tot_nscol / nsuper << std::endl;
      std::cout << "   * numer of leaves: min = " << min_nleave << "\t max = " << max_nleave
                << "\t avg = " << tot_nleave / level << std::endl;
      std::cout << "   * level = " << level << std::endl;
#endif
#ifdef TRISOLVE_SYMB_TIMERS
      std::cout << "   + scheduling time = " << timer.seconds() << std::endl;
#endif
      // Set number of level equal to be the number of supernodal columns
      thandle.set_num_levels(level);
    }
#ifdef TRISOLVE_SYMB_TIMERS
    timer.reset();
#endif
    // workspace size
    if (thandle.get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_SPMV ||
        thandle.get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_SPMV_DAG) {
      max_lwork = thandle.get_nrows();
    }
    thandle.set_workspace_size(max_lwork);
    // workspace offset initialized to be zero
    integer_view_t work_offset = thandle.get_work_offset();
    Kokkos::deep_copy(space, work_offset, work_offset_host);

    // kernel types
    // > off-diagonal
    integer_view_t dkernel_type_by_level = thandle.get_kernel_type();
    Kokkos::deep_copy(space, dkernel_type_by_level, kernel_type_by_level);
    // > diagonal
    integer_view_t ddiag_kernel_type_by_level = thandle.get_diag_kernel_type();
    Kokkos::deep_copy(space, ddiag_kernel_type_by_level, diag_kernel_type_by_level);

    // deep copy to device (of scheduling info)
    Kokkos::deep_copy(space, dnodes_grouped_by_level, nodes_grouped_by_level);
    Kokkos::deep_copy(space, dnodes_per_level, nodes_per_level);
    Kokkos::deep_copy(space, dlevel_list, level_list);

#ifdef TRISOLVE_SYMB_TIMERS
    std::cout << "   + workspace  time = " << timer.seconds() << std::endl;
#endif

    thandle.set_symbolic_complete();
  }
#endif

#ifdef TRISOLVE_SYMB_TIMERS
  std::cout << "  Symbolic (lower tri) Total Time: " << timer_sym_lowertri_total.seconds() << std::endl;
  ;
#endif
}  // end lower_tri_symbolic

template <class ExecutionSpace, class TriSolveHandle, class RowMapType, class EntriesType>
void upper_tri_symbolic(ExecutionSpace& space, TriSolveHandle& thandle, const RowMapType drow_map,
                        const EntriesType dentries) {
#ifdef TRISOLVE_SYMB_TIMERS
  Kokkos::Timer timer_sym_uppertri_total;
  Kokkos::Timer timer;
#endif

  using namespace KokkosSparse::Experimental;
  if (thandle.get_algorithm() == SPTRSVAlgorithm::SEQLVLSCHD_RP ||
      thandle.get_algorithm() == SPTRSVAlgorithm::SEQLVLSCHD_TP1 ||
      /*thandle.get_algorithm () == SPTRSVAlgorithm::SEQLVLSCHED_TP2*/
      thandle.get_algorithm() == SPTRSVAlgorithm::SEQLVLSCHD_TP1CHAIN) {
    // Scheduling currently compute on host - need host copy of all views

    typedef typename TriSolveHandle::size_type size_type;
    typedef typename TriSolveHandle::nnz_lno_view_t DeviceEntriesType;
    typedef typename TriSolveHandle::signed_nnz_lno_view_t DeviceSignedEntriesType;
    typedef typename TriSolveHandle::signed_nnz_lno_view_t::HostMirror HostSignedEntriesType;
    typedef typename TriSolveHandle::signed_integral_t signed_integral_t;

    //  size_type nrows = thandle.get_nrows();
    // Necessary for partitioned persisting sparse matrix
    size_type nrows = drow_map.extent(0) - 1;

    auto row_map = Kokkos::create_mirror_view(drow_map);
    Kokkos::deep_copy(space, row_map, drow_map);

    auto entries = Kokkos::create_mirror_view(dentries);
    Kokkos::deep_copy(space, entries, dentries);

    // get device view - will deep_copy to it at end of this host routine
    DeviceEntriesType dnodes_per_level = thandle.get_nodes_per_level();
    auto nodes_per_level               = thandle.get_host_nodes_per_level();

    // get device view - will deep_copy to it at end of this host routine
    DeviceEntriesType dnodes_grouped_by_level = thandle.get_nodes_grouped_by_level();
    auto nodes_grouped_by_level               = thandle.get_host_nodes_grouped_by_level();

    DeviceSignedEntriesType dlevel_list = thandle.get_level_list();
    HostSignedEntriesType level_list    = Kokkos::create_mirror_view(dlevel_list);
    Kokkos::deep_copy(space, level_list, dlevel_list);

    signed_integral_t level = 0;
    size_type node_count    = 0;

    space.fence();                                                      // Wait for deep copy writes to land
    typename DeviceEntriesType::HostMirror level_ptr("lp", nrows + 1);  // temp View used for index bookkeeping
    level_ptr(0) = 0;
    for (size_type ii = nrows; ii > 0; ii--) {
      size_type i = ii - 1;  // Avoid >= 0 comparison in for-loop to prevent
                             // wraparound errors with unsigned types
      signed_integral_t l = 0;
      size_type rowstart  = row_map(i) + 1;  // skip diag
      size_type rowend    = row_map(i + 1);
      for (size_type j = rowstart; j < rowend; ++j) {
        size_type col = entries(j);
        l             = std::max(l, level_list(col));
      }
      level_list(i) = l + 1;
      nodes_per_level(l) += 1;  // 0-based indexing
      level_ptr(l + 1) += 1;
      level = std::max(level, l + 1);
      node_count++;
    }
    for (signed_integral_t i = 1; i <= level; ++i) {
      level_ptr(i) += level_ptr(i - 1);
    }
    for (size_type i = 0; i < nrows; i++) {
      nodes_grouped_by_level(level_ptr(level_list(i) - 1)) = i;
      level_ptr(level_list(i) - 1) += 1;
    }

    thandle.set_num_levels(level);

    // Create the chain now
    if (thandle.algm_requires_symb_chain()) {
      symbolic_chain_phase(thandle, nodes_per_level);
    }

    thandle.set_symbolic_complete();

    // Output check
#ifdef LVL_OUTPUT_INFO
    std::cout << "  set symbolic complete: " << thandle.is_symbolic_complete() << std::endl;
    std::cout << "  set num levels: " << thandle.get_num_levels() << std::endl;

    std::cout << "  upper_tri_symbolic result: " << std::endl;
    for (size_type i = 0; i < node_count; ++i) {
      std::cout << "node: " << i << "  level_list = " << level_list(i) << std::endl;
    }

    for (size_type i = 0; i < level; ++i) {
      std::cout << "level: " << i << "  nodes_per_level = " << nodes_per_level(i) << std::endl;
    }

    for (size_type i = 0; i < node_count; ++i) {
      std::cout << "i: " << i << "  nodes_grouped_by_level = " << nodes_grouped_by_level(i) << std::endl;
    }
#endif

    // Deep copy to device views
    Kokkos::deep_copy(space, dnodes_grouped_by_level, nodes_grouped_by_level);
    Kokkos::deep_copy(space, dnodes_per_level, nodes_per_level);
    Kokkos::deep_copy(space, dlevel_list, level_list);

    // Extra check:
#ifdef LVL_OUTPUT_INFO
    {
      std::cout << "  End symb - extra checks" << std::endl;
      std::cout << "  node_count = " << node_count << std::endl;
      std::cout << "  nlevel = " << level << std::endl;
      std::cout << "  npl.extent = " << nodes_per_level.extent(0) << std::endl;
      long check_count = 0;
      Kokkos::parallel_reduce(
          "check_count host", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, nodes_per_level.extent(0)),
          KOKKOS_LAMBDA(const long i, long& update) { update += nodes_per_level(i); }, check_count);
      std::cout << "  host check_count= " << check_count << std::endl;

      space.fence();    // wait for deep copy writes to land
      check_count = 0;  // reset
      Kokkos::parallel_reduce(
          "check_count device", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, dnodes_per_level.extent(0)),
          KOKKOS_LAMBDA(const long i, long& update) { update += dnodes_per_level(i); }, check_count);
      std::cout << "  devicecheck_count= " << check_count << std::endl;
    }
#endif
  }
#ifdef KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV
  else if (thandle.get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_NAIVE ||
           thandle.get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_ETREE ||
           thandle.get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_DAG ||
           thandle.get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_SPMV ||
           thandle.get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_SPMV_DAG) {

    using size_type = typename TriSolveHandle::size_type;

    using signed_integral_t   = typename TriSolveHandle::signed_integral_t;
    using integer_view_t      = typename TriSolveHandle::integer_view_t;
    using integer_view_host_t = typename integer_view_t::HostMirror;

    // rowptr: pointer to begining of each row (CRS)
    auto row_map = Kokkos::create_mirror_view(drow_map);
    Kokkos::deep_copy(row_map, drow_map);

    // # of nodes per level
    auto dnodes_per_level = thandle.get_nodes_per_level();
    auto nodes_per_level  = thandle.get_host_nodes_per_level();

    // node ids in each level
    auto dnodes_grouped_by_level = thandle.get_nodes_grouped_by_level();
    auto nodes_grouped_by_level  = thandle.get_host_nodes_grouped_by_level();

    // type of kernels used at each level
    int size_unblocked             = thandle.get_supernode_size_unblocked();
    auto kernel_type_by_level      = thandle.get_kernel_type_host();
    auto diag_kernel_type_by_level = thandle.get_diag_kernel_type_host();

    // map node id to level that this node belongs to
    auto dlevel_list = thandle.get_level_list();
    auto level_list  = Kokkos::create_mirror_view(dlevel_list);

    // # of supernodal columns
    size_type nsuper     = thandle.get_num_supernodes();
    const int* supercols = thandle.get_supercols_host();

    // workspace
    signed_integral_t max_lwork = 0;
    auto work_offset_host       = thandle.get_work_offset_host();
    if (thandle.get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_NAIVE) {
      // >> Naive (sequential) version: going through supernodal column one at a
      // time from 1 to nsuper Set number of level equal to be the number of
      // supernodal columns
      thandle.set_num_levels(nsuper);

      // Set up level sets: going through supernodal column one at a time from 1
      // to nsuper
      for (size_type s = 0; s < nsuper; s++) {
        nodes_per_level(s)         = 1;               // # of nodes per level
        nodes_grouped_by_level(s)  = nsuper - 1 - s;  // only one task per level (task id)
        level_list(nsuper - 1 - s) = s;               // map task id to level

        size_type row           = supercols[s];
        signed_integral_t lwork = row_map(row + 1) - row_map(row);
        if (max_lwork < lwork) {
          max_lwork = lwork;
        }
        work_offset_host(s) = 0;

        if (lwork < size_unblocked) {
          // batched unblocked
          kernel_type_by_level(s)      = 0;
          diag_kernel_type_by_level(s) = 0;
        } else {
          // device
          kernel_type_by_level(s)      = 3;
          diag_kernel_type_by_level(s) = 3;
        }
      }
    } else {
/* schduling from bottom to top (as for L-solve) *
 * then reverse it for U-solve                   */
#ifdef profile_supernodal_etree
      // min, max, tot size of supernodes
      signed_integral_t max_nsrow = 0;
      signed_integral_t min_nsrow = 0;
      signed_integral_t tot_nsrow = 0;

      signed_integral_t max_nscol = 0;
      signed_integral_t min_nscol = 0;
      signed_integral_t tot_nscol = 0;

      // min, max, tot num of leaves
      signed_integral_t max_nleave = 0;
      signed_integral_t min_nleave = 0;
      signed_integral_t tot_nleave = 0;
#endif

      /* initialize the ready tasks with leaves */
      const int* parents = thandle.get_etree_parents();
      integer_view_host_t check("check", nsuper);
      Kokkos::deep_copy(check, 0);

      auto dag         = thandle.get_supernodal_dag();
      auto dag_row_map = dag.row_map;
      auto dag_entries = dag.entries;
      bool use_dag     = (thandle.get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_DAG ||
                      thandle.get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_SPMV_DAG);
      if (use_dag) {
        for (size_type s = 0; s < nsuper; s++) {
          for (size_type e = dag_row_map(s); e < dag_row_map(s + 1); e++) {
            check(dag_entries(e))++;
          }
        }
      } else {
        for (size_type s = 0; s < nsuper; s++) {
          if (parents[s] >= 0) {
            check(parents[s])++;
          }
        }
      }

      // printf( " Init:\n" );
      // for (size_type s = 0; s <nsuper; s++) printf( " check[%d] =
      // %d\n",s,check (s) );

      size_type nrows = thandle.get_nrows();
      integer_view_host_t inverse_nodes_per_level("nodes_per_level", nrows);
      integer_view_host_t inverse_nodes_grouped_by_level("nodes_grouped_by_level", nrows);

      size_type num_done = 0;
      size_type level    = 0;
      while (num_done < nsuper) {
        nodes_per_level(level) = 0;
        // look for ready-tasks
        signed_integral_t lwork     = 0;
        signed_integral_t num_leave = 0;
        for (size_type s = 0; s < nsuper; s++) {
          if (check(s) == 0) {
            inverse_nodes_per_level(level)++;
            inverse_nodes_grouped_by_level(num_done + num_leave) = s;
            // printf( " level=%d: %d/%d: s=%d\n",level,
            // num_done+num_leave,nsuper, s );

            // work offset
            work_offset_host(s) = lwork;

            // update workspace size
            size_type row           = supercols[s];
            signed_integral_t nsrow = row_map(row + 1) - row_map(row);
            // printf( " %d %d %d %d %d\n",num_done+num_leave, level, nsrow,
            // supercols[s+1]-supercols[s],s ); for (int i = supercols[s]; i <
            // supercols[s+1]; i++) printf("%d %d %d\n",i,s,level );  // permute
            // matrix based on scheduling
            lwork += nsrow;

#ifdef profile_supernodal_etree
            // gather static if requested
            signed_integral_t nscol = supercols[s + 1] - supercols[s];
            if (tot_nscol == 0) {
              max_nscol = nscol;
              min_nscol = nscol;

              max_nsrow = nsrow;
              min_nsrow = nsrow;
            } else {
              if (max_nscol < nscol) {
                max_nscol = nscol;
              }
              if (min_nscol > nscol) {
                min_nscol = nscol;
              }

              if (max_nsrow < nsrow) {
                max_nsrow = nsrow;
              }
              if (min_nsrow > nsrow) {
                min_nsrow = nsrow;
              }
            }
            tot_nsrow += nsrow;
            tot_nscol += nscol;
#endif

            num_leave++;
          }
        }
        // printf( " lwork = %d\n",lwork );
        if (lwork > max_lwork) {
          max_lwork = lwork;
        }
#ifdef profile_supernodal_etree
        if (level == 0) {
          max_nleave = num_leave;
          min_nleave = num_leave;
        } else {
          if (max_nleave < num_leave) {
            max_nleave = num_leave;
          }
          if (min_nleave > num_leave) {
            min_nleave = num_leave;
          }
        }
        tot_nleave += num_leave;
#endif

        // free the dependency
        for (signed_integral_t task = 0; task < num_leave; task++) {
          size_type s = inverse_nodes_grouped_by_level(num_done + task);
          check(s)    = -1;
          // printf( " %d: check[%d]=%d ",level,s,check (s));
          if (use_dag) {
            for (size_type e = dag_row_map(s); e < dag_row_map(s + 1); e++) {
              check(dag_entries(e))--;
            }
          } else {
            if (parents[s] >= 0) {
              check(parents[s])--;
              // printf( " -> check[%d]=%d",parents[s],check (parents[s]));
            }
          }
          // printf( "\n" );
        }
        num_done += num_leave;
        // printf( " level=%d: num_done=%d / %d\n",level,num_done,nsuper );
        level++;
      }
#ifdef profile_supernodal_etree
      std::cout << "   * number of supernodes = " << nsuper << std::endl;
      std::cout << "   * supernodal rows: min = " << min_nsrow << "\t max = " << max_nsrow
                << "\t avg = " << tot_nsrow / nsuper << std::endl;
      std::cout << "   * supernodal cols: min = " << min_nscol << "\t max = " << max_nscol
                << "\t avg = " << tot_nscol / nsuper << std::endl;
      std::cout << "   * numer of leaves: min = " << min_nleave << "\t max = " << max_nleave
                << "\t avg = " << tot_nleave / level << std::endl;
      std::cout << "   * level = " << level << std::endl;
#endif

      // now invert the lists
      num_done            = 0;
      size_type num_level = level;
      for (level = 0; level < num_level; level++) {
        signed_integral_t num_leave = inverse_nodes_per_level(num_level - level - 1);
        nodes_per_level(level)      = num_leave;
        // printf( " -> nodes_per_level(%d -> %d) = %d\n",num_level-level-1,
        // level, num_leave );

        signed_integral_t avg_nscol = 0;
        signed_integral_t avg_nsrow = 0;
        for (signed_integral_t task = 0; task < num_leave; task++) {
          // signed_integral_t s = inverse_nodes_grouped_by_level (nsuper -
          // (num_done+task) - 1);
          signed_integral_t s = inverse_nodes_grouped_by_level(nsuper - (num_done + num_leave - 1 - task) - 1);

          nodes_grouped_by_level(num_done + task) = s;
          level_list(s)                           = level;
          // printf( " -> level=%d: %d->%d: s=%d\n",level,
          // nsuper-(num_done+task)-1, num_done+task, s );

          size_type row = supercols[s];
          avg_nsrow += row_map(row + 1) - row_map(row);
          avg_nscol += supercols[s + 1] - supercols[s];
        }
        num_done += num_leave;

        // average supernodal size at this level
        avg_nscol /= num_leave;
        avg_nsrow /= num_leave;
        // kernel type
        if (avg_nscol < size_unblocked) {
          // batched unblocked
          kernel_type_by_level(level)      = 0;
          diag_kernel_type_by_level(level) = 0;
        } else {
          // device
          kernel_type_by_level(level)      = 3;
          diag_kernel_type_by_level(level) = 3;
        }
      }
#ifdef TRISOLVE_SYMB_TIMERS
      std::cout << "   + scheduling time = " << timer.seconds() << std::endl;
#endif

      // Set number of levels
      thandle.set_num_levels(num_level);
    }
#ifdef TRISOLVE_SYMB_TIMERS
    timer.reset();
#endif
    // workspace size
    if (thandle.get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_SPMV ||
        thandle.get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_SPMV_DAG) {
      max_lwork = thandle.get_nrows();
    }
    thandle.set_workspace_size(max_lwork);
    // workspace offset initialized to be zero
    integer_view_t work_offset = thandle.get_work_offset();
    Kokkos::deep_copy(work_offset, work_offset_host);

    // kernel type
    // > off-diagonal
    integer_view_t dkernel_type_by_level = thandle.get_kernel_type();
    Kokkos::deep_copy(dkernel_type_by_level, kernel_type_by_level);
    // > diagonal
    integer_view_t ddiag_kernel_type_by_level = thandle.get_diag_kernel_type();
    Kokkos::deep_copy(ddiag_kernel_type_by_level, diag_kernel_type_by_level);

    // deep copy to device (info about scheduling)
    Kokkos::deep_copy(dnodes_grouped_by_level, nodes_grouped_by_level);
    Kokkos::deep_copy(dnodes_per_level, nodes_per_level);
    Kokkos::deep_copy(dlevel_list, level_list);
#ifdef TRISOLVE_SYMB_TIMERS
    std::cout << "   + workspace  time = " << timer.seconds() << std::endl;
#endif

    thandle.set_symbolic_complete();
  }
#endif

#ifdef TRISOLVE_SYMB_TIMERS
  std::cout << "  Symbolic (upper tri) Total Time: " << timer_sym_uppertri_total.seconds() << std::endl;
  ;
#endif
}  // end upper_tri_symbolic

}  // namespace Experimental
}  // namespace Impl
}  // namespace KokkosSparse

#ifdef LVL_OUTPUT_INFO
#undef LVL_OUTPUT_INFO
#endif

#ifdef CHAIN_LVL_OUTPUT_INFO
#undef CHAIN_LVL_OUTPUT_INFO
#endif

#endif
