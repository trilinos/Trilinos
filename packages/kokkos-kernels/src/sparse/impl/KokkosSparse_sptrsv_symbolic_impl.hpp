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

#ifndef KOKKOSSPARSE_IMPL_SPTRSV_SYMBOLIC_HPP_
#define KOKKOSSPARSE_IMPL_SPTRSV_SYMBOLIC_HPP_

/// \file Kokkos_Sparse_impl_sptrsv_symbolic.hpp
/// \brief Implementation(s) of sparse triangular solve.

#include <KokkosKernels_config.h>
#include <Kokkos_ArithTraits.hpp>
#include <KokkosSparse_sptrsv_handle.hpp>

//#define TRISOLVE_SYMB_TIMERS
//#define LVL_OUTPUT_INFO
//#define CHAIN_LVL_OUTPUT_INFO

// TODO Pass values array and store diagonal entries - should this always be done or optional?

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

template < class TriSolveHandle, class NPLViewType >
void symbolic_chain_phase(TriSolveHandle &thandle, const NPLViewType &nodes_per_level) {

#ifdef TRISOLVE_SYMB_TIMERS
  Kokkos::Timer timer_sym_chain_total;
#endif
  typedef typename TriSolveHandle::size_type size_type;

  size_type nlevels = thandle.get_num_levels();

  // Create the chain now
  // FIXME Implementations will need to be templated on exec space it seems...
 auto cutoff_threshold = thandle.get_chain_threshold();
 if ( thandle.algm_requires_symb_chain() ) {
  auto h_chain_ptr = thandle.get_host_chain_ptr();
  h_chain_ptr(0) = 0;
  size_type chainlinks_length = 0;
  size_type num_chain_entries = 0;
  int chain_state = 0;
  const int cutoff = cutoff_threshold;
  for ( size_type i = 0; i < nlevels; ++i ) {
    auto cnpl = nodes_per_level(i);
    if (cnpl <= cutoff) {
      // this nlevels may be part of a chain passed to the "single_block" solver to reduce kernel launches
      chainlinks_length += 1;
    }
    else {
      // Too many levels to run on single block...
      // If first lvl <= cutoff but next nlevels isn't, the two aren't separately updated and info is lost...
      // if chainlinks_length > 0, take path so that chain-links updated, then current too large chain updated (i.e. 2 updates); if chainlinks_length == 0, then no previous chains and only one update required (npl too large for single-block
      chain_state = chainlinks_length > 0 ? 2 : 1;
    }

    // if we hit final nlevels before a trigger to update the chain, than override it 
    // in this case, there was not a larger value to miss cutoff and reset the update
    if ( chain_state == 0 && i == nlevels-1 ) { chain_state = 1; }

    if (chain_state == 1) {
      num_chain_entries += 1;
      if (chainlinks_length == 0) {
        h_chain_ptr(num_chain_entries) = h_chain_ptr(num_chain_entries-1) + 1;
      }
      else {
        h_chain_ptr(num_chain_entries) = h_chain_ptr(num_chain_entries-1) + chainlinks_length;
      }
      chainlinks_length = 0; //reset
      chain_state = 0; //reset
    }
    // Two updates required - should only occur if chainlinks_length > 0
    // We have found two things: a non-one length chain, and a subsequent one length chain
    if (chain_state == 2) {
      if (chainlinks_length == 0) { std::runtime_error("MAJOR LOGIC ERROR! TERMINATE!"); }

      num_chain_entries += 1;
      h_chain_ptr(num_chain_entries) = h_chain_ptr(num_chain_entries-1) + chainlinks_length;

      num_chain_entries += 1;
      h_chain_ptr(num_chain_entries) = h_chain_ptr(num_chain_entries-1) + 1;

      chainlinks_length = 0; //reset
      chain_state = 0; //reset
    }
  }
  thandle.set_num_chain_entries(num_chain_entries);

#ifdef CHAIN_LVL_OUTPUT_INFO
  std::cout << "  num_chain_entries = " << thandle.get_num_chain_entries() << std::endl;
  for ( size_type i = 0; i < num_chain_entries+1; ++i )
  {
    std::cout << "chain_ptr(" << i << "): " << h_chain_ptr(i) << std::endl;
  }
#endif
 }

#ifdef TRISOLVE_SYMB_TIMERS
 std::cout << "  Symbolic Chain Phase Total Time: " << timer_sym_chain_total.seconds() << std::endl;;
#endif
} // end symbolic_chain_phase


template < class TriSolveHandle, class RowMapType, class EntriesType >
void lower_tri_symbolic (TriSolveHandle &thandle, const RowMapType drow_map, const EntriesType dentries) {
#ifdef TRISOLVE_SYMB_TIMERS
  Kokkos::Timer timer_sym_lowertri_total;
#endif
 if ( thandle.algm_requires_symb_lvlsched() )
 {
  // Scheduling currently computes on host - need host copy of all views

  typedef typename TriSolveHandle::size_type size_type;

  typedef typename TriSolveHandle::nnz_lno_view_t  DeviceEntriesType;

  typedef typename TriSolveHandle::signed_nnz_lno_view_t DeviceSignedEntriesType;
  typedef typename TriSolveHandle::signed_nnz_lno_view_t::HostMirror HostSignedEntriesType;

  typedef typename TriSolveHandle::signed_integral_t signed_integral_t;

  // Necessary for partitioned persisting sparse matrix
  size_type nrows = drow_map.extent(0)-1;

  auto row_map = Kokkos::create_mirror_view(drow_map);
  Kokkos::deep_copy(row_map, drow_map);

  auto entries = Kokkos::create_mirror_view(dentries);
  Kokkos::deep_copy(entries, dentries);
  
  // get device view - will deep_copy to it at end of this host routine
  DeviceEntriesType dnodes_per_level = thandle.get_nodes_per_level();
  auto nodes_per_level = thandle.get_host_nodes_per_level();

  // get device view - will deep_copy to it at end of this host routine
  DeviceEntriesType dnodes_grouped_by_level = thandle.get_nodes_grouped_by_level();
  auto nodes_grouped_by_level = thandle.get_host_nodes_grouped_by_level();

  DeviceSignedEntriesType dlevel_list = thandle.get_level_list();
  HostSignedEntriesType level_list = Kokkos::create_mirror_view(dlevel_list);
  Kokkos::deep_copy(level_list, dlevel_list);

  HostSignedEntriesType previous_level_list( Kokkos::ViewAllocateWithoutInitializing("previous_level_list"), nrows );
  Kokkos::deep_copy( previous_level_list, signed_integral_t(-1) );

  const bool stored_diagonal = thandle.is_stored_diagonal();
  // diagonal_offsets is uninitialized - deep_copy unnecessary at the beginning, only needed at the end
  auto diagonal_offsets = thandle.get_diagonal_offsets();
  auto hdiagonal_offsets = thandle.get_host_diagonal_offsets();

  size_type level = 0;
  auto starting_node = 0;
  auto ending_node = nrows;

  size_type node_count = 0;

  while (node_count < nrows) {

    for ( size_type row = starting_node; row < ending_node; ++row )
    {
      if ( level_list(row) == -1 ) { // unmarked
        bool is_root = true;
        signed_integral_t ptrstart = row_map(row);
        signed_integral_t ptrend   = row_map(row+1);

        for (signed_integral_t offset = ptrstart; offset < ptrend; ++offset) {
          size_type col = entries(offset);
          if ( previous_level_list(col) == -1 && col != row ) { // unmarked
            if ( col < row ) {
              is_root = false;
              break;
            }
          }
          else if ( col == row ) {
            if (stored_diagonal)
              hdiagonal_offsets(row) = offset;
          }
          else if ( col > row ) {
            std::cout << "\nrow = " << row << "  col = " << col << "  offset = " << offset << std::endl;
            std::runtime_error("SYMB ERROR: Lower tri with colid > rowid - SHOULD NOT HAPPEN!!!");
          }
        } // end for offset , i.e. cols of this row

        if ( is_root == true ) {
          level_list(row) = level;
          nodes_per_level(level) += 1;
          nodes_grouped_by_level(node_count) = row;
          node_count += 1;
        }

      } // end if
    } // end for row

    //Kokkos::deep_copy(previous_level_list, level_list);
    for ( size_type i = 0; i < nrows; ++i ) {
      previous_level_list(i) = level_list(i);
    }

    level += 1;
  } // end while

  thandle.set_num_levels(level);

  // Create the chain now
  if ( thandle.algm_requires_symb_chain() ) {
    symbolic_chain_phase(thandle, nodes_per_level);
  }

  thandle.set_symbolic_complete();

  // Output check
#ifdef LVL_OUTPUT_INFO
  std::cout << "  set symbolic complete: " << thandle.is_symbolic_complete() << std::endl;
  std::cout << "  set num levels: " << thandle.get_num_levels() << std::endl;

  std::cout << "  lower_tri_symbolic result: " << std::endl;
  for ( size_type i = 0; i < node_count; ++i )
  { std::cout << "node: " << i << "  level_list = " << level_list(i) << std::endl; }

  for ( size_type i = 0; i < level; ++i )
  { std::cout << "level: " << i << "  nodes_per_level = " << nodes_per_level(i) << std::endl; }

  for ( size_type i = 0; i < node_count; ++i )
  { std::cout << "i: " << i << "  nodes_grouped_by_level = " << nodes_grouped_by_level(i) << std::endl; }
#endif

  // Deep copy to device views
  Kokkos::deep_copy(dnodes_grouped_by_level, nodes_grouped_by_level);
  Kokkos::deep_copy(dnodes_per_level, nodes_per_level);
  Kokkos::deep_copy(dlevel_list, level_list);
  if (stored_diagonal)
    Kokkos::deep_copy(diagonal_offsets, hdiagonal_offsets);

  // Extra check:
#ifdef LVL_OUTPUT_INFO
  {
  std::cout << "  End symb - extra checks" << std::endl;
  std::cout << "  node_count = " << node_count << std::endl;
  std::cout << "  nlevel = " << level << std::endl;
  std::cout << "  npl.extent = " << nodes_per_level.extent(0) << std::endl;
  long check_count = 0;
  Kokkos::parallel_reduce("check_count host", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, nodes_per_level.extent(0)),
    KOKKOS_LAMBDA (const long i, long& update) {
      update+=nodes_per_level(i);
    }, check_count);
  std::cout << "  host check_count= " << check_count << std::endl;

  check_count = 0; // reset
  Kokkos::parallel_reduce("check_count device", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, dnodes_per_level.extent(0)),
    KOKKOS_LAMBDA (const long i, long& update) {
      update+=dnodes_per_level(i);
    }, check_count);
  std::cout << "  devicecheck_count= " << check_count << std::endl;
  }
#endif
 }

#ifdef TRISOLVE_SYMB_TIMERS
 std::cout << "  Symbolic (lower tri) Total Time: " << timer_sym_lowertri_total.seconds() << std::endl;;
#endif
} // end lower_tri_symbolic


template < class TriSolveHandle, class RowMapType, class EntriesType >
void upper_tri_symbolic ( TriSolveHandle &thandle, const RowMapType drow_map, const EntriesType dentries ) {
#ifdef TRISOLVE_SYMB_TIMERS
  Kokkos::Timer timer_sym_uppertri_total;
#endif

 if ( thandle.algm_requires_symb_lvlsched() )
 {
  // Scheduling currently compute on host - need host copy of all views

  typedef typename TriSolveHandle::size_type size_type;
  typedef typename TriSolveHandle::nnz_lno_view_t  DeviceEntriesType;
  typedef typename TriSolveHandle::signed_nnz_lno_view_t DeviceSignedEntriesType;
  typedef typename TriSolveHandle::signed_nnz_lno_view_t::HostMirror HostSignedEntriesType;
  typedef typename TriSolveHandle::signed_integral_t signed_integral_t;

//  size_type nrows = thandle.get_nrows();
  // Necessary for partitioned persisting sparse matrix
  size_type nrows = drow_map.extent(0)-1;

  auto row_map = Kokkos::create_mirror_view(drow_map);
  Kokkos::deep_copy(row_map, drow_map);

  auto entries = Kokkos::create_mirror_view(dentries);
  Kokkos::deep_copy(entries, dentries);
  
  // get device view - will deep_copy to it at end of this host routine
  DeviceEntriesType dnodes_per_level = thandle.get_nodes_per_level();
  auto nodes_per_level = thandle.get_host_nodes_per_level();

  // get device view - will deep_copy to it at end of this host routine
  DeviceEntriesType dnodes_grouped_by_level = thandle.get_nodes_grouped_by_level();
  auto nodes_grouped_by_level = thandle.get_host_nodes_grouped_by_level();

  DeviceSignedEntriesType dlevel_list = thandle.get_level_list();
  HostSignedEntriesType level_list = Kokkos::create_mirror_view(dlevel_list);
  Kokkos::deep_copy(level_list, dlevel_list);

  HostSignedEntriesType previous_level_list( Kokkos::ViewAllocateWithoutInitializing("previous_level_list"), nrows);
  Kokkos::deep_copy( previous_level_list, signed_integral_t(-1) );

  const bool stored_diagonal = thandle.is_stored_diagonal();
  // diagonal_offsets is uninitialized - deep_copy unnecessary at the beginning, only needed at the end
  auto diagonal_offsets = thandle.get_diagonal_offsets();
  auto hdiagonal_offsets = thandle.get_host_diagonal_offsets();

  size_type level = 0;
  auto starting_node = nrows - 1;
  auto ending_node = 0;

  size_type node_count = 0;

  while (node_count < nrows) {

    for ( signed_integral_t row = starting_node; row >= ending_node; --row )
    {
      if ( level_list(row) == -1 ) { // unmarked
        bool is_root = true;
        signed_integral_t ptrstart = row_map(row);
        signed_integral_t ptrend   = row_map(row+1);

        for (signed_integral_t offset = ptrend-1; offset >= ptrstart; --offset) {
          signed_integral_t col = entries(offset);

          if (previous_level_list(col) == -1 && col != row) { // unmarked
            if ( col > row ) {
              is_root = false;
              break;
            }
          }
          else if ( col == row ) {
            if (stored_diagonal)
              hdiagonal_offsets(row) = offset;
          }
        } // end for offset , i.e. cols of this row

        if ( is_root == true ) {
          level_list(row) = level;
          nodes_per_level(level) += 1;
          nodes_grouped_by_level(node_count) = row;
          node_count += 1;
        }

      } // end if
    } // end for row

    //Kokkos::deep_copy(previous_level_list, level_list);
    for ( size_type i = 0; i < nrows; ++i ) {
      previous_level_list(i) = level_list(i);
    }

    level += 1;
  } // end while

  thandle.set_num_levels(level);

  // Create the chain now
  if ( thandle.algm_requires_symb_chain() ) {
    symbolic_chain_phase(thandle, nodes_per_level);
  }

  thandle.set_symbolic_complete();

  // Output check
#ifdef LVL_OUTPUT_INFO
  std::cout << "  set symbolic complete: " << thandle.is_symbolic_complete() << std::endl;
  std::cout << "  set num levels: " << thandle.get_num_levels() << std::endl;

  std::cout << "  upper_tri_symbolic result: " << std::endl;
  for ( size_type i = 0; i < node_count; ++i )
  { std::cout << "node: " << i << "  level_list = " << level_list(i) << std::endl; }

  for ( size_type i = 0; i < level; ++i )
  { std::cout << "level: " << i << "  nodes_per_level = " << nodes_per_level(i) << std::endl; }

  for ( size_type i = 0; i < node_count; ++i )
  { std::cout << "i: " << i << "  nodes_grouped_by_level = " << nodes_grouped_by_level(i) << std::endl; }
#endif

  // Deep copy to device views
  Kokkos::deep_copy(dnodes_grouped_by_level, nodes_grouped_by_level);
  Kokkos::deep_copy(dnodes_per_level, nodes_per_level);
  Kokkos::deep_copy(dlevel_list, level_list);
  if (stored_diagonal)
    Kokkos::deep_copy(diagonal_offsets, hdiagonal_offsets);

  // Extra check:
#ifdef LVL_OUTPUT_INFO
  {
  std::cout << "  End symb - extra checks" << std::endl;
  std::cout << "  node_count = " << node_count << std::endl;
  std::cout << "  nlevel = " << level << std::endl;
  std::cout << "  npl.extent = " << nodes_per_level.extent(0) << std::endl;
  long check_count = 0;
  Kokkos::parallel_reduce("check_count host", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, nodes_per_level.extent(0)),
    KOKKOS_LAMBDA (const long i, long& update) {
      update+=nodes_per_level(i);
    }, check_count);
  std::cout << "  host check_count= " << check_count << std::endl;

  check_count = 0; // reset
  Kokkos::parallel_reduce("check_count device", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, dnodes_per_level.extent(0)),
    KOKKOS_LAMBDA (const long i, long& update) {
      update+=dnodes_per_level(i);
    }, check_count);
  std::cout << "  devicecheck_count= " << check_count << std::endl;
  }
#endif

 }

#ifdef TRISOLVE_SYMB_TIMERS
 std::cout << "  Symbolic (upper tri) Total Time: " << timer_sym_uppertri_total.seconds() << std::endl;;
#endif
} // end upper_tri_symbolic


} // namespace Experimental
} // namespace Impl
} // namespace KokkosSparse

#ifdef LVL_OUTPUT_INFO
#undef LVL_OUTPUT_INFO
#endif

#ifdef CHAIN_LVL_OUTPUT_INFO
#undef CHAIN_LVL_OUTPUT_INFO
#endif

#endif
