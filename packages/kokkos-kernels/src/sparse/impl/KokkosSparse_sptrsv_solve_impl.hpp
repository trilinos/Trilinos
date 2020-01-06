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

#ifndef KOKKOSSPARSE_IMPL_SPTRSV_SOLVE_HPP_
#define KOKKOSSPARSE_IMPL_SPTRSV_SOLVE_HPP_

/// \file KokkosSparse_impl_sptrsv.hpp
/// \brief Implementation(s) of sparse triangular solve.

#include <KokkosKernels_config.h>
#include <Kokkos_ArithTraits.hpp>
#include <KokkosSparse_sptrsv_handle.hpp>
#include <KokkosSparse_spmv.hpp>
#include <KokkosSparse_CrsMatrix.hpp>

#include <KokkosBatched_Trsv_Decl.hpp>
#include <KokkosBatched_Trsv_Serial_Impl.hpp>

//#define TRISOLVE_TIMERS
//#define SERIAL_FOR_LOOP

#define KOKKOSKERNELS_SPTRSV_TRILVLSCHED

//#define KOKKOSPSTRSV_SOLVE_IMPL_PROFILE 1
#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOSPSTRSV_SOLVE_IMPL_PROFILE)
#include "cuda_profiler_api.h"
#endif

namespace KokkosSparse {
namespace Impl {
namespace Experimental {

#if defined(KOKKOS_ENABLE_CUDA) && 10000 < CUDA_VERSION && defined(KOKKOSKERNELS_ENABLE_EXP_CUDAGRAPH)
  #define KOKKOSKERNELS_SPTRSV_CUDAGRAPHSUPPORT
#endif

struct UnsortedTag {};

struct LargerCutoffTag {};

struct UnsortedLargerCutoffTag {};

template <class ViewType>
void print_view1d_solve(const ViewType dv, size_t range = 0) {
  auto v = Kokkos::create_mirror_view(dv);
  Kokkos::deep_copy(v, dv);
  std::cout << "Output for view " << v.label() << std::endl;
  range = range == 0 ? dv.extent(0) : range;
  for (size_t i = 0; i < range; ++i) {
    std::cout << "v(" << i << ") = " << v(i) << " , ";
  }
  std::cout << std::endl;
}

// Needed for cudagraphs
struct EmptyFunctor {
  KOKKOS_INLINE_FUNCTION
  void operator()(const int) const {}
};


// This functor unifies the lower and upper implementations, the hope is the "is_lowertri" check does not add noticable time on larger problems
template <class RowMapType, class EntriesType, class ValuesType, class LHSType, class RHSType, class NGBLType>
struct TriLvlSchedTP1SolverFunctor
{
  typedef typename RowMapType::execution_space execution_space;
  typedef Kokkos::TeamPolicy<execution_space> policy_type;
  typedef typename policy_type::member_type member_type;
  typedef typename EntriesType::non_const_value_type lno_t;
  typedef typename ValuesType::non_const_value_type scalar_t;

  RowMapType row_map;
  EntriesType entries;
  ValuesType values;
  LHSType lhs;
  RHSType rhs;
  NGBLType nodes_grouped_by_level;

  const bool is_lowertri;

  long node_count; // like "block" offset into ngbl, my_league is the "local" offset
  long dense_nrows;


  TriLvlSchedTP1SolverFunctor(const RowMapType &row_map_, const EntriesType &entries_, const ValuesType &values_, LHSType &lhs_, const RHSType &rhs_, const NGBLType &nodes_grouped_by_level_, const bool is_lowertri_, long node_count_, long dense_nrows_ = 0) :
    row_map(row_map_), entries(entries_), values(values_), lhs(lhs_), rhs(rhs_), nodes_grouped_by_level(nodes_grouped_by_level_), is_lowertri(is_lowertri_), node_count(node_count_), dense_nrows(dense_nrows_) {}


  KOKKOS_INLINE_FUNCTION
  void operator()( const member_type & team ) const {
        auto my_league = team.league_rank(); // map to rowid
        auto rowid = nodes_grouped_by_level(my_league + node_count);
        auto my_rank = team.team_rank();

        auto soffset = row_map(rowid);
        auto eoffset = row_map(rowid+1);
        auto rhs_rowid = rhs(rowid);
        scalar_t diff = scalar_t(0.0);

      Kokkos::parallel_reduce( Kokkos::TeamThreadRange( team, soffset, eoffset ), [&] ( const long ptr, scalar_t &tdiff ) {
          auto colid = entries(ptr);

          auto val   = values(ptr);
          if ( colid != rowid ) {
            tdiff = tdiff - val*lhs(colid);
          }
      }, diff );

        team.team_barrier();

        // At end, finalize rowid == colid
        // only one thread should do this; can also use Kokkos::single
        if ( my_rank == 0 )
        {
        // ASSUMPTION: sorted diagonal value located at eoffset - 1
          lhs(rowid) = is_lowertri ? (rhs_rowid+diff)/values(eoffset-1) : (rhs_rowid+diff)/values(soffset);
        }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const UnsortedTag&, const member_type & team) const {
        auto my_league = team.league_rank(); // map to rowid
        auto rowid = nodes_grouped_by_level(my_league + node_count);
        auto my_rank = team.team_rank();

        auto soffset = row_map(rowid);
        auto eoffset = row_map(rowid+1);
        auto rhs_rowid = rhs(rowid);
        scalar_t diff = scalar_t(0.0);

        auto diag = -1;

        Kokkos::parallel_reduce( Kokkos::TeamThreadRange( team, soffset, eoffset ), [&] ( const long ptr, scalar_t &tdiff ) {
          auto colid = entries(ptr);
          auto val   = values(ptr);
          if ( colid != rowid ) {
            tdiff = tdiff - val*lhs(colid);
          }
          else {
            diag = ptr;
          }
        }, diff );
        team.team_barrier();

        // At end, finalize rowid == colid
        // only one thread should do this; can also use Kokkos::single
        if ( my_rank == 0 )
        {
          lhs(rowid) = (rhs_rowid+diff)/values(diag);
        }
  }
};


template <class RowMapType, class EntriesType, class ValuesType, class LHSType, class RHSType, class NGBLType>
struct TriLvlSchedTP1SolverFunctorDiagValues
{
  typedef typename RowMapType::execution_space execution_space;
  typedef Kokkos::TeamPolicy<execution_space> policy_type;
  typedef typename policy_type::member_type member_type;
  typedef typename EntriesType::non_const_value_type lno_t;
  typedef typename ValuesType::non_const_value_type scalar_t;

  RowMapType row_map;
  EntriesType entries;
  ValuesType values;
  LHSType lhs;
  RHSType rhs;
  NGBLType nodes_grouped_by_level;
  ValuesType diagonal_values; // inserted according to rowid

  const bool is_lowertri;

  long node_count; // like "block" offset into ngbl, my_league is the "local" offset
  long dense_nrows;


  TriLvlSchedTP1SolverFunctorDiagValues( const RowMapType &row_map_, const EntriesType &entries_, const ValuesType &values_, LHSType &lhs_, const RHSType &rhs_, const NGBLType &nodes_grouped_by_level_, const ValuesType &diagonal_values_, const bool is_lowertri_, long node_count_, long dense_nrows_ = 0) :
    row_map(row_map_), entries(entries_), values(values_), lhs(lhs_), rhs(rhs_), nodes_grouped_by_level(nodes_grouped_by_level_), diagonal_values(diagonal_values_), is_lowertri(is_lowertri_), node_count(node_count_), dense_nrows(dense_nrows_) {}


  KOKKOS_INLINE_FUNCTION
  void operator()( const member_type & team ) const {
        auto my_league = team.league_rank(); // map to rowid
        auto rowid = nodes_grouped_by_level(my_league + node_count);
        auto my_rank = team.team_rank();

        auto soffset = row_map(rowid);
        auto eoffset = row_map(rowid+1);
        auto rhs_rowid = rhs(rowid);
        scalar_t diff = scalar_t(0.0);

      Kokkos::parallel_reduce( Kokkos::TeamThreadRange( team, soffset, eoffset ), [&] ( const long ptr, scalar_t &tdiff ) {
          auto colid = entries(ptr);
          auto val   = values(ptr);
          if ( colid != rowid ) {
            tdiff = tdiff - val*lhs(colid);
          }
      }, diff );

        team.team_barrier();

        // At end, finalize rowid == colid
        // only one thread should do this; can also use Kokkos::single
        if ( my_rank == 0 )
        {
          //lhs(rowid) = is_lowertri ? (rhs_rowid+diff)/values(eoffset-1) : (rhs_rowid+diff)/values(soffset);
          lhs(rowid) = (rhs_rowid+diff)/diagonal_values(rowid);
        }
  }

};


template <class RowMapType, class EntriesType, class ValuesType, class LHSType, class RHSType, class NGBLType>
struct TriLvlSchedTP2SolverFunctor
{
  typedef typename RowMapType::execution_space execution_space;
  typedef Kokkos::TeamPolicy<execution_space> policy_type;
  typedef typename policy_type::member_type member_type;
  typedef typename EntriesType::non_const_value_type lno_t;
  typedef typename ValuesType::non_const_value_type scalar_t;

  RowMapType row_map;
  EntriesType entries;
  ValuesType values;
  LHSType lhs;
  RHSType rhs;
  NGBLType nodes_grouped_by_level;

  const bool is_lowertri;
  long node_count; // like "block" offset into ngbl, my_league is the "local" offset
  long node_groups;
  long dense_nrows;


  TriLvlSchedTP2SolverFunctor(const RowMapType &row_map_, const EntriesType &entries_, const ValuesType &values_, LHSType &lhs_, const RHSType &rhs_, const NGBLType &nodes_grouped_by_level_, const bool is_lowertri_, long node_count_, long node_groups_ = 0, long dense_nrows_ = 0) :
    row_map(row_map_), entries(entries_), values(values_), lhs(lhs_), rhs(rhs_), nodes_grouped_by_level(nodes_grouped_by_level_), is_lowertri(is_lowertri_), node_count(node_count_), node_groups(node_groups_), dense_nrows(dense_nrows_) {}


  KOKKOS_INLINE_FUNCTION
  void operator()(const member_type & team) const {
        auto my_league = team.league_rank(); // map to rowid

        size_t nrows = row_map.extent(0) - 1;

        Kokkos::parallel_for( Kokkos::TeamThreadRange( team, 0, node_groups ), [&] ( const long ng ) {
          auto rowid = nodes_grouped_by_level(node_count + my_league*node_groups + ng);
          if ( size_t(rowid) < nrows ) {

            auto soffset = row_map(rowid);
            auto eoffset = row_map(rowid+1);
            auto rhs_rowid = rhs(rowid);
            scalar_t diff = scalar_t(0.0);

            Kokkos::parallel_reduce( Kokkos::ThreadVectorRange( team, soffset, eoffset ), [&] ( const long ptr, scalar_t &tdiff ) {
              auto colid = entries(ptr);
              auto val   = values(ptr);
              if ( colid != rowid ) {
                tdiff = tdiff - val*lhs(colid);
              }
            }, diff );

            // ASSUMPTION: sorted diagonal value located at eoffset - 1
            lhs(rowid) = is_lowertri ? (rhs_rowid+diff)/values(eoffset-1) : (rhs_rowid+diff)/values(soffset);
          } // end if
        }); // end TeamThreadRange

        team.team_barrier();
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const UnsortedTag&, const member_type & team) const {
        auto my_league = team.league_rank(); // map to rowid

        size_t nrows = row_map.extent(0) - 1;

        Kokkos::parallel_for( Kokkos::TeamThreadRange( team, 0, node_groups ), [&] ( const long ng ) {
          auto rowid = nodes_grouped_by_level(node_count + my_league*node_groups + ng);
          if ( size_t(rowid) < nrows ) {
            auto soffset = row_map(rowid);
            auto eoffset = row_map(rowid+1);
            auto rhs_rowid = rhs(rowid);
            scalar_t diff = scalar_t(0.0);

            auto diag = -1;
            Kokkos::parallel_reduce( Kokkos::ThreadVectorRange( team, soffset, eoffset ), [&] ( const long ptr, scalar_t &tdiff ) {
              auto colid = entries(ptr);
              auto val   = values(ptr);
              if ( colid != rowid ) {
                tdiff = tdiff - val*lhs(colid);
              }
              else {
                diag = ptr;
              }
            }, diff );

            lhs(rowid) = (rhs_rowid+diff)/values(diag);
          } // end if
        }); // end TeamThreadRange

        team.team_barrier();
  }
};


// Lower vs Upper Multi-block Functors

template <class RowMapType, class EntriesType, class ValuesType, class LHSType, class RHSType, class NGBLType>
struct LowerTriLvlSchedRPSolverFunctor
{
  typedef typename EntriesType::non_const_value_type lno_t;
  RowMapType row_map;
  EntriesType entries;
  ValuesType values;
  LHSType lhs;
  RHSType rhs;
  NGBLType nodes_grouped_by_level;

  LowerTriLvlSchedRPSolverFunctor( const RowMapType &row_map_, const EntriesType &entries_, const ValuesType &values_, LHSType &lhs_, const RHSType &rhs_, NGBLType nodes_grouped_by_level_ ) :
    row_map(row_map_), entries(entries_), values(values_), lhs(lhs_), rhs(rhs_), nodes_grouped_by_level(nodes_grouped_by_level_) {}


  KOKKOS_INLINE_FUNCTION
  void operator()(const lno_t i) const {
    auto rowid = nodes_grouped_by_level(i);
    // Assuming indices are sorted per row, diag entry is final index in the list

    auto soffset = row_map(rowid);
    auto eoffset = row_map(rowid+1);
    auto rhs_rowid = rhs(rowid);

    for ( auto ptr = soffset; ptr < eoffset; ++ptr ) {
      auto colid = entries(ptr);
      auto val   = values(ptr);
      if ( colid != rowid ) {
        rhs_rowid = rhs_rowid - val*lhs(colid);
      }
      else {
        lhs(rowid) = rhs_rowid/val;
      }
    } // end for ptr
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const UnsortedTag&, const lno_t i) const {
    auto rowid = nodes_grouped_by_level(i);
    auto soffset = row_map(rowid);
    auto eoffset = row_map(rowid+1);
    auto rhs_rowid = rhs(rowid);
    auto diag = -1;

    for ( auto ptr = soffset; ptr < eoffset; ++ptr ) {
      auto colid = entries(ptr);
      auto val   = values(ptr);
      if ( colid != rowid ) {
        rhs_rowid = rhs_rowid - val*lhs(colid);
      }
      else {
        diag = ptr;
      }
    } // end for ptr
    lhs(rowid) = rhs_rowid/values(diag);
  }
};



template <class RowMapType, class EntriesType, class ValuesType, class LHSType, class RHSType, class NGBLType>
struct LowerTriLvlSchedTP1SolverFunctor
{
  typedef typename RowMapType::execution_space execution_space;
  typedef Kokkos::TeamPolicy<execution_space> policy_type;
  typedef typename policy_type::member_type member_type;
  typedef typename EntriesType::non_const_value_type lno_t;
  typedef typename ValuesType::non_const_value_type scalar_t;

  RowMapType row_map;
  EntriesType entries;
  ValuesType values;
  LHSType lhs;
  RHSType rhs;
  NGBLType nodes_grouped_by_level;

  long node_count; // like "block" offset into ngbl, my_league is the "local" offset
  long node_groups;


  LowerTriLvlSchedTP1SolverFunctor( const RowMapType &row_map_, const EntriesType &entries_, const ValuesType &values_, LHSType &lhs_, const RHSType &rhs_, const NGBLType &nodes_grouped_by_level_, long node_count_, long node_groups_ = 0) :
    row_map(row_map_), entries(entries_), values(values_), lhs(lhs_), rhs(rhs_), nodes_grouped_by_level(nodes_grouped_by_level_), node_count(node_count_), node_groups(node_groups_) {}


  KOKKOS_INLINE_FUNCTION
  void operator()( const member_type & team ) const {
        auto my_league = team.league_rank(); // map to rowid
        auto rowid = nodes_grouped_by_level(my_league + node_count);
        auto my_rank = team.team_rank();

        auto soffset = row_map(rowid);
        auto eoffset = row_map(rowid+1);
        auto rhs_rowid = rhs(rowid);
        scalar_t diff = scalar_t(0.0);

      Kokkos::parallel_reduce( Kokkos::TeamThreadRange( team, soffset, eoffset ), [&] ( const long ptr, scalar_t &tdiff ) {
          auto colid = entries(ptr);
          auto val   = values(ptr);
          if ( colid != rowid ) {
            tdiff = tdiff - val*lhs(colid);
          }
      }, diff );

        team.team_barrier();

        // At end, finalize rowid == colid
        // only one thread should do this; can also use Kokkos::single
        if ( my_rank == 0 )
        {
        // ASSUMPTION: sorted diagonal value located at eoffset - 1
          lhs(rowid) = (rhs_rowid+diff)/values(eoffset-1);
        }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const UnsortedTag&, const member_type & team) const {
        auto my_league = team.league_rank(); // map to rowid
        auto rowid = nodes_grouped_by_level(my_league + node_count);
        auto my_rank = team.team_rank();

        auto soffset = row_map(rowid);
        auto eoffset = row_map(rowid+1);
        auto rhs_rowid = rhs(rowid);
        scalar_t diff = scalar_t(0.0);

        auto diag = -1;

        Kokkos::parallel_reduce( Kokkos::TeamThreadRange( team, soffset, eoffset ), [&] ( const long ptr, scalar_t &tdiff ) {
          auto colid = entries(ptr);
          auto val   = values(ptr);
          if ( colid != rowid ) {
            tdiff = tdiff - val*lhs(colid);
          }
          else {
            diag = ptr;
          }
        }, diff );
        team.team_barrier();

        // At end, finalize rowid == colid
        // only one thread should do this; can also use Kokkos::single
        if ( my_rank == 0 )
        {
          lhs(rowid) = (rhs_rowid+diff)/values(diag);
        }
  }
};


// FIXME CUDA: This algorithm not working with all integral type combos
// In any case, this serves as a skeleton for 3-level hierarchical parallelism for alg dev
template <class RowMapType, class EntriesType, class ValuesType, class LHSType, class RHSType, class NGBLType>
struct LowerTriLvlSchedTP2SolverFunctor
{
  typedef typename RowMapType::execution_space execution_space;
  typedef Kokkos::TeamPolicy<execution_space> policy_type;
  typedef typename policy_type::member_type member_type;
  typedef typename EntriesType::non_const_value_type lno_t;
  typedef typename ValuesType::non_const_value_type scalar_t;

  RowMapType row_map;
  EntriesType entries;
  ValuesType values;
  LHSType lhs;
  RHSType rhs;
  NGBLType nodes_grouped_by_level;

  long node_count; // like "block" offset into ngbl, my_league is the "local" offset
  long node_groups;


  LowerTriLvlSchedTP2SolverFunctor(const RowMapType &row_map_, const EntriesType &entries_, const ValuesType &values_, LHSType &lhs_, const RHSType &rhs_, const NGBLType &nodes_grouped_by_level_, long node_count_, long node_groups_ = 0) :
    row_map(row_map_), entries(entries_), values(values_), lhs(lhs_), rhs(rhs_), nodes_grouped_by_level(nodes_grouped_by_level_), node_count(node_count_), node_groups(node_groups_) {}


  KOKKOS_INLINE_FUNCTION
  void operator()(const member_type & team) const {
        auto my_league = team.league_rank(); // map to rowid

        size_t nrows = row_map.extent(0) - 1;

        Kokkos::parallel_for( Kokkos::TeamThreadRange( team, 0, node_groups ), [&] ( const long ng ) {
          auto rowid = nodes_grouped_by_level(node_count + my_league*node_groups + ng);
          if ( size_t(rowid) < nrows ) {

            auto soffset = row_map(rowid);
            auto eoffset = row_map(rowid+1);
            auto rhs_rowid = rhs(rowid);
            scalar_t diff = scalar_t(0.0);

            Kokkos::parallel_reduce( Kokkos::ThreadVectorRange( team, soffset, eoffset ), [&] ( const long ptr, scalar_t &tdiff ) {
              auto colid = entries(ptr);
              auto val   = values(ptr);
              if ( colid != rowid ) {
                tdiff = tdiff - val*lhs(colid);
              }
            }, diff );

            // ASSUMPTION: sorted diagonal value located at eoffset - 1
            lhs(rowid) = (rhs_rowid+diff)/values(eoffset-1);
          } // end if
        }); // end TeamThreadRange

        team.team_barrier();
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const UnsortedTag&, const member_type & team) const {
        auto my_league = team.league_rank(); // map to rowid

        size_t nrows = row_map.extent(0) - 1;

        Kokkos::parallel_for( Kokkos::TeamThreadRange( team, 0, node_groups ), [&] ( const long ng ) {
          auto rowid = nodes_grouped_by_level(node_count + my_league*node_groups + ng);
          if ( size_t(rowid) < nrows ) {
            auto soffset = row_map(rowid);
            auto eoffset = row_map(rowid+1);
            auto rhs_rowid = rhs(rowid);
            scalar_t diff = scalar_t(0.0);

            auto diag = -1;
            Kokkos::parallel_reduce( Kokkos::ThreadVectorRange( team, soffset, eoffset ), [&] ( const long ptr, scalar_t &tdiff ) {
              auto colid = entries(ptr);
              auto val   = values(ptr);
              if ( colid != rowid ) {
                tdiff = tdiff - val*lhs(colid);
              }
              else {
                diag = ptr;
              }
            }, diff );

            // ASSUMPTION: sorted diagonal value located at eoffset - 1
            lhs(rowid) = (rhs_rowid+diff)/values(diag);
          } // end if
        }); // end TeamThreadRange

        team.team_barrier();
  }
};

template <class RowMapType, class EntriesType, class ValuesType, class LHSType, class RHSType, class NGBLType>
struct UpperTriLvlSchedRPSolverFunctor
{
  typedef typename EntriesType::non_const_value_type lno_t;
  RowMapType row_map;
  EntriesType entries;
  ValuesType values;
  LHSType lhs;
  RHSType rhs;
  NGBLType nodes_grouped_by_level;


  UpperTriLvlSchedRPSolverFunctor( const RowMapType &row_map_, const EntriesType &entries_, const ValuesType &values_, LHSType &lhs_, const RHSType &rhs_, const NGBLType &nodes_grouped_by_level_ ) :
    row_map(row_map_), entries(entries_), values(values_), lhs(lhs_), rhs(rhs_), nodes_grouped_by_level(nodes_grouped_by_level_) {}


  KOKKOS_INLINE_FUNCTION
  void operator()(const lno_t i) const {
    auto rowid = nodes_grouped_by_level(i);
    // Assuming indices are sorted per row, diag entry is final index in the list
    long soffset = row_map(rowid);
    long eoffset = row_map(rowid+1);
    auto rhs_rowid = rhs(rowid);
    for ( long ptr = eoffset-1; ptr >= soffset; --ptr ) {
      auto colid = entries(ptr);
      auto val   = values(ptr);
      if ( colid != rowid ) {
        rhs_rowid = rhs_rowid - val*lhs(colid);
      }
      else {
        lhs(rowid) = rhs_rowid/val;
      }
    } // end for ptr
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const UnsortedTag&, const lno_t i) const {
    auto rowid = nodes_grouped_by_level(i);
    long soffset = row_map(rowid);
    long eoffset = row_map(rowid+1);
    auto rhs_rowid = rhs(rowid);
    auto diag = -1;
    for ( long ptr = eoffset-1; ptr >= soffset; --ptr ) {
      auto colid = entries(ptr);
      auto val   = values(ptr);
      if ( colid != rowid ) {
        rhs_rowid = rhs_rowid - val*lhs(colid);
      }
      else {
        diag = ptr;
      }
    } // end for ptr
    lhs(rowid) = rhs_rowid/values(diag);
  }

};


template <class RowMapType, class EntriesType, class ValuesType, class LHSType, class RHSType, class NGBLType>
struct UpperTriLvlSchedTP1SolverFunctor
{
  typedef typename RowMapType::execution_space execution_space;
  typedef Kokkos::TeamPolicy<execution_space> policy_type;
  typedef typename policy_type::member_type member_type;
  typedef typename EntriesType::non_const_value_type lno_t;
  typedef typename ValuesType::non_const_value_type scalar_t;

  RowMapType row_map;
  EntriesType entries;
  ValuesType values;
  LHSType lhs;
  RHSType rhs;
  NGBLType nodes_grouped_by_level;

  long node_count; // like "block" offset into ngbl, my_league is the "local" offset
  long node_groups;


  UpperTriLvlSchedTP1SolverFunctor( const RowMapType &row_map_, const EntriesType &entries_, const ValuesType &values_, LHSType &lhs_, const RHSType &rhs_, const NGBLType &nodes_grouped_by_level_, long node_count_, long node_groups_ = 0 ) :
    row_map(row_map_), entries(entries_), values(values_), lhs(lhs_), rhs(rhs_), nodes_grouped_by_level(nodes_grouped_by_level_), node_count(node_count_), node_groups(node_groups_) {}


  KOKKOS_INLINE_FUNCTION
  void operator()(const member_type & team) const {
        auto my_league = team.league_rank(); // map to rowid
        auto rowid = nodes_grouped_by_level(my_league + node_count);
        auto my_rank = team.team_rank();

        auto soffset = row_map(rowid);
        auto eoffset = row_map(rowid+1);
        auto rhs_rowid = rhs(rowid);
        scalar_t diff = scalar_t(0.0);

        Kokkos::parallel_reduce( Kokkos::TeamThreadRange( team, soffset, eoffset ), [&] ( const long ptr, scalar_t &tdiff ) {
          auto colid = entries(ptr);
          auto val   = values(ptr);
          if ( colid != rowid ) {
            tdiff = tdiff - val*lhs(colid);
          }
        }, diff );

        team.team_barrier();

        // At end, finalize rowid == colid
        // only one thread should do this, also can use Kokkos::single
        if ( my_rank == 0 )
        {
        // ASSUMPTION: sorted diagonal value located at start offset
          lhs(rowid) = (rhs_rowid+diff)/values(soffset);
        }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const UnsortedTag&, const member_type & team) const {
        auto my_league = team.league_rank(); // map to rowid
        auto rowid = nodes_grouped_by_level(my_league + node_count);
        auto my_rank = team.team_rank();

        auto soffset = row_map(rowid);
        auto eoffset = row_map(rowid+1);
        auto rhs_rowid = rhs(rowid);
        scalar_t diff = scalar_t(0.0);

        auto diag = -1;

        Kokkos::parallel_reduce( Kokkos::TeamThreadRange( team, soffset, eoffset ), [&] ( const long ptr, scalar_t &tdiff ) {
          auto colid = entries(ptr);
          auto val   = values(ptr);
          if ( colid != rowid ) {
            tdiff = tdiff - val*lhs(colid);
          }
          else {
            diag = ptr;
          }
        }, diff );
        team.team_barrier();

        // At end, finalize rowid == colid
        // only one thread should do this, also can use Kokkos::single
        if ( my_rank == 0 )
        {
          lhs(rowid) = (rhs_rowid+diff)/values(diag);
        }
  }
};


// FIXME CUDA: This algorithm not working with all integral type combos
// In any case, this serves as a skeleton for 3-level hierarchical parallelism for alg dev
template <class RowMapType, class EntriesType, class ValuesType, class LHSType, class RHSType, class NGBLType>
struct UpperTriLvlSchedTP2SolverFunctor
{
  typedef typename RowMapType::execution_space execution_space;
  typedef Kokkos::TeamPolicy<execution_space> policy_type;
  typedef typename policy_type::member_type member_type;
  typedef typename EntriesType::non_const_value_type lno_t;
  typedef typename ValuesType::non_const_value_type scalar_t;

  RowMapType row_map;
  EntriesType entries;
  ValuesType values;
  LHSType lhs;
  RHSType rhs;
  NGBLType nodes_grouped_by_level;

  long node_count; // like "block" offset into ngbl, my_league is the "local" offset
  long node_groups;


  UpperTriLvlSchedTP2SolverFunctor(const RowMapType &row_map_, const EntriesType &entries_, const ValuesType &values_, LHSType &lhs_, const RHSType &rhs_, const NGBLType &nodes_grouped_by_level_, long node_count_, long node_groups_ = 0) :
    row_map(row_map_), entries(entries_), values(values_), lhs(lhs_), rhs(rhs_), nodes_grouped_by_level(nodes_grouped_by_level_), node_count(node_count_), node_groups(node_groups_) {}


  KOKKOS_INLINE_FUNCTION
  void operator()(const member_type & team) const {
        auto my_league = team.league_rank(); // map to rowid

        size_t nrows = row_map.extent(0) - 1;

        Kokkos::parallel_for( Kokkos::TeamThreadRange( team, 0, node_groups ), [&] ( const long ng ) {
          auto rowid = nodes_grouped_by_level(node_count + my_league*node_groups + ng);
          if ( size_t(rowid) < nrows ) {

            auto soffset = row_map(rowid);
            auto eoffset = row_map(rowid+1);
            auto rhs_rowid = rhs(rowid);
            scalar_t diff = scalar_t(0.0);

            Kokkos::parallel_reduce( Kokkos::ThreadVectorRange( team, soffset, eoffset ), [&] ( const long ptr, scalar_t &tdiff ) {
              auto colid = entries(ptr);
              auto val   = values(ptr);
              if ( colid != rowid ) {
                tdiff = tdiff - val*lhs(colid);
              }
            }, diff );

            // ASSUMPTION: sorted diagonal value located at start offset
            lhs(rowid) = (rhs_rowid+diff)/values(soffset);
          } // end if
        }); // end TeamThreadRange

        team.team_barrier();
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const UnsortedTag&, const member_type & team ) const {
        auto my_league = team.league_rank(); // map to rowid

        size_t nrows = row_map.extent(0) - 1;

        Kokkos::parallel_for( Kokkos::TeamThreadRange( team, 0, node_groups ), [&] ( const long ng ) {
          auto rowid = nodes_grouped_by_level(node_count + my_league*node_groups + ng);
          if ( size_t(rowid) < nrows ) {
            auto soffset = row_map(rowid);
            auto eoffset = row_map(rowid+1);
            auto rhs_rowid = rhs(rowid);
            scalar_t diff = scalar_t(0.0);

            auto diag = -1;
            Kokkos::parallel_reduce( Kokkos::ThreadVectorRange( team, soffset, eoffset ), [&] ( const long ptr, scalar_t &tdiff ) {
              auto colid = entries(ptr);
              auto val   = values(ptr);
              if ( colid != rowid ) {
                tdiff = tdiff - val*lhs(colid);
              }
              else {
                diag = ptr;
              }
            }, diff );

            lhs(rowid) = (rhs_rowid+diff)/values(diag);
          } // end if
        }); // end TeamThreadRange

        team.team_barrier();
  }

};


// --------------------------------
// Single-block functors
// --------------------------------

template <class RowMapType, class EntriesType, class ValuesType, class LHSType, class RHSType, class NGBLType>
struct LowerTriLvlSchedTP1SingleBlockFunctor
{
  typedef typename RowMapType::execution_space execution_space;
  typedef Kokkos::TeamPolicy<execution_space> policy_type;
  typedef typename policy_type::member_type member_type;
  typedef typename EntriesType::non_const_value_type lno_t;
  typedef typename ValuesType::non_const_value_type scalar_t;

  RowMapType row_map;
  EntriesType entries;
  ValuesType values;
  LHSType lhs;
  RHSType rhs;
  NGBLType nodes_grouped_by_level;
  NGBLType nodes_per_level;

  long node_count; // like "block" offset into ngbl, my_league is the "local" offset
  long lvl_start;
  long lvl_end;
  long cutoff;
  // team_size: each team can be assigned a row, if there are enough rows...


  LowerTriLvlSchedTP1SingleBlockFunctor( const RowMapType &row_map_, const EntriesType &entries_, const ValuesType &values_, LHSType &lhs_, const RHSType &rhs_, const NGBLType &nodes_grouped_by_level_, NGBLType &nodes_per_level_, long node_count_, long lvl_start_, long lvl_end_, long cutoff_ = 0 ) :
    row_map(row_map_), entries(entries_), values(values_), lhs(lhs_), rhs(rhs_), nodes_grouped_by_level(nodes_grouped_by_level_), nodes_per_level(nodes_per_level_), node_count(node_count_), lvl_start(lvl_start_), lvl_end(lvl_end_), cutoff(cutoff_) {}

  // SingleBlock: Only one block (or league) executing; team_rank used to map thread to row

  KOKKOS_INLINE_FUNCTION
  void operator()( const member_type & team ) const {
    long mut_node_count = node_count;
    typename NGBLType::non_const_value_type rowid {0};
    typename RowMapType::non_const_value_type soffset {0};
    typename RowMapType::non_const_value_type eoffset {0};
    typename RHSType::non_const_value_type rhs_val {0};
    scalar_t diff = scalar_t(0.0);
    for ( auto lvl = lvl_start; lvl < lvl_end; ++lvl ) {
      auto nodes_this_lvl = nodes_per_level(lvl);
      int my_rank = team.team_rank();
      diff = scalar_t(0.0);

      if (my_rank < nodes_this_lvl) {
        // THIS is where the mapping of threadid to rowid happens
        rowid = nodes_grouped_by_level(my_rank + mut_node_count);

        soffset = row_map(rowid);
        eoffset = row_map(rowid+1);
        rhs_val = rhs(rowid);

#ifdef SERIAL_FOR_LOOP
        for (auto ptr = soffset; ptr < eoffset; ++ptr) {
          auto colid = entries(ptr);
          auto val   = values(ptr);
          if ( colid != rowid ) {
            diff -= val*lhs(colid);
          }
        }
#else
        auto trange = eoffset - soffset;
        Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team, trange), [&] (const int loffset, scalar_t& tdiff)
        {
          auto ptr = soffset + loffset;
          auto colid = entries(ptr);
          auto val   = values(ptr);
          if ( colid != rowid ) {
            tdiff -= val*lhs(colid);
          }
        }, diff);
#endif
        // ASSUMPTION: sorted diagonal value located at eoffset - 1
        lhs(rowid) = (rhs_val+diff)/values(eoffset-1);
      } // end if team.team_rank() < nodes_this_lvl
      {
        // Update mut_node_count from nodes_per_level(lvl) each iteration of lvl per thread
        mut_node_count += nodes_this_lvl;
      }
      team.team_barrier();
    } // end for lvl
  } // end operator

  KOKKOS_INLINE_FUNCTION
  void operator()( const UnsortedTag&, const member_type & team ) const {
    long mut_node_count = node_count;
    typename NGBLType::non_const_value_type rowid {0};
    typename RowMapType::non_const_value_type soffset {0};
    typename RowMapType::non_const_value_type eoffset {0};
    typename RHSType::non_const_value_type rhs_val {0};
    scalar_t diff = scalar_t(0.0);
    for ( auto lvl = lvl_start; lvl < lvl_end; ++lvl ) {
      auto nodes_this_lvl = nodes_per_level(lvl);
      int my_rank = team.team_rank();
      diff = scalar_t(0.0);

      if (my_rank < nodes_this_lvl) {
        // THIS is where the mapping of threadid to rowid happens
        rowid = nodes_grouped_by_level(my_rank + mut_node_count);
        soffset = row_map(rowid);
        eoffset = row_map(rowid+1);
        rhs_val = rhs(rowid);

#ifdef SERIAL_FOR_LOOP
        for (auto ptr = soffset; ptr < eoffset; ++ptr) {
          auto colid = entries(ptr);
          auto val   = values(ptr);
          if ( colid != rowid ) {
            diff -= val*lhs(colid);
          }
        }
#else
        auto trange = eoffset - soffset;
        auto diag = -1;

        Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team, trange), [&] (const int loffset, scalar_t& tdiff)
        {
          auto ptr = soffset + loffset;

          auto colid = entries(ptr);
          auto val   = values(ptr);
          if ( colid != rowid ) {
            tdiff -= val*lhs(colid);
          }
          else {
            diag = ptr;
          }
        }, diff);
#endif
        lhs(rowid) = (rhs_val+diff)/values(diag);
      } // end if team.team_rank() < nodes_this_lvl
      {
        // Update mut_node_count from nodes_per_level(lvl) each iteration of lvl per thread
        mut_node_count += nodes_this_lvl;
      }
      team.team_barrier();
    } // end for lvl
  } // end operator


  KOKKOS_INLINE_FUNCTION
  void operator()( const LargerCutoffTag&, const member_type & team ) const {
    long mut_node_count = node_count;
    typename NGBLType::non_const_value_type rowid {0};
    typename RowMapType::non_const_value_type soffset {0};
    typename RowMapType::non_const_value_type eoffset {0};
    typename RHSType::non_const_value_type rhs_val {0};
    scalar_t diff = scalar_t(0.0);
    for ( auto lvl = lvl_start; lvl < lvl_end; ++lvl ) {
      auto nodes_this_lvl = nodes_per_level(lvl);
      int my_team_rank = team.team_rank();
      // If cutoff > team_size, then a thread will be responsible for multiple rows - this may be a helpful scenario depending on occupancy etc.
      for (int my_rank = my_team_rank; my_rank < cutoff; my_rank+=team.team_size() ) {
      diff = scalar_t(0.0);
      if (my_rank < nodes_this_lvl) {
        // THIS is where the mapping of threadid to rowid happens
        rowid = nodes_grouped_by_level(my_rank + mut_node_count);
        soffset = row_map(rowid);
        eoffset = row_map(rowid+1);
        rhs_val = rhs(rowid);

#ifdef SERIAL_FOR_LOOP
        for (auto ptr = soffset; ptr < eoffset; ++ptr) {
          auto colid = entries(ptr);
          auto val   = values(ptr);
          if ( colid != rowid ) {
            diff -= val*lhs(colid);
          }
        }
#else
        auto trange = eoffset - soffset;
        Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team, trange), [&] (const int loffset, scalar_t& tdiff)
        {
          auto ptr = soffset + loffset;
          auto colid = entries(ptr);
          auto val   = values(ptr);
          if ( colid != rowid ) {
            tdiff -= val*lhs(colid);
          }
        },diff);
#endif
         // ASSUMPTION: sorted diagonal value located at eoffset - 1 for lower tri, soffset for upper tri
         lhs(rowid) = (rhs_val+diff)/values(eoffset-1);
       } // end if team.team_rank() < nodes_this_lvl
      } // end for my_rank loop
      {
        // Update mut_node_count from nodes_per_level(lvl) each iteration of lvl per thread
        mut_node_count += nodes_this_lvl;
      }
      team.team_barrier();
    } // end for lvl
  } // end tagged operator

  KOKKOS_INLINE_FUNCTION
  void operator()( const UnsortedLargerCutoffTag&, const member_type & team ) const {
    long mut_node_count = node_count;
    typename NGBLType::non_const_value_type rowid {0};
    typename RowMapType::non_const_value_type soffset {0};
    typename RowMapType::non_const_value_type eoffset {0};
    typename RHSType::non_const_value_type rhs_val {0};
    scalar_t diff = scalar_t(0.0);

    for ( auto lvl = lvl_start; lvl < lvl_end; ++lvl ) {
      auto nodes_this_lvl = nodes_per_level(lvl);
      int my_team_rank = team.team_rank();
      // If cutoff > team_size, then a thread will be responsible for multiple rows - this may be a helpful scenario depending on occupancy etc.
      for (int my_rank = my_team_rank; my_rank < cutoff; my_rank+=team.team_size() ) {
      diff = scalar_t(0.0);
       if (my_rank < nodes_this_lvl) {
        // THIS is where the mapping of threadid to rowid happens
        rowid = nodes_grouped_by_level(my_rank + mut_node_count);
        soffset = row_map(rowid);
        eoffset = row_map(rowid+1);
        rhs_val = rhs(rowid);

#ifdef SERIAL_FOR_LOOP
        for (auto ptr = soffset; ptr < eoffset; ++ptr) {
          auto colid = entries(ptr);
          auto val   = values(ptr);
          if ( colid != rowid ) {
            diff -= val*lhs(colid);
          }
        }
#else
        auto trange = eoffset - soffset;
        auto diag = -1;

        Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team, trange), [&] (const int loffset, scalar_t& tdiff)
        {
          auto ptr = soffset + loffset;
          auto colid = entries(ptr);
          auto val   = values(ptr);
          if ( colid != rowid ) {
            tdiff -= val*lhs(colid);
          }
          else {
            diag = ptr;
          }
        },diff);
#endif
        lhs(rowid) = (rhs_val+diff)/values(diag);
       } // end if team.team_rank() < nodes_this_lvl
      } // end for my_rank loop
      {
        // Update mut_node_count from nodes_per_level(lvl) each iteration of lvl per thread
        mut_node_count += nodes_this_lvl;
      }
      team.team_barrier();
    } // end for lvl
  } // end tagged operator

};


template <class RowMapType, class EntriesType, class ValuesType, class LHSType, class RHSType, class NGBLType>
struct UpperTriLvlSchedTP1SingleBlockFunctor
{
  typedef typename RowMapType::execution_space execution_space;
  typedef Kokkos::TeamPolicy<execution_space> policy_type;
  typedef typename policy_type::member_type member_type;
  typedef typename EntriesType::non_const_value_type lno_t;
  typedef typename ValuesType::non_const_value_type scalar_t;

  RowMapType row_map;
  EntriesType entries;
  ValuesType values;
  LHSType lhs;
  RHSType rhs;
  NGBLType nodes_grouped_by_level;
  NGBLType nodes_per_level;

  long node_count; // like "block" offset into ngbl, my_league is the "local" offset
  long lvl_start;
  long lvl_end;
  long cutoff;
  // team_size: each team can be assigned a row, if there are enough rows...


  UpperTriLvlSchedTP1SingleBlockFunctor( const RowMapType &row_map_, const EntriesType &entries_, const ValuesType &values_, LHSType &lhs_, const RHSType &rhs_, const NGBLType &nodes_grouped_by_level_, NGBLType &nodes_per_level_, long node_count_, long lvl_start_, long lvl_end_, long cutoff_ = 0 ) :
    row_map(row_map_), entries(entries_), values(values_), lhs(lhs_), rhs(rhs_), nodes_grouped_by_level(nodes_grouped_by_level_), nodes_per_level(nodes_per_level_), node_count(node_count_), lvl_start(lvl_start_), lvl_end(lvl_end_), cutoff(cutoff_) {}

  // SingleBlock: Only one block (or league) executing; team_rank used to map thread to row

  KOKKOS_INLINE_FUNCTION
  void operator()( const member_type & team ) const {
    long mut_node_count = node_count;
    typename NGBLType::non_const_value_type rowid {0};
    typename RowMapType::non_const_value_type soffset {0};
    typename RowMapType::non_const_value_type eoffset {0};
    typename RHSType::non_const_value_type rhs_val {0};
    scalar_t diff = scalar_t(0.0);

    for ( auto lvl = lvl_start; lvl < lvl_end; ++lvl ) {
      auto nodes_this_lvl = nodes_per_level(lvl);
      int my_rank = team.team_rank();
      diff = scalar_t(0.0);

      if (my_rank < nodes_this_lvl) {
        // THIS is where the mapping of threadid to rowid happens
        rowid = nodes_grouped_by_level(my_rank + mut_node_count);
        soffset = row_map(rowid);
        eoffset = row_map(rowid+1);
        rhs_val = rhs(rowid);

#ifdef SERIAL_FOR_LOOP
        for (auto ptr = soffset; ptr < eoffset; ++ptr) {
          auto colid = entries(ptr);
          auto val   = values(ptr);
          if ( colid != rowid ) {
            diff -= val*lhs(colid);
          }
        }
#else
        auto trange = eoffset - soffset;
        Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team, trange), [&] (const int loffset, scalar_t& tdiff)
        {
          auto ptr = soffset + loffset;
          auto colid = entries(ptr);
          auto val   = values(ptr);
          if ( colid != rowid ) {
            tdiff -= val*lhs(colid);
          }
        }, diff);
#endif
        // ASSUMPTION: sorted diagonal value located at soffset
        lhs(rowid) = (rhs_val+diff)/values(soffset);
      } // end if
      {
        // Update mut_node_count from nodes_per_level(lvl) each iteration of lvl each thread
        mut_node_count += nodes_this_lvl;
      }
      team.team_barrier();
    } // end for lvl
  } // end operator

  KOKKOS_INLINE_FUNCTION
  void operator()( const UnsortedTag&, const member_type & team ) const {
    long mut_node_count = node_count;
    typename NGBLType::non_const_value_type rowid {0};
    typename RowMapType::non_const_value_type soffset {0};
    typename RowMapType::non_const_value_type eoffset {0};
    typename RHSType::non_const_value_type rhs_val {0};
    scalar_t diff = scalar_t(0.0);

    for ( auto lvl = lvl_start; lvl < lvl_end; ++lvl ) {
      auto nodes_this_lvl = nodes_per_level(lvl);
      int my_rank = team.team_rank();
      diff = scalar_t(0.0);

      if (my_rank < nodes_this_lvl) {
        // THIS is where the mapping of threadid to rowid happens
        rowid = nodes_grouped_by_level(my_rank + mut_node_count);
        soffset = row_map(rowid);
        eoffset = row_map(rowid+1);
        rhs_val = rhs(rowid);

#ifdef SERIAL_FOR_LOOP
        auto diag = -1;
        for (auto ptr = soffset; ptr < eoffset; ++ptr) {
          auto colid = entries(ptr);
          auto val   = values(ptr);
          if ( colid != rowid ) {
            diff -= val*lhs(colid);
          }
          else {
            diag = ptr;
          }
        }
#else
        auto trange = eoffset - soffset;
        auto diag = -1;

        Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team, trange), [&] (const int loffset, scalar_t& tdiff)
        {
          auto ptr = soffset + loffset;
          auto colid = entries(ptr);
          auto val   = values(ptr);
          if ( colid != rowid ) {
            tdiff -= val*lhs(colid);
          }
          else {
            diag = ptr;
          }
        }, diff);
#endif
        lhs(rowid) = (rhs_val+diff)/values(diag);
      } // end if
      {
        // Update mut_node_count from nodes_per_level(lvl) each iteration of lvl each thread
        mut_node_count += nodes_this_lvl;
      }
      team.team_barrier();
    } // end for lvl
  } // end operator


  KOKKOS_INLINE_FUNCTION
  void operator()( const LargerCutoffTag&, const member_type & team ) const {
    long mut_node_count = node_count;
    typename NGBLType::non_const_value_type rowid {0};
    typename RowMapType::non_const_value_type soffset {0};
    typename RowMapType::non_const_value_type eoffset {0};
    typename RHSType::non_const_value_type rhs_val {0};
    scalar_t diff = scalar_t(0.0);

    for ( auto lvl = lvl_start; lvl < lvl_end; ++lvl ) {
      auto nodes_this_lvl = nodes_per_level(lvl);
      int my_team_rank = team.team_rank();
      // If cutoff > team_size, then a thread will be responsible for multiple rows - this may be a helpful scenario depending on occupancy etc.
      for (int my_rank = my_team_rank; my_rank < cutoff; my_rank+=team.team_size() ) {
      diff = scalar_t(0.0);
       if (my_rank < nodes_this_lvl) {
        // THIS is where the mapping of threadid to rowid happens
        rowid = nodes_grouped_by_level(my_rank + mut_node_count);
        soffset = row_map(rowid);
        eoffset = row_map(rowid+1);
        rhs_val = rhs(rowid);

#ifdef SERIAL_FOR_LOOP
        for (auto ptr = soffset; ptr < eoffset; ++ptr) {
          auto colid = entries(ptr);
          auto val   = values(ptr);
          if ( colid != rowid ) {
            diff -= val*lhs(colid);
          }
        }
#else
        auto trange = eoffset - soffset;
        Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team, trange), [&] (const int loffset, scalar_t& tdiff)
        {
          auto ptr = soffset + loffset;
          auto colid = entries(ptr);
          auto val   = values(ptr);
          if ( colid != rowid ) {
            tdiff -= val*lhs(colid);
          }
        }, diff);
#endif
        // ASSUMPTION: sorted diagonal value located at eoffset - 1 for lower tri, soffset for upper tri
          lhs(rowid) = (rhs_val+diff)/values(soffset);
       } // end if team.team_rank() < nodes_this_lvl
      } // end for my_rank loop
      {
        // Update mut_node_count from nodes_per_level(lvl) each iteration of lvl per thread
        mut_node_count += nodes_this_lvl;
      }
      team.team_barrier();
    } // end for lvl
  } // end tagged operator

  KOKKOS_INLINE_FUNCTION
  void operator()( const UnsortedLargerCutoffTag&, const member_type & team ) const {
    long mut_node_count = node_count;
    typename NGBLType::non_const_value_type rowid {0};
    typename RowMapType::non_const_value_type soffset {0};
    typename RowMapType::non_const_value_type eoffset {0};
    typename RHSType::non_const_value_type rhs_val {0};
    scalar_t diff = scalar_t(0.0);

    for ( auto lvl = lvl_start; lvl < lvl_end; ++lvl ) {
      auto nodes_this_lvl = nodes_per_level(lvl);
      int my_team_rank = team.team_rank();
      // If cutoff > team_size, then a thread will be responsible for multiple rows - this may be a helpful scenario depending on occupancy etc.
      for (int my_rank = my_team_rank; my_rank < cutoff; my_rank+=team.team_size() ) {
       diff = scalar_t(0.0);
       if (my_rank < nodes_this_lvl) {
        // THIS is where the mapping of threadid to rowid happens
        rowid = nodes_grouped_by_level(my_rank + mut_node_count);
        soffset = row_map(rowid);
        eoffset = row_map(rowid+1);
        rhs_val = rhs(rowid);

#ifdef SERIAL_FOR_LOOP
        auto diag = -1;
        for (auto ptr = soffset; ptr < eoffset; ++ptr) {
          auto colid = entries(ptr);
          auto val   = values(ptr);
          if ( colid != rowid ) {
            diff -= val*lhs(colid);
          }
          else {
            diag = ptr;
          }
        }
#else
        auto trange = eoffset - soffset;
        auto diag = -1;
        Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team, trange), [&] (const int loffset, scalar_t& tdiff)
        {
          auto ptr = soffset + loffset;
          auto colid = entries(ptr);
          auto val   = values(ptr);
          if ( colid != rowid ) {
            tdiff -= val*lhs(colid);
          }
          else {
            diag = ptr;
          }
        }, diff);
#endif
        lhs(rowid) = (rhs_val+diff)/values(diag);
       } // end if team.team_rank() < nodes_this_lvl
      } // end for my_rank loop
      {
        // Update mut_node_count from nodes_per_level(lvl) each iteration of lvl per thread
        mut_node_count += nodes_this_lvl;
      }
      team.team_barrier();
    } // end for lvl
  } // end tagged operator
};


template <class RowMapType, class EntriesType, class ValuesType, class LHSType, class RHSType, class NGBLType>
struct TriLvlSchedTP1SingleBlockFunctor
{
  typedef typename RowMapType::execution_space execution_space;
  typedef Kokkos::TeamPolicy<execution_space> policy_type;
  typedef typename policy_type::member_type member_type;
  typedef typename EntriesType::non_const_value_type lno_t;
  typedef typename ValuesType::non_const_value_type scalar_t;

  RowMapType row_map;
  EntriesType entries;
  ValuesType values;
  LHSType lhs;
  RHSType rhs;
  NGBLType nodes_grouped_by_level;
  NGBLType nodes_per_level;

  long node_count; // like "block" offset into ngbl, my_league is the "local" offset
  long lvl_start;
  long lvl_end;
  const bool is_lowertri;
  const int dense_nrows;
  const int  cutoff;
  // team_size: each team can be assigned a row, if there are enough rows...


  TriLvlSchedTP1SingleBlockFunctor( const RowMapType &row_map_, const EntriesType &entries_, const ValuesType &values_, LHSType &lhs_, const RHSType &rhs_, const NGBLType &nodes_grouped_by_level_, NGBLType &nodes_per_level_, long node_count_, long lvl_start_, long lvl_end_, const bool is_lower_, const int dense_nrows_ = 0, const int cutoff_ = 0 ) :
    row_map(row_map_), entries(entries_), values(values_), lhs(lhs_), rhs(rhs_), nodes_grouped_by_level(nodes_grouped_by_level_), nodes_per_level(nodes_per_level_), node_count(node_count_), lvl_start(lvl_start_), lvl_end(lvl_end_), is_lowertri(is_lower_), dense_nrows(dense_nrows_), cutoff(cutoff_) {}

  // SingleBlock: Only one block (or league) executing; team_rank used to map thread to row

  KOKKOS_INLINE_FUNCTION
  void operator()( const member_type & team ) const {
    long mut_node_count = node_count;
    typename NGBLType::non_const_value_type rowid {0};
    typename RowMapType::non_const_value_type soffset {0};
    typename RowMapType::non_const_value_type eoffset {0};
    typename RHSType::non_const_value_type rhs_val {0};
    scalar_t diff = scalar_t(0.0);

    for ( auto lvl = lvl_start; lvl < lvl_end; ++lvl ) {
      auto nodes_this_lvl = nodes_per_level(lvl);
      int my_rank = team.team_rank();
      diff = scalar_t(0.0);

      if (my_rank < nodes_this_lvl) {
        // THIS is where the mapping of threadid to rowid happens
        rowid = nodes_grouped_by_level(my_rank + mut_node_count);
        soffset = row_map(rowid);
        eoffset = row_map(rowid+1);
        rhs_val = rhs(rowid);

#ifdef SERIAL_FOR_LOOP
        for (auto ptr = soffset; ptr < eoffset; ++ptr) {
          auto colid = entries(ptr);
          auto val   = values(ptr);
          if ( colid != rowid ) {
            diff -= val*lhs(colid);
          }
        }
#else
        auto trange = eoffset - soffset;
        Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team, trange), [&] (const int loffset, scalar_t& tdiff)
        {
          auto ptr = soffset + loffset;
          auto colid = entries(ptr);
          auto val   = values(ptr);
          if ( colid != rowid ) {
            tdiff -= val*lhs(colid);
          }
        }, diff);
#endif

        // ASSUMPTION: sorted diagonal value located at eoffset - 1 for lower tri, soffset for upper tri
        if (is_lowertri)
          lhs(rowid) = (rhs_val+diff)/values(eoffset-1);
        else
          lhs(rowid) = (rhs_val+diff)/values(soffset);
      } // end if team.team_rank() < nodes_this_lvl
      {
        // Update mut_node_count from nodes_per_level(lvl) each iteration of lvl per thread
        mut_node_count += nodes_this_lvl;
      }
      team.team_barrier();
    } // end for lvl
  } // end operator

  KOKKOS_INLINE_FUNCTION
  void operator()( const UnsortedTag&, const member_type & team ) const {
    long mut_node_count = node_count;
    typename NGBLType::non_const_value_type rowid {0};
    typename RowMapType::non_const_value_type soffset {0};
    typename RowMapType::non_const_value_type eoffset {0};
    typename RHSType::non_const_value_type rhs_val {0};
    scalar_t diff = scalar_t(0.0);

    for ( auto lvl = lvl_start; lvl < lvl_end; ++lvl ) {
      auto nodes_this_lvl = nodes_per_level(lvl);
      int my_rank = team.team_rank();
      diff = scalar_t(0.0);

      if (my_rank < nodes_this_lvl) {
        // THIS is where the mapping of threadid to rowid happens
        rowid = nodes_grouped_by_level(my_rank + mut_node_count);
        soffset = row_map(rowid);
        eoffset = row_map(rowid+1);
        rhs_val = rhs(rowid);

#ifdef SERIAL_FOR_LOOP
        auto diag = -1;
        for (auto ptr = soffset; ptr < eoffset; ++ptr) {
          auto colid = entries(ptr);
          auto val   = values(ptr);
          if ( colid != rowid ) {
            diff -= val*lhs(colid);
          }
          else {
            diag = ptr;
          }
        }
#else
        auto trange = eoffset - soffset;
        auto diag = -1;
        Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team, trange), [&] (const int loffset, scalar_t& tdiff)
        {
          auto ptr = soffset + loffset;
          auto colid = entries(ptr);
          auto val   = values(ptr);
          if ( colid != rowid ) {
            tdiff -= val*lhs(colid);
          }
          else {
            diag = ptr;
          }
        }, diff);
#endif
        lhs(rowid) = (rhs_val+diff)/values(diag);
      } // end if team.team_rank() < nodes_this_lvl
      {
        // Update mut_node_count from nodes_per_level(lvl) each iteration of lvl per thread
        mut_node_count += nodes_this_lvl;
      }
      team.team_barrier();
    } // end for lvl
  } // end operator


  KOKKOS_INLINE_FUNCTION
  void operator()( const LargerCutoffTag&, const member_type & team ) const {
    long mut_node_count = node_count;
    typename NGBLType::non_const_value_type rowid {0};
    typename RowMapType::non_const_value_type soffset {0};
    typename RowMapType::non_const_value_type eoffset {0};
    typename RHSType::non_const_value_type rhs_val {0};
    scalar_t diff = scalar_t(0.0);

    for ( auto lvl = lvl_start; lvl < lvl_end; ++lvl ) {
      auto nodes_this_lvl = nodes_per_level(lvl);
      int my_team_rank = team.team_rank();
      // If cutoff > team_size, then a thread will be responsible for multiple rows - this may be a helpful scenario depending on occupancy etc.
      for (int my_rank = my_team_rank; my_rank < cutoff; my_rank+=team.team_size() ) {
       diff = scalar_t(0.0);
       if (my_rank < nodes_this_lvl) {
        // THIS is where the mapping of threadid to rowid happens
        rowid = nodes_grouped_by_level(my_rank + mut_node_count);
        soffset = row_map(rowid);
        eoffset = row_map(rowid+1);
        rhs_val = rhs(rowid);

#ifdef SERIAL_FOR_LOOP
        for (auto ptr = soffset; ptr < eoffset; ++ptr) {
          auto colid = entries(ptr);
          auto val   = values(ptr);
          if ( colid != rowid ) {
            diff -= val*lhs(colid);
          }
        }
#else
        auto trange = eoffset - soffset;
        Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team, trange), [&] (const int loffset, scalar_t& tdiff)
        {
          auto ptr = soffset + loffset;
          auto colid = entries(ptr);
          auto val   = values(ptr);
          if ( colid != rowid ) {
            tdiff -= val*lhs(colid);
          }
        }, diff);
#endif

        // ASSUMPTION: sorted diagonal value located at eoffset - 1 for lower tri, soffset for upper tri
        if (is_lowertri)
          lhs(rowid) = (rhs_val+diff)/values(eoffset-1);
        else
          lhs(rowid) = (rhs_val+diff)/values(soffset);
       } // end if team.team_rank() < nodes_this_lvl
      } // end for my_rank loop
      {
        // Update mut_node_count from nodes_per_level(lvl) each iteration of lvl per thread
        mut_node_count += nodes_this_lvl;
      }
      team.team_barrier();
    } // end for lvl
  } // end tagged operator

  KOKKOS_INLINE_FUNCTION
  void operator()( const UnsortedLargerCutoffTag&, const member_type & team ) const {
    long mut_node_count = node_count;
    typename NGBLType::non_const_value_type rowid {0};
    typename RowMapType::non_const_value_type soffset {0};
    typename RowMapType::non_const_value_type eoffset {0};
    typename RHSType::non_const_value_type rhs_val {0};
    scalar_t diff = scalar_t(0.0);

    for ( auto lvl = lvl_start; lvl < lvl_end; ++lvl ) {
      auto nodes_this_lvl = nodes_per_level(lvl);
      int my_team_rank = team.team_rank();
      // If cutoff > team_size, then a thread will be responsible for multiple rows - this may be a helpful scenario depending on occupancy etc.
      for (int my_rank = my_team_rank; my_rank < cutoff; my_rank+=team.team_size() ) {
       diff = scalar_t(0.0);
       if (my_rank < nodes_this_lvl) {
        // THIS is where the mapping of threadid to rowid happens
        rowid = nodes_grouped_by_level(my_rank + mut_node_count);
        soffset = row_map(rowid);
        eoffset = row_map(rowid+1);
        rhs_val = rhs(rowid);

#ifdef SERIAL_FOR_LOOP
        auto diag = -1;
        for (auto ptr = soffset; ptr < eoffset; ++ptr) {
          auto colid = entries(ptr);
          auto val   = values(ptr);
          if ( colid != rowid ) {
            diff -= val*lhs(colid);
          }
          else {
            diag = ptr;
          }
        }
#else
        auto trange = eoffset - soffset;
        auto diag = -1;
        Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team, trange), [&] (const int loffset, scalar_t& tdiff)
        {
          auto ptr = soffset + loffset;
          auto colid = entries(ptr);
          auto val   = values(ptr);
          if ( colid != rowid ) {
            tdiff -= val*lhs(colid);
          }
          else {
            diag = ptr;
          }
        }, diff);
#endif
          lhs(rowid) = (rhs_val+diff)/values(diag);
       } // end if team.team_rank() < nodes_this_lvl
      } // end for my_rank loop
      {
        // Update mut_node_count from nodes_per_level(lvl) each iteration of lvl per thread
        mut_node_count += nodes_this_lvl;
      }
      team.team_barrier();
    } // end for lvl
  } // end tagged operator

};


template <class RowMapType, class EntriesType, class ValuesType, class LHSType, class RHSType, class NGBLType>
struct TriLvlSchedTP1SingleBlockFunctorDiagValues
{
  typedef typename RowMapType::execution_space execution_space;
  typedef Kokkos::TeamPolicy<execution_space> policy_type;
  typedef typename policy_type::member_type member_type;
  typedef typename EntriesType::non_const_value_type lno_t;
  typedef typename ValuesType::non_const_value_type scalar_t;

  RowMapType row_map;
  EntriesType entries;
  ValuesType values;
  LHSType lhs;
  RHSType rhs;
  NGBLType nodes_grouped_by_level;
  NGBLType nodes_per_level;
  ValuesType diagonal_values;

  long node_count; // like "block" offset into ngbl, my_league is the "local" offset
  long lvl_start;
  long lvl_end;
  const bool is_lowertri;
  const int dense_nrows;
  const int  cutoff;
  // team_size: each team can be assigned a row, if there are enough rows...


  TriLvlSchedTP1SingleBlockFunctorDiagValues( const RowMapType &row_map_, const EntriesType &entries_, const ValuesType &values_, LHSType &lhs_, const RHSType &rhs_, const NGBLType &nodes_grouped_by_level_, const NGBLType &nodes_per_level_, const ValuesType &diagonal_values_, long node_count_, const long lvl_start_, const long lvl_end_, const bool is_lower_, const int dense_nrows_ = 0, const int cutoff_ = 0 ) :
    row_map(row_map_), entries(entries_), values(values_), lhs(lhs_), rhs(rhs_), nodes_grouped_by_level(nodes_grouped_by_level_), nodes_per_level(nodes_per_level_), diagonal_values(diagonal_values_), node_count(node_count_), lvl_start(lvl_start_), lvl_end(lvl_end_), is_lowertri(is_lower_), dense_nrows(dense_nrows_), cutoff(cutoff_) {}

  // SingleBlock: Only one block (or league) executing; team_rank used to map thread to row

  KOKKOS_INLINE_FUNCTION
  void operator()( const member_type & team ) const {
    long mut_node_count = node_count;
    typename NGBLType::non_const_value_type rowid {0};
    typename RowMapType::non_const_value_type soffset {0};
    typename RowMapType::non_const_value_type eoffset {0};
    typename RHSType::non_const_value_type rhs_val {0};
    scalar_t diff = scalar_t(0.0);

    for ( auto lvl = lvl_start; lvl < lvl_end; ++lvl ) {
      auto nodes_this_lvl = nodes_per_level(lvl);
      int my_rank = team.team_rank();
      diff = scalar_t(0.0);

      if (my_rank < nodes_this_lvl) {
        // THIS is where the mapping of threadid to rowid happens
        rowid = nodes_grouped_by_level(my_rank + mut_node_count);
        soffset = row_map(rowid);
        eoffset = row_map(rowid+1);
        rhs_val = rhs(rowid);

#ifdef SERIAL_FOR_LOOP
        for (auto ptr = soffset; ptr < eoffset; ++ptr) {
          auto colid = entries(ptr);
          auto val   = values(ptr);
          if ( colid != rowid ) {
            diff -= val*lhs(colid);
          }
        }
#else
      auto trange = eoffset - soffset;
      Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team, trange), [&] (const int loffset, scalar_t& tdiff)
      {
        auto ptr = soffset + loffset;
        auto colid = entries(ptr);
        auto val   = values(ptr);

        if ( colid != rowid ) {
          tdiff -= val*lhs(colid);
        }
      }, diff);
#endif
        // ASSUMPTION: sorted diagonal value located at eoffset - 1 for lower tri, soffset for upper tri
        lhs(rowid) = (rhs_val+diff)/diagonal_values(rowid);
      } // end if team.team_rank() < nodes_this_lvl
      {
        // Update mut_node_count from nodes_per_level(lvl) each iteration of lvl per thread
        mut_node_count += nodes_this_lvl;
      }
      team.team_barrier();
    } // end for lvl
  } // end operator


  KOKKOS_INLINE_FUNCTION
  void operator()( const LargerCutoffTag&, const member_type & team ) const {
    long mut_node_count = node_count;
    typename NGBLType::non_const_value_type rowid {0};
    typename RowMapType::non_const_value_type soffset {0};
    typename RowMapType::non_const_value_type eoffset {0};
    typename RHSType::non_const_value_type rhs_val {0};
    scalar_t diff = scalar_t(0.0);

    for ( auto lvl = lvl_start; lvl < lvl_end; ++lvl ) {
      auto nodes_this_lvl = nodes_per_level(lvl);
      int my_team_rank = team.team_rank();
      // If cutoff > team_size, then a thread will be responsible for multiple rows - this may be a helpful scenario depending on occupancy etc.
      for (int my_rank = my_team_rank; my_rank < cutoff; my_rank+=team.team_size() ) {
       diff = scalar_t(0.0);
       if (my_rank < nodes_this_lvl) {
        // THIS is where the mapping of threadid to rowid happens
        rowid = nodes_grouped_by_level(my_rank + mut_node_count);
        soffset = row_map(rowid);
        eoffset = row_map(rowid+1);
        rhs_val = rhs(rowid);

#ifdef SERIAL_FOR_LOOP
        for (auto ptr = soffset; ptr < eoffset; ++ptr) {
          auto colid = entries(ptr);
          auto val   = values(ptr);
          if ( colid != rowid ) {
            diff -= val*lhs(colid);
          }
        }
#else
      auto trange = eoffset - soffset;
      Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team, trange), [&] (const int loffset, scalar_t& tdiff)
      {
        auto ptr = soffset + loffset;
          auto colid = entries(ptr);
          auto val   = values(ptr);
          if ( colid != rowid ) {
            tdiff -= val*lhs(colid);
          }
        }, diff);
#endif
        lhs(rowid) = (rhs_val+diff)/diagonal_values(rowid);
       } // end if team.team_rank() < nodes_this_lvl
      } // end for my_rank loop
      {
        // Update mut_node_count from nodes_per_level(lvl) each iteration of lvl per thread
        mut_node_count += nodes_this_lvl;
      }
      team.team_barrier();
    } // end for lvl
  } // end tagged operator

};


#ifdef KOKKOSKERNELS_SPTRSV_CUDAGRAPHSUPPORT
template <class SpaceType>
struct ReturnTeamPolicyType;

#ifdef KOKKOS_ENABLE_SERIAL
template <>
struct ReturnTeamPolicyType<Kokkos::Serial> {
  using PolicyType = Kokkos::TeamPolicy<Kokkos::Serial>;

  static inline
  PolicyType get_policy(int nt, int ts) {
    return PolicyType(nt,ts);
  }

  template <class ExecInstanceType>
  static inline
  PolicyType get_policy(int nt, int ts, ExecInstanceType ) {
    return PolicyType(nt,ts);
    //return PolicyType(ExecInstanceType(),nt,ts);
  }
};
#endif
#ifdef KOKKOS_ENABLE_OPENMP
template <>
struct ReturnTeamPolicyType<Kokkos::OpenMP> {
  using PolicyType = Kokkos::TeamPolicy<Kokkos::OpenMP>;

  static inline
  PolicyType get_policy(int nt, int ts) {
    return PolicyType(nt,ts);
  }

  template <class ExecInstanceType>
  static inline
  PolicyType get_policy(int nt, int ts, ExecInstanceType ) {
    return PolicyType(nt,ts);
    //return PolicyType(ExecInstanceType(),nt,ts);
  }
};
#endif
#ifdef KOKKOS_ENABLE_CUDA
template <>
struct ReturnTeamPolicyType<Kokkos::Cuda> {
  using PolicyType = Kokkos::TeamPolicy<Kokkos::Cuda>;

  static inline
  PolicyType get_policy(int nt, int ts) {
    return PolicyType(nt,ts);
  }

  template <class ExecInstanceType>
  static inline
  PolicyType get_policy(int nt, int ts, ExecInstanceType stream) {
    return PolicyType(stream,nt,ts);
  }
};
#endif

template <class SpaceType>
struct ReturnRangePolicyType;

#ifdef KOKKOS_ENABLE_SERIAL
template <>
struct ReturnRangePolicyType<Kokkos::Serial> {
  using PolicyType = Kokkos::RangePolicy<Kokkos::Serial>;

  static inline
  PolicyType get_policy(int nt, int ts) {
    return PolicyType(nt,ts);
  }

  template <class ExecInstanceType>
  static inline
  PolicyType get_policy(int nt, int ts, ExecInstanceType ) {
    return PolicyType(nt,ts);
    //return PolicyType(ExecInstanceType(),nt,ts);
  }
};
#endif
#ifdef KOKKOS_ENABLE_OPENMP
template <>
struct ReturnRangePolicyType<Kokkos::OpenMP> {
  using PolicyType = Kokkos::RangePolicy<Kokkos::OpenMP>;

  static inline
  PolicyType get_policy(int nt, int ts) {
    return PolicyType(nt,ts);
  }

  template <class ExecInstanceType>
  static inline
  PolicyType get_policy(int nt, int ts, ExecInstanceType ) {
    return PolicyType(nt,ts);
    //return PolicyType(ExecInstanceType(),nt,ts);
  }
};
#endif
#ifdef KOKKOS_ENABLE_CUDA
template <>
struct ReturnRangePolicyType<Kokkos::Cuda> {
  using PolicyType = Kokkos::RangePolicy<Kokkos::Cuda>;

  static inline
  PolicyType get_policy(int nt, int ts) {
    return PolicyType(nt,ts);
  }

  template <class ExecInstanceType>
  static inline
  PolicyType get_policy(int nt, int ts, ExecInstanceType stream) {
    return PolicyType(stream,nt,ts);
  }
};
#endif

template < class TriSolveHandle, class RowMapType, class EntriesType, class ValuesType, class RHSType, class LHSType >
void lower_tri_solve_cg( TriSolveHandle & thandle, const RowMapType row_map, const EntriesType entries, const ValuesType values, const RHSType & rhs, LHSType &lhs) {

    typedef typename TriSolveHandle::nnz_lno_view_t NGBLType;
    typedef typename TriSolveHandle::execution_space execution_space;
    typedef typename TriSolveHandle::size_type size_type;
    typename TriSolveHandle::SPTRSVcudaGraphWrapperType* lcl_cudagraph = thandle.get_sptrsvCudaGraph();

    auto nlevels = thandle.get_num_levels();

    auto stream1 = lcl_cudagraph->stream;
    Kokkos::Cuda cuda1(stream1);
    auto graph = lcl_cudagraph->cudagraph;

    Kokkos::parallel_for("Init", Kokkos::RangePolicy<execution_space>(0,1), EmptyFunctor());
    Kokkos::Cuda().fence();
    cudaStreamSynchronize(stream1);
    //Kokkos::fence();

    auto hnodes_per_level = thandle.get_host_nodes_per_level();
    auto nodes_grouped_by_level = thandle.get_nodes_grouped_by_level();

    size_type node_count = 0;

    int team_size = thandle.get_team_size();
    team_size = team_size == -1 ? 64 : team_size;

    // Start capturing stream
    if(thandle.cudagraphCreated == false) {
    Kokkos::fence();
    cudaStreamBeginCapture(stream1, cudaStreamCaptureModeGlobal);
    {
      for (int iter = 0; iter < nlevels; ++iter) {
        size_type lvl_nodes = hnodes_per_level(iter);

        using policy_type = ReturnTeamPolicyType<execution_space>;

        Kokkos::parallel_for("parfor_l_team_cudagraph",  Kokkos::Experimental::require(ReturnTeamPolicyType<execution_space>::get_policy(lvl_nodes,team_size,cuda1), Kokkos::Experimental::WorkItemProperty::HintLightWeight), LowerTriLvlSchedTP1SolverFunctor<RowMapType, EntriesType, ValuesType, LHSType, RHSType, NGBLType>(row_map, entries, values, lhs, rhs, nodes_grouped_by_level, node_count));

        node_count += hnodes_per_level(iter);
      }
    }
    cudaStreamEndCapture(stream1, &graph);

    // Create graphExec
    cudaGraphInstantiate(&(lcl_cudagraph->cudagraphinstance), graph, NULL, NULL, 0);
      thandle.cudagraphCreated = true;
    }
    // Run graph
    Kokkos::fence();
    cudaGraphLaunch(lcl_cudagraph->cudagraphinstance, stream1);

    cudaStreamSynchronize(stream1);
    Kokkos::fence();
} // end lower_tri_solve_cg


template < class TriSolveHandle, class RowMapType, class EntriesType, class ValuesType, class RHSType, class LHSType >
void upper_tri_solve_cg( TriSolveHandle & thandle, const RowMapType row_map, const EntriesType entries, const ValuesType values, const RHSType & rhs, LHSType &lhs) {

    typedef typename TriSolveHandle::nnz_lno_view_t NGBLType;
    typedef typename TriSolveHandle::execution_space execution_space;
    typedef typename TriSolveHandle::size_type size_type;
    typename TriSolveHandle::SPTRSVcudaGraphWrapperType* lcl_cudagraph = thandle.get_sptrsvCudaGraph();

    auto nlevels = thandle.get_num_levels();

    auto stream1 = lcl_cudagraph->stream;
    Kokkos::Cuda cuda1(stream1);
    auto graph = lcl_cudagraph->cudagraph;

    Kokkos::parallel_for("Init", Kokkos::RangePolicy<execution_space>(0,1), EmptyFunctor());
    Kokkos::Cuda().fence();
    cudaStreamSynchronize(stream1);

    auto hnodes_per_level = thandle.get_host_nodes_per_level();
    auto nodes_grouped_by_level = thandle.get_nodes_grouped_by_level();

    size_type node_count = 0;

    int team_size = thandle.get_team_size();
    team_size = team_size == -1 ? 64 : team_size;

    // Start capturing stream
    if(thandle.cudagraphCreated == false) {
    Kokkos::fence();
    cudaStreamBeginCapture(stream1, cudaStreamCaptureModeGlobal);
    {
      for (int iter = 0; iter < nlevels; ++iter) {
        size_type lvl_nodes = hnodes_per_level(iter);

        using policy_type = ReturnTeamPolicyType<execution_space>;

        Kokkos::parallel_for("parfor_u_team_cudagraph",  Kokkos::Experimental::require(ReturnTeamPolicyType<execution_space>::get_policy(lvl_nodes,team_size,cuda1), Kokkos::Experimental::WorkItemProperty::HintLightWeight), UpperTriLvlSchedTP1SolverFunctor<RowMapType, EntriesType, ValuesType, LHSType, RHSType, NGBLType>(row_map, entries, values, lhs, rhs, nodes_grouped_by_level, node_count));

        node_count += hnodes_per_level(iter);
      }
    }
    cudaStreamEndCapture(stream1, &graph);

    // Create graphExec
    cudaGraphInstantiate(&(lcl_cudagraph->cudagraphinstance), graph, NULL, NULL, 0);
      thandle.cudagraphCreated = true;
    }
    // Run graph
    Kokkos::fence();
    cudaGraphLaunch(lcl_cudagraph->cudagraphinstance, stream1);

    cudaStreamSynchronize(stream1);
    Kokkos::fence();
} // end upper_tri_solve_cg

#endif


template < class TriSolveHandle, class RowMapType, class EntriesType, class ValuesType, class RHSType, class LHSType >
void lower_tri_solve(TriSolveHandle & thandle, const RowMapType row_map, const EntriesType entries, const ValuesType values, const RHSType & rhs, LHSType &lhs) {

#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOSPSTRSV_SOLVE_IMPL_PROFILE)
cudaProfilerStop();
#endif

  typedef typename TriSolveHandle::execution_space execution_space;
  typedef typename TriSolveHandle::size_type size_type;
  typedef typename TriSolveHandle::nnz_lno_view_t NGBLType;

  auto nlevels = thandle.get_num_levels();
  // Keep this a host View, create device version and copy to back to host during scheduling
  // This requires making sure the host view in the handle is properly updated after the symbolic phase
  auto nodes_per_level = thandle.get_nodes_per_level();
  auto hnodes_per_level = thandle.get_host_nodes_per_level();
  auto nodes_grouped_by_level = thandle.get_nodes_grouped_by_level();

  size_type node_count = 0;

  for ( size_type lvl = 0; lvl < nlevels; ++lvl ) {
   {
    size_type lvl_nodes = hnodes_per_level(lvl);

    if ( lvl_nodes != 0 ) {

#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOSPSTRSV_SOLVE_IMPL_PROFILE)
cudaProfilerStart();
#endif
      if ( thandle.get_algorithm() == KokkosSparse::Experimental::SPTRSVAlgorithm::SEQLVLSCHD_RP ) {
        Kokkos::parallel_for( "parfor_fixed_lvl", Kokkos::RangePolicy<execution_space>( node_count, node_count+lvl_nodes ), LowerTriLvlSchedRPSolverFunctor<RowMapType, EntriesType, ValuesType, LHSType, RHSType, NGBLType> (row_map, entries, values, lhs, rhs, nodes_grouped_by_level) );
      }
      else if ( thandle.get_algorithm() == KokkosSparse::Experimental::SPTRSVAlgorithm::SEQLVLSCHD_TP1 ) {
        typedef Kokkos::TeamPolicy<execution_space> policy_type;
        int team_size = thandle.get_team_size();

#ifdef KOKKOSKERNELS_SPTRSV_TRILVLSCHED
        TriLvlSchedTP1SolverFunctor<RowMapType, EntriesType, ValuesType, LHSType, RHSType, NGBLType> tstf(row_map, entries, values, lhs, rhs, nodes_grouped_by_level, true, node_count);
#else
        LowerTriLvlSchedTP1SolverFunctor<RowMapType, EntriesType, ValuesType, LHSType, RHSType, NGBLType> tstf(row_map, entries, values, lhs, rhs, nodes_grouped_by_level, node_count);
#endif
        if ( team_size == -1 )
          Kokkos::parallel_for("parfor_l_team", policy_type( lvl_nodes , Kokkos::AUTO ), tstf);
        else
          Kokkos::parallel_for("parfor_l_team", policy_type( lvl_nodes , team_size ), tstf);
      }
      // TP2 algorithm has issues with some offset-ordinal combo to be addressed
      /*
      else if ( thandle.get_algorithm() == KokkosSparse::Experimental::SPTRSVAlgorithm::SEQLVLSCHED_TP2 ) {
        typedef Kokkos::TeamPolicy<execution_space> tvt_policy_type;

        int team_size = thandle.get_team_size();
        if ( team_size == -1 ) {
          team_size = std::is_same< typename Kokkos::DefaultExecutionSpace::memory_space, Kokkos::HostSpace >::value ? 1 : 64;
        }
        int vector_size = thandle.get_team_size();
        if ( vector_size == -1 ) {
          vector_size = std::is_same< typename Kokkos::DefaultExecutionSpace::memory_space, Kokkos::HostSpace >::value ? 1 : 4;
        }

        // This impl: "chunk" lvl_nodes into node_groups; a league_rank is responsible for processing team_size # nodes
        //       TeamThreadRange over number nodes of node_groups
        //       To avoid masking threads, 1 thread (team) per node in node_group (thread has full ownership of a node)
        //       ThreadVectorRange responsible for the actual solve computation
        //const int node_groups = team_size;
        const int node_groups = vector_size;

#ifdef KOKKOSKERNELS_SPTRSV_TRILVLSCHED
        TriLvlSchedTP2SolverFunctor<RowMapType, EntriesType, ValuesType, LHSType, RHSType, NGBLType> tstf(row_map, entries, values, lhs, rhs, nodes_grouped_by_level, true, node_count, vector_size, 0);
#else
        LowerTriLvlSchedTP2SolverFunctor<RowMapType, EntriesType, ValuesType, LHSType, RHSType, NGBLType> tstf(row_map, entries, values, lhs, rhs, nodes_grouped_by_level, node_count, node_groups);
#endif
        Kokkos::parallel_for("parfor_u_team_vector", tvt_policy_type( (int)std::ceil((float)lvl_nodes/(float)node_groups) , team_size, vector_size ), tstf);
      } // end elseif
      */

      node_count += lvl_nodes;

#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOSPSTRSV_SOLVE_IMPL_PROFILE)
cudaProfilerStop();
#endif
    } // end if
   } // scope for if-block

  } // end for lvl

} // end lower_tri_solve



template < class TriSolveHandle, class RowMapType, class EntriesType, class ValuesType, class RHSType, class LHSType >
void upper_tri_solve(TriSolveHandle & thandle, const RowMapType row_map, const EntriesType entries, const ValuesType values, const RHSType & rhs, LHSType &lhs) {

#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOSPSTRSV_SOLVE_IMPL_PROFILE)
cudaProfilerStop();
#endif
  typedef typename TriSolveHandle::execution_space execution_space;
  typedef typename TriSolveHandle::size_type size_type;
  typedef typename TriSolveHandle::nnz_lno_view_t NGBLType;


  auto nlevels = thandle.get_num_levels();
  // Keep this a host View, create device version and copy to back to host during scheduling
  // This requires making sure the host view in the handle is properly updated after the symbolic phase
  auto nodes_per_level = thandle.get_nodes_per_level();
  auto hnodes_per_level = thandle.get_host_nodes_per_level();
  //auto hnodes_per_level = Kokkos::create_mirror_view(nodes_per_level);
  //Kokkos::deep_copy(hnodes_per_level, nodes_per_level);

  auto nodes_grouped_by_level = thandle.get_nodes_grouped_by_level();

  size_type node_count = 0;

  // This must stay serial; would be nice to try out Cuda's graph stuff to reduce kernel launch overhead
  for ( size_type lvl = 0; lvl < nlevels; ++lvl ) {
    size_type lvl_nodes = hnodes_per_level(lvl);

    if ( lvl_nodes != 0 ) {
#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOSPSTRSV_SOLVE_IMPL_PROFILE)
cudaProfilerStart();
#endif

      if ( thandle.get_algorithm() == KokkosSparse::Experimental::SPTRSVAlgorithm::SEQLVLSCHD_RP ) {
        Kokkos::parallel_for( "parfor_fixed_lvl", Kokkos::RangePolicy<execution_space>( node_count, node_count+lvl_nodes ), UpperTriLvlSchedRPSolverFunctor<RowMapType, EntriesType, ValuesType, LHSType, RHSType, NGBLType> (row_map, entries, values, lhs, rhs, nodes_grouped_by_level) );
      }
      else if ( thandle.get_algorithm() == KokkosSparse::Experimental::SPTRSVAlgorithm::SEQLVLSCHD_TP1 ) {
        typedef Kokkos::TeamPolicy<execution_space> policy_type;

        int team_size = thandle.get_team_size();

#ifdef KOKKOSKERNELS_SPTRSV_TRILVLSCHED
        TriLvlSchedTP1SolverFunctor<RowMapType, EntriesType, ValuesType, LHSType, RHSType, NGBLType> tstf(row_map, entries, values, lhs, rhs, nodes_grouped_by_level, false, node_count);
#else
        UpperTriLvlSchedTP1SolverFunctor<RowMapType, EntriesType, ValuesType, LHSType, RHSType, NGBLType> tstf(row_map, entries, values, lhs, rhs, nodes_grouped_by_level, node_count);
#endif
        if ( team_size == -1 )
          Kokkos::parallel_for("parfor_u_team", policy_type( lvl_nodes , Kokkos::AUTO ), tstf);
        else
          Kokkos::parallel_for("parfor_u_team", policy_type( lvl_nodes , team_size ), tstf);
      }
      // TP2 algorithm has issues with some offset-ordinal combo to be addressed
      /*
      else if ( thandle.get_algorithm() == KokkosSparse::Experimental::SPTRSVAlgorithm::SEQLVLSCHED_TP2 ) {
        typedef Kokkos::TeamPolicy<execution_space> tvt_policy_type;

        int team_size = thandle.get_team_size();
        if ( team_size == -1 ) {
          team_size = std::is_same< typename Kokkos::DefaultExecutionSpace::memory_space, Kokkos::HostSpace >::value ? 1 : 64;
        }
        int vector_size = thandle.get_team_size();
        if ( vector_size == -1 ) {
          vector_size = std::is_same< typename Kokkos::DefaultExecutionSpace::memory_space, Kokkos::HostSpace >::value ? 1 : 4;
        }

        // This impl: "chunk" lvl_nodes into node_groups; a league_rank is responsible for processing that many nodes
        //       TeamThreadRange over number nodes of node_groups
        //       To avoid masking threads, 1 thread (team) per node in node_group (thread has full ownership of a node)
        //       ThreadVectorRange responsible for the actual solve computation
        //const int node_groups = team_size;
        const int node_groups = vector_size;

#ifdef KOKKOSKERNELS_SPTRSV_TRILVLSCHED
        TriLvlSchedTP2SolverFunctor<RowMapType, EntriesType, ValuesType, LHSType, RHSType, NGBLType> tstf(row_map, entries, values, lhs, rhs, nodes_grouped_by_level, false, node_count, vector_size, 0);
#else
        UpperTriLvlSchedTP2SolverFunctor<RowMapType, EntriesType, ValuesType, LHSType, RHSType, NGBLType> tstf(row_map, entries, values, lhs, rhs, nodes_grouped_by_level, node_count, node_groups);
#endif

        Kokkos::parallel_for("parfor_u_team_vector", tvt_policy_type( (int)std::ceil((float)lvl_nodes/(float)node_groups) , team_size, vector_size ), tstf);
      } // end elseif
      */

      node_count += lvl_nodes;

#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOSPSTRSV_SOLVE_IMPL_PROFILE)
cudaProfilerStop();
#endif
    } // end if
  } // end for lvl

} // end upper_tri_solve


template < class TriSolveHandle, class RowMapType, class EntriesType, class ValuesType, class RHSType, class LHSType >
void tri_solve_chain(TriSolveHandle & thandle, const RowMapType row_map, const EntriesType entries, const ValuesType values, const RHSType & rhs, LHSType &lhs, const bool is_lowertri_) {

#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOSPSTRSV_SOLVE_IMPL_PROFILE)
cudaProfilerStop();
#endif
  typedef typename TriSolveHandle::execution_space execution_space;
  typedef typename TriSolveHandle::size_type size_type;
  typedef typename TriSolveHandle::nnz_lno_view_t NGBLType;

  // Algorithm is checked before this function is called
  auto h_chain_ptr = thandle.get_host_chain_ptr();
  size_type num_chain_entries = thandle.get_num_chain_entries();

  // Keep this a host View, create device version and copy to back to host during scheduling
  // This requires making sure the host view in the handle is properly updated after the symbolic phase
  auto nodes_per_level = thandle.get_nodes_per_level();
  auto hnodes_per_level = thandle.get_host_nodes_per_level();

  auto nodes_grouped_by_level = thandle.get_nodes_grouped_by_level();

  const bool is_lowertri =  thandle.is_lower_tri();

  size_type node_count = 0;

// REFACTORED to cleanup; next, need debug and timer routines
  using policy_type = Kokkos::TeamPolicy<execution_space>;
  using large_cutoff_policy_type = Kokkos::TeamPolicy<LargerCutoffTag, execution_space>;
/*
  using TP1Functor = TriLvlSchedTP1SolverFunctor<RowMapType, EntriesType, ValuesType, LHSType, RHSType, NGBLType>;
  using LTP1Functor = LowerTriLvlSchedTP1SolverFunctor<RowMapType, EntriesType, ValuesType, LHSType, RHSType, NGBLType>;
  using UTP1Functor = UpperTriLvlSchedTP1SolverFunctor<RowMapType, EntriesType, ValuesType, LHSType, RHSType, NGBLType>;
  using LSingleBlockFunctor = LowerTriLvlSchedTP1SingleBlockFunctor<RowMapType, EntriesType, ValuesType, LHSType, RHSType, NGBLType>;
  using USingleBlockFunctor = UpperTriLvlSchedTP1SingleBlockFunctor<RowMapType, EntriesType, ValuesType, LHSType, RHSType, NGBLType>;
*/
  using SingleBlockFunctor = TriLvlSchedTP1SingleBlockFunctor<RowMapType, EntriesType, ValuesType, LHSType, RHSType, NGBLType>;

  int team_size = thandle.get_team_size();
  int vector_size = thandle.get_vector_size() > 0 ? thandle.get_vector_size() : 1;

  auto cutoff = thandle.get_chain_threshold();
  int team_size_singleblock = team_size;

  // Enumerate options
  // ts -1,0 | cu 0 - select default ts == 1
  // ts -1,0 | cu > 0 - select default ts; restriction: ts <= tsmax (auto)
  // ts > 0 | cu 0 - set
  // ts > 0 | cu > 0 - set
  // Controls ts,cu > 0
  // co > ts  - not all rows can be mapped to a thread - must call largercutoff impl
  // co <= ts - okay, kernel must be careful not to access out-of-bounds; some threads idol
  if (team_size_singleblock <= 0 && cutoff == 0) {
    team_size_singleblock = 1;
    // If cutoff == 0, no single-block calls will be made, team_size_singleblock is unimportant
  }

  // This is only necessary for Lower,UpperTri functor versions; else, is_lowertri can be passed as arg to the generic Tri functor...
  if (is_lowertri) {

    for ( size_type chainlink = 0; chainlink < num_chain_entries; ++chainlink ) {
      size_type schain = h_chain_ptr(chainlink);
      size_type echain = h_chain_ptr(chainlink+1);

      if ( echain - schain == 1 ) {

        // if team_size is -1 (unset), get recommended size from Kokkos
#ifdef KOKKOSKERNELS_SPTRSV_TRILVLSCHED
        TriLvlSchedTP1SolverFunctor<RowMapType, EntriesType, ValuesType, LHSType, RHSType, NGBLType> tstf(row_map, entries, values, lhs, rhs, nodes_grouped_by_level, true, node_count);
#else
        LowerTriLvlSchedTP1SolverFunctor<RowMapType, EntriesType, ValuesType, LHSType, RHSType, NGBLType> tstf(row_map, entries, values, lhs, rhs, nodes_grouped_by_level, node_count);
#endif
        if (team_size == - 1) {
          team_size = policy_type(1, 1, vector_size).team_size_recommended(tstf, Kokkos::ParallelForTag());
        }

        size_type lvl_nodes = hnodes_per_level(schain); //lvl == echain????
        Kokkos::parallel_for("parfor_l_team_chain1", policy_type(lvl_nodes , team_size, vector_size), tstf);
        node_count += lvl_nodes;

      }
      else {
        size_type lvl_nodes = 0;

        for (size_type i = schain; i < echain; ++i) {
          lvl_nodes += hnodes_per_level(i);
        }

        if (team_size_singleblock <= 0) {
          team_size_singleblock = policy_type(1, 1, vector_size).team_size_recommended(SingleBlockFunctor(row_map, entries, values, lhs, rhs, nodes_grouped_by_level, nodes_per_level, node_count, schain, echain, is_lowertri), Kokkos::ParallelForTag());
        }

        if (cutoff <= team_size_singleblock) {
#ifdef KOKKOSKERNELS_SPTRSV_TRILVLSCHED
          TriLvlSchedTP1SingleBlockFunctor<RowMapType, EntriesType, ValuesType, LHSType, RHSType, NGBLType> tstf(row_map, entries, values, lhs, rhs, nodes_grouped_by_level, nodes_per_level, node_count, schain, echain, true);
#else
          LowerTriLvlSchedTP1SingleBlockFunctor<RowMapType, EntriesType, ValuesType, LHSType, RHSType, NGBLType> tstf(row_map, entries, values, lhs, rhs, nodes_grouped_by_level, nodes_per_level, node_count, schain, echain);
#endif
          Kokkos::parallel_for("parfor_l_team_chainmulti", policy_type(1, team_size_singleblock, vector_size), tstf);
        }
        else {
          // team_size_singleblock < cutoff => kernel must allow for a block-stride internally
#ifdef KOKKOSKERNELS_SPTRSV_TRILVLSCHED
          TriLvlSchedTP1SingleBlockFunctor<RowMapType, EntriesType, ValuesType, LHSType, RHSType, NGBLType> tstf(row_map, entries, values, lhs, rhs, nodes_grouped_by_level, nodes_per_level, node_count, schain, echain, true, 0, cutoff);
#else
          LowerTriLvlSchedTP1SingleBlockFunctor<RowMapType, EntriesType, ValuesType, LHSType, RHSType, NGBLType> tstf(row_map, entries, values, lhs, rhs, nodes_grouped_by_level, nodes_per_level, node_count, schain, echain, cutoff);
#endif
          Kokkos::parallel_for("parfor_l_team_chainmulti_cutoff", large_cutoff_policy_type(1, team_size_singleblock, vector_size), tstf);
        }
        node_count += lvl_nodes;
      }
      Kokkos::fence(); // TODO - is this necessary? that is, can the parallel_for launch before the s/echain values have been updated?
    }

  }
  else {

    for ( size_type chainlink = 0; chainlink < num_chain_entries; ++chainlink ) {
      size_type schain = h_chain_ptr(chainlink);
      size_type echain = h_chain_ptr(chainlink+1);

      if ( echain - schain == 1 ) {

        // if team_size is -1 (unset), get recommended size from Kokkos
#ifdef KOKKOSKERNELS_SPTRSV_TRILVLSCHED
        TriLvlSchedTP1SolverFunctor<RowMapType, EntriesType, ValuesType, LHSType, RHSType, NGBLType> tstf(row_map, entries, values, lhs, rhs, nodes_grouped_by_level, is_lowertri, node_count);
#else
        UpperTriLvlSchedTP1SolverFunctor<RowMapType, EntriesType, ValuesType, LHSType, RHSType, NGBLType> tstf(row_map, entries, values, lhs, rhs, nodes_grouped_by_level, node_count);
#endif
        if (team_size == - 1) {
          team_size = policy_type(1, 1, vector_size).team_size_recommended(tstf, Kokkos::ParallelForTag());
        }

        // TODO To use cudagraph here, need to know how many non-unit chains there are, create a graph for each and launch accordingly
        size_type lvl_nodes = hnodes_per_level(schain); //lvl == echain????
        Kokkos::parallel_for("parfor_u_team_chain1", policy_type(lvl_nodes , team_size, vector_size), tstf);
        node_count += lvl_nodes;

      }
      else {
        size_type lvl_nodes = 0;

        for (size_type i = schain; i < echain; ++i) {
          lvl_nodes += hnodes_per_level(i);
        }

        if (team_size_singleblock <= 0) {
          //team_size_singleblock = policy_type(1, 1, 1).team_size_recommended(SingleBlockFunctor(row_map, entries, values, lhs, rhs, nodes_grouped_by_level, is_lowertri, node_count), Kokkos::ParallelForTag());
          team_size_singleblock = policy_type(1, 1, vector_size).team_size_recommended(SingleBlockFunctor(row_map, entries, values, lhs, rhs, nodes_grouped_by_level, nodes_per_level, node_count, schain, echain, is_lowertri), Kokkos::ParallelForTag());
        }

        if (cutoff <= team_size_singleblock) {
#ifdef KOKKOSKERNELS_SPTRSV_TRILVLSCHED
          TriLvlSchedTP1SingleBlockFunctor<RowMapType, EntriesType, ValuesType, LHSType, RHSType, NGBLType> tstf(row_map, entries, values, lhs, rhs, nodes_grouped_by_level, nodes_per_level, node_count, schain, echain, is_lowertri);
#else
          UpperTriLvlSchedTP1SingleBlockFunctor<RowMapType, EntriesType, ValuesType, LHSType, RHSType, NGBLType> tstf(row_map, entries, values, lhs, rhs, nodes_grouped_by_level, nodes_per_level, node_count, schain, echain);
#endif
          Kokkos::parallel_for("parfor_u_team_chainmulti", policy_type(1, team_size_singleblock, vector_size), tstf);
        }
        else {
          // team_size_singleblock < cutoff => kernel must allow for a block-stride internally
#ifdef KOKKOSKERNELS_SPTRSV_TRILVLSCHED
          TriLvlSchedTP1SingleBlockFunctor<RowMapType, EntriesType, ValuesType, LHSType, RHSType, NGBLType> tstf(row_map, entries, values, lhs, rhs, nodes_grouped_by_level, nodes_per_level, node_count, schain, echain, is_lowertri, 0, cutoff);
#else
          UpperTriLvlSchedTP1SingleBlockFunctor<RowMapType, EntriesType, ValuesType, LHSType, RHSType, NGBLType> tstf(row_map, entries, values, lhs, rhs, nodes_grouped_by_level, nodes_per_level, node_count, schain, echain, cutoff);
#endif
          Kokkos::parallel_for("parfor_u_team_chainmulti_cutoff", large_cutoff_policy_type(1, team_size_singleblock, vector_size), tstf);
        }
        node_count += lvl_nodes;
      }
      Kokkos::fence(); // TODO - is this necessary? that is, can the parallel_for launch before the s/echain values have been updated?
    }

  }

} // end tri_solve_chain

} // namespace Experimental
} // namespace Impl
} // namespace KokkosSparse

#endif
