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

#ifndef KOKKOSSPARSE_IMPL_SPILUK_NUMERIC_HPP_
#define KOKKOSSPARSE_IMPL_SPILUK_NUMERIC_HPP_

/// \file KokkosSparse_spiluk_numeric_impl.hpp
/// \brief Implementation(s) of the numeric phase of sparse ILU(k).

#include <KokkosKernels_config.h>
#include <Kokkos_ArithTraits.hpp>
#include <KokkosSparse_spiluk_handle.hpp>

//#define NUMERIC_OUTPUT_INFO

namespace KokkosSparse {
namespace Impl {
namespace Experimental {


//struct UnsortedTag {};

template <class ARowMapType,
          class AEntriesType,
          class AValuesType,
          class LRowMapType,
          class LEntriesType,
          class LValuesType,
          class URowMapType,
          class UEntriesType,
          class UValuesType,
          class LevelViewType,
          class WorkViewType,
          class nnz_lno_t>
struct ILUKLvlSchedRPNumericFunctor
{
  using lno_t    = typename AEntriesType::non_const_value_type;
  using scalar_t = typename AValuesType::non_const_value_type;
  ARowMapType   A_row_map;
  AEntriesType  A_entries;
  AValuesType   A_values;
  LRowMapType   L_row_map;
  LEntriesType  L_entries;
  LValuesType   L_values;
  URowMapType   U_row_map;
  UEntriesType  U_entries;
  UValuesType   U_values;
  LevelViewType level_idx;
  WorkViewType  iw;
  nnz_lno_t     lev_start;

  ILUKLvlSchedRPNumericFunctor( const ARowMapType &A_row_map_, const AEntriesType &A_entries_, const AValuesType &A_values_, const LRowMapType &L_row_map_, const LEntriesType &L_entries_, LValuesType &L_values_, const URowMapType &U_row_map_, const UEntriesType &U_entries_, UValuesType &U_values_, const LevelViewType &level_idx_, WorkViewType &iw_, const nnz_lno_t &lev_start_ ) :
    A_row_map(A_row_map_), A_entries(A_entries_), A_values(A_values_), L_row_map(L_row_map_), L_entries(L_entries_), L_values(L_values_), U_row_map(U_row_map_), U_entries(U_entries_), U_values(U_values_), level_idx(level_idx_), iw(iw_), lev_start(lev_start_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const lno_t i) const {
    auto rowid  = level_idx(i);
    auto tid    = i-lev_start;
    auto k1 = L_row_map(rowid); 
    auto k2 = L_row_map(rowid+1);
#ifdef KEEP_DIAG
    for (auto k = k1; k < k2-1; ++k) {
#else
    for (auto k = k1; k < k2; ++k) {
#endif
      auto col    = L_entries(k);
      L_values(k) = 0.0;
      iw(tid,col) = k;
    }
#ifdef KEEP_DIAG
    L_values(k2-1) = scalar_t(1.0);
#endif

    k1 = U_row_map(rowid); 
    k2 = U_row_map(rowid+1);
    for (auto k = k1; k < k2; ++k) {
      auto col    = U_entries(k);
      U_values(k) = 0.0;
      iw(tid,col) = k;
    }

    //Unpack the ith row of A
    k1 = A_row_map(rowid);
    k2 = A_row_map(rowid+1);
    for (auto k = k1; k < k2; ++k) {
      auto col  = A_entries(k);
      auto ipos = iw(tid,col);
      if (col < rowid)
        L_values(ipos) = A_values(k);
      else
        U_values(ipos) = A_values(k);
    }

    //Eliminate prev rows
    k1 = L_row_map(rowid); 
    k2 = L_row_map(rowid+1);
#ifdef KEEP_DIAG
    for (auto k = k1; k < k2-1; ++k) {
#else
    for (auto k = k1; k < k2; ++k) {
#endif
      auto prev_row = L_entries(k);
#ifdef KEEP_DIAG
      auto fact = L_values(k) / U_values(U_row_map(prev_row));
#else
      auto fact = L_values(k) * U_values(U_row_map(prev_row));
#endif
      L_values(k) = fact;
      for (auto kk = U_row_map(prev_row)+1; kk < U_row_map(prev_row+1); ++kk) {
        auto col  = U_entries(kk);
        auto ipos = iw(tid,col);
        if (ipos == -1) continue;
        auto lxu = -U_values(kk) * fact;
        if (col < rowid)
          L_values(ipos) += lxu;
        else
          U_values(ipos) += lxu;
      }// end for kk
    }// end for k

#ifdef KEEP_DIAG
    if (U_values(iw(tid,rowid)) == 0.0) {
      U_values(iw(tid,rowid)) = 1e6;
    }
#else
    if (U_values(iw(tid,rowid)) == 0.0) {
      U_values(iw(tid,rowid)) = 1e6;
    }
    else {
      U_values(iw(tid,rowid)) = 1.0 / U_values(iw(tid,rowid));
    }
#endif

    //Reset
    k1 = L_row_map(rowid); 
    k2 = L_row_map(rowid+1);
#ifdef KEEP_DIAG
    for (auto k = k1; k < k2-1; ++k)
#else
    for (auto k = k1; k < k2; ++k)
#endif
      iw(tid,L_entries(k)) = -1;

    k1 = U_row_map(rowid); 
    k2 = U_row_map(rowid+1);
    for (auto k = k1; k < k2; ++k)
      iw(tid,U_entries(k)) = -1;
  }
};

template <class ARowMapType,
          class AEntriesType,
          class AValuesType,
          class LRowMapType,
          class LEntriesType,
          class LValuesType,
          class URowMapType,
          class UEntriesType,
          class UValuesType,
          class LevelViewType,
          class WorkViewType,
          class nnz_lno_t>
struct ILUKLvlSchedTP1NumericFunctor
{
  using execution_space = typename ARowMapType::execution_space;
  using policy_type     = Kokkos::TeamPolicy<execution_space>;
  using member_type     = typename policy_type::member_type;
  using size_type       = typename ARowMapType::non_const_value_type;
  using lno_t           = typename AEntriesType::non_const_value_type;
  using scalar_t        = typename AValuesType::non_const_value_type ;

  ARowMapType   A_row_map;
  AEntriesType  A_entries;
  AValuesType   A_values;
  LRowMapType   L_row_map;
  LEntriesType  L_entries;
  LValuesType   L_values;
  URowMapType   U_row_map;
  UEntriesType  U_entries;
  UValuesType   U_values;
  LevelViewType level_idx;
  WorkViewType  iw;
  nnz_lno_t     lev_start;

  ILUKLvlSchedTP1NumericFunctor( const ARowMapType &A_row_map_, const AEntriesType &A_entries_, const AValuesType &A_values_, const LRowMapType &L_row_map_, const LEntriesType &L_entries_, LValuesType &L_values_, const URowMapType &U_row_map_, const UEntriesType &U_entries_, UValuesType &U_values_, const LevelViewType &level_idx_, WorkViewType &iw_, const nnz_lno_t &lev_start_ ) :
    A_row_map(A_row_map_), A_entries(A_entries_), A_values(A_values_), L_row_map(L_row_map_), L_entries(L_entries_), L_values(L_values_), U_row_map(U_row_map_), U_entries(U_entries_), U_values(U_values_), level_idx(level_idx_), iw(iw_), lev_start(lev_start_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const member_type & team ) const {
    auto my_league = team.league_rank(); // map to rowid
    auto rowid     = level_idx(my_league + lev_start);
    auto my_team   = team.team_rank();

    auto k1 = L_row_map(rowid); 
    auto k2 = L_row_map(rowid+1);
#ifdef KEEP_DIAG
    Kokkos::parallel_for( Kokkos::TeamThreadRange( team, k1, k2-1 ), [&] ( const size_type k ) { 
      auto col          = L_entries(k);
      L_values(k)       = 0.0;
      iw(my_league,col) = k;
    });
#else
    Kokkos::parallel_for( Kokkos::TeamThreadRange( team, k1, k2 ), [&] ( const size_type k ) { 
      auto col          = L_entries(k);
      L_values(k)       = 0.0;
      iw(my_league,col) = k;
    });
#endif

#ifdef KEEP_DIAG
    if ( my_team == 0 ) L_values(k2-1) = scalar_t(1.0);
#endif

    team.team_barrier();

    k1 = U_row_map(rowid); 
    k2 = U_row_map(rowid+1);
    Kokkos::parallel_for( Kokkos::TeamThreadRange( team, k1, k2 ), [&] ( const size_type k ) { 
      auto col          = U_entries(k);
      U_values(k)       = 0.0;
      iw(my_league,col) = k;
    });

    team.team_barrier();
	
    //Unpack the ith row of A
    k1 = A_row_map(rowid);
    k2 = A_row_map(rowid+1);
    Kokkos::parallel_for( Kokkos::TeamThreadRange( team, k1, k2 ), [&] ( const size_type k ) {
      auto col  = A_entries(k);
      auto ipos = iw(my_league,col);
      if (col < rowid)
        L_values(ipos) = A_values(k);
      else
        U_values(ipos) = A_values(k);	  
    });

    team.team_barrier();
	
    //Eliminate prev rows
    k1 = L_row_map(rowid); 
    k2 = L_row_map(rowid+1);
#ifdef KEEP_DIAG
    for (auto k = k1; k < k2-1; ++k) {
#else
    for (auto k = k1; k < k2; ++k) {
#endif
      auto prev_row = L_entries(k);
#ifdef KEEP_DIAG
      auto fact = L_values(k) / U_values(U_row_map(prev_row));
#else
      auto fact = L_values(k) * U_values(U_row_map(prev_row));
#endif
      if ( my_team == 0 ) L_values(k) = fact; 

      team.team_barrier();

      Kokkos::parallel_for( Kokkos::TeamThreadRange( team, U_row_map(prev_row)+1, U_row_map(prev_row+1) ), [&] ( const size_type kk ) {
        auto col  = U_entries(kk);
        auto ipos = iw(my_league,col);
        if (ipos != -1) {
          auto lxu = -U_values(kk) * fact;
          if (col < rowid)
            L_values(ipos) += lxu;
          else
            U_values(ipos) += lxu;
        }
      });// end for kk

      team.team_barrier();
    }// end for k

    if ( my_team == 0 ) {
#ifdef KEEP_DIAG
      if (U_values(iw(my_league,rowid)) == 0.0) {
        U_values(iw(my_league,rowid)) = 1e6;
      }
#else
      if (U_values(iw(my_league,rowid)) == 0.0) {
        U_values(iw(my_league,rowid)) = 1e6;
      }
      else {
        U_values(iw(my_league,rowid)) = 1.0 / U_values(iw(my_league,rowid));
      }
#endif
    }

    team.team_barrier();

    //Reset
    k1 = L_row_map(rowid); 
    k2 = L_row_map(rowid+1);
#ifdef KEEP_DIAG
    Kokkos::parallel_for( Kokkos::TeamThreadRange( team, k1, k2-1 ), [&] ( const size_type k ) {
      iw(my_league,L_entries(k)) = -1;
    });
#else
    Kokkos::parallel_for( Kokkos::TeamThreadRange( team, k1, k2 ), [&] ( const size_type k ) {
      iw(my_league,L_entries(k)) = -1;
    });
#endif

    k1 = U_row_map(rowid); 
    k2 = U_row_map(rowid+1);
    Kokkos::parallel_for( Kokkos::TeamThreadRange( team, k1, k2 ), [&] ( const size_type k ) {
      iw(my_league,U_entries(k)) = -1;
    });
  }
};

template <class IlukHandle,
          class ARowMapType,
          class AEntriesType,
          class AValuesType,
          class LRowMapType,
          class LEntriesType,
          class LValuesType,
          class URowMapType,
          class UEntriesType,
          class UValuesType>
void iluk_numeric ( IlukHandle& thandle,
                    const ARowMapType&  A_row_map,
                    const AEntriesType& A_entries,
                    const AValuesType&  A_values,
                    const LRowMapType&  L_row_map,
                    const LEntriesType& L_entries,
                          LValuesType&  L_values,
                    const URowMapType&  U_row_map,
                    const UEntriesType& U_entries,
                          UValuesType&  U_values ) {

  using execution_space = typename IlukHandle::execution_space;
  using memory_space    = typename IlukHandle::memory_space;
  using size_type       = typename IlukHandle::size_type;
  using nnz_lno_t       = typename IlukHandle::nnz_lno_t;
  using HandleDeviceEntriesType = typename IlukHandle::nnz_lno_view_t;
  using HandleHostEntriesType   = typename IlukHandle::nnz_lno_view_t::HostMirror;

  size_type nlevels = thandle.get_num_levels();
  size_type nrows   = thandle.get_nrows();

  // Keep this as host View, create device version and copy to back to host
  HandleDeviceEntriesType level_ptr = thandle.get_level_ptr();
  HandleHostEntriesType level_ptr_h = Kokkos::create_mirror_view(level_ptr);
  Kokkos::deep_copy(level_ptr_h, level_ptr);

  HandleDeviceEntriesType level_idx = thandle.get_level_idx();

  using WorkViewType = Kokkos::View<nnz_lno_t**, Kokkos::Device<execution_space,memory_space>>;
  
  WorkViewType iw ( "iw", thandle.get_level_maxrows(), nrows );
  Kokkos::deep_copy(iw, nnz_lno_t(-1));

  // Main loop must be performed sequential. Question: Try out Cuda's graph stuff to reduce kernel launch overhead
  for ( size_type lvl = 0; lvl < nlevels; ++lvl ) {
    nnz_lno_t lev_start = level_ptr_h(lvl);
    nnz_lno_t lev_end   = level_ptr_h(lvl+1);

    if ( (lev_end - lev_start) != 0 ) {

      if ( thandle.get_algorithm() == KokkosSparse::Experimental::SPILUKAlgorithm::SEQLVLSCHD_RP ) {
        Kokkos::parallel_for( "parfor_fixed_lvl", Kokkos::RangePolicy<execution_space>( lev_start, lev_end ), 
                                                  ILUKLvlSchedRPNumericFunctor<ARowMapType,
                                                                               AEntriesType,
                                                                               AValuesType,
                                                                               LRowMapType,
                                                                               LEntriesType,
                                                                               LValuesType,
                                                                               URowMapType,
                                                                               UEntriesType,
                                                                               UValuesType,
                                                                               HandleDeviceEntriesType,
                                                                               WorkViewType,
                                                                               nnz_lno_t> (A_row_map, A_entries, A_values, L_row_map, L_entries, L_values, U_row_map, U_entries, U_values, level_idx, iw, lev_start) );
      }
      else if ( thandle.get_algorithm() == KokkosSparse::Experimental::SPILUKAlgorithm::SEQLVLSCHD_TP1 ) {
        using policy_type = Kokkos::TeamPolicy<execution_space>;
        int team_size = thandle.get_team_size();

        ILUKLvlSchedTP1NumericFunctor<ARowMapType,
                                      AEntriesType,
                                      AValuesType,
                                      LRowMapType,
                                      LEntriesType,
                                      LValuesType,
                                      URowMapType,
                                      UEntriesType,
                                      UValuesType,
                                      HandleDeviceEntriesType,
                                      WorkViewType,
                                      nnz_lno_t> tstf(A_row_map, A_entries, A_values, L_row_map, L_entries, L_values, U_row_map, U_entries, U_values, level_idx, iw, lev_start);
        if ( team_size == -1 )
          Kokkos::parallel_for("parfor_l_team", policy_type( lev_end - lev_start , Kokkos::AUTO ), tstf);
        else
          Kokkos::parallel_for("parfor_l_team", policy_type( lev_end - lev_start , team_size ), tstf);
      }
//      /*
//      // TP2 algorithm has issues with some offset-ordinal combo to be addressed
//      else if ( thandle.get_algorithm() == KokkosSparse::Experimental::SPTRSVAlgorithm::SEQLVLSCHED_TP2 ) {
//        typedef Kokkos::TeamPolicy<execution_space> tvt_policy_type;
//
//        int team_size = thandle.get_team_size();
//        if ( team_size == -1 ) {
//          team_size = std::is_same< typename Kokkos::DefaultExecutionSpace::memory_space, Kokkos::HostSpace >::value ? 1 : 128;
//        }
//        int vector_size = thandle.get_team_size();
//        if ( vector_size == -1 ) {
//          vector_size = std::is_same< typename Kokkos::DefaultExecutionSpace::memory_space, Kokkos::HostSpace >::value ? 1 : 4;
//        }
//
//        // This impl: "chunk" lvl_nodes into node_groups; a league_rank is responsible for processing that many nodes
//        //       TeamThreadRange over number of node_groups
//        //       To avoid masking threads, 1 thread (team) per node in node_group
//        //       ThreadVectorRange responsible for the actual solve computation
//        const int node_groups = team_size;
//
//        LowerTriLvlSchedTP2SolverFunctor<RowMapType, EntriesType, ValuesType, LHSType, RHSType, HandleDeviceEntriesType> tstf(row_map, entries, values, lhs, rhs, nodes_grouped_by_level, row_count, node_groups);
//        Kokkos::parallel_for("parfor_u_team_vector", tvt_policy_type( (int)std::ceil((float)lvl_nodes/(float)node_groups) , team_size, vector_size ), tstf);
//      } // end elseif
//      */

    } // end if
  } // end for lvl

// Output check
#ifdef NUMERIC_OUTPUT_INFO
  std::cout << "  iluk_numeric result: " << std::endl;

  std::cout << "  nnzL: " << thandle.get_nnzL() << std::endl;
  std::cout << "  L_row_map = ";
  for ( size_type i = 0; i < nrows+1; ++i )
  { std::cout << L_row_map(i) << " "; }
  std::cout << std::endl;

  std::cout << "  L_entries = ";
  for ( size_type i = 0; i < thandle.get_nnzL(); ++i )
  { std::cout << L_entries(i) << " "; }
  std::cout << std::endl;

  std::cout << "  L_values = ";
  for ( size_type i = 0; i < thandle.get_nnzL(); ++i )
  { std::cout << L_values(i) << " "; }
  std::cout << std::endl;

  std::cout << "  nnzU: " << thandle.get_nnzU() << std::endl;
  std::cout << "  U_row_map = ";
  for ( size_type i = 0; i < nrows+1; ++i )
  { std::cout << U_row_map(i) << " "; }
  std::cout << std::endl;

  std::cout << "  U_entries = ";
  for ( size_type i = 0; i < thandle.get_nnzU(); ++i )
  { std::cout << U_entries(i) << " "; }
  std::cout << std::endl;

  std::cout << "  U_values = ";
  for ( size_type i = 0; i < thandle.get_nnzU(); ++i )
  { std::cout << U_values(i) << " "; }
  std::cout << std::endl;
#endif

} // end iluk_numeric


} // namespace Experimental
} // namespace Impl
} // namespace KokkosSparse

#endif
