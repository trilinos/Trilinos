// clang-format off
// @HEADER
// *****************************************************************************
//                            Tacho package
//
// Copyright 2022 NTESS and the Tacho contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
// clang-format on
#ifndef __TACHO_TEAMFUNCTOR_EXTRACT_CRS_HPP__
#define __TACHO_TEAMFUNCTOR_EXTRACT_CRS_HPP__

/// \file Tacho_TeamFunctor_FactorizeChol.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

#include "Tacho_SupernodeInfo.hpp"

namespace Tacho {


struct rowptr_sum {
  int* _rowptr;

  rowptr_sum(int* rowptr)
      : _rowptr(rowptr) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const ordinal_type i, ordinal_type& update,
                  const bool is_final) const {
    const ordinal_type val_i = _rowptr[i];
    update += val_i;
    if (is_final) {
      _rowptr[i] = update;
    }
  }
};

template <typename SupernodeInfoType> struct TeamFunctor_ExtractCrs {
public:
  typedef Kokkos::pair<ordinal_type, ordinal_type> range_type;

  typedef SupernodeInfoType supernode_info_type;
  typedef typename supernode_info_type::supernode_type supernode_type;

  typedef typename supernode_info_type::value_type value_type;
  typedef typename supernode_info_type::value_type_matrix value_type_matrix;
  typedef typename supernode_info_type::ordinal_type_array ordinal_type_array;

private:
  supernode_info_type _info;
  ordinal_type_array _compute_mode, _level_sids;
  ordinal_type _pbeg, _pend;
  ordinal_type _m;

  // in CRS format
  int* _rowptr;
  int* _colind;
  value_type* _nzvals;
  // in CRS format, transpose
  int* _rowptrT;
  int* _colindT;
  value_type* _nzvalsT;
  // pivot
  ordinal_type_array _piv;

public:
  KOKKOS_INLINE_FUNCTION
  TeamFunctor_ExtractCrs() = delete;

  KOKKOS_INLINE_FUNCTION
  TeamFunctor_ExtractCrs(const supernode_info_type &info, const ordinal_type_array &compute_mode,
                         const ordinal_type_array &level_sids)
      : _info(info), _compute_mode(compute_mode), _level_sids(level_sids) {}

  inline void setGlobalSize(const ordinal_type m) {
    _m = m;
  }

  inline void setRange(const ordinal_type pbeg, const ordinal_type pend) {
    _pbeg = pbeg;
    _pend = pend;
  }

  inline void setRowPtr(int* rowptr) { _rowptr = rowptr; }
  inline void setCrsView(int *colind, value_type *nzvals) {
    _colind = colind;
    _nzvals = nzvals;
  }
  inline void setRowPtrT(int* rowptrT) { _rowptrT = rowptrT; }
  inline void setCrsViewT(int *colindT, value_type *nzvalsT) {
    _colindT = colindT;
    _nzvalsT = nzvalsT;
  }
  inline void setPivPtr(ordinal_type_array &piv) { _piv = piv; }

  struct ExtractPtrTag {};
  struct ExtractValTag {};
  struct ExtractPtrColTag {};
  struct ExtractValColTag {};
  struct TransPtrTag {};
  struct TransMatTag {};


  // ---------------------------------------
  // Functors to convert to CRS format
  //  from row-major
  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const ExtractPtrTag &, const MemberType &member) const {
    const ordinal_type id = member.league_rank();
    const ordinal_type p = _pbeg + id;

    const value_type zero(0);
    const ordinal_type sid  = (p == _pend ? 0 : _level_sids(p));
    const ordinal_type mode = (p == _pend ? 0 : _compute_mode(sid));
    if (mode == 0) {
      // extract this supernode
      const auto &s  = _info.supernodes(sid);
      const ordinal_type offm = (p == _pend ? _m : s.row_begin);
      #define TACHO_INSERT_DIAGONALS
      #ifdef TACHO_INSERT_DIAGONALS
      // last row of previous supernode
      ordinal_type row_id = 0;
      if (p > _pbeg) {
        const ordinal_type prev_sid = _level_sids(p-1);
        const auto &prev_s = _info.supernodes(prev_sid);
        row_id = prev_s.row_begin + prev_s.m;
      }
      // add diagonal entry
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member, offm-row_id),
                          [&](const int& i) { _rowptr[row_id+i+1] = 1; });
      #endif
      if (p < _pend) {
        if (s.m > 0) {
          // extract this supernode
          //  stored by row, but checking for nonzereo (instead of just taking all s.n nonzereos)
          value_type *aptr = s.u_buf;
          UnmanagedViewType<value_type_matrix> AT(aptr, s.m, s.n);
          Kokkos::parallel_for(Kokkos::TeamThreadRange(member, s.m),
                               [&](const int& i) { 
                                 _rowptr[1+i+offm] = 0;
                                 for (ordinal_type j = 0; j < s.n; j++) {
                                   if (AT(i,j) != zero) {
                                     _rowptr[1+i+offm] ++;
                                   }
                                 }
                               });
        }
      }
    }
  }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const ExtractValTag &, const MemberType &member) const {
    const ordinal_type id = member.league_rank();
    const ordinal_type p = _pbeg + id;

    const value_type one (1);
    const value_type zero(0);
    const ordinal_type sid  = (p == _pend ? 0 : _level_sids(p));
    const ordinal_type mode = (p == _pend ? 0 : _compute_mode(sid));
    if (mode == 0) {
      // extract this supernode
      const auto &s  = _info.supernodes(sid);
      const ordinal_type offm = (p == _pend ? _m : s.row_begin);
      const ordinal_type offn = (p == _pend ?  0 : s.gid_col_begin);
      #ifdef TACHO_INSERT_DIAGONALS
      // last row of previous supernode
      ordinal_type row_id = 0;
      if (p > _pbeg) {
        const ordinal_type prev_sid = _level_sids(p-1);
        const auto &prev_s = _info.supernodes(prev_sid);
        row_id = prev_s.row_begin + prev_s.m;
      }
      // insert diagonals for the missing rows between previous and this block
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member, offm-row_id),
                          [&](const int& i) {
                            int nnz = _rowptr[row_id+i];
                            _colind[nnz] = row_id+i;
                            _nzvals[nnz] = one;
                            _rowptr[row_id+i]++;
                          });
      #endif
      if (p < _pend) {
        if (s.m > 0) {
          // extract this supernode
          value_type *aptr = s.u_buf;
          UnmanagedViewType<value_type_matrix> AT(aptr, s.m, s.n);
          Kokkos::parallel_for(Kokkos::TeamThreadRange(member, s.m),
                               [&](const int& i) {
                                 // diagonal block
                                 ordinal_type j;
                                 for (ordinal_type j = i; j < s.m; j++) {
                                   if (AT(i,j) != zero) {
                                     int nnz = _rowptr[i+offm];
                                     _colind[nnz] = j+offm;
                                     _nzvals[nnz] = AT(i,j);
                                     _rowptr[i+offm] ++;
                                   }
                                 }
                                 // off-diagonal blocksa
                                 j = s.m;
                                 for (ordinal_type id = s.sid_col_begin + 1; id < s.sid_col_end - 1; id++) {
                                   for (ordinal_type k = _info.sid_block_colidx(id).second; k < _info.sid_block_colidx(id + 1).second; k++) {
                                     if (AT(i,j) != zero) {
                                       int nnz = _rowptr[i+offm];
                                       _colind[nnz] = _info.gid_colidx(k+offn);
                                       _nzvals[nnz] = AT(i,j);
                                       _rowptr[i+offm] ++;
                                     }
                                     j++;
                                   }
                                 }
                               });
        }
      }
    }
  }

  // ---------------------------------------
  // Functors to convert to CRS format
  //  from col-major
  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const ExtractPtrColTag &, const MemberType &member) const {
    const ordinal_type id = member.league_rank();
    const ordinal_type p = _pbeg + id;

    const value_type zero(0);
    const ordinal_type sid  = (p == _pend ? 0 : _level_sids(p));
    const ordinal_type mode = (p == _pend ? 0 : _compute_mode(sid));
    if (mode == 0) {
      // extract this supernode
      const auto &s  = _info.supernodes(sid);
      const ordinal_type offm = (p == _pend ? _m : s.row_begin);
      const ordinal_type offn = (p == _pend ?  0 : s.gid_col_begin);
      #ifdef TACHO_INSERT_DIAGONALS
      // last row of previous supernode
      ordinal_type row_id = 0;
      if (p > _pbeg) {
        const ordinal_type prev_sid = _level_sids(p-1);
        const auto &prev_s = _info.supernodes(prev_sid);
        row_id = prev_s.row_begin + prev_s.m;
      }
      // add diagonal entry
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member, offm-row_id),
                          [&](const int& i) { Kokkos::atomic_add(&(_rowptr[row_id+i+1]), 1); });
      #endif
      if (p < _pend) {
        // extract this supernode (AL is stored by col)
        if (s.m > 0) {
          value_type *aptr = s.l_buf;
          UnmanagedViewType<value_type_matrix> AL(aptr, s.n, s.m);
          Kokkos::parallel_for(Kokkos::TeamThreadRange(member, s.m),
                               [&](const int& j) {
                                 // first extract diagnal block (each thread extract row in parallel)
                                 for (ordinal_type i = 0; i < s.m; i++) {
                                   if (AL(i,j) != zero) {
                                     Kokkos::atomic_add(&(_rowptr[1+i+offm]), 1);
                                   }
                                 }
                                 // off-diagonals (each thread extract col, needing atomic-add)
                                 ordinal_type i = s.m;
                                 for (ordinal_type id = s.sid_col_begin + 1; id < s.sid_col_end - 1; id++) {
                                   for (ordinal_type k = _info.sid_block_colidx(id).second; k < _info.sid_block_colidx(id + 1).second; k++) {
                                     if (AL(i, j) != zero) {
                                       ordinal_type gid_i = _info.gid_colidx(k+offn);
                                       Kokkos::atomic_add(&(_rowptr[1+gid_i]), 1);
                                     }
                                     i++;
                                   }
                                 }
                               });
        }
      }
    }
  }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const ExtractValColTag &, const MemberType &member) const {
    const ordinal_type id = member.league_rank();
    const ordinal_type p = _pbeg + id;

    const value_type zero(0);
    const value_type one (1);
    const ordinal_type sid  = (p == _pend ? 0 : _level_sids(p));
    const ordinal_type mode = (p == _pend ? 0 : _compute_mode(sid));
    if (mode == 0) {
      // extract this supernode
      const auto &s  = _info.supernodes(sid);
      const ordinal_type offm = (p == _pend ? _m : s.row_begin);
      const ordinal_type offn = (p == _pend ?  0 : s.gid_col_begin);
      #ifdef TACHO_INSERT_DIAGONALS
      // last row of previous supernode
      ordinal_type row_id = 0;
      if (p > _pbeg) {
        const ordinal_type prev_sid = _level_sids(p-1);
        const auto &prev_s = _info.supernodes(prev_sid);
        row_id = prev_s.row_begin + prev_s.m;
      }
      // add diagonal entry
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member, offm-row_id),
                          [&](const int& i) {
                            ordinal_type nnz = Kokkos::atomic_fetch_add(&(_rowptr[row_id+i]), 1);
                            _colind[nnz] = row_id+i;
                            _nzvals[nnz] = one;
                          });
      #endif
      if (p < _pend) {
        // extract this supernode
        //  stored by col
        if (s.m > 0) {
          bool no_perm = s.do_not_apply_pivots;
          ConstUnmanagedViewType<ordinal_type_array> perm(_piv.data() + 4 * offm + 2 * s.m, s.m);

          value_type *aptr = s.l_buf;
          UnmanagedViewType<value_type_matrix> AL(aptr, s.n, s.m);
          Kokkos::parallel_for(Kokkos::TeamThreadRange(member, s.m),
                               [&](const int& j) {
                                 ordinal_type gid_j = (no_perm ? j+offm : perm(j)+offm);
                                 // diagnal block 
                                 for (ordinal_type i = 0; i < s.m; i++) {
                                   if (AL(i,j) != zero) {
                                     ordinal_type nnz = Kokkos::atomic_fetch_add(&(_rowptr[offm+i]), 1);
                                     _colind[nnz] = gid_j;
                                     _nzvals[nnz] = AL(i,j);
                                   }
                                 }
                                 // off-diagonals (each thread extract col, needing atomic-add)
                                 ordinal_type i = s.m;
                                 for (ordinal_type id = s.sid_col_begin + 1; id < s.sid_col_end - 1; id++) {
                                   for (ordinal_type k = _info.sid_block_colidx(id).second; k < _info.sid_block_colidx(id + 1).second; k++) {
                                     if (AL(i, j) != zero) {
                                       ordinal_type gid_i = _info.gid_colidx(k+offn);
                                       ordinal_type nnz = Kokkos::atomic_fetch_add(&(_rowptr[gid_i]), 1);
                                       _colind[nnz] = gid_j;
                                       _nzvals[nnz] = AL(i,j);
                                     }
                                     i++;
                                   }
                                 }
                               });
        }
      }
    }
  }


  // ---------------------------------------
  // Functors to transpose
  KOKKOS_INLINE_FUNCTION void operator()(const TransPtrTag &, const int i) const {
    // count offset rowptrT
    for (ordinal_type k = _rowptr[i]; k < _rowptr[i+1]; k++) {
      Kokkos::atomic_add(&(_rowptrT[_colind[k]+1]), 1);
    }
  }

  KOKKOS_INLINE_FUNCTION void operator()(const TransMatTag &, const int i) const {
    // count offset rowptrT
    for (ordinal_type k = _rowptr[i]; k < _rowptr[i+1]; k++) {
      int nnz = Kokkos::atomic_fetch_add(&(_rowptrT[_colind[k]]), 1);
      _colindT[nnz] = i;
      _nzvalsT[nnz] = _nzvals[k];
    }
  }
};
} // namespace Tacho

#endif
