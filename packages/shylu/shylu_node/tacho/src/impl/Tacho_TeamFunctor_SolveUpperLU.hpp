#ifndef __TACHO_TEAMFUNCTOR_SOLVE_UPPER_LU_HPP__
#define __TACHO_TEAMFUNCTOR_SOLVE_UPPER_LU_HPP__

/// \file Tacho_TeamFunctor_SolveUpperLU.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

#include "Tacho_SupernodeInfo.hpp"

namespace Tacho {

template <typename SupernodeInfoType> struct TeamFunctor_SolveUpperLU {
public:
  using range_type = Kokkos::pair<ordinal_type, ordinal_type>;

  using supernode_info_type = SupernodeInfoType;
  using supernode_type = typename supernode_info_type::supernode_type;
  using supernode_type_array = typename supernode_info_type::supernode_type_array;

  using ordinal_type_array = typename supernode_info_type::ordinal_type_array;
  using size_type_array = typename supernode_info_type::size_type_array;

  using ordinal_pair_type_array = typename supernode_info_type::ordinal_pair_type_array;

  using value_type = typename supernode_info_type::value_type;
  using value_type_array = typename supernode_info_type::value_type_array;
  using value_type_matrix = typename supernode_info_type::value_type_matrix;

private:
  ConstUnmanagedViewType<supernode_type_array> _supernodes;
  ConstUnmanagedViewType<ordinal_type_array> _gid_colidx;

  ConstUnmanagedViewType<ordinal_type_array> _compute_mode, _level_sids;
  ordinal_type _pbeg, _pend;

  UnmanagedViewType<value_type_matrix> _t;
  ordinal_type _nrhs;

  UnmanagedViewType<size_type_array> _buf_ptr;
  UnmanagedViewType<value_type_array> _buf;

public:
  KOKKOS_INLINE_FUNCTION
  TeamFunctor_SolveUpperLU() = delete;

  KOKKOS_INLINE_FUNCTION
  TeamFunctor_SolveUpperLU(const supernode_info_type &info, const ordinal_type_array &compute_mode,
                           const ordinal_type_array &level_sids, const value_type_matrix &t,
                           const value_type_array &buf)
      : _supernodes(info.supernodes), _gid_colidx(info.gid_colidx), _compute_mode(compute_mode),
        _level_sids(level_sids), _t(t), _nrhs(t.extent(1)), _buf(buf) {}

  inline void setRange(const ordinal_type pbeg, const ordinal_type pend) {
    _pbeg = pbeg;
    _pend = pend;
  }

  inline void setBufferPtr(const size_type_array &buf_ptr) { _buf_ptr = buf_ptr; }

  ///
  /// Algorithm Variant 0: gemv - trsv
  ///
  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void solve(MemberType &member, const supernode_type &s, value_type *bptr) const {
    using TrsvAlgoType = typename TrsvAlgorithm::type;
    using GemvAlgoType = typename GemvAlgorithm::type;

    const value_type minus_one(-1), one(1);
    {
      const ordinal_type m = s.m, n = s.n, n_m = n - m;
      if (m > 0) {
        value_type *uptr = s.u_buf;
        // solve
        const UnmanagedViewType<value_type_matrix> AL(uptr, m, m);
        uptr += m * m;
        const ordinal_type offm = s.row_begin;
        const auto tT = Kokkos::subview(_t, range_type(offm, offm + m), Kokkos::ALL());

        if (n_m > 0) {
          // update
          const UnmanagedViewType<value_type_matrix> AR(uptr, m, n_m); // uptr += m*n;
          const UnmanagedViewType<value_type_matrix> bB(bptr, n_m, _nrhs);
          Gemv<Trans::NoTranspose, GemvAlgoType>::invoke(member, minus_one, AR, bB, one, tT);
          member.team_barrier();
        }
        Trsv<Uplo::Upper, Trans::NoTranspose, TrsvAlgoType>::invoke(member, Diag::NonUnit(), AL, tT);
      }
    }
  }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void update(MemberType &member, const supernode_type &s, value_type *bptr) const {
    {
      const ordinal_type m = s.m, n = s.n, n_m = n - m;
      if (n_m > 0) {
        // update
        const UnmanagedViewType<value_type_matrix> bB(bptr, n_m, _nrhs);
        const ordinal_type goffset = s.gid_col_begin + s.m;
        Kokkos::parallel_for(
            Kokkos::TeamVectorRange(member, n_m),
            [&,
             goffset](const ordinal_type &i) { // Value capture is a workaround for cuda + gcc-7.2 compiler bug w/c++14
              // for (ordinal_type i=0;i<n;++i) {
              const ordinal_type row = _gid_colidx(i + goffset);
              for (ordinal_type j = 0; j < _nrhs; ++j)
                bB(i, j) = _t(row, j);
            });
      }
    }
  }

  struct SolveTag {};
  struct UpdateTag {};
  struct DummyTag {};

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const SolveTag &, const MemberType &member) const {
    const ordinal_type p = _pbeg + member.league_rank();
    const ordinal_type sid = _level_sids(p);
    const ordinal_type mode = _compute_mode(sid);
    if (p < _pend && mode == 1) {
      const supernode_type &s = _supernodes(sid);
      value_type *bptr = _buf.data() + _buf_ptr(member.league_rank());
      solve(member, s, bptr);
    } else if (mode == -1) {
      printf("Error: TeamFunctorSolveUpperChol::SolveTag, computing mode is not determined\n");
    } else {
      // skip
    }
  }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const UpdateTag &, const MemberType &member) const {
    const ordinal_type p = _pbeg + member.league_rank();
    if (p < _pend) {
      const ordinal_type sid = _level_sids(p);
      const supernode_type &s = _supernodes(sid);
      value_type *bptr = _buf.data() + _buf_ptr(member.league_rank());
      update(member, s, bptr);
    } else {
      // skip
    }
  }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const DummyTag &, const MemberType &member) const {}
};

} // namespace Tacho

#endif
