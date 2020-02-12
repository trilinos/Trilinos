#ifndef __TACHO_TEAMFUNCTOR_SOLVE_UPPER_CHOL_HPP__
#define __TACHO_TEAMFUNCTOR_SOLVE_UPPER_CHOL_HPP__

/// \file Tacho_TeamFunctor_SolveUpperChol.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

#include "Tacho_SupernodeInfo.hpp"

namespace Tacho {

  template<typename SupernodeInfoType>
  struct TeamFunctor_SolveUpperChol {
  public:
    typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;

    typedef SupernodeInfoType supernode_info_type;
    typedef typename supernode_info_type::supernode_type supernode_type;

    typedef typename supernode_info_type::ordinal_type_array ordinal_type_array;
    typedef typename supernode_info_type::size_type_array size_type_array;

    typedef typename supernode_info_type::value_type value_type;
    typedef typename supernode_info_type::value_type_array value_type_array;
    typedef typename supernode_info_type::value_type_matrix value_type_matrix;

    typedef typename std::conditional
    <std::is_same<Kokkos::Impl::ActiveExecutionMemorySpace,Kokkos::HostSpace>::value,
     Algo::External,Algo::Internal>::type TrsvAlgoType;

    typedef typename std::conditional
    <std::is_same<Kokkos::Impl::ActiveExecutionMemorySpace,Kokkos::HostSpace>::value,
     Algo::External,Algo::Internal>::type GemvAlgoType;

  private:
    supernode_info_type _info;
    ordinal_type_array _compute_mode, _level_sids;
    ordinal_type _pbeg, _pend;

    value_type_matrix _t;
    ordinal_type _nrhs;

    size_type_array _buf_ptr;
    value_type_array _buf;

  public:
    KOKKOS_INLINE_FUNCTION
    TeamFunctor_SolveUpperChol() = delete;

    KOKKOS_INLINE_FUNCTION
    TeamFunctor_SolveUpperChol(const supernode_info_type &info,
                               const ordinal_type_array &compute_mode,
                               const ordinal_type_array &level_sids,
                               const value_type_matrix t,
                               const size_type_array buf_ptr,
                               const value_type_array buf)
      :
      _info(info),
      _compute_mode(compute_mode),
      _level_sids(level_sids),
      _t(t),
      _nrhs(t.extent(1)),
      _buf_ptr(buf_ptr),
      _buf(buf)
    {}

    inline
    void setRange(const ordinal_type pbeg,
                  const ordinal_type pend) {
      _pbeg = pbeg; _pend = pend;
    }

    ///
    /// Algorithm Variant 0: gemv - trsv
    ///
    template<typename MemberType>
    KOKKOS_INLINE_FUNCTION
    void solve_var0(MemberType &member, const ordinal_type sid) const {
      const value_type minus_one(-1), one(1);
      const auto &s = _info.supernodes(sid);
      {
        const ordinal_type m = s.m, n = s.n, n_m = n-m;
        if (m > 0) {
          value_type *aptr = s.buf;
          // solve
          const UnmanagedViewType<value_type_matrix> AL(aptr, m, m); aptr += m*m;
          const ordinal_type offm = s.row_begin;
          const auto tT = Kokkos::subview(_t, range_type(offm, offm+m), Kokkos::ALL());

          if (n_m > 0) {
            value_type *bptr = _buf.data() + _buf_ptr(sid);
            // update
            const UnmanagedViewType<value_type_matrix> AR(aptr, m, n_m); // aptr += m*n;
            const UnmanagedViewType<value_type_matrix> bB(bptr, n_m, _nrhs);
            Gemv<Trans::NoTranspose,GemvAlgoType>
              ::invoke(member, minus_one, AR, bB, one, tT);
          }
          Trsv<Uplo::Upper,Trans::NoTranspose,TrsvAlgoType>
            ::invoke(member, Diag::NonUnit(), AL, tT);
        }
      }
    }

    template<typename MemberType>
    KOKKOS_INLINE_FUNCTION
    void update_var0(MemberType &member, const ordinal_type sid) const {
      const auto &s = _info.supernodes(sid);
      {
        const ordinal_type m = s.m, n = s.n, n_m = n-m;
        if (n_m > 0) {
          value_type *bptr = _buf.data() + _buf_ptr(sid);
          // update
          const UnmanagedViewType<value_type_matrix> bB(bptr, n_m, _nrhs);
          const ordinal_type goffset = s.gid_col_begin + s.m;
          Kokkos::parallel_for
            (Kokkos::TeamVectorRange(member, n_m),
             [&](const ordinal_type &i) {
              //for (ordinal_type i=0;i<n;++i) {
              const ordinal_type row = _info.gid_colidx(i+goffset);
              for (ordinal_type j=0;j<_nrhs;++j)
                bB(i,j) = _t(row,j);
            });
        }
      }
    }

    ///
    /// Algorithm Variant 1: gemv - gemv
    ///
    template<typename MemberType>
    KOKKOS_INLINE_FUNCTION
    void solve_var1(MemberType &member, const ordinal_type sid) const {
      const value_type minus_one(-1), one(1), zero(0);
      const auto &s = _info.supernodes(sid);
      {
        const ordinal_type m = s.m, n = s.n, n_m = n-m;
        if (m > 0) {
          value_type *aptr = s.buf, *bptr = _buf.data() + _buf_ptr(sid);

          const UnmanagedViewType<value_type_matrix> AL(aptr, m, m); aptr += m*m;
          UnmanagedViewType<value_type_matrix> b(bptr, n, _nrhs);
          auto bT = Kokkos::subview(b, range_type(0, m), Kokkos::ALL());

          const ordinal_type offm = s.row_begin;
          const auto tT = Kokkos::subview(_t, range_type(offm, offm+m), Kokkos::ALL());

          if (n_m > 0) {
            // update
            const UnmanagedViewType<value_type_matrix> AR(aptr, m, n_m); // aptr += m*n;
            auto bB = Kokkos::subview(b, range_type(m, n), Kokkos::ALL());
            Gemv<Trans::NoTranspose,GemvAlgoType>
              ::invoke(member, minus_one, AR, bB, one, tT);
          }
          Gemv<Trans::NoTranspose,GemvAlgoType>
            ::invoke(member, one, AL, tT, zero, bT);
          member.team_barrier();
          // copy to t
          Kokkos::parallel_for
            (Kokkos::TeamVectorRange(member, m*_nrhs),
             [&](const ordinal_type &k) {
              const ordinal_type i = k%m, j = k/m;
              tT(i,j) = bT(i,j);
            });
        }
      }
    }

    template<typename MemberType>
    KOKKOS_INLINE_FUNCTION
    void update_var1(MemberType &member, const ordinal_type sid) const {
      const auto &s = _info.supernodes(sid);
      {
        const ordinal_type m = s.m, n = s.n, n_m = n-m;
        if (n_m > 0) {
          value_type *bptr = _buf.data() + _buf_ptr(sid);
          UnmanagedViewType<value_type_matrix> b(bptr, n, _nrhs);
          auto bB = Kokkos::subview(b, range_type(m, n), Kokkos::ALL());

          // update
          const ordinal_type goffset = s.gid_col_begin + s.m;
          Kokkos::parallel_for
            (Kokkos::TeamVectorRange(member, n_m),
             [&](const ordinal_type &i) {
              //for (ordinal_type i=0;i<n;++i) {
              const ordinal_type row = _info.gid_colidx(i+goffset);
              for (ordinal_type j=0;j<_nrhs;++j)
                bB(i,j) = _t(row,j);
            });
        }
      }
    }

    ///
    /// Algorithm Variant 2: gemv
    ///
    template<typename MemberType>
    KOKKOS_INLINE_FUNCTION
    void solve_var2(MemberType &member, const ordinal_type sid) const {
      const value_type one(1), zero(0);
      const auto &s = _info.supernodes(sid);
      {
        const ordinal_type m = s.m, n = s.n; //, n_m = n-m;
        if (m > 0 && n > 0) {
          value_type *aptr = s.buf, *bptr = _buf.data() + _buf_ptr(sid);
          UnmanagedViewType<value_type_matrix> A(aptr, m, n); // aptr += m*n;
          UnmanagedViewType<value_type_matrix> t(bptr, n, _nrhs); bptr += n*_nrhs;
          UnmanagedViewType<value_type_matrix> b(bptr, n, _nrhs);
          Gemv<Trans::NoTranspose,GemvAlgoType>
            ::invoke(member, one, A, t, zero, b);
          member.team_barrier();

          const ordinal_type offm = s.row_begin;
          auto tT = Kokkos::subview(_t, range_type(offm, offm+m), Kokkos::ALL());

          // copy to t
          Kokkos::parallel_for
            (Kokkos::TeamVectorRange(member, m*_nrhs),
             [&](const ordinal_type &k) {
              const ordinal_type i = k%m, j = k/m;
              tT(i,j) = b(i,j);
            });
        }
      }
    }

    template<typename MemberType>
    KOKKOS_INLINE_FUNCTION
    void update_var2(MemberType &member, const ordinal_type sid) const {
      const auto &s = _info.supernodes(sid);
      {
        const ordinal_type m = s.m, n = s.n, n_m = n-m;
        if (m > 0) {
          value_type *bptr = _buf.data() + _buf_ptr(sid);
          UnmanagedViewType<value_type_matrix> t(bptr, n, _nrhs); bptr += n*_nrhs;

          auto tT = Kokkos::subview(t, range_type(0, m), Kokkos::ALL());

          const ordinal_type offm = s.row_begin;
          auto tT_src = Kokkos::subview(_t, range_type(offm, offm+m), Kokkos::ALL());

          // copy to t
          Kokkos::parallel_for
            (Kokkos::TeamVectorRange(member, m*_nrhs),
             [&](const ordinal_type &k) {
              const ordinal_type i = k%m, j = k/m;
              tT(i,j) = tT_src(i,j);
            });

          if (n_m > 0) {
            auto tB = Kokkos::subview(t, range_type(m, n), Kokkos::ALL());
            // update
            const ordinal_type goffset = s.gid_col_begin + s.m;
            Kokkos::parallel_for
              (Kokkos::TeamVectorRange(member, n_m),
               [&](const ordinal_type &i) {
                //for (ordinal_type i=0;i<n;++i) {
                const ordinal_type row = _info.gid_colidx(i+goffset);
                for (ordinal_type j=0;j<_nrhs;++j)
                  tB(i,j) = _t(row,j);
              });
          }
        }
      }
    }

    template<int Var> struct SolveTag  { enum { variant = Var }; };
    template<int Var> struct UpdateTag { enum { variant = Var }; };
    struct DummyTag {};

    template<typename MemberType, int Var>
    KOKKOS_INLINE_FUNCTION
    void operator()(const SolveTag<Var> &, const MemberType &member) const {
      const ordinal_type p = _pbeg + member.league_rank();
      const ordinal_type sid = _level_sids(p);
      const ordinal_type mode = _compute_mode(sid);
      if (p < _pend && mode == 1) {
        typedef SolveTag<Var> solve_tag_type;
        if      (solve_tag_type::variant == 0) solve_var0(member, sid);
        else if (solve_tag_type::variant == 1) solve_var1(member, sid);
        else if (solve_tag_type::variant == 2) solve_var2(member, sid);
        else printf("Error: abort\n");
      } else if (mode == -1) {
        printf("Error: abort\n");
      } else {
        // skip
      }
    }

    template<typename MemberType, int Var>
    KOKKOS_INLINE_FUNCTION
    void operator()(const UpdateTag<Var> &, const MemberType &member) const {
      const ordinal_type p = _pbeg + member.league_rank();
      if (p < _pend) {
        const ordinal_type sid = _level_sids(p);
        typedef UpdateTag<Var> update_tag_type;
        if      (update_tag_type::variant == 0) update_var0(member, sid);
        else if (update_tag_type::variant == 1) update_var1(member, sid);
        else if (update_tag_type::variant == 2) update_var2(member, sid);
        else printf("Error: abort\n");
      } else {
        // skip
      }
    }

    template<typename MemberType>
    KOKKOS_INLINE_FUNCTION
    void operator()(const DummyTag &, const MemberType &member) const {

    }

  };

}

#endif




// template<typename MemberType>
// KOKKOS_INLINE_FUNCTION
// void solve_and_update_recursive(MemberType &member, const ordinal_type sid) const {
//   update(member, sid);
//   solve(member, sid);

//   const auto &s = _info.supernodes(sid);
//   for (ordinal_type i=0;i<s.nchildren;++i)
//     solve_and_update_recursive(member, s.children[i]);
// }

// template<typename MemberType>
// KOKKOS_INLINE_FUNCTION
// void operator()(const SolveUpdateTag &, const MemberType &member) const {
//   const ordinal_type p = _pbeg + member.league_rank();
//   const ordinal_type sid = _level_sids(p);
//   const ordinal_type mode = _compute_mode(sid);

//   // supernodes below this level are all solved and updated
//   if (p < _pend && mode == 2) {
//     solve_and_update_recursive(member, sid);
//   } else if (mode == -1) {
//     printf("Error: abort\n");
//   } else {
//     // skip
//   }
// }
