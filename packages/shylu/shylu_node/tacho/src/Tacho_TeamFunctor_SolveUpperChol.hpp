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

    template<typename MemberType>
    KOKKOS_INLINE_FUNCTION
    void solve(MemberType &member, const ordinal_type sid) const {
      const value_type minus_one(-1), one(1);
      const auto &s = _info.supernodes(sid);
      value_type *ptr = s.buf; 
      {
        const ordinal_type m = s.m, n = s.n - s.m;
        {
          const UnmanagedViewType<value_type_matrix> tB(&_buf(_buf_ptr(sid)), n, _nrhs); 
          if (m > 0) {
            // solve
            const UnmanagedViewType<value_type_matrix> AL(ptr, m, m); ptr += m*m;
            const ordinal_type offm = s.row_begin;
            const auto tT = Kokkos::subview(_t, range_type(offm, offm+m), Kokkos::ALL());
              
            if (n > 0) {
              const UnmanagedViewType<value_type_matrix> AR(ptr, m, n); // ptr += m*n;
              Gemv<Trans::NoTranspose,GemvAlgoType>
                ::invoke(member, minus_one, AR, tB, one, tT);
            }
            Trsv<Uplo::Upper,Trans::NoTranspose,TrsvAlgoType>
              ::invoke(member, Diag::NonUnit(), AL, tT);
          }
        }
      }
    }

    template<typename MemberType>
    KOKKOS_INLINE_FUNCTION
    void update(MemberType &member, const ordinal_type sid) const {
      const auto &s = _info.supernodes(sid);
      //value_type *ptr = s.buf; 
      {
        const ordinal_type n = s.n - s.m;
        {
          if (n > 0) {
            // update
            const UnmanagedViewType<value_type_matrix> tB(&_buf(_buf_ptr(sid)), n, _nrhs); 
            const ordinal_type goffset = s.gid_col_begin + s.m;
            Kokkos::parallel_for
              (Kokkos::TeamVectorRange(member, n),
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

    template<typename MemberType>
    KOKKOS_INLINE_FUNCTION
    void solve_and_update_recursive(MemberType &member, const ordinal_type sid) const {
      update(member, sid);
      solve(member, sid);

      const auto &s = _info.supernodes(sid);
      for (ordinal_type i=0;i<s.nchildren;++i)
        solve_and_update_recursive(member, s.children[i]);
    }

    struct SolveTag {};
    struct UpdateTag {};
    struct SolveUpdateTag {};

    template<typename MemberType>
    KOKKOS_INLINE_FUNCTION
    void operator()(const SolveTag &, const MemberType &member) const {
      const ordinal_type p = _pbeg + member.league_rank();
      const ordinal_type sid = _level_sids(p);
      const ordinal_type mode = _compute_mode(sid);
      if (p < _pend && mode == 1) {
        solve(member, sid);
      } else if (mode == -1) {
        printf("Error: abort\n"); 
      } else {
        // skip
      }         
    }

    template<typename MemberType>
    KOKKOS_INLINE_FUNCTION
    void operator()(const UpdateTag &, const MemberType &member) const {
      const ordinal_type p = _pbeg + member.league_rank();
      if (p < _pend) { 
        const ordinal_type sid = _level_sids(p);
        update(member, sid);
      } else {
        // skip
      }         
    }
    
    template<typename MemberType>
    KOKKOS_INLINE_FUNCTION
    void operator()(const SolveUpdateTag &, const MemberType &member) const {
      const ordinal_type p = _pbeg + member.league_rank();
      const ordinal_type sid = _level_sids(p);
      const ordinal_type mode = _compute_mode(sid);

      // supernodes below this level are all solved and updated
      if (p < _pend && mode == 2) {
        solve_and_update_recursive(member, sid);
      } else if (mode == -1) {
        printf("Error: abort\n");
      } else {
        // skip
      }
        
    }

  };

}

#endif

            
            

