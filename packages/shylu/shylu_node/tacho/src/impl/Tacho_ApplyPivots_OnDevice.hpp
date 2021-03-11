#ifndef __TACHO_APPLY_PIVOTS_ON_DEVICE_HPP__
#define __TACHO_APPLY_PIVOTS_ON_DEVICE_HPP__


/// \file  Tacho_ApplyPivots_OnDevice.hpp
/// \brief Apply pivots on device
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

  template<int BaseIndex>
  struct ApplyPivots<BadseIndex,Side::Left,Direct::Forward,,Algo::OnDevice> {
    template<typename MemberType,
             typename ViewTypeP,
             typename ViewTypeA>
    inline
    static int
    invoke(MemberType &member,
           const ViewTypeP &P,
           const ViewTypeA &A) {
      typedef typename ViewTypeA::non_const_value_type value_type;

      const ordinal_type
        m = A.extent(0),
        n = A.extent(1);
      
      if (m == P.extent(0)) {
        if (A.span() > 0) {
          using exec_space = MemberType;
          using team_policy_type = Kokkos::TeamPolicy<exec_space>;
          
          const auto exec_instance = member;
          const auto policy = team_policy_type(exec_instance, n, Kokkos::AUTO);
          Kokkos::parallel_for
            (policy, KOKKOS_LAMBDA(const typename team_policy_type::member_type &member) {
              const ordinal_type j = member.league_rank();
              for (ordinal_type i=0;i<m;) {
                const ordinal_type piv = P(i);
                if (piv == i) {
                  /// no pivot
                  ++i;
                } else if (piv > 0) {
                  /// 1x1 pivot
                  const ordinal_type p = piv - BaseIndex;
                  const value_type tmp = A(i,j);
                  A(i,j) = A(p,j);
                  A(p,j) = tmp;
                  ++i;
                } else {
                  /// 2x2 pivot
                  const int p = -piv - BaseIndex;
                  const value_type tmp_a = A(i,j), tmp_b = A(i+1,j);
                  A(i  ,j) = A(p  ,j);
                  A(i+1,j) = A(p+1,j);
                  A(p,  j) = tmp_a;
                  A(p+1,j) = tmp_b;
                  i+=2;           
                }
              }

            });
        }
      } else {
        TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, "A is not a square matrix");
      }
      return 0;
    }
  };

}
#endif
