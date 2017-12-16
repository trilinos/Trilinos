/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "gtest/gtest.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"

//#include "KokkosBatched_Vector.hpp"

#include "KokkosBatched_Trsv_Decl.hpp"
#include "KokkosBatched_Trsv_Serial_Impl.hpp"
#include "KokkosBatched_Trsv_Team_Impl.hpp"

#include "KokkosKernels_TestUtils.hpp"

using namespace KokkosBatched::Experimental;

namespace Test {

  template<typename U, typename T, typename D>
  struct ParamTag {
    typedef U uplo;
    typedef T trans;
    typedef D diag;
  };

  template<typename DeviceType,
           typename ViewType,
           typename ScalarType,
           typename ParamTagType,
           typename AlgoTagType>
  struct Functor_TestBatchedTeamTrsv {
    ViewType _a, _b;
    
    ScalarType _alpha;

    KOKKOS_INLINE_FUNCTION
    Functor_TestBatchedTeamTrsv(const ScalarType alpha, 
            const ViewType &a,
            const ViewType &b) 
      : _a(a), _b(b), _alpha(alpha) {}

    template<typename MemberType>
    KOKKOS_INLINE_FUNCTION
    void operator()(const ParamTagType &, const MemberType &member) const {
      const int k = member.league_rank();

      auto aa = Kokkos::subview(_a, k, Kokkos::ALL(), Kokkos::ALL());
      auto bb = Kokkos::subview(_b, k, Kokkos::ALL(), 0);
      
      TeamTrsv<MemberType,
        typename ParamTagType::uplo,
        typename ParamTagType::trans,
        typename ParamTagType::diag,
        AlgoTagType>::
        invoke(member, _alpha, aa, bb);
    }

    inline
    void run() {
      const int league_size = _b.dimension_0();
      Kokkos::TeamPolicy<DeviceType,ParamTagType> policy(league_size, Kokkos::AUTO);
      Kokkos::parallel_for(policy, *this);
    }
  };

  template<typename DeviceType,
           typename ViewType,
           typename ScalarType,
           typename ParamTagType,
           typename AlgoTagType>
  void impl_test_batched_trsv(const int N, const int BlkSize) {
    typedef typename ViewType::value_type value_type;
    typedef Kokkos::Details::ArithTraits<value_type> ats;

    /// randomized input testing views
    ScalarType alpha(1.5);

    ViewType
      a0("a0", N, BlkSize, BlkSize), a1("a1", N, BlkSize, BlkSize),
      b0("b0", N, BlkSize, 1),       b1("b1", N, BlkSize, 1);

    Kokkos::Random_XorShift64_Pool<typename DeviceType::execution_space> random(13718);
    Kokkos::fill_random(a0, random, value_type(1.0));
    Kokkos::fill_random(b0, random, value_type(1.0));

    Kokkos::fence();

    Kokkos::deep_copy(b0, 1.0);

    Kokkos::deep_copy(a1, a0);
    Kokkos::deep_copy(b1, b0);

    Functor_TestBatchedTeamTrsv<DeviceType,ViewType,ScalarType,
      ParamTagType,Algo::Trsv::Unblocked>(alpha, a0, b0).run();
    Functor_TestBatchedTeamTrsv<DeviceType,ViewType,ScalarType,
      ParamTagType,AlgoTagType>(alpha, a1, b1).run();

    Kokkos::fence();

    /// for comparison send it to host
    typename ViewType::HostMirror b0_host = Kokkos::create_mirror_view(b0);
    typename ViewType::HostMirror b1_host = Kokkos::create_mirror_view(b1);

    Kokkos::deep_copy(b0_host, b0);
    Kokkos::deep_copy(b1_host, b1);

    /// check b0 = b1 ; this eps is about 10^-14
    typedef typename ats::mag_type mag_type;
    mag_type sum(1), diff(0);
    const mag_type eps = 1.0e3 * ats::epsilon();

    for (int k=0;k<N;++k)
      for (int i=0;i<BlkSize;++i)
        for (int j=0;j<1;++j) {
          sum  += ats::abs(b0_host(k,i,j));
          diff += ats::abs(b0_host(k,i,j)-b1_host(k,i,j));
        }
    EXPECT_NEAR_KK( diff/sum, 0.0, eps);
  }
}


template<typename DeviceType,
         typename ValueType,
         typename ScalarType,
         typename ParamTagType,
         typename AlgoTagType>
int test_batched_trsv() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    typedef Kokkos::View<ValueType***,Kokkos::LayoutLeft,DeviceType> ViewType;
    Test::impl_test_batched_trsv<DeviceType,ViewType,ScalarType,ParamTagType,AlgoTagType>(     0, 10);
    for (int i=0;i<10;++i) {
      printf("Testing: LayoutLeft,  Blksize %d\n", i);
      Test::impl_test_batched_trsv<DeviceType,ViewType,ScalarType,ParamTagType,AlgoTagType>(1024,  i);
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    typedef Kokkos::View<ValueType***,Kokkos::LayoutRight,DeviceType> ViewType;
    Test::impl_test_batched_trsv<DeviceType,ViewType,ScalarType,ParamTagType,AlgoTagType>(     0, 10);
    for (int i=0;i<10;++i) {
      printf("Testing: LayoutRight, Blksize %d\n", i);
      Test::impl_test_batched_trsv<DeviceType,ViewType,ScalarType,ParamTagType,AlgoTagType>(1024,  i);
    }
  }
#endif

  return 0;
}

