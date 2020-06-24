#ifndef USE_STK_SIMD_NONE  // somehow these tests are broken for scalar SIMD

#include <gtest/gtest.h>                // for AssertHelper, etc
#include <iomanip>                      // for operator<<
#include <iostream>                     // for basic_ostream::operator<<, etc
#include <vector>                       // for vector
#include <stk_simd_view/simd_view.hpp>
#include <stk_simd_view/simd_parallel.hpp>
#include <stk_simd_view/simd_while.hpp>
#include <cmath>
#include <algorithm>
#include "StkSimdViewFixture.hpp"


stk::simd::View<double*[3]> create_velocity_field(int numNodes) {
  stk::simd::View<double*[3]> veloField("velo_field", numNodes);
  stk::simd::parallel_for("fill_velo", numNodes, STK_LAMBDA(stk::simd::DeviceIndex i) {
    // Use DeviceIndex and DeviceDouble so that it will work on both CPU and GPU.  
    // For CPU, these have a simd length, on GPU they have length 1.
    stk::simd::DeviceDouble one = 1.0;
    veloField(i, 0) = one;
    veloField(i, 1) = 2*one;
    veloField(i, 2) = one+2;
  });
  return veloField;
}

TEST_F(StkSimdViewFixture, SimdParallelFor) {
  int numNodes = 123;
  auto veloField = create_velocity_field(numNodes);
  auto veloFieldHost = copy_from_device(veloField);
  for (int i=0; i < numNodes; ++i) {
    EXPECT_EQ(1, veloFieldHost(i,0));
    EXPECT_EQ(2, veloFieldHost(i,1));
    EXPECT_EQ(3, veloFieldHost(i,2));
  }
}


double reduce_sum_velocity_field(stk::simd::View<double*[3]> veloField, int direction) {
  double reducedVelo=0.0;
  stk::simd::parallel_reduce_sum("reduce_velo", veloField.extent(0),
                                 STK_LAMBDA(stk::simd::DeviceIndex i, stk::simd::DeviceDouble& sum) {
    sum += veloField(i, direction);
  }, reducedVelo);
  return reducedVelo;
}

TEST_F(StkSimdViewFixture, SimdParallelReduce) {
  int numNodes = 123;
  auto veloField = create_velocity_field(numNodes);
  EXPECT_EQ(1*numNodes, reduce_sum_velocity_field(veloField, 0));
  EXPECT_EQ(2*numNodes, reduce_sum_velocity_field(veloField, 1));
  EXPECT_EQ(3*numNodes, reduce_sum_velocity_field(veloField, 2));
}

#ifndef KOKKOS_ENABLE_CUDA

// simpler for_each and reduce_sum_each only work on host for now.

double double_and_reduce_sum_each(stk::simd::View<double*[3]>::HostMirror veloField, int direction) {
  const int numNodes = veloField.extent(0);
  stk::simd::for_each(numNodes, [&](stk::simd::Index i) {
    veloField(i, direction) *= 2;
  });
  return stk::simd::reduce_sum_each(numNodes, [&](stk::simd::Index i) {
    stk::simd::Double val = veloField(i, direction);
    return val;
  });
}

TEST_F(StkSimdViewFixture, ForEach_ReduceSumEach) {
  int numNodes = 123;
  auto veloField = create_velocity_field(numNodes);
  EXPECT_EQ(2*numNodes, double_and_reduce_sum_each(veloField, 0));
  EXPECT_EQ(4*numNodes, double_and_reduce_sum_each(veloField, 1));
  EXPECT_EQ(6*numNodes, double_and_reduce_sum_each(veloField, 2));
}

#endif

template <typename Real, typename Layout=void>
class StkSimdView1dTester {
 public:

  StkSimdView1dTester(int size) : N(size) {}

  typedef stk::simd::DeviceIndex DeviceIndex;
  typedef typename stk::simd::DeviceTraits<Real>::simd_type DeviceReal;

  typedef Kokkos::View<Real*> View;
  typedef Kokkos::View<const Real*> ConstView;

  typedef stk::simd::View<Real*, Layout> SimdView;
  typedef stk::simd::View<const Real*, Layout> ConstSimdView;
  typedef typename stk::simd::View<Real*, Layout>::HostMirror HostSimdView;
  typedef typename stk::simd::View<const Real*, Layout>::HostMirror ConstHostSimdView;

  void const_cast_test() const {
    View a = get_view();
    ConstView b = a;
    SimdView aSimd = get_simd_view();
    ConstSimdView bSimd = aSimd;
  }

  void host_mirror_test() const {
    SimdView a = get_simd_view();
    auto aHost = stk::simd::create_mirror_view(a);
    stk::simd::deep_copy(aHost, a);
    for (int i=0; i < N; ++i) {
      EXPECT_EQ(i, aHost(i));
    }
  }

  void parallel_for_test() const {
    SimdView a = get_simd_view();
    SimdView b = get_simd_view();
    
    set_b_to_a_plus_b_lambda(a, b);
    test_view_equal_to_index_multiple(stk::simd::copy_from_device(b), 2);
    
    set_b_to_a_plus_b_functor(a, b);
    test_view_equal_to_index_multiple(stk::simd::copy_from_device(b), 3);

    set_b_to_2a_plus_b_functor_with_tag(a, b);
    test_view_equal_to_index_multiple(stk::simd::copy_from_device(b), 5);
    
    auto aHost = stk::simd::copy_from_device(a);
    auto bHost = stk::simd::copy_from_device(b);
    set_b_to_3a_plus_b_for_each_lambda(aHost, bHost);
    test_view_equal_to_index_multiple(bHost, 8);
  }
  
  void parallel_reduce_test(int loopSize) const {

    ASSERT_LE(loopSize, N);
    SimdView a = get_simd_view();

    EXPECT_EQ( running_sum(loopSize-1), reduce_sum_lambda(a, loopSize) ); 
    EXPECT_EQ( running_sum(loopSize-1), reduce_sum_functor(a, loopSize) );
    EXPECT_EQ( running_sum(loopSize-1), reduce_sum_functor_with_tag(a, loopSize) );
    auto aHost = stk::simd::copy_from_device(a);
    EXPECT_EQ( running_sum(loopSize-1), reduce_sum_each_lambda(aHost, loopSize) );
  }
  
  View get_view() const {
    View data("data", N);
    Kokkos::parallel_for("fill_1d_view", N, STK_LAMBDA(const int i) {
      data(i) = i;
    });
    return data;
  }

  SimdView get_simd_view() const {
    SimdView simdData("simd_data", N);
    Kokkos::parallel_for("fill_1d_simd_view", N, STK_LAMBDA(const int i) {
      simdData(i) = i;
    });
    return simdData;
  }

  struct SomeTag {};
    
  struct AddFunctor {

    AddFunctor(ConstSimdView a_, SimdView b_) : a(a_), b(b_) {}

    STK_INLINE
    void operator() (const DeviceIndex& i) const {
      DeviceReal tmp = a(i)+b(i);
      b(i) = tmp;
    }

    STK_INLINE
    void operator() (SomeTag, const DeviceIndex& i) const {
      DeviceReal tmp = 2*a(i)+b(i);
      b(i) = tmp;
    }

   private:
    ConstSimdView a;
    SimdView b;
  };

  void set_b_to_a_plus_b_lambda(ConstSimdView a, SimdView b) const {
    stk::simd::parallel_for<Real>("b=a+b", N, STK_LAMBDA(const DeviceIndex& i) {
      DeviceReal tmp = a(i)+b(i);
      b(i) = tmp;
    });
  }

  void set_b_to_a_plus_b_functor(ConstSimdView a, SimdView b) const {
    stk::simd::parallel_for<Real>("b=a+b", N, AddFunctor(a,b));
  }

  void set_b_to_2a_plus_b_functor_with_tag(ConstSimdView a, SimdView b) const {
    stk::simd::parallel_for<Real>("b=2a+b",
                                  Kokkos::RangePolicy<SomeTag>(0, N), 
                                  AddFunctor(a,b));
  }

  // figure out a float version
  void set_b_to_3a_plus_b_for_each_lambda(ConstHostSimdView a,
                                          HostSimdView b) const {
    stk::simd::for_each(N, [&](const DeviceIndex& i) {
      b(i) = 3*a(i)+b(i);
    });
  }
    
  struct ReduceSumFunctor {
    ReduceSumFunctor(SimdView a_) : a(a_) {}
    STK_INLINE void operator() (const DeviceIndex& i, DeviceReal& v) const {
      v += a(i);
    }
   private:
    SimdView a;
  };

  struct ReduceSumFunctorWithTag {
    ReduceSumFunctorWithTag(SimdView a_) : a(a_) {}
    STK_INLINE void operator() (SomeTag, const DeviceIndex& i, DeviceReal& v) const {
      v += a(i);
    }
   private:
    SimdView a;
  };

  Real reduce_sum_lambda(SimdView a, int loopSize) const {
    Real reducedRunningSum=0;
    stk::simd::parallel_reduce_sum("running_sum", loopSize, 
                                   STK_LAMBDA(const DeviceIndex& i, DeviceReal& v) {
      v += a(i);
    }, reducedRunningSum);
    return reducedRunningSum;
  }

  Real reduce_sum_functor(SimdView a, int loopSize) const {
    Real reducedRunningSum=0.0;
    stk::simd::parallel_reduce_sum("running_sum", loopSize,
                                   ReduceSumFunctor(a), reducedRunningSum);
    return reducedRunningSum;
  }

  Real reduce_sum_functor_with_tag(SimdView a, int loopSize) const {
    Real reducedRunningSum=0.0;
    stk::simd::parallel_reduce_sum("running_sum", 
                                   Kokkos::RangePolicy<SomeTag>(0, loopSize),
                                   ReduceSumFunctorWithTag(a),
                                   reducedRunningSum);
    return reducedRunningSum;
  }

  Real reduce_sum_each_lambda(HostSimdView a, int loopSize) const {
    return stk::simd::reduce_sum_each(loopSize, [&](const DeviceIndex& i) {
      return a(i);
    });
  }

 private:

  void test_view_equal_to_index_multiple(HostSimdView view, int n) const {
    for (int i=0; i < N; ++i) {
      EXPECT_EQ(n*i, view(i));
    }
  }
  int N;
};



TEST_F(StkSimdViewFixture, SimdParallelFor_DefaultLayout) {
  StkSimdView1dTester<double> tester(237);
  tester.const_cast_test();
  tester.host_mirror_test();
  tester.parallel_for_test();
  tester.parallel_reduce_test(237);
  tester.parallel_reduce_test(236);
  tester.parallel_reduce_test(3);
  tester.parallel_reduce_test(64);
}

TEST_F(StkSimdViewFixture, SimdParallelFor_LayoutRight) {
  StkSimdView1dTester<double, stk::simd::LayoutRight<double> > tester(321);
  tester.const_cast_test();
  tester.host_mirror_test();
  tester.parallel_for_test();
  tester.parallel_reduce_test(321);
  tester.parallel_reduce_test(36);
  tester.parallel_reduce_test(64);
  tester.parallel_reduce_test(121);
}

TEST_F(StkSimdViewFixture, SimdParallelFor_LayoutLeft) {
  StkSimdView1dTester<double, stk::simd::LayoutLeft<double> > tester(139);
  tester.const_cast_test();
  tester.host_mirror_test();
  tester.parallel_for_test();
  tester.parallel_reduce_test(139);
  tester.parallel_reduce_test(89);
  tester.parallel_reduce_test(100);
  tester.parallel_reduce_test(64);
}


// permutations:
// * View type: View ranks 1->3 
// * View type: Const, non-const Views ^
// * Casting: from non-const to const ^
// * Mirror view ^
// * Data type: floats, doubles
// * Layout: default, layout left, layout right ^
// * simd_for: for, reduce ^
// * parallel_for: for, reduce ^
// * timing vs. standard view

#endif
