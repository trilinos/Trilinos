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


template <typename ArrayStyle, typename Real, typename Layout>
struct ViewMaker {
  stk::simd::View<ArrayStyle, Layout> make_view(int N, int M) {
    return stk::simd::View<ArrayStyle, Layout>("name", N);
  }
};

template <typename Real, typename Layout>
struct ViewMaker<Real**, Real, Layout> {
  stk::simd::View<Real**, Layout> make_view(int N, int M) {
    return stk::simd::View<Real**, Layout>("name", N, M);
  }
};


// need to allocate the view outside the class for some reason
template <typename ArrayStyle, typename Real, int dim, typename Layout>
stk::simd::View<ArrayStyle, Layout> get_2d_simd_view(int N, int scaling) {
  ViewMaker<ArrayStyle, Real, Layout> maker;
  stk::simd::View<ArrayStyle, Layout> a = maker.make_view(N, dim);
  Kokkos::parallel_for("2d_view_fill", N, KOKKOS_LAMBDA(const int i) {
    for (int j=0; j < dim; ++j) {
      a(i,j) = i+scaling*j;
    }
  });
  return a;
}


template <typename ArrayStyle, typename Real, int dim, typename Layout=void>
class StkSimdView2dTester {
 public:

  StkSimdView2dTester(int firstDim) : N(firstDim), scaling(100) {}

  typedef stk::simd::DeviceIndex DeviceIndex;
  typedef typename stk::simd::DeviceTraits<Real>::simd_type DeviceReal;

  typedef stk::simd::View<ArrayStyle, Layout> SimdView;
  typedef stk::simd::View<const Real*[dim], Layout> ConstSimdView;
  typedef typename stk::simd::View<ArrayStyle, Layout>::HostMirror HostSimdView;
  typedef typename stk::simd::View<const Real*[dim], Layout>::HostMirror ConstHostSimdView;

  void const_cast_test() const {
    SimdView aSimd = get_2d_simd_view<ArrayStyle, Real, dim, Layout>(N, scaling);
    ConstSimdView bSimd = aSimd;
  }

  void host_mirror_test() const {
    const SimdView a = get_2d_simd_view<ArrayStyle, Real, dim, Layout>(N, scaling);
    auto aHost = stk::simd::create_mirror_view(a);
    stk::simd::deep_copy(aHost, a);
    for (int i=0; i < N; ++i) {
      for (int j=0; j < dim; ++j) {
        EXPECT_EQ(i+scaling*j, aHost(i,j));
      }
    }
  }

  void parallel_for_test() const {
    SimdView a = get_2d_simd_view<ArrayStyle, Real,dim, Layout>(N, scaling);
    SimdView b = get_2d_simd_view<ArrayStyle, Real,dim, Layout>(N, scaling);
    
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
  
  void parallel_reduce_test(int loopSize, int j) const {
    ASSERT_LE(loopSize, N);
    ASSERT_LT(j, dim);
    SimdView a = get_2d_simd_view<ArrayStyle, Real, dim, Layout>(N, scaling);

    int expectedSum = running_sum(loopSize-1) + loopSize*scaling*j;
    EXPECT_EQ( expectedSum, reduce_sum_lambda(a, loopSize, j) ); 
    EXPECT_EQ( expectedSum, reduce_sum_functor(a, loopSize, j) );
    EXPECT_EQ( expectedSum, reduce_sum_functor_with_tag(a, loopSize, j) );
    auto aHost = stk::simd::copy_from_device(a);
    EXPECT_EQ( expectedSum, reduce_sum_each_lambda(aHost, loopSize, j) );
  }
    
  struct SomeTag {};
    
  struct AddFunctor {
    AddFunctor(ConstSimdView a_, SimdView b_) : a(a_), b(b_) {}

    STK_INLINE
    void operator() (const DeviceIndex& i) const {
      for (int j=0; j < dim; ++j) {
        DeviceReal tmp = a(i,j)+b(i,j);
        b(i,j) = tmp;
      }
    }

    STK_INLINE
    void operator() (SomeTag, const DeviceIndex& i) const {
      for (int j=0; j < dim; ++j) {
        DeviceReal tmp = 2*a(i,j)+b(i,j);
        b(i,j) = tmp;
      }
    }

   private:
    ConstSimdView a;
    SimdView b;
  };

  void set_b_to_a_plus_b_lambda(ConstSimdView a, SimdView b) const {
    stk::simd::parallel_for<Real>("b=a+b", N, STK_LAMBDA(const DeviceIndex& i) {
      for (int j=0; j < dim; ++j) {
        DeviceReal tmp = a(i,j)+b(i,j);
        b(i,j) = tmp;
      }
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
      for (int j=0; j < dim; ++j) {
        b(i,j) = 3*a(i,j)+b(i,j);
      }
    });
  }
    
  struct ReduceSumFunctor {
    ReduceSumFunctor(SimdView a_, int j_) : a(a_), j(j_) {}
    STK_INLINE void operator() (const DeviceIndex& i, DeviceReal& v) const {
      v += a(i,j);
    }
   private:
    SimdView a;
    int j;
  };

  struct ReduceSumFunctorWithTag {
    ReduceSumFunctorWithTag(SimdView a_, int j_) : a(a_), j(j_) {}
    STK_INLINE void operator() (SomeTag, const DeviceIndex& i, DeviceReal& v) const {
      v += a(i,j);
    }
   private:
    SimdView a;
    int j;
  };

  Real reduce_sum_lambda(SimdView a, int loopSize, int j) const {
    Real reducedRunningSum=0;
    stk::simd::parallel_reduce_sum("running_sum", loopSize, STK_LAMBDA(const DeviceIndex& i, DeviceReal& v) {
      v += a(i,j);
    }, reducedRunningSum);
    return reducedRunningSum;
  }

  Real reduce_sum_functor(SimdView a, int loopSize, int j) const {
    Real reducedRunningSum=0;
    stk::simd::parallel_reduce_sum("running_sum", loopSize,
                                   ReduceSumFunctor(a,j), reducedRunningSum);
    return reducedRunningSum;
  }

  Real reduce_sum_functor_with_tag(SimdView a, int loopSize, int j) const {
    Real reducedRunningSum=0;
    stk::simd::parallel_reduce_sum("running_sum", 
                                   Kokkos::RangePolicy<SomeTag>(0, loopSize),
                                   ReduceSumFunctorWithTag(a,j),
                                   reducedRunningSum);
    return reducedRunningSum;
  }

  Real reduce_sum_each_lambda(HostSimdView a, int loopSize, int j) const {
    return stk::simd::reduce_sum_each(loopSize, [&](const DeviceIndex& i) {
      return a(i,j);
    });
  }

 private:

  void test_view_equal_to_index_multiple(HostSimdView view, int n) const {
    for (int i=0; i < N; ++i) {
      for (int j=0; j < dim; ++j) {
        EXPECT_EQ( n*(i+scaling*j), view(i,j) );
      }
    }    
  }

  int N;
  int scaling;
};


TEST_F(StkSimdViewFixture, SimdParallelFor2d_DefaultLayout) {
  StkSimdView2dTester<double*[3], double, 3> tester(237);
  tester.const_cast_test();
  tester.host_mirror_test();
  tester.parallel_for_test();
  tester.parallel_reduce_test(237, 0);
  tester.parallel_reduce_test(3, 0);
  tester.parallel_reduce_test(64, 1);
  tester.parallel_reduce_test(101, 2);
}

// Specifying a simd layout other than the default doesn't work on the GPU
#ifndef KOKKOS_ENABLE_CUDA
TEST_F(StkSimdViewFixture, SimdParallelFor2d_LayoutRight) {
  StkSimdView2dTester<double*[4], double, 4, stk::simd::LayoutRight<double> > tester(237);
  tester.const_cast_test();
  tester.host_mirror_test();
  tester.parallel_for_test();
  tester.parallel_reduce_test(144, 0);
  tester.parallel_reduce_test(7, 3);
  tester.parallel_reduce_test(64, 1);
  tester.parallel_reduce_test(43, 2);
}

TEST_F(StkSimdViewFixture, SimdParallelFor2d_LayoutLeft) {
  StkSimdView2dTester<double*[2], double, 2, stk::simd::LayoutLeft<double> > tester(154);
  tester.const_cast_test();
  tester.host_mirror_test();
  tester.parallel_for_test();
  tester.parallel_reduce_test(123, 0);
  tester.parallel_reduce_test(98, 0);
  tester.parallel_reduce_test(3, 1);
  tester.parallel_reduce_test(64, 1);
}
#endif

// Two runtime dimensions

TEST_F(StkSimdViewFixture, SimdParallelFor2d_DefaultLayout_PtrPtr) {
  StkSimdView2dTester<double**, double, 3> tester(237);
  tester.const_cast_test();
  tester.host_mirror_test();
  tester.parallel_for_test();
  tester.parallel_reduce_test(237, 0);
  tester.parallel_reduce_test(3, 0);
  tester.parallel_reduce_test(64, 1);
  tester.parallel_reduce_test(101, 2);
}

// Specifying a simd layout other than the default doesn't work on the GPU
#ifndef KOKKOS_ENABLE_CUDA
TEST_F(StkSimdViewFixture, SimdParallelFor2d_LayoutRight_PtrPtr) {
  StkSimdView2dTester<double**, double, 4, stk::simd::LayoutRight<double> > tester(237);
  tester.const_cast_test();
  tester.host_mirror_test();
  tester.parallel_for_test();
  tester.parallel_reduce_test(144, 0);
  tester.parallel_reduce_test(7, 3);
  tester.parallel_reduce_test(64, 1);
  tester.parallel_reduce_test(43, 2);
}

TEST_F(StkSimdViewFixture, SimdParallelFor2d_LayoutLeft_PtrPtr) {
  StkSimdView2dTester<double**, double, 2, stk::simd::LayoutLeft<double> > tester(154);
  tester.const_cast_test();
  tester.host_mirror_test();
  tester.parallel_for_test();
  tester.parallel_reduce_test(123, 0);
  tester.parallel_reduce_test(98, 0);
  tester.parallel_reduce_test(3, 1);
  tester.parallel_reduce_test(64, 1);
}
#endif

#endif
