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


template <typename ArrayStyle, typename Real, typename Layout, int dim>
struct ViewMaker {
  stk::simd::View<ArrayStyle, Layout> make_view(int N, int M, int L) {
    return stk::simd::View<ArrayStyle, Layout>("name", N);
  }
};

template <typename Real, typename Layout, int dim>
struct ViewMaker<Real**[dim], Real, Layout, dim> {
  stk::simd::View<Real**[dim], Layout> make_view(int N, int M, int L) {
    return stk::simd::View<Real**[dim], Layout>("name", N, M);
  }
};

template <typename Real, typename Layout, int dim>
struct ViewMaker<Real***, Real, Layout, dim> {
  stk::simd::View<Real***, Layout> make_view(int N, int M, int L) {
    return stk::simd::View<Real***, Layout>("name", N, M, L);
  }
};


// need to allocate the view outside the class for some reason
template <typename ArrayStyle, typename Real, int dim1, int dim2, typename Layout>
stk::simd::View<ArrayStyle,Layout> get_3d_simd_view(int N, int scaling1, int scaling2) {
  ViewMaker<ArrayStyle, Real, Layout, dim2> maker;
  stk::simd::View<ArrayStyle, Layout> a = maker.make_view(N, dim1, dim2);
  Kokkos::parallel_for("3d_view_fill", N, KOKKOS_LAMBDA(const int i) {
    for (int j=0; j < dim1; ++j) {
      for (int k=0; k < dim2; ++k) {
        a(i,j,k) = i+scaling1*j+scaling2*k;
      }
    }
  });
  return a;
}

template <typename ArrayStyle, typename Real, int dim1, int dim2, typename Layout=void>
class StkSimdView3dTester {
 public:

  StkSimdView3dTester(int firstDim) : N(firstDim), scaling1(100), scaling2(20000) {}

  typedef stk::simd::DeviceIndex DeviceIndex;
  typedef typename stk::simd::DeviceTraits<Real>::simd_type DeviceReal;

  typedef stk::simd::View<ArrayStyle, Layout> SimdView;
  typedef stk::simd::View<const Real*[dim1][dim2], Layout> ConstSimdView;
  typedef typename stk::simd::View<ArrayStyle, Layout>::HostMirror HostSimdView;
  typedef typename stk::simd::View<const Real*[dim1][dim2], Layout>::HostMirror ConstHostSimdView;

  void const_cast_test() const {
    SimdView aSimd = get_3d_simd_view<ArrayStyle, Real, dim1, dim2, Layout>(N, scaling1, scaling2);
    ConstSimdView bSimd = aSimd;
  }

  void host_mirror_test() const {
    const SimdView a = get_3d_simd_view<ArrayStyle, Real, dim1, dim2, Layout>(N, scaling1, scaling2);
    auto aHost = stk::simd::create_mirror_view(a);
    stk::simd::deep_copy(aHost, a);
    for (int i=0; i < N; ++i) {
      for (int j=0; j < dim1; ++j) {
        for (int k=0; k < dim2; ++k) {
          EXPECT_EQ(i+scaling1*j+scaling2*k, aHost(i,j,k));
        }
      }
    }
  }

  void parallel_for_test() const {
    SimdView a = get_3d_simd_view<ArrayStyle, Real, dim1, dim2, Layout>(N, scaling1, scaling2);
    SimdView b = get_3d_simd_view<ArrayStyle, Real, dim1, dim2, Layout>(N, scaling1, scaling2);
    
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
  
  void parallel_reduce_test(int loopSize, int j, int k) const {
    ASSERT_LE(loopSize, N);
    ASSERT_LT(j, dim1);
    ASSERT_LT(k, dim2);
    SimdView a = get_3d_simd_view<ArrayStyle, Real, dim1, dim2, Layout>(N, scaling1, scaling2);

    int expectedSum = running_sum(loopSize-1) + loopSize*scaling1*j + loopSize*scaling2*k;
    EXPECT_EQ( expectedSum, reduce_sum_lambda(a, loopSize, j, k) ); 
    EXPECT_EQ( expectedSum, reduce_sum_functor(a, loopSize, j, k) );
    EXPECT_EQ( expectedSum, reduce_sum_functor_with_tag(a, loopSize, j, k) );
    auto aHost = stk::simd::copy_from_device(a);
    EXPECT_EQ( expectedSum, reduce_sum_each_lambda(aHost, loopSize, j, k) );
  }
    
  struct SomeTag {};
    
  struct AddFunctor {
    AddFunctor(ConstSimdView a_, SimdView b_) : a(a_), b(b_) {}

    STK_INLINE
    void operator() (const DeviceIndex& i) const {
      for (int j=0; j < dim1; ++j) {
        for (int k=0; k < dim2; ++k) {
          DeviceReal tmp = a(i,j,k)+b(i,j,k);
          b(i,j,k) = tmp;
        }
      }
    }

    STK_INLINE
    void operator() (SomeTag, const DeviceIndex& i) const {
      for (int j=0; j < dim1; ++j) {
        for (int k=0; k < dim2; ++k) {
          DeviceReal tmp = 2*a(i,j,k)+b(i,j,k);
          b(i,j,k) = tmp;
        }
      }
    }

   private:
    ConstSimdView a;
    SimdView b;
  };

  void set_b_to_a_plus_b_lambda(ConstSimdView a, SimdView b) const {
    stk::simd::parallel_for<Real>("b=a+b", N, STK_LAMBDA(const DeviceIndex& i) {
      for (int j=0; j < dim1; ++j) {
        for (int k=0; k < dim2; ++k) {
          DeviceReal tmp = a(i,j,k)+b(i,j,k);
          b(i,j,k) = tmp;
        }
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
      for (int j=0; j < dim1; ++j) {
        for (int k=0; k < dim2; ++k) {
          b(i,j,k) = 3*a(i,j,k)+b(i,j,k);
        }
      }
    });
  }
    
  struct ReduceSumFunctor {
    ReduceSumFunctor(SimdView a_, int j_, int k_) : a(a_), j(j_), k(k_) {}
    STK_INLINE void operator() (const DeviceIndex& i, DeviceReal& v) const {
      v += a(i,j,k);
    }
   private:
    SimdView a;
    int j,k;
  };

  struct ReduceSumFunctorWithTag {
    ReduceSumFunctorWithTag(SimdView a_, int j_, int k_) : a(a_), j(j_), k(k_) {}
    STK_INLINE void operator() (SomeTag, const DeviceIndex& i, DeviceReal& v) const {
      v += a(i,j,k);
    }
   private:
    SimdView a;
    int j,k;
  };

  Real reduce_sum_lambda(SimdView a, int loopSize, int j, int k) const {
    Real reducedRunningSum=0;
    stk::simd::parallel_reduce_sum("running_sum", loopSize, STK_LAMBDA(const DeviceIndex& i, DeviceReal& v) {
      v += a(i,j,k);
    }, reducedRunningSum);
    return reducedRunningSum;
  }

  Real reduce_sum_functor(SimdView a, int loopSize, int j, int k) const {
    Real reducedRunningSum=0;
    stk::simd::parallel_reduce_sum("running_sum", loopSize,
                                   ReduceSumFunctor(a,j,k), reducedRunningSum);
    return reducedRunningSum;
  }

  Real reduce_sum_functor_with_tag(SimdView a, int loopSize, int j, int k) const {
    Real reducedRunningSum=0;
    stk::simd::parallel_reduce_sum("running_sum", 
                                   Kokkos::RangePolicy<SomeTag>(0, loopSize),
                                   ReduceSumFunctorWithTag(a,j,k),
                                   reducedRunningSum);
    return reducedRunningSum;
  }

  Real reduce_sum_each_lambda(HostSimdView a, int loopSize, int j, int k) const {
    return stk::simd::reduce_sum_each(loopSize, [&](const DeviceIndex& i) {
      return a(i,j,k);
    });
  }

 private:

  void test_view_equal_to_index_multiple(HostSimdView view, int n) const {
    for (int i=0; i < N; ++i) {
      for (int j=0; j < dim1; ++j) {
        for (int k=0; k < dim2; ++k) {
          EXPECT_EQ( n*(i+scaling1*j+scaling2*k), view(i,j,k) );
        }
      }
    }    
  }

  int N;
  int scaling1;
  int scaling2;
};

TEST_F(StkSimdViewFixture, SimdParallelFor3d_DefaultLayout) {
  StkSimdView3dTester<double*[3][4], double, 3, 4> tester(237);
  tester.const_cast_test();
  tester.host_mirror_test();
  tester.parallel_for_test();
  tester.parallel_reduce_test(237, 0, 1);
  tester.parallel_reduce_test(3, 0, 3);
  tester.parallel_reduce_test(64, 1, 0);
  tester.parallel_reduce_test(101, 2, 2);
}

// Specifying a simd layout other than the default doesn't work on the GPU
#ifndef KOKKOS_ENABLE_CUDA
TEST_F(StkSimdViewFixture, SimdParallelFor3d_LayoutRight) {
  StkSimdView3dTester<double*[4][2], double, 4, 2, stk::simd::LayoutRight<double> > tester(237);
  tester.const_cast_test();
  tester.host_mirror_test();
  tester.parallel_for_test();
  tester.parallel_reduce_test(144, 0, 1);
  tester.parallel_reduce_test(7, 3, 0);
  tester.parallel_reduce_test(64, 1, 1);
  tester.parallel_reduce_test(43, 2, 0);
}

TEST_F(StkSimdViewFixture, SimdParallelFor3d_LayoutLeft) {
  StkSimdView3dTester<double*[2][4], double, 2, 4, stk::simd::LayoutLeft<double> > tester(154);
  tester.const_cast_test();
  tester.host_mirror_test();
  tester.parallel_for_test();
  tester.parallel_reduce_test(123, 0, 1);
  tester.parallel_reduce_test(98, 0, 0);
  tester.parallel_reduce_test(3, 1, 3);
  tester.parallel_reduce_test(64, 1, 2);
}
#endif

// Two runtime dimensions

TEST_F(StkSimdViewFixture, SimdParallelFor3d_DefaultLayout_PtrPtr) {
  StkSimdView3dTester<double**[4], double, 3, 4> tester(237);
  tester.const_cast_test();
  tester.host_mirror_test();
  tester.parallel_for_test();
  tester.parallel_reduce_test(237, 0, 1);
  tester.parallel_reduce_test(3, 0, 3);
  tester.parallel_reduce_test(64, 1, 0);
  tester.parallel_reduce_test(101, 2, 2);
}

// Specifying a simd layout other than the default doesn't work on the GPU
#ifndef KOKKOS_ENABLE_CUDA
TEST_F(StkSimdViewFixture, SimdParallelFor3d_LayoutRight_PtrPtr) {
  StkSimdView3dTester<double**[2], double, 4, 2, stk::simd::LayoutRight<double> > tester(237);
  tester.const_cast_test();
  tester.host_mirror_test();
  tester.parallel_for_test();
  tester.parallel_reduce_test(144, 0, 1);
  tester.parallel_reduce_test(7, 3, 0);
  tester.parallel_reduce_test(64, 1, 1);
  tester.parallel_reduce_test(43, 2, 0);
}

TEST_F(StkSimdViewFixture, SimdParallelFor3d_LayoutLeft_PtrPtr) {
  StkSimdView3dTester<double**[4], double, 2, 4, stk::simd::LayoutLeft<double> > tester(154);
  tester.const_cast_test();
  tester.host_mirror_test();
  tester.parallel_for_test();
  tester.parallel_reduce_test(123, 0, 1);
  tester.parallel_reduce_test(98, 0, 0);
  tester.parallel_reduce_test(3, 1, 3);
  tester.parallel_reduce_test(64, 1, 2);
}
#endif

// Three runtime dimensions

TEST_F(StkSimdViewFixture, SimdParallelFor3d_DefaultLayout_PtrPtrPtr) {
  StkSimdView3dTester<double***, double, 3, 4> tester(237);
  tester.const_cast_test();
  tester.host_mirror_test();
  tester.parallel_for_test();
  tester.parallel_reduce_test(237, 0, 1);
  tester.parallel_reduce_test(3, 0, 3);
  tester.parallel_reduce_test(64, 1, 0);
  tester.parallel_reduce_test(101, 2, 2);
}

// Specifying a simd layout other than the default doesn't work on the GPU
#ifndef KOKKOS_ENABLE_CUDA
TEST_F(StkSimdViewFixture, SimdParallelFor3d_LayoutRight_PtrPtrPtr) {
  StkSimdView3dTester<double***, double, 4, 2, stk::simd::LayoutRight<double> > tester(237);
  tester.const_cast_test();
  tester.host_mirror_test();
  tester.parallel_for_test();
  tester.parallel_reduce_test(144, 0, 1);
  tester.parallel_reduce_test(7, 3, 0);
  tester.parallel_reduce_test(64, 1, 1);
  tester.parallel_reduce_test(43, 2, 0);
}

TEST_F(StkSimdViewFixture, SimdParallelFor3d_LayoutLeft_PtrPtrPtr) {
  StkSimdView3dTester<double***, double, 2, 4, stk::simd::LayoutLeft<double> > tester(154);
  tester.const_cast_test();
  tester.host_mirror_test();
  tester.parallel_for_test();
  tester.parallel_reduce_test(123, 0, 1);
  tester.parallel_reduce_test(98, 0, 0);
  tester.parallel_reduce_test(3, 1, 3);
  tester.parallel_reduce_test(64, 1, 2);
}
#endif

#endif
