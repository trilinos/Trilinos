//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

// **********************************************************************
// The tests executed by the code below cover many combinations for
// the operation y += a * x + b * y:
// 01) Type of 'x' and 'a' components: float, double, complex, ...
// 02) Type of 'y' and 'b' components: float, double, complex, ...
// 03) Execution space: serial, threads, OpenMP, Cuda, ...
// 04) Layout of 'x' and 'a'
// 05) Layout of 'y' and 'b'
// 06) Ranks of 'x' and 'y': rank-1 or rank-2
// 07) Ranks of 'a' and 'b': scalars or rank-0 or rank-1
//
// Choices (01)-(03) are selected in the routines TEST_F() at the very
// bottom of the file, when calling:
// - either test_axpby_unification<...>(),
// - or test_axpby_mv_unification<...>().
//
// Choices (04)-(05) are selected in routines:
// - test_axpby_unification<...>(), when calling
//   Test::impl_test_axpby_unification<...>(), and
// - test_axpby_mv_unification<...>(), when calling
//   Test::impl_test_axpby_mv_unification<...>().
//
// Choices (06)-(07) are selected in routines:
// - Test::impl_test_axpby_unification<...>(), through
//   16 different combinations and calls to
//   Test::impl_test_axpby_unification_compare<...>(), and
// - Test::impl_test_axpby_mv_unification<...>(), through
//   36 different combinations and calls to
//   Test::impl_test_axpby_mv_unification_compare<...>().
//
// The constexpr integer value 15 for 'numVecsAxpbyTest' was chosen to
// force the test of the three unrolling values 8, 4, and 1, in routine
// Axpby_MV_Invoke_Left<...>(...) in file KokkosBlas1_axpby_mv_impl.hpp
// **********************************************************************

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <KokkosBlas1_axpby.hpp>
#include <KokkosKernels_TestUtils.hpp>

static constexpr int numVecsAxpbyTest = 15;

namespace Test {

template <typename T, bool Enable = false>
struct getScalarTypeFromT {
  using type = T;
};

template <typename T>
struct getScalarTypeFromT<T, true> {
  using type = typename T::value_type;
};

template <class T>
constexpr bool isRank0() {
  if constexpr (Kokkos::is_view_v<T>) {
    return (T::rank == 0);
  }
  return false;
}

template <class tScalarA, class tA, class tX, class tScalarB, class tB, class tY, class Device>
void impl_test_axpby_unification_compare(tA const& a, tX const& x, tB const& b, tY const& y, int N, bool testWithNanY,
                                         typename Kokkos::ArithTraits<tScalarB>::mag_type const max_val,
                                         typename Kokkos::ArithTraits<tScalarB>::mag_type const max_error,
                                         tScalarA const inputValueA = Kokkos::ArithTraits<tScalarA>::zero(),
                                         tScalarB const inputValueB = Kokkos::ArithTraits<tScalarB>::zero()) {
  using ScalarTypeX = typename std::remove_const<typename tX::DView::value_type>::type;
  using ScalarTypeY = typename std::remove_const<typename tY::DView::value_type>::type;

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);

  {
    ScalarTypeX randStart, randEnd;
    Test::getRandomBounds(max_val, randStart, randEnd);
    Kokkos::fill_random(x.d_view, rand_pool, randStart, randEnd);
  }
  Kokkos::deep_copy(x.h_base, x.d_base);

  {
    ScalarTypeY randStart, randEnd;
    Test::getRandomBounds(max_val, randStart, randEnd);
    if (testWithNanY) {
      Kokkos::deep_copy(y.d_view, Kokkos::ArithTraits<ScalarTypeY>::nan());
    } else {
      Kokkos::fill_random(y.d_view, rand_pool, randStart, randEnd);
    }
  }
  tY org_y("Org_Y", N);
  Kokkos::deep_copy(org_y.h_base, y.d_base);

  tScalarA valueA(Kokkos::ArithTraits<tScalarA>::zero());
  tScalarB valueB(Kokkos::ArithTraits<tScalarB>::zero());

  if constexpr (std::is_same_v<tA, tScalarA>) {
    valueA = a;
    if constexpr (std::is_same_v<tB, tScalarB>) {
      valueB = b;
      KokkosBlas::axpby(a, x.d_view, b, y.d_view);
    } else if constexpr (isRank0<tB>()) {
      if constexpr (std::is_same_v<typename tB::array_layout, Kokkos::LayoutStride>) {
        valueB = inputValueB;
      } else {
        typename tB::HostMirror h_b("h_B");
        Kokkos::deep_copy(h_b, b);
        valueB = h_b();
      }
      KokkosBlas::axpby(a, x.d_view, b, y.d_view);
    } else {
      Kokkos::deep_copy(b.h_base, b.d_base);
      valueB = b.h_view(0);
      KokkosBlas::axpby(a, x.d_view, b.d_view, y.d_view);
    }
  } else if constexpr (isRank0<tA>()) {
    if constexpr (std::is_same_v<typename tA::array_layout, Kokkos::LayoutStride>) {
      valueA = inputValueA;
    } else {
      typename tA::HostMirror h_a("h_A");
      Kokkos::deep_copy(h_a, a);
      valueA = h_a();
    }
    if constexpr (std::is_same_v<tB, tScalarB>) {
      valueB = b;
      KokkosBlas::axpby(a, x.d_view, b, y.d_view);
    } else if constexpr (isRank0<tB>()) {
      if constexpr (std::is_same_v<typename tB::array_layout, Kokkos::LayoutStride>) {
        valueB = inputValueB;
      } else {
        typename tB::HostMirror h_b("h_B");
        Kokkos::deep_copy(h_b, b);
        valueB = h_b();
      }
      KokkosBlas::axpby(a, x.d_view, b, y.d_view);
    } else {
      Kokkos::deep_copy(b.h_base, b.d_base);
      valueB = b.h_view(0);
      KokkosBlas::axpby(a, x.d_view, b.d_view, y.d_view);
    }
  } else {
    Kokkos::deep_copy(a.h_base, a.d_base);
    valueA = a.h_view(0);
    if constexpr (std::is_same_v<tB, tScalarB>) {
      valueB = b;
      KokkosBlas::axpby(a.d_view, x.d_view, b, y.d_view);
    } else if constexpr (isRank0<tB>()) {
      if constexpr (std::is_same_v<typename tB::array_layout, Kokkos::LayoutStride>) {
        valueB = inputValueB;
      } else {
        typename tB::HostMirror h_b("h_B");
        Kokkos::deep_copy(h_b, b);
        valueB = h_b();
      }
      KokkosBlas::axpby(a.d_view, x.d_view, b, y.d_view);
    } else {
      Kokkos::deep_copy(b.h_base, b.d_base);
      valueB = b.h_view(0);
      KokkosBlas::axpby(a.d_view, x.d_view, b.d_view, y.d_view);
    }
  }

  Kokkos::deep_copy(y.h_base, y.d_base);

  if (testWithNanY == false) {
    for (int i(0); i < N; ++i) {
      EXPECT_NEAR_KK(static_cast<ScalarTypeY>(valueA * x.h_view(i) + valueB * org_y.h_view(i)), y.h_view(i),
                     4. * max_error);
    }
  } else {
    // ********************************************************
    // Tests with 'Y == nan()' are called only for cases where
    // b == Kokkos::ArithTraits<tScalarB>::zero()
    // ********************************************************
    for (int i(0); i < N; ++i) {
#if 0
      ScalarTypeY tmp = static_cast<ScalarTypeY>(valueA * x.h_view(i) + valueB * org_y.h_view(i));
      std::cout << "i = "                 << i
                << ", valueA = "          << valueA
                << ", x.h_view(i) = "     << x.h_view(i)
                << ", valueB = "          << valueB
                << ", org_y.h_view(i) = " << org_y.h_view(i)
                << ", tmp = "             << tmp
                << ", y.h_view(i) = "     << y.h_view(i)
                << std::endl;
#endif
      if constexpr (std::is_same_v<ScalarTypeY, int>) {
        // ****************************************************************
        // 'nan()' converts to '-1' in case of 'int' => no need to compare
        // ****************************************************************
        if (y.h_view(i) != -1) {
          EXPECT_NE(y.h_view(i), Kokkos::ArithTraits<ScalarTypeY>::nan());
        }
      } else {
        EXPECT_NE(y.h_view(i), Kokkos::ArithTraits<ScalarTypeY>::nan());
      }
      EXPECT_NEAR_KK(static_cast<ScalarTypeY>(valueA * x.h_view(i)), y.h_view(i), 4. * max_error);
    }
  }
}

template <class tScalarA, class tA, class tX, class tScalarB, class tB, class tY, class Device>
void impl_test_axpby_mv_unification_compare(tA const& a, tX const& x, tB const& b, tY const& y, int N, int K,
                                            bool testWithNanY,
                                            typename Kokkos::ArithTraits<tScalarB>::mag_type const max_val,
                                            typename Kokkos::ArithTraits<tScalarB>::mag_type const max_error,
                                            tScalarA const inputValueA = Kokkos::ArithTraits<tScalarA>::zero(),
                                            tScalarB const inputValueB = Kokkos::ArithTraits<tScalarB>::zero()) {
  using ScalarTypeX = typename std::remove_const<typename tY::DView::value_type>::type;
  using ScalarTypeY = typename std::remove_const<typename tY::DView::value_type>::type;

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);

  {
    ScalarTypeX randStart, randEnd;
    Test::getRandomBounds(max_val, randStart, randEnd);
    Kokkos::fill_random(x.d_view, rand_pool, randStart, randEnd);
  }
  Kokkos::deep_copy(x.h_base, x.d_base);

  {
    ScalarTypeY randStart, randEnd;
    Test::getRandomBounds(max_val, randStart, randEnd);
    if (testWithNanY) {
      Kokkos::deep_copy(y.d_view, Kokkos::ArithTraits<ScalarTypeY>::nan());
    } else {
      Kokkos::fill_random(y.d_view, rand_pool, randStart, randEnd);
    }
  }
  tY org_y("Org_Y", N, K);
  Kokkos::deep_copy(org_y.h_base, y.d_base);

  // Cannot use "if constexpr (isRank1<tA>()) {" because rank-1 variables
  // are passed to current routine with view_stride_adapter<...>
  bool constexpr aIsRank1 = !std::is_same_v<tA, tScalarA> && !isRank0<tA>();
  if constexpr (aIsRank1) {
    Kokkos::deep_copy(a.h_base, a.d_base);
  }

  // Cannot use "if constexpr (isRank1<tB>()) {" because rank-1 variables
  // are passed to current routine with view_stride_adapter<...>
  bool constexpr bIsRank1 = !std::is_same_v<tB, tScalarB> && !isRank0<tB>();
  if constexpr (bIsRank1) {
    Kokkos::deep_copy(b.h_base, b.d_base);
  }

  tScalarA valueA(Kokkos::ArithTraits<tScalarA>::zero());
  tScalarB valueB(Kokkos::ArithTraits<tScalarB>::zero());
  if constexpr (std::is_same_v<tA, tScalarA>) {
    valueA = a;
    if constexpr (std::is_same_v<tB, tScalarB>) {
      valueB = b;
      KokkosBlas::axpby(a, x.d_view, b, y.d_view);
    } else if constexpr (isRank0<tB>()) {
      if constexpr (std::is_same_v<typename tB::array_layout, Kokkos::LayoutStride>) {
        valueB = inputValueB;
      } else {
        typename tB::HostMirror h_b("h_B");
        Kokkos::deep_copy(h_b, b);
        valueB = h_b();
      }
      KokkosBlas::axpby(a, x.d_view, b, y.d_view);
    } else {
      valueB = b.h_view(0);
      KokkosBlas::axpby(a, x.d_view, b.d_view, y.d_view);
    }
  } else if constexpr (isRank0<tA>()) {
    if constexpr (std::is_same_v<typename tA::array_layout, Kokkos::LayoutStride>) {
      valueA = inputValueA;
    } else {
      typename tA::HostMirror h_a("h_A");
      Kokkos::deep_copy(h_a, a);
      valueA = h_a();
    }
    if constexpr (std::is_same_v<tB, tScalarB>) {
      valueB = b;
      KokkosBlas::axpby(a, x.d_view, b, y.d_view);
    } else if constexpr (isRank0<tB>()) {
      if constexpr (std::is_same_v<typename tB::array_layout, Kokkos::LayoutStride>) {
        valueB = inputValueB;
      } else {
        typename tB::HostMirror h_b("h_B");
        Kokkos::deep_copy(h_b, b);
        valueB = h_b();
      }
      KokkosBlas::axpby(a, x.d_view, b, y.d_view);
    } else {
      valueB = b.h_view(0);
      KokkosBlas::axpby(a, x.d_view, b.d_view, y.d_view);
    }
  } else {
    valueA = a.h_view(0);
    if constexpr (std::is_same_v<tB, tScalarB>) {
      valueB = b;
      KokkosBlas::axpby(a.d_view, x.d_view, b, y.d_view);
    } else if constexpr (isRank0<tB>()) {
      if constexpr (std::is_same_v<typename tB::array_layout, Kokkos::LayoutStride>) {
        valueB = inputValueB;
      } else {
        typename tB::HostMirror h_b("h_B");
        Kokkos::deep_copy(h_b, b);
        valueB = h_b();
      }
      KokkosBlas::axpby(a.d_view, x.d_view, b, y.d_view);
    } else {
      valueB = b.h_view(0);
      KokkosBlas::axpby(a.d_view, x.d_view, b.d_view, y.d_view);
    }
  }

  Kokkos::deep_copy(y.h_base, y.d_base);

  if (testWithNanY == false) {
    for (int i(0); i < N; ++i) {
      for (int k(0); k < K; ++k) {
        ScalarTypeY vanillaValue(Kokkos::ArithTraits<ScalarTypeY>::zero());
        if constexpr (aIsRank1) {
          (void)valueA;  // Avoid "set but not used" error
          if constexpr (bIsRank1) {
            (void)valueB;  // Avoid "set but not used" error
            int a_k(a.h_view.extent(0) == 1 ? 0 : k);
            int b_k(b.h_view.extent(0) == 1 ? 0 : k);
#if 0
            std::cout << "In impl_test_axpby_mv_unification_compare()"
                      << ": i = " << i
                      << ", k = " << k
                      << ", a.h_view.extent(0) = " << a.h_view.extent(0)
                      << ", a_k = "                << a_k
                      << ", b.h_view.extent(0) = " << b.h_view.extent(0)
                      << ", b_k = "                << b_k
                      << ", a.h_view(a_k) = "      << a.h_view(a_k)
                      << ", x.h_view(i, k) = "     << x.h_view(i, k)
                      << ", b.h_view(b_k) = "      << b.h_view(b_k)
                      << ", org_y.h_view(i, k) = " << org_y.h_view(i, k)
                      << std::endl;
#endif
            vanillaValue =
                static_cast<ScalarTypeY>(a.h_view(a_k) * x.h_view(i, k) + b.h_view(b_k) * org_y.h_view(i, k));
          } else {
            int a_k(a.h_view.extent(0) == 1 ? 0 : k);
            vanillaValue = static_cast<ScalarTypeY>(a.h_view(a_k) * x.h_view(i, k) + valueB * org_y.h_view(i, k));
          }
        } else {
          if constexpr (bIsRank1) {
            (void)valueB;  // Avoid "set but not used" error
            int b_k(b.h_view.extent(0) == 1 ? 0 : k);
            vanillaValue = static_cast<ScalarTypeY>(valueA * x.h_view(i, k) + b.h_view(b_k) * org_y.h_view(i, k));
          } else {
            vanillaValue = static_cast<ScalarTypeY>(valueA * x.h_view(i, k) + valueB * org_y.h_view(i, k));
          }
        }
#if 0
        std::cout << "In impl_test_axpby_mv_unification_compare(1)"
                  << ": i = " << i
                  << ", k = " << k
                  << ", y.h_view(i, k) = " << y.h_view(i, k)
                  << ", vanillaValue = "   << vanillaValue
                  << std::endl;
#endif
        EXPECT_NEAR_KK(vanillaValue, y.h_view(i, k), 4. * max_error);
      }
    }
  } else {
    // ********************************************************
    // Tests with 'Y == nan()' are called only for cases where
    // b == Kokkos::ArithTraits<tScalarB>::zero()
    // ********************************************************
    for (int i(0); i < N; ++i) {
      for (int k(0); k < K; ++k) {
        ScalarTypeY vanillaValue(Kokkos::ArithTraits<ScalarTypeY>::zero());
        if constexpr (aIsRank1) {
          (void)valueA;  // Avoid "set but not used" error
          int a_k(a.h_view.extent(0) == 1 ? 0 : k);
          vanillaValue = static_cast<ScalarTypeY>(a.h_view(a_k) * x.h_view(i, k));
#if 0
          ScalarTypeY tmp = static_cast<ScalarTypeY>(a.h_view(a_k) * x.h_view(i, k) + valueB * org_y.h_view(i, k));
          std::cout << "i = "                    << i
                    << ", k = "                  << k
                    << ", a_k = "                << a_k
                    << ", a.h_view(a_k) = "      << a.h_view(a_k)
                    << ", x.h_view(i, k) = "     << x.h_view(i, k)
                    << ", valueB = "             << valueB
                    << ", org_y.h_view(i, k) = " << org_y.h_view(i, k)
                    << ", tmp = "                << tmp
                    << ", vanillaValue = "       << vanillaValue
                    << ", y.h_view(i, k) = "     << y.h_view(i, k)
                    << std::endl;
#endif
        } else {
          vanillaValue = static_cast<ScalarTypeY>(valueA * x.h_view(i, k));
#if 0
          ScalarTypeY tmp = static_cast<ScalarTypeY>(valueA * x.h_view(i, k) + valueB * org_y.h_view(i, k));
          std::cout << "i = "                    << i
                    << ", k = "                  << k
                    << ", valueA = "             << valueA
                    << ", x.h_view(i, k) = "     << x.h_view(i, k)
                    << ", valueB = "             << valueB
                    << ", org_y.h_view(i, k) = " << org_y.h_view(i, k)
                    << ", tmp = "                << tmp
                    << ", vanillaValue = "       << vanillaValue
                    << ", y.h_view(i, k) = "     << y.h_view(i, k)
                    << std::endl;
#endif
        }

        if constexpr (std::is_same_v<ScalarTypeY, int>) {
          // ****************************************************************
          // 'nan()' converts to '-1' in case of 'int' => no need to compare
          // ****************************************************************
          if (y.h_view(i, k) != -1) {
            EXPECT_NE(y.h_view(i, k), Kokkos::ArithTraits<ScalarTypeY>::nan());
          }
        } else {
          EXPECT_NE(y.h_view(i, k), Kokkos::ArithTraits<ScalarTypeY>::nan());
        }
#if 0
        std::cout << "In impl_test_axpby_mv_unification_compare(2)"
                  << ": i = " << i
                  << ", k = " << k
                  << ", y.h_view(i, k) = " << y.h_view(i, k)
                  << ", vanillaValue = "   << vanillaValue
                  << std::endl;
#endif
        EXPECT_NEAR_KK(vanillaValue, y.h_view(i, k), 4. * max_error);
      }
    }
  }
}

template <class tScalarA, class tLayoutA, class tScalarX, class tLayoutX, class tScalarB, class tLayoutB,
          class tScalarY, class tLayoutY, class Device>
void impl_test_axpby_unification(int const N) {
  using ViewTypeAr0    = Kokkos::View<tScalarA, tLayoutA, Device>;
  using ViewTypeAr1s_1 = Kokkos::View<tScalarA[1], tLayoutA, Device>;
  using ViewTypeAr1d   = Kokkos::View<tScalarA*, tLayoutA, Device>;

  using ViewTypeX = Kokkos::View<tScalarX*, tLayoutX, Device>;

  using ViewTypeBr0    = Kokkos::View<tScalarB, tLayoutB, Device>;
  using ViewTypeBr1s_1 = Kokkos::View<tScalarB[1], tLayoutB, Device>;
  using ViewTypeBr1d   = Kokkos::View<tScalarB*, tLayoutB, Device>;

  using ViewTypeY = Kokkos::View<tScalarY*, tLayoutY, Device>;

  std::array<tScalarA, 4> const valuesA{-1, Kokkos::ArithTraits<tScalarA>::zero(), 1, 3};
  std::array<tScalarB, 4> const valuesB{-1, Kokkos::ArithTraits<tScalarB>::zero(), 1, 5};

  // eps should probably be based on tScalarB since that is the type
  // in which the result is computed.
  using MagnitudeB         = typename Kokkos::ArithTraits<tScalarB>::mag_type;
  MagnitudeB const eps     = Kokkos::ArithTraits<tScalarB>::epsilon();
  MagnitudeB const max_val = 10;
  MagnitudeB const max_error =
      static_cast<MagnitudeB>(Kokkos::ArithTraits<tScalarA>::abs(valuesA[valuesA.size() - 1]) +
                              Kokkos::ArithTraits<tScalarB>::abs(valuesB[valuesB.size() - 1])) *
      max_val * eps;

  // ************************************************************
  // Case 01/16: Ascalar + Bscalar
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 01/16" << std::endl;
#endif
  for (size_t i(0); i < valuesA.size(); ++i) {
    tScalarA const valueA(valuesA[i]);
    for (size_t j(0); j < valuesB.size(); ++j) {
      tScalarB const valueB(valuesB[j]);
      {
        tScalarA a;
        view_stride_adapter<ViewTypeX> x("X", N);
        tScalarB b;
        view_stride_adapter<ViewTypeY> y("Y", N);

        a = valueA;
        b = valueB;
        impl_test_axpby_unification_compare<tScalarA, tScalarA, view_stride_adapter<ViewTypeX>, tScalarB, tScalarB,
                                            view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, false, max_val,
                                                                                    max_error);
        if (valueB == Kokkos::ArithTraits<tScalarB>::zero()) {
          impl_test_axpby_unification_compare<tScalarA, tScalarA, view_stride_adapter<ViewTypeX>, tScalarB, tScalarB,
                                              view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, true, max_val,
                                                                                      max_error);
        }
      }
    }
  }

  // ************************************************************
  // Case 02/16: Ascalar + Br0
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 02/16" << std::endl;
#endif
  if constexpr (std::is_same_v<tLayoutB, Kokkos::LayoutStride>) {
    // Avoid the test, due to compilation errors
    // ViewTypeBr0 b;
    // Kokkos::deep_copy(b, valueB);
    // //std::cout << "b() = " << b() << std::endl;
  } else {
    for (size_t i(0); i < valuesA.size(); ++i) {
      tScalarA const valueA(valuesA[i]);
      for (size_t j(0); j < valuesB.size(); ++j) {
        tScalarB const valueB(valuesB[j]);
        {
          tScalarA a;
          view_stride_adapter<ViewTypeX> x("X", N);
          ViewTypeBr0 b("B");
          view_stride_adapter<ViewTypeY> y("Y", N);

          a = valueA;
          Kokkos::deep_copy(b, valueB);
          impl_test_axpby_unification_compare<tScalarA, tScalarA, view_stride_adapter<ViewTypeX>, tScalarB, ViewTypeBr0,
                                              view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, false, max_val,
                                                                                      max_error);
          if (valueB == Kokkos::ArithTraits<tScalarB>::zero()) {
            impl_test_axpby_unification_compare<tScalarA, tScalarA, view_stride_adapter<ViewTypeX>, tScalarB,
                                                ViewTypeBr0, view_stride_adapter<ViewTypeY>, Device>(
                a, x, b, y, N, true, max_val, max_error);
          }
        }
      }
    }
  }

  // ************************************************************
  // Case 03/16: Ascalar + Br1s_1
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 03/16" << std::endl;
#endif
  for (size_t i(0); i < valuesA.size(); ++i) {
    tScalarA const valueA(valuesA[i]);
    for (size_t j(0); j < valuesB.size(); ++j) {
      tScalarB const valueB(valuesB[j]);
      {
        tScalarA a;
        view_stride_adapter<ViewTypeX> x("X", N);
        view_stride_adapter<ViewTypeBr1s_1> b("B", 1);
        view_stride_adapter<ViewTypeY> y("Y", N);

        a = valueA;
        Kokkos::deep_copy(b.d_base, valueB);
        impl_test_axpby_unification_compare<tScalarA, tScalarA, view_stride_adapter<ViewTypeX>, tScalarB,
                                            view_stride_adapter<ViewTypeBr1s_1>, view_stride_adapter<ViewTypeY>,
                                            Device>(a, x, b, y, N, false, max_val, max_error);
        if (valueB == Kokkos::ArithTraits<tScalarB>::zero()) {
          impl_test_axpby_unification_compare<tScalarA, tScalarA, view_stride_adapter<ViewTypeX>, tScalarB,
                                              view_stride_adapter<ViewTypeBr1s_1>, view_stride_adapter<ViewTypeY>,
                                              Device>(a, x, b, y, N, true, max_val, max_error);
        }
      }
    }
  }

  // ************************************************************
  // Case 04/16: Ascalar + Br1d
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 04/16" << std::endl;
#endif
  for (size_t i(0); i < valuesA.size(); ++i) {
    tScalarA const valueA(valuesA[i]);
    for (size_t j(0); j < valuesB.size(); ++j) {
      tScalarB const valueB(valuesB[j]);
      {
        tScalarA a;
        view_stride_adapter<ViewTypeX> x("X", N);
        view_stride_adapter<ViewTypeBr1d> b("B", 1);
        view_stride_adapter<ViewTypeY> y("Y", N);

        a = valueA;
        Kokkos::deep_copy(b.d_base, valueB);
        impl_test_axpby_unification_compare<tScalarA, tScalarA, view_stride_adapter<ViewTypeX>, tScalarB,
                                            view_stride_adapter<ViewTypeBr1d>, view_stride_adapter<ViewTypeY>, Device>(
            a, x, b, y, N, false, max_val, max_error);
        if (valueB == Kokkos::ArithTraits<tScalarB>::zero()) {
          impl_test_axpby_unification_compare<tScalarA, tScalarA, view_stride_adapter<ViewTypeX>, tScalarB,
                                              view_stride_adapter<ViewTypeBr1d>, view_stride_adapter<ViewTypeY>,
                                              Device>(a, x, b, y, N, true, max_val, max_error);
        }
      }
    }
  }

  // ************************************************************
  // Case 05/16: Ar0 + Bscalar
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 05/16" << std::endl;
#endif
  if constexpr (std::is_same_v<tLayoutA, Kokkos::LayoutStride>) {
    // Avoid the test, due to compilation errors
  } else {
    for (size_t i(0); i < valuesA.size(); ++i) {
      tScalarA const valueA(valuesA[i]);
      for (size_t j(0); j < valuesB.size(); ++j) {
        tScalarB const valueB(valuesB[j]);
        {
          ViewTypeAr0 a("A");
          view_stride_adapter<ViewTypeX> x("X", N);
          tScalarB b;
          view_stride_adapter<ViewTypeY> y("Y", N);

          Kokkos::deep_copy(a, valueA);
          b = valueB;
          impl_test_axpby_unification_compare<tScalarA, ViewTypeAr0, view_stride_adapter<ViewTypeX>, tScalarB, tScalarB,
                                              view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, false, max_val,
                                                                                      max_error);
          if (valueB == Kokkos::ArithTraits<tScalarB>::zero()) {
            impl_test_axpby_unification_compare<tScalarA, ViewTypeAr0, view_stride_adapter<ViewTypeX>, tScalarB,
                                                tScalarB, view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, true,
                                                                                                  max_val, max_error);
          }
        }
      }
    }
  }

  // ************************************************************
  // Case 06/16: Ar0 + Br0
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 06/16" << std::endl;
#endif
  if constexpr ((std::is_same_v<tLayoutA, Kokkos::LayoutStride>) || (std::is_same_v<tLayoutB, Kokkos::LayoutStride>)) {
    // Avoid the test, due to compilation errors
  } else {
    for (size_t i(0); i < valuesA.size(); ++i) {
      tScalarA const valueA(valuesA[i]);
      for (size_t j(0); j < valuesB.size(); ++j) {
        tScalarB const valueB(valuesB[j]);
        {
          ViewTypeAr0 a("A");
          view_stride_adapter<ViewTypeX> x("X", N);
          ViewTypeBr0 b("B");
          view_stride_adapter<ViewTypeY> y("Y", N);

          Kokkos::deep_copy(a, valueA);
          Kokkos::deep_copy(b, valueB);
          impl_test_axpby_unification_compare<tScalarA, ViewTypeAr0, view_stride_adapter<ViewTypeX>, tScalarB,
                                              ViewTypeBr0, view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, false,
                                                                                                   max_val, max_error);
          if (valueB == Kokkos::ArithTraits<tScalarB>::zero()) {
            impl_test_axpby_unification_compare<tScalarA, ViewTypeAr0, view_stride_adapter<ViewTypeX>, tScalarB,
                                                ViewTypeBr0, view_stride_adapter<ViewTypeY>, Device>(
                a, x, b, y, N, true, max_val, max_error);
          }
        }
      }
    }
  }

  // ************************************************************
  // Case 07/16: Ar0 + Br1s_1
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 07/16" << std::endl;
#endif
  if constexpr (std::is_same_v<tLayoutA, Kokkos::LayoutStride>) {
    // Avoid the test, due to compilation errors
  } else {
    for (size_t i(0); i < valuesA.size(); ++i) {
      tScalarA const valueA(valuesA[i]);
      for (size_t j(0); j < valuesB.size(); ++j) {
        tScalarB const valueB(valuesB[j]);
        {
          ViewTypeAr0 a("A");
          view_stride_adapter<ViewTypeX> x("X", N);
          view_stride_adapter<ViewTypeBr1s_1> b("B", 1);
          view_stride_adapter<ViewTypeY> y("Y", N);

          Kokkos::deep_copy(a, valueA);
          Kokkos::deep_copy(b.d_base, valueB);
          impl_test_axpby_unification_compare<tScalarA, ViewTypeAr0, view_stride_adapter<ViewTypeX>, tScalarB,
                                              view_stride_adapter<ViewTypeBr1s_1>, view_stride_adapter<ViewTypeY>,
                                              Device>(a, x, b, y, N, false, max_val, max_error);
          if (valueB == Kokkos::ArithTraits<tScalarB>::zero()) {
            impl_test_axpby_unification_compare<tScalarA, ViewTypeAr0, view_stride_adapter<ViewTypeX>, tScalarB,
                                                view_stride_adapter<ViewTypeBr1s_1>, view_stride_adapter<ViewTypeY>,
                                                Device>(a, x, b, y, N, true, max_val, max_error);
          }
        }
      }
    }
  }

  // ************************************************************
  // Case 08/16: Ar0 + Br1d
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 08/16" << std::endl;
#endif
  if constexpr (std::is_same_v<tLayoutA, Kokkos::LayoutStride>) {
    // Avoid the test, due to compilation errors
  } else {
    for (size_t i(0); i < valuesA.size(); ++i) {
      tScalarA const valueA(valuesA[i]);
      for (size_t j(0); j < valuesB.size(); ++j) {
        tScalarB const valueB(valuesB[j]);
        {
          ViewTypeAr0 a("A");
          view_stride_adapter<ViewTypeX> x("X", N);
          view_stride_adapter<ViewTypeBr1d> b("B", 1);
          view_stride_adapter<ViewTypeY> y("Y", N);

          Kokkos::deep_copy(a, valueA);
          Kokkos::deep_copy(b.d_base, valueB);
          impl_test_axpby_unification_compare<tScalarA, ViewTypeAr0, view_stride_adapter<ViewTypeX>, tScalarB,
                                              view_stride_adapter<ViewTypeBr1d>, view_stride_adapter<ViewTypeY>,
                                              Device>(a, x, b, y, N, false, max_val, max_error);
          if (valueB == Kokkos::ArithTraits<tScalarB>::zero()) {
            impl_test_axpby_unification_compare<tScalarA, ViewTypeAr0, view_stride_adapter<ViewTypeX>, tScalarB,
                                                view_stride_adapter<ViewTypeBr1d>, view_stride_adapter<ViewTypeY>,
                                                Device>(a, x, b, y, N, true, max_val, max_error);
          }
        }
      }
    }
  }

  // ************************************************************
  // Case 09/16: Ar1s_1 + Bscalar
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 09/16" << std::endl;
#endif
  for (size_t i(0); i < valuesA.size(); ++i) {
    tScalarA const valueA(valuesA[i]);
    for (size_t j(0); j < valuesB.size(); ++j) {
      tScalarB const valueB(valuesB[j]);
      {
        view_stride_adapter<ViewTypeAr1s_1> a("A", 1);
        view_stride_adapter<ViewTypeX> x("X", N);
        tScalarB b;
        view_stride_adapter<ViewTypeY> y("Y", N);

        Kokkos::deep_copy(a.d_base, valueA);
        b = valueB;
        impl_test_axpby_unification_compare<tScalarA, view_stride_adapter<ViewTypeAr1s_1>,
                                            view_stride_adapter<ViewTypeX>, tScalarB, tScalarB,
                                            view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, false, max_val,
                                                                                    max_error);
        if (valueB == Kokkos::ArithTraits<tScalarB>::zero()) {
          impl_test_axpby_unification_compare<tScalarA, view_stride_adapter<ViewTypeAr1s_1>,
                                              view_stride_adapter<ViewTypeX>, tScalarB, tScalarB,
                                              view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, true, max_val,
                                                                                      max_error);
        }
      }
    }
  }

  // ************************************************************
  // Case 10/16: Ar1s_1 + Br0
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 10/16" << std::endl;
#endif
  if constexpr (std::is_same_v<tLayoutB, Kokkos::LayoutStride>) {
    // Avoid the test, due to compilation errors
  } else {
    for (size_t i(0); i < valuesA.size(); ++i) {
      tScalarA const valueA(valuesA[i]);
      for (size_t j(0); j < valuesB.size(); ++j) {
        tScalarB const valueB(valuesB[j]);
        {
          view_stride_adapter<ViewTypeAr1s_1> a("A", 1);
          view_stride_adapter<ViewTypeX> x("X", N);
          ViewTypeBr0 b("B");
          view_stride_adapter<ViewTypeY> y("Y", N);

          Kokkos::deep_copy(a.d_base, valueA);
          Kokkos::deep_copy(b, valueB);
          impl_test_axpby_unification_compare<tScalarA, view_stride_adapter<ViewTypeAr1s_1>,
                                              view_stride_adapter<ViewTypeX>, tScalarB, ViewTypeBr0,
                                              view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, false, max_val,
                                                                                      max_error);
          if (valueB == Kokkos::ArithTraits<tScalarB>::zero()) {
            impl_test_axpby_unification_compare<tScalarA, view_stride_adapter<ViewTypeAr1s_1>,
                                                view_stride_adapter<ViewTypeX>, tScalarB, ViewTypeBr0,
                                                view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, true, max_val,
                                                                                        max_error);
          }
        }
      }
    }
  }

  // ************************************************************
  // Case 11/16: Ar1s_1 + Br1s_1
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 11/16" << std::endl;
#endif
  for (size_t i(0); i < valuesA.size(); ++i) {
    tScalarA const valueA(valuesA[i]);
    for (size_t j(0); j < valuesB.size(); ++j) {
      tScalarB const valueB(valuesB[j]);
      {
        view_stride_adapter<ViewTypeAr1s_1> a("A", 1);
        view_stride_adapter<ViewTypeX> x("X", N);
        view_stride_adapter<ViewTypeBr1s_1> b("B", 1);
        view_stride_adapter<ViewTypeY> y("Y", N);

        Kokkos::deep_copy(a.d_base, valueA);
        Kokkos::deep_copy(b.d_base, valueB);
        impl_test_axpby_unification_compare<
            tScalarA, view_stride_adapter<ViewTypeAr1s_1>, view_stride_adapter<ViewTypeX>, tScalarB,
            view_stride_adapter<ViewTypeBr1s_1>, view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, false, max_val,
                                                                                         max_error);
        if (valueB == Kokkos::ArithTraits<tScalarB>::zero()) {
          impl_test_axpby_unification_compare<
              tScalarA, view_stride_adapter<ViewTypeAr1s_1>, view_stride_adapter<ViewTypeX>, tScalarB,
              view_stride_adapter<ViewTypeBr1s_1>, view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, true, max_val,
                                                                                           max_error);
        }
      }
    }
  }

  // ************************************************************
  // Case 12/16: Ar1s_1 + Br1d
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 12/16" << std::endl;
#endif
  for (size_t i(0); i < valuesA.size(); ++i) {
    tScalarA const valueA(valuesA[i]);
    for (size_t j(0); j < valuesB.size(); ++j) {
      tScalarB const valueB(valuesB[j]);
      {
        view_stride_adapter<ViewTypeAr1s_1> a("A", 1);
        view_stride_adapter<ViewTypeX> x("X", N);
        view_stride_adapter<ViewTypeBr1d> b("B", 1);
        view_stride_adapter<ViewTypeY> y("Y", N);

        Kokkos::deep_copy(a.d_base, valueA);
        Kokkos::deep_copy(b.d_base, valueB);
        impl_test_axpby_unification_compare<tScalarA, view_stride_adapter<ViewTypeAr1s_1>,
                                            view_stride_adapter<ViewTypeX>, tScalarB, view_stride_adapter<ViewTypeBr1d>,
                                            view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, false, max_val,
                                                                                    max_error);
        if (valueB == Kokkos::ArithTraits<tScalarB>::zero()) {
          impl_test_axpby_unification_compare<
              tScalarA, view_stride_adapter<ViewTypeAr1s_1>, view_stride_adapter<ViewTypeX>, tScalarB,
              view_stride_adapter<ViewTypeBr1d>, view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, true, max_val,
                                                                                         max_error);
        }
      }
    }
  }

  // ************************************************************
  // Case 13/16: Ar1d + Bscalar
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 13/16" << std::endl;
#endif
  for (size_t i(0); i < valuesA.size(); ++i) {
    tScalarA const valueA(valuesA[i]);
    for (size_t j(0); j < valuesB.size(); ++j) {
      tScalarB const valueB(valuesB[j]);
      {
        view_stride_adapter<ViewTypeAr1d> a("A", 1);
        view_stride_adapter<ViewTypeX> x("X", N);
        tScalarB b;
        view_stride_adapter<ViewTypeY> y("Y", N);

        Kokkos::deep_copy(a.d_base, valueA);
        b = valueB;
        impl_test_axpby_unification_compare<tScalarA, view_stride_adapter<ViewTypeAr1d>, view_stride_adapter<ViewTypeX>,
                                            tScalarB, tScalarB, view_stride_adapter<ViewTypeY>, Device>(
            a, x, b, y, N, false, max_val, max_error);
        if (valueB == Kokkos::ArithTraits<tScalarB>::zero()) {
          impl_test_axpby_unification_compare<tScalarA, view_stride_adapter<ViewTypeAr1d>,
                                              view_stride_adapter<ViewTypeX>, tScalarB, tScalarB,
                                              view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, true, max_val,
                                                                                      max_error);
        }
      }
    }
  }

  // ************************************************************
  // Case 14/16: Ar1d + Br0
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 14/16" << std::endl;
#endif
  if constexpr (std::is_same_v<tLayoutB, Kokkos::LayoutStride>) {
    // Avoid the test, due to compilation errors
  } else {
    for (size_t i(0); i < valuesA.size(); ++i) {
      tScalarA const valueA(valuesA[i]);
      for (size_t j(0); j < valuesB.size(); ++j) {
        tScalarB const valueB(valuesB[j]);
        {
          view_stride_adapter<ViewTypeAr1d> a("A", 1);
          view_stride_adapter<ViewTypeX> x("X", N);
          ViewTypeBr0 b("B");
          view_stride_adapter<ViewTypeY> y("Y", N);

          Kokkos::deep_copy(a.d_base, valueA);
          Kokkos::deep_copy(b, valueB);
          impl_test_axpby_unification_compare<tScalarA, view_stride_adapter<ViewTypeAr1d>,
                                              view_stride_adapter<ViewTypeX>, tScalarB, ViewTypeBr0,
                                              view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, false, max_val,
                                                                                      max_error);
          if (valueB == Kokkos::ArithTraits<tScalarB>::zero()) {
            impl_test_axpby_unification_compare<tScalarA, view_stride_adapter<ViewTypeAr1d>,
                                                view_stride_adapter<ViewTypeX>, tScalarB, ViewTypeBr0,
                                                view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, true, max_val,
                                                                                        max_error);
          }
        }
      }
    }
  }

  // ************************************************************
  // Case 15/16: Ar1d + Br1s_1
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 15/16" << std::endl;
#endif
  for (size_t i(0); i < valuesA.size(); ++i) {
    tScalarA const valueA(valuesA[i]);
    for (size_t j(0); j < valuesB.size(); ++j) {
      tScalarB const valueB(valuesB[j]);
      {
        view_stride_adapter<ViewTypeAr1d> a("A", 1);
        view_stride_adapter<ViewTypeX> x("X", N);
        view_stride_adapter<ViewTypeBr1s_1> b("B", 1);
        view_stride_adapter<ViewTypeY> y("Y", N);

        Kokkos::deep_copy(a.d_base, valueA);
        Kokkos::deep_copy(b.d_base, valueB);
        impl_test_axpby_unification_compare<tScalarA, view_stride_adapter<ViewTypeAr1d>, view_stride_adapter<ViewTypeX>,
                                            tScalarB, view_stride_adapter<ViewTypeBr1s_1>,
                                            view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, false, max_val,
                                                                                    max_error);
        if (valueB == Kokkos::ArithTraits<tScalarB>::zero()) {
          impl_test_axpby_unification_compare<
              tScalarA, view_stride_adapter<ViewTypeAr1d>, view_stride_adapter<ViewTypeX>, tScalarB,
              view_stride_adapter<ViewTypeBr1s_1>, view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, true, max_val,
                                                                                           max_error);
        }
      }
    }
  }

  // ************************************************************
  // Case 16/16: Ar1d + Br1d
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 16/16" << std::endl;
#endif
  for (size_t i(0); i < valuesA.size(); ++i) {
    tScalarA const valueA(valuesA[i]);
    for (size_t j(0); j < valuesB.size(); ++j) {
      tScalarB const valueB(valuesB[j]);
      {
        view_stride_adapter<ViewTypeAr1d> a("A", 1);
        view_stride_adapter<ViewTypeX> x("X", N);
        view_stride_adapter<ViewTypeBr1d> b("B", 1);
        view_stride_adapter<ViewTypeY> y("Y", N);

        Kokkos::deep_copy(a.d_base, valueA);
        Kokkos::deep_copy(b.d_base, valueB);
        impl_test_axpby_unification_compare<tScalarA, view_stride_adapter<ViewTypeAr1d>, view_stride_adapter<ViewTypeX>,
                                            tScalarB, view_stride_adapter<ViewTypeBr1d>, view_stride_adapter<ViewTypeY>,
                                            Device>(a, x, b, y, N, false, max_val, max_error);
        if (valueB == Kokkos::ArithTraits<tScalarB>::zero()) {
          impl_test_axpby_unification_compare<
              tScalarA, view_stride_adapter<ViewTypeAr1d>, view_stride_adapter<ViewTypeX>, tScalarB,
              view_stride_adapter<ViewTypeBr1d>, view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, true, max_val,
                                                                                         max_error);
        }
      }
    }
  }
}

template <class tScalarA, class tLayoutA, class tScalarX, class tLayoutX, class tScalarB, class tLayoutB,
          class tScalarY, class tLayoutY, class Device>
void impl_test_axpby_mv_unification(int const N, int const K) {
  // std::cout << "=========================================" << std::endl;
  // std::cout << "Entering impl_test_axpby_mv_unification()"
  //          << ": tLayoutA = " << typeid(tLayoutA).name()
  //          << ": tLayoutX = " << typeid(tLayoutX).name()
  //          << ", tLayoutB = " << typeid(tLayoutB).name()
  //          << ": tLayoutY = " << typeid(tLayoutY).name()
  //          << std::endl;
  using ViewTypeAr0    = Kokkos::View<tScalarA, tLayoutA, Device>;
  using ViewTypeAr1s_1 = Kokkos::View<tScalarA[1], tLayoutA, Device>;
  using ViewTypeAr1s_k = Kokkos::View<tScalarA[numVecsAxpbyTest], tLayoutA,
                                      Device>;  // Yes, hard coded
  using ViewTypeAr1d   = Kokkos::View<tScalarA*, tLayoutA, Device>;

  using ViewTypeX = Kokkos::View<tScalarX**, tLayoutX, Device>;

  using ViewTypeBr0    = Kokkos::View<tScalarB, tLayoutB, Device>;
  using ViewTypeBr1s_1 = Kokkos::View<tScalarB[1], tLayoutB, Device>;
  using ViewTypeBr1s_k = Kokkos::View<tScalarB[numVecsAxpbyTest], tLayoutB,
                                      Device>;  // Yes, hard coded
  using ViewTypeBr1d   = Kokkos::View<tScalarB*, tLayoutB, Device>;

  using ViewTypeY = Kokkos::View<tScalarY**, tLayoutY, Device>;

  std::array<tScalarA, 4> const valuesA{-1, Kokkos::ArithTraits<tScalarA>::zero(), 1, 3};
  std::array<tScalarB, 4> const valuesB{-1, Kokkos::ArithTraits<tScalarB>::zero(), 1, 5};

  // eps should probably be based on tScalarB since that is the type
  // in which the result is computed.
  using MagnitudeB         = typename Kokkos::ArithTraits<tScalarB>::mag_type;
  MagnitudeB const eps     = Kokkos::ArithTraits<tScalarB>::epsilon();
  MagnitudeB const max_val = 10;
  MagnitudeB const max_error =
      static_cast<MagnitudeB>(Kokkos::ArithTraits<tScalarA>::abs(valuesA[valuesA.size() - 1]) +
                              Kokkos::ArithTraits<tScalarB>::abs(valuesB[valuesB.size() - 1])) *
      max_val * eps;

  // ************************************************************
  // Case 01/36: Ascalar + Bscalar
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 01/36" << std::endl;
#endif
  for (size_t i(0); i < valuesA.size(); ++i) {
    tScalarA const valueA(valuesA[i]);
    for (size_t j(0); j < valuesB.size(); ++j) {
      tScalarB const valueB(valuesB[j]);
      {
        tScalarA a;
        view_stride_adapter<ViewTypeX> x("X", N, K);
        tScalarB b;
        view_stride_adapter<ViewTypeY> y("Y", N, K);

        a = valueA;
        b = valueB;
        impl_test_axpby_mv_unification_compare<tScalarA, tScalarA, view_stride_adapter<ViewTypeX>, tScalarB, tScalarB,
                                               view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, K, false, max_val,
                                                                                       max_error);
        if (valueB == Kokkos::ArithTraits<tScalarB>::zero()) {
          impl_test_axpby_mv_unification_compare<tScalarA, tScalarA, view_stride_adapter<ViewTypeX>, tScalarB, tScalarB,
                                                 view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, K, true,
                                                                                         max_val, max_error);
        }
      }
    }
  }

  // ************************************************************
  // Case 02/36: Ascalar + Br0
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 02/36" << std::endl;
#endif
  if constexpr (std::is_same_v<tLayoutB, Kokkos::LayoutStride>) {
    // Avoid the test, due to compilation errors
  } else {
    for (size_t i(0); i < valuesA.size(); ++i) {
      tScalarA const valueA(valuesA[i]);
      for (size_t j(0); j < valuesB.size(); ++j) {
        tScalarB const valueB(valuesB[j]);
        {
          tScalarA a;
          view_stride_adapter<ViewTypeX> x("X", N, K);
          ViewTypeBr0 b("B");
          view_stride_adapter<ViewTypeY> y("Y", N, K);

          a = valueA;
          Kokkos::deep_copy(b, valueB);
          impl_test_axpby_mv_unification_compare<tScalarA, tScalarA, view_stride_adapter<ViewTypeX>, tScalarB,
                                                 ViewTypeBr0, view_stride_adapter<ViewTypeY>, Device>(
              a, x, b, y, N, K, false, max_val, max_error);
          if (valueB == Kokkos::ArithTraits<tScalarB>::zero()) {
            impl_test_axpby_mv_unification_compare<tScalarA, tScalarA, view_stride_adapter<ViewTypeX>, tScalarB,
                                                   ViewTypeBr0, view_stride_adapter<ViewTypeY>, Device>(
                a, x, b, y, N, K, true, max_val, max_error);
          }
        }
      }
    }
  }

  // ************************************************************
  // Case 03/36: Ascalar + Br1s_1
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 03/36" << std::endl;
#endif
  for (size_t i(0); i < valuesA.size(); ++i) {
    tScalarA const valueA(valuesA[i]);
    for (size_t j(0); j < valuesB.size(); ++j) {
      tScalarB const valueB(valuesB[j]);
      {
        tScalarA a;
        view_stride_adapter<ViewTypeX> x("X", N, K);
        view_stride_adapter<ViewTypeBr1s_1> b("B", 1);
        view_stride_adapter<ViewTypeY> y("Y", N, K);

        a = valueA;
        Kokkos::deep_copy(b.d_base, valueB);
        impl_test_axpby_mv_unification_compare<tScalarA, tScalarA, view_stride_adapter<ViewTypeX>, tScalarB,
                                               view_stride_adapter<ViewTypeBr1s_1>, view_stride_adapter<ViewTypeY>,
                                               Device>(a, x, b, y, N, K, false, max_val, max_error);
        if (valueB == Kokkos::ArithTraits<tScalarB>::zero()) {
          impl_test_axpby_mv_unification_compare<tScalarA, tScalarA, view_stride_adapter<ViewTypeX>, tScalarB,
                                                 view_stride_adapter<ViewTypeBr1s_1>, view_stride_adapter<ViewTypeY>,
                                                 Device>(a, x, b, y, N, K, true, max_val, max_error);
        }
      }
    }
  }

  // ************************************************************
  // Case 04/36: Ascalar + Br1s_k
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 04/36" << std::endl;
#endif
  for (size_t i(0); i < valuesA.size(); ++i) {
    tScalarA const valueA(valuesA[i]);
    for (size_t j(0); j < valuesB.size(); ++j) {
      tScalarB const valueB(valuesB[j]);
      if (K == numVecsAxpbyTest) {
        tScalarA a;
        view_stride_adapter<ViewTypeX> x("X", N, K);
        view_stride_adapter<ViewTypeBr1s_k> b("B", K);
        view_stride_adapter<ViewTypeY> y("Y", N, K);

        a = valueA;
        if constexpr (std::is_same_v<tLayoutB, Kokkos::LayoutStride>) {
          for (int k(0); k < K; ++k) {
            b.h_view[k] = valueB + k;
          }
          Kokkos::deep_copy(b.d_base, b.h_base);
        } else {
          for (int k(0); k < K; ++k) {
            b.h_base[k] = valueB + k;
          }
          Kokkos::deep_copy(b.d_base, b.h_base);
        }
        impl_test_axpby_mv_unification_compare<tScalarA, tScalarA, view_stride_adapter<ViewTypeX>, tScalarB,
                                               view_stride_adapter<ViewTypeBr1s_k>, view_stride_adapter<ViewTypeY>,
                                               Device>(a, x, b, y, N, K, false, max_val, max_error);
      }
    }
  }

  // ************************************************************
  // Case 05/36: Ascalar + Br1d,1
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 05/36" << std::endl;
#endif
  for (size_t i(0); i < valuesA.size(); ++i) {
    tScalarA const valueA(valuesA[i]);
    for (size_t j(0); j < valuesB.size(); ++j) {
      tScalarB const valueB(valuesB[j]);
      {
        tScalarA a;
        view_stride_adapter<ViewTypeX> x("X", N, K);
        view_stride_adapter<ViewTypeBr1d> b("B", 1);
        view_stride_adapter<ViewTypeY> y("Y", N, K);

        a = valueA;
        Kokkos::deep_copy(b.d_base, valueB);
        impl_test_axpby_mv_unification_compare<tScalarA, tScalarA, view_stride_adapter<ViewTypeX>, tScalarB,
                                               view_stride_adapter<ViewTypeBr1d>, view_stride_adapter<ViewTypeY>,
                                               Device>(a, x, b, y, N, K, false, max_val, max_error);
        if (valueB == Kokkos::ArithTraits<tScalarB>::zero()) {
          impl_test_axpby_mv_unification_compare<tScalarA, tScalarA, view_stride_adapter<ViewTypeX>, tScalarB,
                                                 view_stride_adapter<ViewTypeBr1d>, view_stride_adapter<ViewTypeY>,
                                                 Device>(a, x, b, y, N, K, true, max_val, max_error);
        }
      }
    }
  }

  // ************************************************************
  // Case 06/36: Ascalar + Br1d,k
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 06/36" << std::endl;
#endif
  for (size_t i(0); i < valuesA.size(); ++i) {
    tScalarA const valueA(valuesA[i]);
    for (size_t j(0); j < valuesB.size(); ++j) {
      tScalarB const valueB(valuesB[j]);
      {
        tScalarA a;
        view_stride_adapter<ViewTypeX> x("X", N, K);
        view_stride_adapter<ViewTypeBr1d> b("B", K);
        view_stride_adapter<ViewTypeY> y("Y", N, K);

        a = valueA;
        if constexpr (std::is_same_v<tLayoutB, Kokkos::LayoutStride>) {
          for (int k(0); k < K; ++k) {
            b.h_view[k] = valueB + k;
          }
          Kokkos::deep_copy(b.d_base, b.h_base);
        } else {
          for (int k(0); k < K; ++k) {
            b.h_base[k] = valueB + k;
          }
          Kokkos::deep_copy(b.d_base, b.h_base);
        }
        impl_test_axpby_mv_unification_compare<tScalarA, tScalarA, view_stride_adapter<ViewTypeX>, tScalarB,
                                               view_stride_adapter<ViewTypeBr1d>, view_stride_adapter<ViewTypeY>,
                                               Device>(a, x, b, y, N, K, false, max_val, max_error);
      }
    }
  }

  // ************************************************************
  // Case 07/36: Ar0 + Bscalar
  // ************************************************************w
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 07/36" << std::endl;
#endif
  if constexpr (std::is_same_v<tLayoutA, Kokkos::LayoutStride>) {
    // Avoid the test, due to compilation errors
  } else {
    for (size_t i(0); i < valuesA.size(); ++i) {
      tScalarA const valueA(valuesA[i]);
      for (size_t j(0); j < valuesB.size(); ++j) {
        tScalarB const valueB(valuesB[j]);
        {
          ViewTypeAr0 a("A");
          view_stride_adapter<ViewTypeX> x("X", N, K);
          tScalarB b;
          view_stride_adapter<ViewTypeY> y("Y", N, K);

          Kokkos::deep_copy(a, valueA);
          b = valueB;
          impl_test_axpby_mv_unification_compare<tScalarA, ViewTypeAr0, view_stride_adapter<ViewTypeX>, tScalarB,
                                                 tScalarB, view_stride_adapter<ViewTypeY>, Device>(
              a, x, b, y, N, K, false, max_val, max_error);
          if (valueB == Kokkos::ArithTraits<tScalarB>::zero()) {
            impl_test_axpby_mv_unification_compare<tScalarA, ViewTypeAr0, view_stride_adapter<ViewTypeX>, tScalarB,
                                                   tScalarB, view_stride_adapter<ViewTypeY>, Device>(
                a, x, b, y, N, K, true, max_val, max_error);
          }
        }
      }
    }
  }

  // ************************************************************
  // Case 08/36: Ar0 + Br0
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 08/36" << std::endl;
#endif
  if constexpr ((std::is_same_v<tLayoutA, Kokkos::LayoutStride>) || (std::is_same_v<tLayoutB, Kokkos::LayoutStride>)) {
    // Avoid the test, due to compilation errors
  } else {
    for (size_t i(0); i < valuesA.size(); ++i) {
      tScalarA const valueA(valuesA[i]);
      for (size_t j(0); j < valuesB.size(); ++j) {
        tScalarB const valueB(valuesB[j]);
        {
          ViewTypeAr0 a("A");
          view_stride_adapter<ViewTypeX> x("X", N, K);
          ViewTypeBr0 b("B");
          view_stride_adapter<ViewTypeY> y("Y", N, K);

          Kokkos::deep_copy(a, valueA);
          Kokkos::deep_copy(b, valueB);
          impl_test_axpby_mv_unification_compare<tScalarA, ViewTypeAr0, view_stride_adapter<ViewTypeX>, tScalarB,
                                                 ViewTypeBr0, view_stride_adapter<ViewTypeY>, Device>(
              a, x, b, y, N, K, false, max_val, max_error);
          if (valueB == Kokkos::ArithTraits<tScalarB>::zero()) {
            impl_test_axpby_mv_unification_compare<tScalarA, ViewTypeAr0, view_stride_adapter<ViewTypeX>, tScalarB,
                                                   ViewTypeBr0, view_stride_adapter<ViewTypeY>, Device>(
                a, x, b, y, N, K, true, max_val, max_error);
          }
        }
      }
    }
  }

  // ************************************************************
  // Case 09/36: Ar0 + Br1s_1
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 09/36" << std::endl;
#endif
  if constexpr (std::is_same_v<tLayoutA, Kokkos::LayoutStride>) {
    // Avoid the test, due to compilation errors
  } else {
    for (size_t i(0); i < valuesA.size(); ++i) {
      tScalarA const valueA(valuesA[i]);
      for (size_t j(0); j < valuesB.size(); ++j) {
        tScalarB const valueB(valuesB[j]);
        {
          ViewTypeAr0 a("A");
          view_stride_adapter<ViewTypeX> x("X", N, K);
          view_stride_adapter<ViewTypeBr1s_1> b("B", 1);
          view_stride_adapter<ViewTypeY> y("Y", N, K);

          Kokkos::deep_copy(a, valueA);
          Kokkos::deep_copy(b.d_base, valueB);
          impl_test_axpby_mv_unification_compare<tScalarA, ViewTypeAr0, view_stride_adapter<ViewTypeX>, tScalarB,
                                                 view_stride_adapter<ViewTypeBr1s_1>, view_stride_adapter<ViewTypeY>,
                                                 Device>(a, x, b, y, N, K, false, max_val, max_error);
          if (valueB == Kokkos::ArithTraits<tScalarB>::zero()) {
            impl_test_axpby_mv_unification_compare<tScalarA, ViewTypeAr0, view_stride_adapter<ViewTypeX>, tScalarB,
                                                   view_stride_adapter<ViewTypeBr1s_1>, view_stride_adapter<ViewTypeY>,
                                                   Device>(a, x, b, y, N, K, true, max_val, max_error);
          }
        }
      }
    }
  }

  // ************************************************************
  // Case 10/36: Ar0 + Br1s_k
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 10/36" << std::endl;
#endif
  if constexpr (std::is_same_v<tLayoutA, Kokkos::LayoutStride>) {
    // Avoid the test, due to compilation errors
  } else {
    for (size_t i(0); i < valuesA.size(); ++i) {
      tScalarA const valueA(valuesA[i]);
      for (size_t j(0); j < valuesB.size(); ++j) {
        tScalarB const valueB(valuesB[j]);
        if (K == numVecsAxpbyTest) {
          ViewTypeAr0 a("A");
          view_stride_adapter<ViewTypeX> x("X", N, K);
          view_stride_adapter<ViewTypeBr1s_k> b("B", K);
          view_stride_adapter<ViewTypeY> y("Y", N, K);

          Kokkos::deep_copy(a, valueA);
          if constexpr (std::is_same_v<tLayoutB, Kokkos::LayoutStride>) {
            for (int k(0); k < K; ++k) {
              b.h_view[k] = valueB + k;
            }
            Kokkos::deep_copy(b.d_base, b.h_base);
          } else {
            for (int k(0); k < K; ++k) {
              b.h_base[k] = valueB + k;
            }
            Kokkos::deep_copy(b.d_base, b.h_base);
          }
          impl_test_axpby_mv_unification_compare<tScalarA, ViewTypeAr0, view_stride_adapter<ViewTypeX>, tScalarB,
                                                 view_stride_adapter<ViewTypeBr1s_k>, view_stride_adapter<ViewTypeY>,
                                                 Device>(a, x, b, y, N, K, false, max_val, max_error);
        }
      }
    }
  }

  // ************************************************************
  // Case 11/36: Ar0 + Br1d,1
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 11/36" << std::endl;
#endif
  if constexpr (std::is_same_v<tLayoutA, Kokkos::LayoutStride>) {
    // Avoid the test, due to compilation errors
  } else {
    for (size_t i(0); i < valuesA.size(); ++i) {
      tScalarA const valueA(valuesA[i]);
      for (size_t j(0); j < valuesB.size(); ++j) {
        tScalarB const valueB(valuesB[j]);
        {
          ViewTypeAr0 a("A");
          view_stride_adapter<ViewTypeX> x("X", N, K);
          view_stride_adapter<ViewTypeBr1d> b("B", 1);
          view_stride_adapter<ViewTypeY> y("Y", N, K);

          Kokkos::deep_copy(a, valueA);
          Kokkos::deep_copy(b.d_base, valueB);
          impl_test_axpby_mv_unification_compare<tScalarA, ViewTypeAr0, view_stride_adapter<ViewTypeX>, tScalarB,
                                                 view_stride_adapter<ViewTypeBr1d>, view_stride_adapter<ViewTypeY>,
                                                 Device>(a, x, b, y, N, K, false, max_val, max_error);
          if (valueB == Kokkos::ArithTraits<tScalarB>::zero()) {
            impl_test_axpby_mv_unification_compare<tScalarA, ViewTypeAr0, view_stride_adapter<ViewTypeX>, tScalarB,
                                                   view_stride_adapter<ViewTypeBr1d>, view_stride_adapter<ViewTypeY>,
                                                   Device>(a, x, b, y, N, K, true, max_val, max_error);
          }
        }
      }
    }
  }

  // ************************************************************
  // Case 12/36: Ar0 + Br1d,k
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 12/36" << std::endl;
#endif
  if constexpr (std::is_same_v<tLayoutA, Kokkos::LayoutStride>) {
    // Avoid the test, due to compilation errors
  } else {
    for (size_t i(0); i < valuesA.size(); ++i) {
      tScalarA const valueA(valuesA[i]);
      for (size_t j(0); j < valuesB.size(); ++j) {
        tScalarB const valueB(valuesB[j]);
        {
          ViewTypeAr0 a("A");
          view_stride_adapter<ViewTypeX> x("X", N, K);
          view_stride_adapter<ViewTypeBr1d> b("B", K);
          view_stride_adapter<ViewTypeY> y("Y", N, K);

          Kokkos::deep_copy(a, valueA);
          if constexpr (std::is_same_v<tLayoutB, Kokkos::LayoutStride>) {
            for (int k(0); k < K; ++k) {
              b.h_view[k] = valueB + k;
            }
            Kokkos::deep_copy(b.d_base, b.h_base);
          } else {
            for (int k(0); k < K; ++k) {
              b.h_base[k] = valueB + k;
            }
            Kokkos::deep_copy(b.d_base, b.h_base);
          }
          impl_test_axpby_mv_unification_compare<tScalarA, ViewTypeAr0, view_stride_adapter<ViewTypeX>, tScalarB,
                                                 view_stride_adapter<ViewTypeBr1d>, view_stride_adapter<ViewTypeY>,
                                                 Device>(a, x, b, y, N, K, false, max_val, max_error);
        }
      }
    }
  }

  // ************************************************************
  // Case 13/36: Ar1s_1 + Bscalar
  // ************************************************************w
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 13/36" << std::endl;
#endif
  for (size_t i(0); i < valuesA.size(); ++i) {
    tScalarA const valueA(valuesA[i]);
    for (size_t j(0); j < valuesB.size(); ++j) {
      tScalarB const valueB(valuesB[j]);
      {
        view_stride_adapter<ViewTypeAr1s_1> a("A", 1);
        view_stride_adapter<ViewTypeX> x("X", N, K);
        tScalarB b;
        view_stride_adapter<ViewTypeY> y("Y", N, K);

        Kokkos::deep_copy(a.d_base, valueA);
        b = valueB;
        impl_test_axpby_mv_unification_compare<tScalarA, view_stride_adapter<ViewTypeAr1s_1>,
                                               view_stride_adapter<ViewTypeX>, tScalarB, tScalarB,
                                               view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, K, false, max_val,
                                                                                       max_error);
        if (valueB == Kokkos::ArithTraits<tScalarB>::zero()) {
          impl_test_axpby_mv_unification_compare<tScalarA, view_stride_adapter<ViewTypeAr1s_1>,
                                                 view_stride_adapter<ViewTypeX>, tScalarB, tScalarB,
                                                 view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, K, true,
                                                                                         max_val, max_error);
        }
      }
    }
  }

  // ************************************************************
  // Case 14/36: Ar1s_1 + Br0
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 14/36" << std::endl;
#endif
  if constexpr (std::is_same_v<tLayoutB, Kokkos::LayoutStride>) {
    // Avoid the test, due to compilation errors
  } else {
    for (size_t i(0); i < valuesA.size(); ++i) {
      tScalarA const valueA(valuesA[i]);
      for (size_t j(0); j < valuesB.size(); ++j) {
        tScalarB const valueB(valuesB[j]);
        {
          view_stride_adapter<ViewTypeAr1s_1> a("A", 1);
          view_stride_adapter<ViewTypeX> x("X", N, K);
          ViewTypeBr0 b("B");
          view_stride_adapter<ViewTypeY> y("Y", N, K);

          Kokkos::deep_copy(a.d_base, valueA);
          Kokkos::deep_copy(b, valueB);
          impl_test_axpby_mv_unification_compare<tScalarA, view_stride_adapter<ViewTypeAr1s_1>,
                                                 view_stride_adapter<ViewTypeX>, tScalarB, ViewTypeBr0,
                                                 view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, K, false,
                                                                                         max_val, max_error);
          if (valueB == Kokkos::ArithTraits<tScalarB>::zero()) {
            impl_test_axpby_mv_unification_compare<tScalarA, view_stride_adapter<ViewTypeAr1s_1>,
                                                   view_stride_adapter<ViewTypeX>, tScalarB, ViewTypeBr0,
                                                   view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, K, true,
                                                                                           max_val, max_error);
          }
        }
      }
    }
  }

  // ************************************************************
  // Case 15/36: Ar1s_1 + Br1s_1
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 15/36" << std::endl;
#endif
  for (size_t i(0); i < valuesA.size(); ++i) {
    tScalarA const valueA(valuesA[i]);
    for (size_t j(0); j < valuesB.size(); ++j) {
      tScalarB const valueB(valuesB[j]);
      {
        view_stride_adapter<ViewTypeAr1s_1> a("A", 1);
        view_stride_adapter<ViewTypeX> x("X", N, K);
        view_stride_adapter<ViewTypeBr1s_1> b("B", 1);
        view_stride_adapter<ViewTypeY> y("Y", N, K);

        Kokkos::deep_copy(a.d_base, valueA);
        Kokkos::deep_copy(b.d_base, valueB);
        impl_test_axpby_mv_unification_compare<
            tScalarA, view_stride_adapter<ViewTypeAr1s_1>, view_stride_adapter<ViewTypeX>, tScalarB,
            view_stride_adapter<ViewTypeBr1s_1>, view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, K, false,
                                                                                         max_val, max_error);
        if (valueB == Kokkos::ArithTraits<tScalarB>::zero()) {
          impl_test_axpby_mv_unification_compare<
              tScalarA, view_stride_adapter<ViewTypeAr1s_1>, view_stride_adapter<ViewTypeX>, tScalarB,
              view_stride_adapter<ViewTypeBr1s_1>, view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, K, true,
                                                                                           max_val, max_error);
        }
      }
    }
  }

  // ************************************************************
  // Case 16/36: Ar1s_1 + Br1s_k
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 16/36" << std::endl;
#endif
  for (size_t i(0); i < valuesA.size(); ++i) {
    tScalarA const valueA(valuesA[i]);
    for (size_t j(0); j < valuesB.size(); ++j) {
      tScalarB const valueB(valuesB[j]);
      if (K == numVecsAxpbyTest) {
        view_stride_adapter<ViewTypeAr1s_1> a("A", 1);
        view_stride_adapter<ViewTypeX> x("X", N, K);
        view_stride_adapter<ViewTypeBr1s_k> b("B", K);
        view_stride_adapter<ViewTypeY> y("Y", N, K);

        Kokkos::deep_copy(a.d_base, valueA);
        if constexpr (std::is_same_v<tLayoutB, Kokkos::LayoutStride>) {
          for (int k(0); k < K; ++k) {
            b.h_view[k] = valueB + k;
          }
          Kokkos::deep_copy(b.d_base, b.h_base);
        } else {
          for (int k(0); k < K; ++k) {
            b.h_base[k] = valueB + k;
          }
          Kokkos::deep_copy(b.d_base, b.h_base);
        }
        impl_test_axpby_mv_unification_compare<
            tScalarA, view_stride_adapter<ViewTypeAr1s_1>, view_stride_adapter<ViewTypeX>, tScalarB,
            view_stride_adapter<ViewTypeBr1s_k>, view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, K, false,
                                                                                         max_val, max_error);
      }
    }
  }

  // ************************************************************
  // Case 17/36: Ar1s_1 + Br1d,1
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 17/36" << std::endl;
#endif
  for (size_t i(0); i < valuesA.size(); ++i) {
    tScalarA const valueA(valuesA[i]);
    for (size_t j(0); j < valuesB.size(); ++j) {
      tScalarB const valueB(valuesB[j]);
      {
        view_stride_adapter<ViewTypeAr1s_1> a("A", 1);
        view_stride_adapter<ViewTypeX> x("X", N, K);
        view_stride_adapter<ViewTypeBr1d> b("B", 1);
        view_stride_adapter<ViewTypeY> y("Y", N, K);

        Kokkos::deep_copy(a.d_base, valueA);
        Kokkos::deep_copy(b.d_base, valueB);
        impl_test_axpby_mv_unification_compare<
            tScalarA, view_stride_adapter<ViewTypeAr1s_1>, view_stride_adapter<ViewTypeX>, tScalarB,
            view_stride_adapter<ViewTypeBr1d>, view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, K, false, max_val,
                                                                                       max_error);
        if (valueB == Kokkos::ArithTraits<tScalarB>::zero()) {
          impl_test_axpby_mv_unification_compare<
              tScalarA, view_stride_adapter<ViewTypeAr1s_1>, view_stride_adapter<ViewTypeX>, tScalarB,
              view_stride_adapter<ViewTypeBr1d>, view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, K, true,
                                                                                         max_val, max_error);
        }
      }
    }
  }

  // ************************************************************
  // Case 18/36: Ar1s_1 + Br1d,k
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 18/36" << std::endl;
#endif
  for (size_t i(0); i < valuesA.size(); ++i) {
    tScalarA const valueA(valuesA[i]);
    for (size_t j(0); j < valuesB.size(); ++j) {
      tScalarB const valueB(valuesB[j]);
      {
        view_stride_adapter<ViewTypeAr1s_1> a("A", 1);
        view_stride_adapter<ViewTypeX> x("X", N, K);
        view_stride_adapter<ViewTypeBr1d> b("B", K);
        view_stride_adapter<ViewTypeY> y("Y", N, K);

        Kokkos::deep_copy(a.d_base, valueA);
        if constexpr (std::is_same_v<tLayoutB, Kokkos::LayoutStride>) {
          for (int k(0); k < K; ++k) {
            b.h_view[k] = valueB + k;
          }
          Kokkos::deep_copy(b.d_base, b.h_base);
        } else {
          for (int k(0); k < K; ++k) {
            b.h_base[k] = valueB + k;
          }
          Kokkos::deep_copy(b.d_base, b.h_base);
        }
        impl_test_axpby_mv_unification_compare<
            tScalarA, view_stride_adapter<ViewTypeAr1s_1>, view_stride_adapter<ViewTypeX>, tScalarB,
            view_stride_adapter<ViewTypeBr1d>, view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, K, false, max_val,
                                                                                       max_error);
      }
    }
  }

  // ************************************************************
  // Case 19/36: Ar1s_k + Bscalar
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 19/36" << std::endl;
#endif
  for (size_t i(0); i < valuesA.size(); ++i) {
    tScalarA const valueA(valuesA[i]);
    for (size_t j(0); j < valuesB.size(); ++j) {
      tScalarB const valueB(valuesB[j]);
      if (K == numVecsAxpbyTest) {
        view_stride_adapter<ViewTypeAr1s_k> a("A", K);
        view_stride_adapter<ViewTypeX> x("X", N, K);
        tScalarB b;
        view_stride_adapter<ViewTypeY> y("Y", N, K);

        if constexpr (std::is_same_v<tLayoutA, Kokkos::LayoutStride>) {
          for (int k(0); k < K; ++k) {
            a.h_view[k] = valueA + k;
          }
          Kokkos::deep_copy(a.d_base, a.h_base);
        } else {
          for (int k(0); k < K; ++k) {
            a.h_base[k] = valueA + k;
          }
          Kokkos::deep_copy(a.d_base, a.h_base);
        }
        b = valueB;
        impl_test_axpby_mv_unification_compare<tScalarA, view_stride_adapter<ViewTypeAr1s_k>,
                                               view_stride_adapter<ViewTypeX>, tScalarB, tScalarB,
                                               view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, K, false, max_val,
                                                                                       max_error);
        if (valueB == Kokkos::ArithTraits<tScalarB>::zero()) {
          impl_test_axpby_mv_unification_compare<tScalarA, view_stride_adapter<ViewTypeAr1s_k>,
                                                 view_stride_adapter<ViewTypeX>, tScalarB, tScalarB,
                                                 view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, K, true,
                                                                                         max_val, max_error);
        }
      }
    }
  }

  // ************************************************************
  // Case 20/36: Ar1s_k + Br0
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 20/36" << std::endl;
#endif
  if constexpr (std::is_same_v<tLayoutB, Kokkos::LayoutStride>) {
    // Avoid the test, due to compilation errors
  } else {
    for (size_t i(0); i < valuesA.size(); ++i) {
      tScalarA const valueA(valuesA[i]);
      for (size_t j(0); j < valuesB.size(); ++j) {
        tScalarB const valueB(valuesB[j]);
        if (K == numVecsAxpbyTest) {
          view_stride_adapter<ViewTypeAr1s_k> a("A", K);
          view_stride_adapter<ViewTypeX> x("X", N, K);
          ViewTypeBr0 b("B");
          view_stride_adapter<ViewTypeY> y("Y", N, K);

          if constexpr (std::is_same_v<tLayoutA, Kokkos::LayoutStride>) {
            for (int k(0); k < K; ++k) {
              a.h_view[k] = valueA + k;
            }
            Kokkos::deep_copy(a.d_base, a.h_base);
          } else {
            for (int k(0); k < K; ++k) {
              a.h_base[k] = valueA + k;
            }
            Kokkos::deep_copy(a.d_base, a.h_base);
          }
          Kokkos::deep_copy(b, valueB);
          impl_test_axpby_mv_unification_compare<tScalarA, view_stride_adapter<ViewTypeAr1s_k>,
                                                 view_stride_adapter<ViewTypeX>, tScalarB, ViewTypeBr0,
                                                 view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, K, false,
                                                                                         max_val, max_error);
          if (valueB == Kokkos::ArithTraits<tScalarB>::zero()) {
            impl_test_axpby_mv_unification_compare<tScalarA, view_stride_adapter<ViewTypeAr1s_k>,
                                                   view_stride_adapter<ViewTypeX>, tScalarB, ViewTypeBr0,
                                                   view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, K, true,
                                                                                           max_val, max_error);
          }
        }
      }
    }
  }

  // ************************************************************
  // Case 21/36: Ar1s_k + Br1s_1
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 21/36" << std::endl;
#endif
  for (size_t i(0); i < valuesA.size(); ++i) {
    tScalarA const valueA(valuesA[i]);
    for (size_t j(0); j < valuesB.size(); ++j) {
      tScalarB const valueB(valuesB[j]);
      if (K == numVecsAxpbyTest) {
        view_stride_adapter<ViewTypeAr1s_k> a("A", K);
        view_stride_adapter<ViewTypeX> x("X", N, K);
        view_stride_adapter<ViewTypeBr1s_1> b("B", 1);
        view_stride_adapter<ViewTypeY> y("Y", N, K);

        if constexpr (std::is_same_v<tLayoutA, Kokkos::LayoutStride>) {
          for (int k(0); k < K; ++k) {
            a.h_view[k] = valueA + k;
          }
          Kokkos::deep_copy(a.d_base, a.h_base);
        } else {
          for (int k(0); k < K; ++k) {
            a.h_base[k] = valueA + k;
          }
          Kokkos::deep_copy(a.d_base, a.h_base);
        }
        Kokkos::deep_copy(b.d_base, valueB);
        impl_test_axpby_mv_unification_compare<
            tScalarA, view_stride_adapter<ViewTypeAr1s_k>, view_stride_adapter<ViewTypeX>, tScalarB,
            view_stride_adapter<ViewTypeBr1s_1>, view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, K, false,
                                                                                         max_val, max_error);
        if (valueB == Kokkos::ArithTraits<tScalarB>::zero()) {
          impl_test_axpby_mv_unification_compare<
              tScalarA, view_stride_adapter<ViewTypeAr1s_k>, view_stride_adapter<ViewTypeX>, tScalarB,
              view_stride_adapter<ViewTypeBr1s_1>, view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, K, true,
                                                                                           max_val, max_error);
        }
      }
    }
  }

  // ************************************************************
  // Case 22/36: Ar1s_k + Br1s_k
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 22/36" << std::endl;
#endif
  for (size_t i(0); i < valuesA.size(); ++i) {
    tScalarA const valueA(valuesA[i]);
    for (size_t j(0); j < valuesB.size(); ++j) {
      tScalarB const valueB(valuesB[j]);
      if (K == numVecsAxpbyTest) {
        view_stride_adapter<ViewTypeAr1s_k> a("A", K);
        view_stride_adapter<ViewTypeX> x("X", N, K);
        view_stride_adapter<ViewTypeBr1s_k> b("B", K);
        view_stride_adapter<ViewTypeY> y("Y", N, K);

        if constexpr (std::is_same_v<tLayoutA, Kokkos::LayoutStride>) {
          for (int k(0); k < K; ++k) {
            a.h_view[k] = valueA + k;
          }
          Kokkos::deep_copy(a.d_base, a.h_base);
        } else {
          for (int k(0); k < K; ++k) {
            a.h_base[k] = valueA + k;
          }
          Kokkos::deep_copy(a.d_base, a.h_base);
        }

        if constexpr (std::is_same_v<tLayoutB, Kokkos::LayoutStride>) {
          for (int k(0); k < K; ++k) {
            b.h_view[k] = valueB + k;
          }
          Kokkos::deep_copy(b.d_base, b.h_base);
        } else {
          for (int k(0); k < K; ++k) {
            b.h_base[k] = valueB + k;
          }
          Kokkos::deep_copy(b.d_base, b.h_base);
        }
        impl_test_axpby_mv_unification_compare<
            tScalarA, view_stride_adapter<ViewTypeAr1s_k>, view_stride_adapter<ViewTypeX>, tScalarB,
            view_stride_adapter<ViewTypeBr1s_k>, view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, K, false,
                                                                                         max_val, max_error);
      }
    }
  }

  // ************************************************************
  // Case 23/36: Ar1s_k + Br1d,1
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 23/36" << std::endl;
#endif
  for (size_t i(0); i < valuesA.size(); ++i) {
    tScalarA const valueA(valuesA[i]);
    for (size_t j(0); j < valuesB.size(); ++j) {
      tScalarB const valueB(valuesB[j]);
      if (K == numVecsAxpbyTest) {
        view_stride_adapter<ViewTypeAr1s_k> a("A", K);
        view_stride_adapter<ViewTypeX> x("X", N, K);
        view_stride_adapter<ViewTypeBr1d> b("B", 1);
        view_stride_adapter<ViewTypeY> y("Y", N, K);

        if constexpr (std::is_same_v<tLayoutA, Kokkos::LayoutStride>) {
          for (int k(0); k < K; ++k) {
            a.h_view[k] = valueA + k;
          }
          Kokkos::deep_copy(a.d_base, a.h_base);
        } else {
          for (int k(0); k < K; ++k) {
            a.h_base[k] = valueA + k;
          }
          Kokkos::deep_copy(a.d_base, a.h_base);
        }
        Kokkos::deep_copy(b.d_base, valueB);
        impl_test_axpby_mv_unification_compare<
            tScalarA, view_stride_adapter<ViewTypeAr1s_k>, view_stride_adapter<ViewTypeX>, tScalarB,
            view_stride_adapter<ViewTypeBr1d>, view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, K, false, max_val,
                                                                                       max_error);
        if (valueB == Kokkos::ArithTraits<tScalarB>::zero()) {
          impl_test_axpby_mv_unification_compare<
              tScalarA, view_stride_adapter<ViewTypeAr1s_k>, view_stride_adapter<ViewTypeX>, tScalarB,
              view_stride_adapter<ViewTypeBr1d>, view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, K, true,
                                                                                         max_val, max_error);
        }
      }
    }
  }

  // ************************************************************
  // Case 24/36: Ar1s_k + Br1d,k
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 24/36" << std::endl;
#endif
  for (size_t i(0); i < valuesA.size(); ++i) {
    tScalarA const valueA(valuesA[i]);
    for (size_t j(0); j < valuesB.size(); ++j) {
      tScalarB const valueB(valuesB[j]);
      if (K == numVecsAxpbyTest) {
        view_stride_adapter<ViewTypeAr1s_k> a("A", K);
        view_stride_adapter<ViewTypeX> x("X", N, K);
        view_stride_adapter<ViewTypeBr1d> b("B", K);
        view_stride_adapter<ViewTypeY> y("Y", N, K);

        if constexpr (std::is_same_v<tLayoutA, Kokkos::LayoutStride>) {
          for (int k(0); k < K; ++k) {
            a.h_view[k] = valueA + k;
          }
          Kokkos::deep_copy(a.d_base, a.h_base);
        } else {
          for (int k(0); k < K; ++k) {
            a.h_base[k] = valueA + k;
          }
          Kokkos::deep_copy(a.d_base, a.h_base);
        }

        if constexpr (std::is_same_v<tLayoutB, Kokkos::LayoutStride>) {
          for (int k(0); k < K; ++k) {
            b.h_view[k] = valueB + k;
          }
          Kokkos::deep_copy(b.d_base, b.h_base);
        } else {
          for (int k(0); k < K; ++k) {
            b.h_base[k] = valueB + k;
          }
          Kokkos::deep_copy(b.d_base, b.h_base);
        }

        impl_test_axpby_mv_unification_compare<
            tScalarA, view_stride_adapter<ViewTypeAr1s_k>, view_stride_adapter<ViewTypeX>, tScalarB,
            view_stride_adapter<ViewTypeBr1d>, view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, K, false, max_val,
                                                                                       max_error);
      }
    }
  }

  // ************************************************************
  // Case 25/36: Ar1d,1 + Bscalar
  // ************************************************************w
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 25/36" << std::endl;
#endif
  for (size_t i(0); i < valuesA.size(); ++i) {
    tScalarA const valueA(valuesA[i]);
    for (size_t j(0); j < valuesB.size(); ++j) {
      tScalarB const valueB(valuesB[j]);
      {
        view_stride_adapter<ViewTypeAr1d> a("A", 1);
        view_stride_adapter<ViewTypeX> x("X", N, K);
        tScalarB b;
        view_stride_adapter<ViewTypeY> y("Y", N, K);

        Kokkos::deep_copy(a.d_base, valueA);
        b = valueB;
        impl_test_axpby_mv_unification_compare<tScalarA, view_stride_adapter<ViewTypeAr1d>,
                                               view_stride_adapter<ViewTypeX>, tScalarB, tScalarB,
                                               view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, K, false, max_val,
                                                                                       max_error);
        if (valueB == Kokkos::ArithTraits<tScalarB>::zero()) {
          impl_test_axpby_mv_unification_compare<tScalarA, view_stride_adapter<ViewTypeAr1d>,
                                                 view_stride_adapter<ViewTypeX>, tScalarB, tScalarB,
                                                 view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, K, true,
                                                                                         max_val, max_error);
        }
      }
    }
  }

  // ************************************************************
  // Case 26/36: Ar1d,1 + Br0
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 26/36" << std::endl;
#endif
  if constexpr (std::is_same_v<tLayoutB, Kokkos::LayoutStride>) {
    // Avoid the test, due to compilation errors
  } else {
    for (size_t i(0); i < valuesA.size(); ++i) {
      tScalarA const valueA(valuesA[i]);
      for (size_t j(0); j < valuesB.size(); ++j) {
        tScalarB const valueB(valuesB[j]);
        {
          view_stride_adapter<ViewTypeAr1d> a("A", 1);
          view_stride_adapter<ViewTypeX> x("X", N, K);
          ViewTypeBr0 b("B");
          view_stride_adapter<ViewTypeY> y("Y", N, K);

          Kokkos::deep_copy(a.d_base, valueA);
          Kokkos::deep_copy(b, valueB);
          impl_test_axpby_mv_unification_compare<tScalarA, view_stride_adapter<ViewTypeAr1d>,
                                                 view_stride_adapter<ViewTypeX>, tScalarB, ViewTypeBr0,
                                                 view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, K, false,
                                                                                         max_val, max_error);
          if (valueB == Kokkos::ArithTraits<tScalarB>::zero()) {
            impl_test_axpby_mv_unification_compare<tScalarA, view_stride_adapter<ViewTypeAr1d>,
                                                   view_stride_adapter<ViewTypeX>, tScalarB, ViewTypeBr0,
                                                   view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, K, true,
                                                                                           max_val, max_error);
          }
        }
      }
    }
  }

  // ************************************************************
  // Case 27/36: Ar1d,1 + Br1s_1
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 27/36" << std::endl;
#endif
  for (size_t i(0); i < valuesA.size(); ++i) {
    tScalarA const valueA(valuesA[i]);
    for (size_t j(0); j < valuesB.size(); ++j) {
      tScalarB const valueB(valuesB[j]);
      {
        view_stride_adapter<ViewTypeAr1d> a("A", 1);
        view_stride_adapter<ViewTypeX> x("X", N, K);
        view_stride_adapter<ViewTypeBr1s_1> b("B", 1);
        view_stride_adapter<ViewTypeY> y("Y", N, K);

        Kokkos::deep_copy(a.d_base, valueA);
        Kokkos::deep_copy(b.d_base, valueB);
        impl_test_axpby_mv_unification_compare<
            tScalarA, view_stride_adapter<ViewTypeAr1d>, view_stride_adapter<ViewTypeX>, tScalarB,
            view_stride_adapter<ViewTypeBr1s_1>, view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, K, false,
                                                                                         max_val, max_error);
        if (valueB == Kokkos::ArithTraits<tScalarB>::zero()) {
          impl_test_axpby_mv_unification_compare<
              tScalarA, view_stride_adapter<ViewTypeAr1d>, view_stride_adapter<ViewTypeX>, tScalarB,
              view_stride_adapter<ViewTypeBr1s_1>, view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, K, true,
                                                                                           max_val, max_error);
        }
      }
    }
  }

  // ************************************************************
  // Case 28/36: Ar1d,1 + Br1s_k
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 28/36" << std::endl;
#endif
  for (size_t i(0); i < valuesA.size(); ++i) {
    tScalarA const valueA(valuesA[i]);
    for (size_t j(0); j < valuesB.size(); ++j) {
      tScalarB const valueB(valuesB[j]);
      if (K == numVecsAxpbyTest) {
        view_stride_adapter<ViewTypeAr1d> a("A", 1);
        view_stride_adapter<ViewTypeX> x("X", N, K);
        view_stride_adapter<ViewTypeBr1s_k> b("B", K);
        view_stride_adapter<ViewTypeY> y("Y", N, K);

        Kokkos::deep_copy(a.d_base, valueA);
        if constexpr (std::is_same_v<tLayoutB, Kokkos::LayoutStride>) {
          for (int k(0); k < K; ++k) {
            b.h_view[k] = valueB + k;
          }
          Kokkos::deep_copy(b.d_base, b.h_base);
        } else {
          for (int k(0); k < K; ++k) {
            b.h_base[k] = valueB + k;
          }
          Kokkos::deep_copy(b.d_base, b.h_base);
        }
        impl_test_axpby_mv_unification_compare<
            tScalarA, view_stride_adapter<ViewTypeAr1d>, view_stride_adapter<ViewTypeX>, tScalarB,
            view_stride_adapter<ViewTypeBr1s_k>, view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, K, false,
                                                                                         max_val, max_error);
      }
    }
  }

  // ************************************************************
  // Case 29/36: Ar1d,1 + Br1d,1
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 29/36" << std::endl;
#endif
  for (size_t i(0); i < valuesA.size(); ++i) {
    tScalarA const valueA(valuesA[i]);
    for (size_t j(0); j < valuesB.size(); ++j) {
      tScalarB const valueB(valuesB[j]);
      {
        view_stride_adapter<ViewTypeAr1d> a("A", 1);
        view_stride_adapter<ViewTypeX> x("X", N, K);
        view_stride_adapter<ViewTypeBr1d> b("B", 1);
        view_stride_adapter<ViewTypeY> y("Y", N, K);

        Kokkos::deep_copy(a.d_base, valueA);
        Kokkos::deep_copy(b.d_base, valueB);
        impl_test_axpby_mv_unification_compare<
            tScalarA, view_stride_adapter<ViewTypeAr1d>, view_stride_adapter<ViewTypeX>, tScalarB,
            view_stride_adapter<ViewTypeBr1d>, view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, K, false, max_val,
                                                                                       max_error);
        if (valueB == Kokkos::ArithTraits<tScalarB>::zero()) {
          impl_test_axpby_mv_unification_compare<
              tScalarA, view_stride_adapter<ViewTypeAr1d>, view_stride_adapter<ViewTypeX>, tScalarB,
              view_stride_adapter<ViewTypeBr1d>, view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, K, true,
                                                                                         max_val, max_error);
        }
      }
    }
  }

  // ************************************************************
  // Case 30/36: Ar1d,1 + Br1d,k
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 30/36" << std::endl;
#endif
  for (size_t i(0); i < valuesA.size(); ++i) {
    tScalarA const valueA(valuesA[i]);
    for (size_t j(0); j < valuesB.size(); ++j) {
      tScalarB const valueB(valuesB[j]);
      {
        view_stride_adapter<ViewTypeAr1d> a("A", 1);
        view_stride_adapter<ViewTypeX> x("X", N, K);
        view_stride_adapter<ViewTypeBr1d> b("B", K);
        view_stride_adapter<ViewTypeY> y("Y", N, K);

        Kokkos::deep_copy(a.d_base, valueA);
        if constexpr (std::is_same_v<tLayoutB, Kokkos::LayoutStride>) {
          for (int k(0); k < K; ++k) {
            b.h_view[k] = valueB + k;
          }
          Kokkos::deep_copy(b.d_base, b.h_base);
        } else {
          for (int k(0); k < K; ++k) {
            b.h_base[k] = valueB + k;
          }
          Kokkos::deep_copy(b.d_base, b.h_base);
        }
        impl_test_axpby_mv_unification_compare<
            tScalarA, view_stride_adapter<ViewTypeAr1d>, view_stride_adapter<ViewTypeX>, tScalarB,
            view_stride_adapter<ViewTypeBr1d>, view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, K, false, max_val,
                                                                                       max_error);
      }
    }
  }

  // ************************************************************
  // Case 31/36: Ar1d,k + Bscalar
  // ************************************************************w
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 31/36" << std::endl;
#endif
  for (size_t i(0); i < valuesA.size(); ++i) {
    tScalarA const valueA(valuesA[i]);
    for (size_t j(0); j < valuesB.size(); ++j) {
      tScalarB const valueB(valuesB[j]);
      {
        view_stride_adapter<ViewTypeAr1d> a("A", K);
        view_stride_adapter<ViewTypeX> x("X", N, K);
        tScalarB b;
        view_stride_adapter<ViewTypeY> y("Y", N, K);

        if constexpr (std::is_same_v<tLayoutA, Kokkos::LayoutStride>) {
          for (int k(0); k < K; ++k) {
            a.h_view[k] = valueA + k;
          }
          Kokkos::deep_copy(a.d_base, a.h_base);
        } else {
          for (int k(0); k < K; ++k) {
            a.h_base[k] = valueA + k;
          }
          Kokkos::deep_copy(a.d_base, a.h_base);
        }
        b = valueB;
        impl_test_axpby_mv_unification_compare<tScalarA, view_stride_adapter<ViewTypeAr1d>,
                                               view_stride_adapter<ViewTypeX>, tScalarB, tScalarB,
                                               view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, K, false, max_val,
                                                                                       max_error);
        if (valueB == Kokkos::ArithTraits<tScalarB>::zero()) {
          impl_test_axpby_mv_unification_compare<tScalarA, view_stride_adapter<ViewTypeAr1d>,
                                                 view_stride_adapter<ViewTypeX>, tScalarB, tScalarB,
                                                 view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, K, true,
                                                                                         max_val, max_error);
        }
      }
    }
  }

  // ************************************************************
  // Case 32/36: Ar1d,k + Br0
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 32/36" << std::endl;
#endif
  if constexpr (std::is_same_v<tLayoutB, Kokkos::LayoutStride>) {
    // Avoid the test, due to compilation errors
  } else {
    for (size_t i(0); i < valuesA.size(); ++i) {
      tScalarA const valueA(valuesA[i]);
      for (size_t j(0); j < valuesB.size(); ++j) {
        tScalarB const valueB(valuesB[j]);
        {
          view_stride_adapter<ViewTypeAr1d> a("A", K);
          view_stride_adapter<ViewTypeX> x("X", N, K);
          ViewTypeBr0 b("B");
          view_stride_adapter<ViewTypeY> y("Y", N, K);

          if constexpr (std::is_same_v<tLayoutA, Kokkos::LayoutStride>) {
            for (int k(0); k < K; ++k) {
              a.h_view[k] = valueA + k;
            }
            Kokkos::deep_copy(a.d_base, a.h_base);
          } else {
            for (int k(0); k < K; ++k) {
              a.h_base[k] = valueA + k;
            }
            Kokkos::deep_copy(a.d_base, a.h_base);
          }
          Kokkos::deep_copy(b, valueB);
          impl_test_axpby_mv_unification_compare<tScalarA, view_stride_adapter<ViewTypeAr1d>,
                                                 view_stride_adapter<ViewTypeX>, tScalarB, ViewTypeBr0,
                                                 view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, K, false,
                                                                                         max_val, max_error);
          if (valueB == Kokkos::ArithTraits<tScalarB>::zero()) {
            impl_test_axpby_mv_unification_compare<tScalarA, view_stride_adapter<ViewTypeAr1d>,
                                                   view_stride_adapter<ViewTypeX>, tScalarB, ViewTypeBr0,
                                                   view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, K, true,
                                                                                           max_val, max_error);
          }
        }
      }
    }
  }

  // ************************************************************
  // Case 33/36: Ar1d,k + Br1s_1
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 33/36" << std::endl;
#endif
  for (size_t i(0); i < valuesA.size(); ++i) {
    tScalarA const valueA(valuesA[i]);
    for (size_t j(0); j < valuesB.size(); ++j) {
      tScalarB const valueB(valuesB[j]);
      {
        view_stride_adapter<ViewTypeAr1d> a("A", K);
        view_stride_adapter<ViewTypeX> x("X", N, K);
        view_stride_adapter<ViewTypeBr1s_1> b("B", 1);
        view_stride_adapter<ViewTypeY> y("Y", N, K);

        if constexpr (std::is_same_v<tLayoutA, Kokkos::LayoutStride>) {
          for (int k(0); k < K; ++k) {
            a.h_view[k] = valueA + k;
          }
          Kokkos::deep_copy(a.d_base, a.h_base);
        } else {
          for (int k(0); k < K; ++k) {
            a.h_base[k] = valueA + k;
          }
          Kokkos::deep_copy(a.d_base, a.h_base);
        }
        Kokkos::deep_copy(b.d_base, valueB);
        impl_test_axpby_mv_unification_compare<
            tScalarA, view_stride_adapter<ViewTypeAr1d>, view_stride_adapter<ViewTypeX>, tScalarB,
            view_stride_adapter<ViewTypeBr1s_1>, view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, K, false,
                                                                                         max_val, max_error);
        if (valueB == Kokkos::ArithTraits<tScalarB>::zero()) {
          impl_test_axpby_mv_unification_compare<
              tScalarA, view_stride_adapter<ViewTypeAr1d>, view_stride_adapter<ViewTypeX>, tScalarB,
              view_stride_adapter<ViewTypeBr1s_1>, view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, K, true,
                                                                                           max_val, max_error);
        }
      }
    }
  }

  // ************************************************************
  // Case 34/36: Ar1d,k + Br1s_k
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 34/36" << std::endl;
#endif
  for (size_t i(0); i < valuesA.size(); ++i) {
    tScalarA const valueA(valuesA[i]);
    for (size_t j(0); j < valuesB.size(); ++j) {
      tScalarB const valueB(valuesB[j]);
      if (K == numVecsAxpbyTest) {
        view_stride_adapter<ViewTypeAr1d> a("A", K);
        view_stride_adapter<ViewTypeX> x("X", N, K);
        view_stride_adapter<ViewTypeBr1s_k> b("B", K);
        view_stride_adapter<ViewTypeY> y("Y", N, K);

        if constexpr (std::is_same_v<tLayoutA, Kokkos::LayoutStride>) {
          for (int k(0); k < K; ++k) {
            a.h_view[k] = valueA + k;
          }
          Kokkos::deep_copy(a.d_base, a.h_base);
        } else {
          for (int k(0); k < K; ++k) {
            a.h_base[k] = valueA + k;
          }
          Kokkos::deep_copy(a.d_base, a.h_base);
        }

        if constexpr (std::is_same_v<tLayoutB, Kokkos::LayoutStride>) {
          for (int k(0); k < K; ++k) {
            b.h_view[k] = valueB + k;
          }
          Kokkos::deep_copy(b.d_base, b.h_base);
        } else {
          for (int k(0); k < K; ++k) {
            b.h_base[k] = valueB + k;
          }
          Kokkos::deep_copy(b.d_base, b.h_base);
        }

        impl_test_axpby_mv_unification_compare<
            tScalarA, view_stride_adapter<ViewTypeAr1d>, view_stride_adapter<ViewTypeX>, tScalarB,
            view_stride_adapter<ViewTypeBr1s_k>, view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, K, false,
                                                                                         max_val, max_error);
      }
    }
  }

  // ************************************************************
  // Case 35/36: Ar1d,k + Br1d,1
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 35/36" << std::endl;
#endif
  for (size_t i(0); i < valuesA.size(); ++i) {
    tScalarA const valueA(valuesA[i]);
    for (size_t j(0); j < valuesB.size(); ++j) {
      tScalarB const valueB(valuesB[j]);
      {
        view_stride_adapter<ViewTypeAr1d> a("A", K);
        view_stride_adapter<ViewTypeX> x("X", N, K);
        view_stride_adapter<ViewTypeBr1d> b("B", 1);
        view_stride_adapter<ViewTypeY> y("Y", N, K);

        if constexpr (std::is_same_v<tLayoutA, Kokkos::LayoutStride>) {
          for (int k(0); k < K; ++k) {
            a.h_view[k] = valueA + k;
          }
          Kokkos::deep_copy(a.d_base, a.h_base);
        } else {
          for (int k(0); k < K; ++k) {
            a.h_base[k] = valueA + k;
          }
          Kokkos::deep_copy(a.d_base, a.h_base);
        }
        Kokkos::deep_copy(b.d_base, valueB);
        impl_test_axpby_mv_unification_compare<
            tScalarA, view_stride_adapter<ViewTypeAr1d>, view_stride_adapter<ViewTypeX>, tScalarB,
            view_stride_adapter<ViewTypeBr1d>, view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, K, false, max_val,
                                                                                       max_error);
        if (valueB == Kokkos::ArithTraits<tScalarB>::zero()) {
          impl_test_axpby_mv_unification_compare<
              tScalarA, view_stride_adapter<ViewTypeAr1d>, view_stride_adapter<ViewTypeX>, tScalarB,
              view_stride_adapter<ViewTypeBr1d>, view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, K, true,
                                                                                         max_val, max_error);
        }
      }
    }
  }

  // ************************************************************
  // Case 36/36: Ar1d,k + Br1d,k
  // ************************************************************
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Starting case 36/36" << std::endl;
#endif
  for (size_t i(0); i < valuesA.size(); ++i) {
    tScalarA const valueA(valuesA[i]);
    for (size_t j(0); j < valuesB.size(); ++j) {
      tScalarB const valueB(valuesB[j]);
      {
        view_stride_adapter<ViewTypeAr1d> a("A", K);
        view_stride_adapter<ViewTypeX> x("X", N, K);
        view_stride_adapter<ViewTypeBr1d> b("B", K);
        view_stride_adapter<ViewTypeY> y("Y", N, K);

        if constexpr (std::is_same_v<tLayoutA, Kokkos::LayoutStride>) {
          for (int k(0); k < K; ++k) {
            a.h_view[k] = valueA + k;
          }
          Kokkos::deep_copy(a.d_base, a.h_base);
        } else {
          for (int k(0); k < K; ++k) {
            a.h_base[k] = valueA + k;
          }
          Kokkos::deep_copy(a.d_base, a.h_base);
        }

        if constexpr (std::is_same_v<tLayoutB, Kokkos::LayoutStride>) {
          for (int k(0); k < K; ++k) {
            b.h_view[k] = valueB + k;
          }
          Kokkos::deep_copy(b.d_base, b.h_base);
        } else {
          for (int k(0); k < K; ++k) {
            b.h_base[k] = valueB + k;
          }
          Kokkos::deep_copy(b.d_base, b.h_base);
        }

        impl_test_axpby_mv_unification_compare<
            tScalarA, view_stride_adapter<ViewTypeAr1d>, view_stride_adapter<ViewTypeX>, tScalarB,
            view_stride_adapter<ViewTypeBr1d>, view_stride_adapter<ViewTypeY>, Device>(a, x, b, y, N, K, false, max_val,
                                                                                       max_error);
      }
    }
  }

  // std::cout << "Leaving impl_test_axpby_mv_unification()" << std::endl;
  // std::cout << "=========================================" << std::endl;
}

}  // namespace Test

template <class tScalarA, class tScalarX, class tScalarB, class tScalarY, class Device>
int test_axpby_unification() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Calling impl_test_axpby_unif(), L-LLL" << std::endl;
#endif
  Test::impl_test_axpby_unification<tScalarA, Kokkos::LayoutLeft, tScalarX, Kokkos::LayoutLeft, tScalarB,
                                    Kokkos::LayoutLeft, tScalarY, Kokkos::LayoutLeft, Device>(14);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Calling impl_test_axpby_unif(), L-RRR" << std::endl;
#endif
  Test::impl_test_axpby_unification<tScalarA, Kokkos::LayoutRight, tScalarX, Kokkos::LayoutRight, tScalarB,
                                    Kokkos::LayoutRight, tScalarY, Kokkos::LayoutRight, Device>(14);
#endif

#if (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Calling impl_test_axpby_unif(), L-SSS" << std::endl;
#endif
  Test::impl_test_axpby_unification<tScalarA, Kokkos::LayoutStride, tScalarX, Kokkos::LayoutStride, tScalarB,
                                    Kokkos::LayoutStride, tScalarY, Kokkos::LayoutStride, Device>(14);
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Calling impl_test_axpby_unif(), L-SLL" << std::endl;
#endif
  Test::impl_test_axpby_unification<tScalarA, Kokkos::LayoutStride, tScalarX, Kokkos::LayoutStride, tScalarB,
                                    Kokkos::LayoutLeft, tScalarY, Kokkos::LayoutLeft, Device>(14);

#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Calling impl_test_axpby_unif(), L-LSS" << std::endl;
#endif
  Test::impl_test_axpby_unification<tScalarA, Kokkos::LayoutLeft, tScalarX, Kokkos::LayoutLeft, tScalarB,
                                    Kokkos::LayoutStride, tScalarY, Kokkos::LayoutStride, Device>(14);

#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Calling impl_test_axpby_unif(), L-SRS" << std::endl;
#endif
  Test::impl_test_axpby_unification<tScalarA, Kokkos::LayoutLeft, tScalarX, Kokkos::LayoutStride, tScalarB,
                                    Kokkos::LayoutRight, tScalarY, Kokkos::LayoutStride, Device>(14);

#ifdef HAVE_KOKKOSKERNELS_DEBUG
  std::cout << "Calling impl_test_axpby_unif(), L-LSR" << std::endl;
#endif
  Test::impl_test_axpby_unification<tScalarA, Kokkos::LayoutStride, tScalarX, Kokkos::LayoutLeft, tScalarB,
                                    Kokkos::LayoutStride, tScalarY, Kokkos::LayoutRight, Device>(14);
#endif
  return 1;
}

template <class tScalarA, class tScalarX, class tScalarB, class tScalarY, class Device>
int test_axpby_mv_unification() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  Test::impl_test_axpby_mv_unification<tScalarA, Kokkos::LayoutLeft, tScalarX, Kokkos::LayoutLeft, tScalarB,
                                       Kokkos::LayoutLeft, tScalarY, Kokkos::LayoutLeft, Device>(14, numVecsAxpbyTest);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  Test::impl_test_axpby_mv_unification<tScalarA, Kokkos::LayoutRight, tScalarX, Kokkos::LayoutRight, tScalarB,
                                       Kokkos::LayoutRight, tScalarY, Kokkos::LayoutRight, Device>(14,
                                                                                                   numVecsAxpbyTest);
#endif

#if (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  Test::impl_test_axpby_mv_unification<tScalarA, Kokkos::LayoutStride, tScalarX, Kokkos::LayoutStride, tScalarB,
                                       Kokkos::LayoutStride, tScalarY, Kokkos::LayoutStride, Device>(14,
                                                                                                     numVecsAxpbyTest);
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
  Test::impl_test_axpby_mv_unification<tScalarA, Kokkos::LayoutStride, tScalarX, Kokkos::LayoutStride, tScalarB,
                                       Kokkos::LayoutLeft, tScalarY, Kokkos::LayoutLeft, Device>(14, numVecsAxpbyTest);
  Test::impl_test_axpby_mv_unification<tScalarA, Kokkos::LayoutLeft, tScalarX, Kokkos::LayoutLeft, tScalarB,
                                       Kokkos::LayoutStride, tScalarY, Kokkos::LayoutStride, Device>(14,
                                                                                                     numVecsAxpbyTest);

  Test::impl_test_axpby_mv_unification<tScalarA, Kokkos::LayoutLeft, tScalarX, Kokkos::LayoutStride, tScalarB,
                                       Kokkos::LayoutRight, tScalarY, Kokkos::LayoutStride, Device>(14,
                                                                                                    numVecsAxpbyTest);

  Test::impl_test_axpby_mv_unification<tScalarA, Kokkos::LayoutStride, tScalarX, Kokkos::LayoutLeft, tScalarB,
                                       Kokkos::LayoutStride, tScalarY, Kokkos::LayoutRight, Device>(14,
                                                                                                    numVecsAxpbyTest);
#endif
  return 1;
}

#if defined(KOKKOSKERNELS_INST_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, axpby_unification_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::axpby_unification_float");
  test_axpby_unification<float, float, float, float, TestDevice>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, axpby_mv_unification_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::axpby_mv_unification_float");
  test_axpby_mv_unification<float, float, float, float, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, axpby_unification_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::axpby_unification_double");
  test_axpby_unification<double, double, double, double, TestDevice>();
}
TEST_F(TestCategory, axpby_mv_unification_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::axpby_mv_unification_double");
  test_axpby_mv_unification<double, double, double, double, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, axpby_unification_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::axpby_unification_complex_double");
  test_axpby_unification<Kokkos::complex<double>, Kokkos::complex<double>, Kokkos::complex<double>,
                         Kokkos::complex<double>, TestDevice>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, axpby_mv_unification_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::axpby_mv_unification_complex_double");
  test_axpby_mv_unification<Kokkos::complex<double>, Kokkos::complex<double>, Kokkos::complex<double>,
                            Kokkos::complex<double>, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_INT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, axpby_unification_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::axpby_unification_int");
  test_axpby_unification<int, int, int, int, TestDevice>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, axpby_mv_unification_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::axpby_mv_unification_int");
  test_axpby_mv_unification<int, int, int, int, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
TEST_F(TestCategory, axpby_unification_double_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::axpby_unification_double_int");
  test_axpby_unification<double, double, int, int, TestDevice>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, axpby_double_mv_unification_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::axpby_mv_unification_double_int");
  test_axpby_mv_unification<double, double, int, int, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif
