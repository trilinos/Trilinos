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
#ifndef KOKKOS_BLAS1_IMPL_DOT_SPEC_HPP_
#define KOKKOS_BLAS1_IMPL_DOT_SPEC_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>
#include <Kokkos_InnerProductSpaceTraits.hpp>

// Include the actual functors
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include <KokkosBlas1_dot_impl.hpp>
#include <KokkosBlas1_dot_mv_impl.hpp>
#endif

namespace KokkosBlas {
namespace Impl {

// Some platforms, such as Mac Clang, seem to get poor accuracy with
// float and complex<float>.  Work around some Trilinos test
// failures by using a higher-precision type for intermediate dot
// product sums.
//
// Note that this is not the same thing as
// InnerProductSpaceTraits<scalar>::dot_type
template <typename scalar_t>
struct DotAccumulatingScalar {
  using type = scalar_t;
};

template <>
struct DotAccumulatingScalar<float> {
  using type = double;
};

template <>
struct DotAccumulatingScalar<Kokkos::complex<float>> {
  using type = Kokkos::complex<double>;
};

template <typename scalar_t>
struct HasSpecialAccumulator {
  enum : bool { value = !std::is_same<scalar_t, typename DotAccumulatingScalar<scalar_t>::type>::value };
};

// Specialization struct which defines whether a specialization exists
template <class execution_space, class AV, class XV, class YV, int Xrank = XV::rank, int Yrank = YV::rank>
struct dot_eti_spec_avail {
  enum : bool { value = false };
};
}  // namespace Impl
}  // namespace KokkosBlas

//
// Macro for declaration of full specialization availability
// KokkosBlas::Impl::Dot for rank == 1.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _INST macro below) across one or
// more .cpp files.
#define KOKKOSBLAS1_DOT_ETI_SPEC_AVAIL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                                         \
  template <>                                                                                                         \
  struct dot_eti_spec_avail<EXEC_SPACE,                                                                               \
                            Kokkos::View<SCALAR, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
                            Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                \
                                         Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                    \
                            Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                \
                                         Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                    \
                            1, 1> {                                                                                   \
    enum : bool { value = true };                                                                                     \
  };                                                                                                                  \
  template <>                                                                                                         \
  struct dot_eti_spec_avail<                                                                                          \
      EXEC_SPACE,                                                                                                     \
      Kokkos::View<SCALAR, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,   \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                          \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                          \
      1, 1> {                                                                                                         \
    enum : bool { value = true };                                                                                     \
  };

//
// Macro for declaration of full specialization availability
// KokkosBlas::Impl::Dot for rank == 2.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _DEF macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_DOT_MV_ETI_SPEC_AVAIL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                          \
  template <>                                                                                             \
  struct dot_eti_spec_avail<                                                                              \
      EXEC_SPACE,                                                                                         \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                              \
      Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                         \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                              \
      Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                         \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                              \
      2, 2> {                                                                                             \
    enum : bool { value = true };                                                                         \
  };                                                                                                      \
  template <>                                                                                             \
  struct dot_eti_spec_avail<                                                                              \
      EXEC_SPACE,                                                                                         \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                              \
      Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                         \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                              \
      Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                         \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                              \
      2, 1> {                                                                                             \
    enum : bool { value = true };                                                                         \
  };                                                                                                      \
  template <>                                                                                             \
  struct dot_eti_spec_avail<                                                                              \
      EXEC_SPACE,                                                                                         \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                              \
      Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                         \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                              \
      Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                         \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                              \
      1, 2> {                                                                                             \
    enum : bool { value = true };                                                                         \
  };

// Include the actual specialization declarations
#include <KokkosBlas1_dot_tpl_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosBlas1_dot_eti_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosBlas1_dot_mv_eti_spec_avail.hpp>

namespace KokkosBlas {
namespace Impl {

// Unification layer
template <class execution_space, class RV, class XV, class YV, int XV_Rank = XV::rank, int YV_Rank = YV::rank,
          bool tpl_spec_avail = dot_tpl_spec_avail<execution_space, RV, XV, YV>::value,
          bool eti_spec_avail = dot_eti_spec_avail<execution_space, RV, XV, YV>::value>
struct Dot {
  static void dot(const execution_space& space, const RV&, const XV& R, const YV& X);
};

// This version never has TPL support, but it does use the same ETI system
template <class execution_space, class RV, class XV, class YV,
          bool eti_spec_avail = dot_eti_spec_avail<execution_space, RV, XV, YV>::value>
struct DotSpecialAccumulator {
  // Note: not doing the static_asserts to validate RV, XV, YV since those
  // errors would have already arisen when building the library.
  using size_type  = typename YV::size_type;
  using dot_type   = typename Kokkos::Details::InnerProductSpaceTraits<typename XV::non_const_value_type>::dot_type;
  using accum_type = typename DotAccumulatingScalar<dot_type>::type;
  // This is the same View type as RV, but using the special accumulator as the
  // value type
  using RV_Result = Kokkos::View<accum_type, typename RV::array_layout, typename RV::device_type,
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  static void dot(const execution_space& space, const RV_Result& R, const XV& X, const YV& Y);
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
//! Full specialization of Dot for single vectors (1-D Views).
//  The rank-1 case is currently the only one that may use a different
//  accumulator type than <tt>InnerProductSpaceTraits::dot_type</tt>.
template <class execution_space, class RV, class XV, class YV>
struct Dot<execution_space, RV, XV, YV, 1, 1, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  // Check some things about the template parameters at compile time to get nice
  // error messages, before using them under the assumption they are valid.
  static_assert(Kokkos::is_view<XV>::value,
                "KokkosBlas::Impl::"
                "Dot<1-D>: XV is not a Kokkos::View.");
  static_assert(Kokkos::is_view<YV>::value,
                "KokkosBlas::Impl::"
                "Dot<1-D>: YV is not a Kokkos::View.");
  static_assert(Kokkos::is_view<RV>::value,
                "KokkosBlas::Impl::"
                "Dot<1-D>: RV is not a Kokkos::View.");
  static_assert(RV::rank == 0,
                "KokkosBlas::Impl::Dot<1-D>: "
                "RV is not rank 0.");
  static_assert(XV::rank == 1,
                "KokkosBlas::Impl::Dot<1-D>: "
                "XV is not rank 1.");
  static_assert(YV::rank == 1,
                "KokkosBlas::Impl::Dot<1-D>: "
                "YV is not rank 1.");
  static_assert(std::is_same<typename RV::value_type, typename RV::non_const_value_type>::value,
                "KokkosBlas::Dot<1D>: R is const.  "
                "It must be nonconst, because it is an output argument "
                "(we have to be able to write to its entries).");

  typedef typename YV::size_type size_type;
  typedef typename RV::non_const_value_type dot_type;
  typedef typename DotAccumulatingScalar<dot_type>::type special_result_type;

  // This is the same View type as RV, but using the special accumulator as the
  // value type
  typedef Kokkos::View<special_result_type, typename RV::array_layout, typename RV::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      RV_Result;

  static void dot(const execution_space& space, const RV& R, const XV& X, const YV& Y) {
    Kokkos::Profiling::pushRegion(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY ? "KokkosBlas::dot[ETI]"
                                                                     : "KokkosBlas::dot[noETI]");
#ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
    if (KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
      printf("KokkosBlas::dot<> ETI specialization for < %s , %s >\n", typeid(XV).name(), typeid(YV).name());
    else {
      printf("KokkosBlas::dot<> non-ETI specialization for < %s , %s >\n", typeid(XV).name(), typeid(YV).name());
    }
#endif
    const size_type numElems = X.extent(0);

    if (numElems < static_cast<size_type>(INT_MAX)) {
      typedef int index_type;
      DotFunctor<execution_space, RV, XV, YV, index_type> f(X, Y);
      f.run("KokkosBlas::dot<1D>", space, R);
    } else {
      typedef int64_t index_type;
      DotFunctor<execution_space, RV, XV, YV, index_type> f(X, Y);
      f.run("KokkosBlas::dot<1D>", space, R);
    }
    Kokkos::Profiling::popRegion();
  }
};

// Implementation that has the same template args as Dot, but which internally
// uses DotAccumulatingScalar for the result view.
//
// Is never supported by TPLs, but uses the same dot_eti_spec_avail::value.
template <class execution_space, class RV, class XV, class YV>
struct DotSpecialAccumulator<execution_space, RV, XV, YV, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  static_assert(Kokkos::is_view<XV>::value,
                "KokkosBlas::Impl::"
                "DotSpecialAccumulator: XV is not a Kokkos::View.");
  static_assert(Kokkos::is_view<YV>::value,
                "KokkosBlas::Impl::"
                "DotSpecialAccumulator: YV is not a Kokkos::View.");
  static_assert(static_cast<int>(XV::rank) == static_cast<int>(YV::rank),
                "KokkosBlas::Impl::"
                "DotSpecialAccumulator: X and Y have different ranks.");
  static_assert(XV::rank == 1,
                "KokkosBlas::Impl::"
                "DotSpecialAccumulator: X and Y are not rank-1 Views.");
  static_assert(Kokkos::is_view<RV>::value,
                "KokkosBlas::Impl::"
                "DotSpecialAccumulator: RV is not a Kokkos::View.");
  static_assert(std::is_same<typename XV::non_const_value_type, typename YV::non_const_value_type>::value,
                "KokkosBlas::Impl::DotSpecialAccumulator: X and Y have "
                "different scalar types.");
  static_assert(std::is_same<typename RV::value_type, typename RV::non_const_value_type>::value,
                "KokkosBlas::Dot<1D>: R is const.  "
                "It must be nonconst, because it is an output argument "
                "(we have to be able to write to its entries).");

  using size_type  = typename YV::size_type;
  using dot_type   = typename Kokkos::Details::InnerProductSpaceTraits<typename XV::non_const_value_type>::dot_type;
  using accum_type = typename DotAccumulatingScalar<dot_type>::type;
  // This is the same View type as RV, but using the special accumulator as the
  // value type
  using RV_Result = Kokkos::View<accum_type, typename RV::array_layout, typename RV::device_type,
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  static void dot(const execution_space& space, const RV_Result& R, const XV& X, const YV& Y) {
    Kokkos::Profiling::pushRegion(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY ? "KokkosBlas::dot[ETI]"
                                                                     : "KokkosBlas::dot[noETI]");
#ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
    if (KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
      printf("KokkosBlas::dot<> ETI specialization for < %s , %s >\n", typeid(XV).name(), typeid(YV).name());
    else {
      printf("KokkosBlas::dot<> non-ETI specialization for < %s , %s >\n", typeid(XV).name(), typeid(YV).name());
    }
#endif
    const size_type numElems = X.extent(0);

    if (numElems < static_cast<size_type>(INT_MAX)) {
      typedef int index_type;
      DotFunctor<execution_space, RV_Result, XV, YV, index_type> f(X, Y);
      f.run("KokkosBlas::dot<1D>", space, R);
    } else {
      typedef int64_t index_type;
      DotFunctor<execution_space, RV_Result, XV, YV, index_type> f(X, Y);
      f.run("KokkosBlas::dot<1D>", space, R);
    }
    Kokkos::Profiling::popRegion();
  }
};

template <class execution_space, class RV, class XV, class YV, int X_Rank, int Y_Rank>
struct Dot<execution_space, RV, XV, YV, X_Rank, Y_Rank, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  static_assert(Kokkos::is_view<XV>::value,
                "KokkosBlas::Impl::"
                "Dot<2-D>: XV is not a Kokkos::View.");
  static_assert(Kokkos::is_view<YV>::value,
                "KokkosBlas::Impl::"
                "Dot<2-D>: YV is not a Kokkos::View.");
  static_assert(RV::rank == 1,
                "KokkosBlas::Impl::Dot<2-D>: "
                "RV is not rank 1.");

  typedef typename YV::size_type size_type;

  // Helper to get the first column of a rank-1 or rank-2 view.
  // This makes it easier to add a path for single-column dot.
  template <typename V>
  static auto getFirstColumn(const V& v, typename std::enable_if<V::rank == 2>::type* = nullptr) {
    return Kokkos::subview(v, Kokkos::ALL(), 0);
  }

  template <typename V>
  static V getFirstColumn(const V& v, typename std::enable_if<V::rank == 1>::type* = nullptr) {
    return v;
  }

  static void dot(const execution_space& space, const RV& R, const XV& X, const YV& Y) {
    Kokkos::Profiling::pushRegion(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY ? "KokkosBlas::dot[ETI]"
                                                                     : "KokkosBlas::dot[noETI]");
#ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
    if (KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
      printf("KokkosBlas1::dot<> ETI specialization for < %s , %s , %s >\n", typeid(RV).name(), typeid(XV).name(),
             typeid(YV).name());
    else {
      printf("KokkosBlas1::dot<> non-ETI specialization for < %s , %s , %s >\n", typeid(RV).name(), typeid(XV).name(),
             typeid(YV).name());
    }
#endif

    const size_type numRows = X.extent(0);
    const size_type numDots = std::max(X.extent(1), Y.extent(1));
    if (numDots == Kokkos::ArithTraits<size_type>::one()) {
      auto R0 = Kokkos::subview(R, 0);
      auto X0 = getFirstColumn(X);
      auto Y0 = getFirstColumn(Y);
      if (numRows < static_cast<size_type>(INT_MAX)) {
        typedef int index_type;
        DotFunctor<execution_space, decltype(R0), decltype(X0), decltype(Y0), index_type> f(X0, Y0);
        f.run("KokkosBlas::dot<1D>", space, R0);
      } else {
        typedef int64_t index_type;
        DotFunctor<execution_space, decltype(R0), decltype(X0), decltype(Y0), index_type> f(X0, Y0);
        f.run("KokkosBlas::dot<1D>", space, R0);
      }
    } else {
      if (numRows < static_cast<size_type>(INT_MAX) && numRows * numDots < static_cast<size_type>(INT_MAX)) {
        typedef int index_type;
        MV_Dot_Invoke<execution_space, RV, XV, YV, index_type>(space, R, X, Y);
      } else {
        typedef std::int64_t index_type;
        MV_Dot_Invoke<execution_space, RV, XV, YV, index_type>(space, R, X, Y);
      }
    }
    Kokkos::Profiling::popRegion();
  }
};
#endif

}  // namespace Impl
}  // namespace KokkosBlas

//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::Dot for rank == 1.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _DEF macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS1_DOT_ETI_SPEC_DECL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                                           \
  extern template struct Dot<EXEC_SPACE,                                                                               \
                             Kokkos::View<SCALAR, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
                             Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                \
                                          Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                    \
                             Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                \
                                          Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                    \
                             1, 1, false, true>;                                                                       \
  extern template struct Dot<                                                                                          \
      EXEC_SPACE,                                                                                                      \
      Kokkos::View<SCALAR, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,    \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      1, 1, false, true>;                                                                                              \
  extern template struct DotSpecialAccumulator<                                                                        \
      EXEC_SPACE,                                                                                                      \
      Kokkos::View<SCALAR, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,    \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      true>;                                                                                                           \
  extern template struct DotSpecialAccumulator<                                                                        \
      EXEC_SPACE, Kokkos::View<SCALAR, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,            \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      true>;

#define KOKKOSBLAS1_DOT_ETI_SPEC_INST(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                                        \
  template struct Dot<EXEC_SPACE,                                                                                   \
                      Kokkos::View<SCALAR, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,     \
                      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                    \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                        \
                      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                    \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                        \
                      1, 1, false, true>;                                                                           \
  template struct Dot<                                                                                              \
      EXEC_SPACE,                                                                                                   \
      Kokkos::View<SCALAR, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                    \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                        \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                    \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                        \
      1, 1, false, true>;                                                                                           \
  template struct DotSpecialAccumulator<                                                                            \
      EXEC_SPACE, Kokkos::View<SCALAR, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,         \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                    \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                        \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                    \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                        \
      true>;                                                                                                        \
  template struct DotSpecialAccumulator<                                                                            \
      EXEC_SPACE,                                                                                                   \
      Kokkos::View<SCALAR, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                    \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                        \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                    \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                        \
      true>;

//
//
// Macro for definition of full specialization of
// KokkosBlas::Impl::Dot for rank == 2.  This is NOT for users!!!  We
// use this macro in one or more .cpp files in this directory.
//
#define KOKKOSBLAS1_DOT_MV_ETI_SPEC_DECL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                           \
  extern template struct Dot<                                                                             \
      EXEC_SPACE,                                                                                         \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                              \
      Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                         \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                              \
      Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                         \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                              \
      2, 2, false, true>;                                                                                 \
  extern template struct Dot<                                                                             \
      EXEC_SPACE,                                                                                         \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                              \
      Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                         \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                              \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                          \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                              \
      2, 1, false, true>;                                                                                 \
  extern template struct Dot<                                                                             \
      EXEC_SPACE,                                                                                         \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                              \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                          \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                              \
      Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                         \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                              \
      1, 2, false, true>;

#define KOKKOSBLAS1_DOT_MV_ETI_SPEC_INST(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                           \
  template struct Dot<                                                                                    \
      EXEC_SPACE,                                                                                         \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                              \
      Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                         \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                              \
      Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                         \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                              \
      2, 2, false, true>;                                                                                 \
  template struct Dot<                                                                                    \
      EXEC_SPACE,                                                                                         \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                              \
      Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                         \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                              \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                          \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                              \
      2, 1, false, true>;                                                                                 \
  template struct Dot<                                                                                    \
      EXEC_SPACE,                                                                                         \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                              \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                          \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                              \
      Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                         \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                              \
      1, 2, false, true>;

#include <KokkosBlas1_dot_tpl_spec_decl.hpp>

#endif  // KOKKOS_BLAS1_MV_IMPL_DOT_HPP_
