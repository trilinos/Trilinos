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
#ifndef KOKKOS_BLAS1_AXPBY_UNIFICATION_ATTEMPT_TRAITS_HPP_
#define KOKKOS_BLAS1_AXPBY_UNIFICATION_ATTEMPT_TRAITS_HPP_

#include <KokkosKernels_helpers.hpp>
#include <KokkosKernels_ExecSpaceUtils.hpp>
#include <sstream>

namespace KokkosBlas {
namespace Impl {

// --------------------------------

template <class T>
constexpr int typeRank() {
  if constexpr (Kokkos::is_view_v<T>) {
    return T::rank;
  }
  return -1;
}

// --------------------------------

template <class T>
constexpr typename std::enable_if<Kokkos::is_view_v<T>, bool>::type Tr0_val() {
  return (T::rank == 0);
}

template <class T>
constexpr typename std::enable_if<!Kokkos::is_view_v<T>, bool>::type Tr0_val() {
  return false;
}

// --------------------------------

template <class T>
constexpr typename std::enable_if<Kokkos::is_view_v<T>, bool>::type Tr1s_val() {
  return (T::rank == 1) && (T::rank_dynamic == 0);
}

template <class T>
constexpr typename std::enable_if<!Kokkos::is_view_v<T>, bool>::type Tr1s_val() {
  return false;
}

// --------------------------------

template <class T>
constexpr typename std::enable_if<Kokkos::is_view_v<T>, bool>::type Tr1d_val() {
  return (T::rank == 1) && (T::rank_dynamic == 1);
}

template <class T>
constexpr typename std::enable_if<!Kokkos::is_view_v<T>, bool>::type Tr1d_val() {
  return false;
}

// --------------------------------

template <typename T, bool Enable = false>
struct getScalarTypeFromView {
  using type = void;
};

template <typename T>
struct getScalarTypeFromView<T, true> {
  using type = typename T::value_type;
};

// --------------------------------

template <typename T, bool Enable = false>
struct getLayoutFromView {
  using type = void;
};

template <typename T>
struct getLayoutFromView<T, true> {
  using type = typename T::array_layout;
};

// --------------------------------

template <class tExecSpace, class AV, class XMV, class BV, class YMV>
struct AxpbyUnificationAttemptTraits {
  // ********************************************************************
  // Terminology:
  // - variable names begin with lower case letters
  // - type names begin with upper case letters
  // ********************************************************************
 public:
  static constexpr bool onDevice = KokkosKernels::Impl::kk_is_gpu_exec_space<tExecSpace>();

 private:
  static constexpr bool onHost = !onDevice;

 public:
  static constexpr bool a_is_scalar = !Kokkos::is_view_v<AV>;

 private:
  static constexpr bool a_is_r0  = Tr0_val<AV>();
  static constexpr bool a_is_r1s = Tr1s_val<AV>();
  static constexpr bool a_is_r1d = Tr1d_val<AV>();

  static constexpr bool x_is_r1 = Kokkos::is_view_v<XMV> && (XMV::rank == 1);
  static constexpr bool x_is_r2 = Kokkos::is_view_v<XMV> && (XMV::rank == 2);

 public:
  static constexpr bool b_is_scalar = !Kokkos::is_view_v<BV>;

 private:
  static constexpr bool b_is_r0  = Tr0_val<BV>();
  static constexpr bool b_is_r1s = Tr1s_val<BV>();
  static constexpr bool b_is_r1d = Tr1d_val<BV>();

  static constexpr bool y_is_r1 = Kokkos::is_view_v<YMV> && (YMV::rank == 1);
  static constexpr bool y_is_r2 = Kokkos::is_view_v<YMV> && (YMV::rank == 2);

  static constexpr bool xyRank1Case = x_is_r1 && y_is_r1;
  static constexpr bool xyRank2Case = x_is_r2 && y_is_r2;

  // ********************************************************************
  // Declare 'AtInputScalarTypeA_nonConst'
  // ********************************************************************
  using ScalarTypeA2_onDevice = typename getScalarTypeFromView<AV, (a_is_r0 || a_is_r1s || a_is_r1d) && onDevice>::type;
  using ScalarTypeA1_onDevice = std::conditional_t<a_is_scalar && onDevice, AV, ScalarTypeA2_onDevice>;

  using ScalarTypeA2_onHost = typename getScalarTypeFromView<AV, (a_is_r0 || a_is_r1s || a_is_r1d) && onHost>::type;
  using ScalarTypeA1_onHost = std::conditional_t<a_is_scalar && onHost, AV, ScalarTypeA2_onHost>;

  using AtInputScalarTypeA = std::conditional_t<onHost, ScalarTypeA1_onHost, ScalarTypeA1_onDevice>;

  using AtInputScalarTypeA_nonConst = typename std::remove_const<AtInputScalarTypeA>::type;

  // ********************************************************************
  // Declare 'AtInputScalarTypeX_nonConst'
  // ********************************************************************
  using AtInputScalarTypeX = typename XMV::value_type;

  using AtInputScalarTypeX_nonConst = typename XMV::non_const_value_type;

  // ********************************************************************
  // Declare 'AtInputScalarTypeB_nonConst'
  // ********************************************************************
  using ScalarTypeB2_onDevice = typename getScalarTypeFromView<BV, (b_is_r0 || b_is_r1s || b_is_r1d) && onDevice>::type;
  using ScalarTypeB1_onDevice = std::conditional_t<b_is_scalar && onDevice, BV, ScalarTypeB2_onDevice>;

  using ScalarTypeB2_onHost = typename getScalarTypeFromView<BV, (b_is_r0 || b_is_r1s || b_is_r1d) && onHost>::type;
  using ScalarTypeB1_onHost = std::conditional_t<b_is_scalar && onHost, BV, ScalarTypeB2_onHost>;

  using AtInputScalarTypeB = std::conditional_t<onHost, ScalarTypeB1_onHost, ScalarTypeB1_onDevice>;

  using AtInputScalarTypeB_nonConst = typename std::remove_const<AtInputScalarTypeB>::type;

  // ********************************************************************
  // Declare 'AtInputScalarTypeY_nonConst'
  // ********************************************************************
  using AtInputScalarTypeY = typename YMV::value_type;

  using AtInputScalarTypeY_nonConst = typename YMV::non_const_value_type;

  // ********************************************************************
  // Declare 'InternalLayoutX' and 'InternalLayoutY'
  // ********************************************************************
  using InternalLayoutX = typename KokkosKernels::Impl::GetUnifiedLayout<XMV>::array_layout;
  using InternalLayoutY = typename KokkosKernels::Impl::GetUnifiedLayoutPreferring<YMV, InternalLayoutX>::array_layout;

  // ********************************************************************
  // Declare 'InternalTypeA_tmp'
  // ********************************************************************
  using AtInputLayoutA = typename getLayoutFromView<AV, (a_is_r0 || a_is_r1s || a_is_r1d)>::type;

 public:
  static constexpr bool atInputLayoutA_isStride = std::is_same_v<AtInputLayoutA, Kokkos::LayoutStride>;

 private:
  using InternalLayoutA =
      std::conditional_t<(a_is_r1d || a_is_r1s) && atInputLayoutA_isStride, AtInputLayoutA, InternalLayoutX>;

  static constexpr bool atInputScalarTypeA_mustRemain = Kokkos::ArithTraits<AtInputScalarTypeA_nonConst>::is_complex &&
                                                        !Kokkos::ArithTraits<AtInputScalarTypeX_nonConst>::is_complex;

  using InternalScalarTypeA =
      std::conditional_t<atInputScalarTypeA_mustRemain || ((a_is_r1d || a_is_r1s) && xyRank2Case),
                         AtInputScalarTypeA_nonConst  // Yes, keep the input scalar type
                         ,
                         AtInputScalarTypeX_nonConst  // Yes, instead of
                                                      // 'AtInputScalarTypeA_nonConst'
                         >;

  using InternalTypeA_onDevice =
      std::conditional_t<a_is_scalar && b_is_scalar && onDevice,  // Keep 'a' as scalar
                         InternalScalarTypeA,
                         Kokkos::View<const InternalScalarTypeA*, InternalLayoutA, typename XMV::device_type,
                                      Kokkos::MemoryTraits<Kokkos::Unmanaged>>>;

  using InternalTypeA_onHost =
      std::conditional_t<(a_is_r1d || a_is_r1s) && xyRank2Case && onHost,
                         Kokkos::View<const InternalScalarTypeA*, InternalLayoutA, typename XMV::device_type,
                                      Kokkos::MemoryTraits<Kokkos::Unmanaged>>,
                         InternalScalarTypeA>;

  using InternalTypeA_tmp = std::conditional_t<onHost, InternalTypeA_onHost, InternalTypeA_onDevice>;

  // ********************************************************************
  // Declare 'InternalTypeX'
  // ********************************************************************
 public:
  using InternalTypeX =
      std::conditional_t<x_is_r2,
                         Kokkos::View<const AtInputScalarTypeX_nonConst**, InternalLayoutX, typename XMV::device_type,
                                      Kokkos::MemoryTraits<Kokkos::Unmanaged>>,
                         Kokkos::View<const AtInputScalarTypeX_nonConst*, InternalLayoutX, typename XMV::device_type,
                                      Kokkos::MemoryTraits<Kokkos::Unmanaged>>>;

  // ********************************************************************
  // Declare 'InternalTypeB_tmp'
  // ********************************************************************
 private:
  using AtInputLayoutB = typename getLayoutFromView<BV, (b_is_r0 || b_is_r1s || b_is_r1d)>::type;

 public:
  static constexpr bool atInputLayoutB_isStride = std::is_same_v<AtInputLayoutB, Kokkos::LayoutStride>;

 private:
  using InternalLayoutB =
      std::conditional_t<(b_is_r1d || b_is_r1s) && atInputLayoutB_isStride, AtInputLayoutB, InternalLayoutY>;

  static constexpr bool atInputScalarTypeB_mustRemain = Kokkos::ArithTraits<AtInputScalarTypeB_nonConst>::is_complex &&
                                                        !Kokkos::ArithTraits<AtInputScalarTypeY_nonConst>::is_complex;

  using InternalScalarTypeB =
      std::conditional_t<atInputScalarTypeB_mustRemain || ((b_is_r1d || b_is_r1s) && xyRank2Case),
                         AtInputScalarTypeB_nonConst  // Yes, keep the input scalar type
                         ,
                         AtInputScalarTypeY_nonConst  // Yes, instead of
                                                      // 'AtInputScalarTypeB_nonConst'
                         >;

  using InternalTypeB_onDevice =
      std::conditional_t<a_is_scalar && b_is_scalar && onDevice,  // Keep 'b' as scalar
                         InternalScalarTypeB,
                         Kokkos::View<const InternalScalarTypeB*, InternalLayoutB, typename YMV::device_type,
                                      Kokkos::MemoryTraits<Kokkos::Unmanaged>>>;

  using InternalTypeB_onHost =
      std::conditional_t<(b_is_r1d || b_is_r1s) && xyRank2Case && onHost,
                         Kokkos::View<const InternalScalarTypeB*, InternalLayoutB, typename YMV::device_type,
                                      Kokkos::MemoryTraits<Kokkos::Unmanaged>>,
                         InternalScalarTypeB>;

  using InternalTypeB_tmp = std::conditional_t<onHost, InternalTypeB_onHost, InternalTypeB_onDevice>;

  // ********************************************************************
  // Declare 'InternalTypeY'
  // ********************************************************************
 public:
  using InternalTypeY =
      std::conditional_t<y_is_r2,
                         Kokkos::View<AtInputScalarTypeY_nonConst**, InternalLayoutY, typename YMV::device_type,
                                      Kokkos::MemoryTraits<Kokkos::Unmanaged>>,
                         Kokkos::View<AtInputScalarTypeY_nonConst*, InternalLayoutY, typename YMV::device_type,
                                      Kokkos::MemoryTraits<Kokkos::Unmanaged>>>;

  // ********************************************************************
  // Declare 'InternalTypeA': if 'InternalTypeB_tmp' is a view then
  // make sure 'InternalTypeA' is a view as well
  // ********************************************************************
  using InternalTypeA =
      std::conditional_t<!Kokkos::is_view_v<InternalTypeA_tmp> && Kokkos::is_view_v<InternalTypeB_tmp>,
                         Kokkos::View<const InternalScalarTypeA*, InternalLayoutA, typename XMV::device_type,
                                      Kokkos::MemoryTraits<Kokkos::Unmanaged>>,
                         InternalTypeA_tmp>;

  // ********************************************************************
  // Declare 'InternalTypeA_managed' with the same scalar type in
  // 'InternalTypeA'
  // ********************************************************************
 private:
  using InternalLayoutA_managed = InternalLayoutA;

 public:
  using InternalTypeA_managed =
      std::conditional_t<Kokkos::is_view_v<InternalTypeA>,
                         Kokkos::View<InternalScalarTypeA*, InternalLayoutA_managed, typename XMV::device_type>, void>;

  // ********************************************************************
  // Declare 'InternalTypeB' if 'InternalTypeA_tmp' is a view then
  // make sure 'InternalTypeB' is a view as well
  // ********************************************************************
  using InternalTypeB =
      std::conditional_t<Kokkos::is_view_v<InternalTypeA_tmp> && !Kokkos::is_view_v<InternalTypeB_tmp>,
                         Kokkos::View<const InternalScalarTypeB*, InternalLayoutB, typename YMV::device_type,
                                      Kokkos::MemoryTraits<Kokkos::Unmanaged>>,
                         InternalTypeB_tmp>;

  // ********************************************************************
  // Declare 'InternalTypeB_managed' with the same scalar type in
  // 'InternalTypeB'
  // ********************************************************************
 private:
  using InternalLayoutB_managed = InternalLayoutB;

 public:
  using InternalTypeB_managed =
      std::conditional_t<Kokkos::is_view_v<InternalTypeB>,
                         Kokkos::View<InternalScalarTypeB*, InternalLayoutB_managed, typename YMV::device_type>, void>;

  // ********************************************************************
  // Auxiliary Boolean results on internal types
  // ********************************************************************
 private:
  static constexpr bool internalTypeA_is_scalar = !Kokkos::is_view_v<InternalTypeA>;
  static constexpr bool internalTypeA_is_r1d    = Tr1d_val<InternalTypeA>();

  static constexpr bool internalTypeB_is_scalar = !Kokkos::is_view_v<InternalTypeB>;
  static constexpr bool internalTypeB_is_r1d    = Tr1d_val<InternalTypeB>();

 public:
  static constexpr bool internalTypesAB_bothScalars = (internalTypeA_is_scalar && internalTypeB_is_scalar);
  static constexpr bool internalTypesAB_bothViews   = (internalTypeA_is_r1d && internalTypeB_is_r1d);

  // ********************************************************************
  // Routine to perform checks (both compile time and run time)
  // ********************************************************************
  static void performChecks(const AV& a, const XMV& X, const BV& b, const YMV& Y) {
    // ******************************************************************
    // Check 1/6: General checks
    // ******************************************************************
    static_assert(Kokkos::is_execution_space_v<tExecSpace>,
                  "KokkosBlas::Impl::AxpbyUnificationAttemptTraits::performChecks()"
                  ": tExecSpace must be a valid Kokkos execution space.");

    static_assert((xyRank1Case && !xyRank2Case) || (!xyRank1Case && xyRank2Case),
                  "KokkosBlas::Impl::AxpbyUnificationAttemptTraits::performChecks()"
                  ": one must have either both X and Y as rank 1, or both X and Y as "
                  "rank 2");

    if constexpr (!Kokkos::ArithTraits<AtInputScalarTypeY_nonConst>::is_complex) {
      static_assert((!Kokkos::ArithTraits<AtInputScalarTypeA_nonConst>::is_complex) &&
                        (!Kokkos::ArithTraits<AtInputScalarTypeX_nonConst>::is_complex) &&
                        (!Kokkos::ArithTraits<AtInputScalarTypeB_nonConst>::is_complex),
                    "KokkosBlas::Impl::AxpbyUnificationAttemptTraits::performChecks()"
                    ": if Y is not complex, then A, X and B cannot be complex");
    }

    // ******************************************************************
    // Check 2/6: YMV is valid
    // ******************************************************************
    static_assert(Kokkos::is_view<YMV>::value,
                  "KokkosBlas::Impl::AxpbyUnificationAttemptTraits::performChecks()"
                  ": Y is not a Kokkos::View.");
    static_assert(std::is_same<typename YMV::value_type, typename YMV::non_const_value_type>::value,
                  "KokkosBlas::Impl::AxpbyUnificationAttemptTraits::performChecks()"
                  ": Y is const.  It must be nonconst, "
                  "because it is an output argument "
                  "(we must be able to write to its entries).");
    static_assert(Kokkos::SpaceAccessibility<tExecSpace, typename YMV::memory_space>::accessible,
                  "KokkosBlas::Impl::AxpbyUnificationAttemptTraits::performChecks()"
                  ": XMV must be accessible from tExecSpace");

    // ******************************************************************
    // Check 3/6: XMV is valid
    // ******************************************************************
    static_assert(Kokkos::is_view<XMV>::value,
                  "KokkosBlas::Impl::AxpbyUnificationAttemptTraits::performChecks()"
                  ": X is not a Kokkos::View.");
    static_assert(Kokkos::SpaceAccessibility<tExecSpace, typename XMV::memory_space>::accessible,
                  "KokkosBlas::Impl::AxpbyUnificationAttemptTraits::performChecks()"
                  ": XMV must be accessible from tExecSpace");

    if constexpr (xyRank1Case) {
      if (X.extent(0) != Y.extent(0)) {
        std::ostringstream msg;
        msg << "KokkosBlas::Impl::AxpbyUnificationAttemptTraits::performChecks("
               ")"
            << ", invalid rank-1 X extent"
            << ": X.extent(0) = " << X.extent(0) << ", Y.extent(0) = " << Y.extent(0);
        KokkosKernels::Impl::throw_runtime_exception(msg.str());
      }
    } else {
      if ((X.extent(0) != Y.extent(0)) || (X.extent(1) != Y.extent(1))) {
        std::ostringstream msg;
        msg << "KokkosBlas::Impl::AxpbyUnificationAttemptTraits::performChecks("
               ")"
            << ", invalid rank-2 X extents"
            << ": X.extent(0) = " << X.extent(0) << ", X.extent(1) = " << X.extent(1)
            << ", Y.extent(0) = " << Y.extent(0) << ", Y.extent(1) = " << Y.extent(1);
        KokkosKernels::Impl::throw_runtime_exception(msg.str());
      }
    }

    // ******************************************************************
    // Check 4/6: AV is valid
    // ******************************************************************
    static_assert(
        (a_is_scalar && !a_is_r0 && !a_is_r1s && !a_is_r1d) || (!a_is_scalar && a_is_r0 && !a_is_r1s && !a_is_r1d) ||
            (!a_is_scalar && !a_is_r0 && a_is_r1s && !a_is_r1d) || (!a_is_scalar && !a_is_r0 && !a_is_r1s && a_is_r1d),
        "KokkosBlas::Impl::AxpbyUnificationAttemptTraits::performChecks()"
        ": 'a' must be either scalar or rank 0 or rank 1 static or rank 1 "
        "dynamic");

    if constexpr (a_is_r1d || a_is_r1s) {
      if constexpr (xyRank1Case) {
        if (a.extent(0) != 1) {
          std::ostringstream msg;
          msg << "KokkosBlas::Impl::AxpbyUnificationAttemptTraits::"
                 "performChecks()"
              << ": view 'a' must have extent(0) == 1 for xyRank1Case"
              << ", a.extent(0) = " << a.extent(0);
          KokkosKernels::Impl::throw_runtime_exception(msg.str());
        }
      } else {
        if ((a.extent(0) == 1) || (a.extent(0) == Y.extent(1))) {  // Yes, 'Y' is the reference
          // Ok
        } else {
          std::ostringstream msg;
          msg << "KokkosBlas::Impl::AxpbyUnificationAttemptTraits::"
                 "performChecks()"
              << ": view 'a' must have extent(0) == 1 or Y.extent(1) for "
                 "xyRank2Case"
              << ", a.extent(0) = " << a.extent(0) << ", Y.extent(0) = " << Y.extent(0)
              << ", Y.extent(1) = " << Y.extent(1);
          KokkosKernels::Impl::throw_runtime_exception(msg.str());
        }
      }  // if (rank1Case) else
    }    // if a_is_r1d

    // ******************************************************************
    // Check 5/6: BV is valid
    // ******************************************************************
    static_assert(
        (b_is_scalar && !b_is_r0 && !b_is_r1s && !b_is_r1d) || (!b_is_scalar && b_is_r0 && !b_is_r1s && !b_is_r1d) ||
            (!b_is_scalar && !b_is_r0 && b_is_r1s && !b_is_r1d) || (!b_is_scalar && !b_is_r0 && !b_is_r1s && b_is_r1d),
        "KokkosBlas::Impl::AxpbyUnificationAttemptTraits::performChecks()"
        ": 'b' must be either scalar or rank 0 or rank 1 static or rank 1 "
        "dynamic");

    if constexpr (b_is_r1d || b_is_r1s) {
      if constexpr (xyRank1Case) {
        if (b.extent(0) != 1) {
          std::ostringstream msg;
          msg << "KokkosBlas::Impl::AxpbyUnificationAttemptTraits::"
                 "performChecks()"
              << ": view 'b' must have extent(0) == 1 for xyRank1Case"
              << ", b.extent(0) = " << b.extent(0);
          KokkosKernels::Impl::throw_runtime_exception(msg.str());
        }
      } else {
        if ((b.extent(0) == 1) || (b.extent(0) == Y.extent(1))) {
          // Ok
        } else {
          std::ostringstream msg;
          msg << "KokkosBlas::Impl::AxpbyUnificationAttemptTraits::"
                 "performChecks()"
              << ": view 'b' must have extent(0) == 1 or Y.extent(1) for "
                 "xyRank2Case"
              << ", b.extent(0) = " << b.extent(0) << ", Y.extent(0) = " << Y.extent(0)
              << ", Y.extent(1) = " << Y.extent(1);
          KokkosKernels::Impl::throw_runtime_exception(msg.str());
        }
      }  // if (rank1Case) else
    }    // if b_is_r1d

    // ******************************************************************
    // Check 6/6: Checks on InternalTypeA, X, B, Y
    // ******************************************************************
    if constexpr (onHost) {
      if constexpr (xyRank1Case) {
        constexpr bool internalTypeA_isOk = (internalTypeA_is_scalar || internalTypeA_is_r1d);
        static_assert(internalTypeA_isOk,
                      "KokkosBlas::Impl::AxpbyUnificationAttemptTraits::performChecks()"
                      ", onHost, xyRank1Case: InternalTypeA is wrong");

        constexpr bool internalTypeX_isOk =
            std::is_same_v<InternalTypeX,
                           Kokkos::View<const AtInputScalarTypeX_nonConst*, InternalLayoutX, typename XMV::device_type,
                                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>>;
        static_assert(internalTypeX_isOk,
                      "KokkosBlas::Impl::AxpbyUnificationAttemptTraits::performChecks()"
                      ", onHost, xyRank1Case: InternalTypeX is wrong");

        constexpr bool internalTypeB_isOk = (internalTypeB_is_scalar || internalTypeB_is_r1d);
        static_assert(internalTypeB_isOk,
                      "KokkosBlas::Impl::AxpbyUnificationAttemptTraits::performChecks()"
                      ", onHost, xyRank1Case: InternalTypeB is wrong");

        constexpr bool internalTypeY_isOk =
            std::is_same_v<InternalTypeY,
                           Kokkos::View<AtInputScalarTypeY_nonConst*, InternalLayoutY, typename YMV::device_type,
                                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>>;
        static_assert(internalTypeY_isOk,
                      "KokkosBlas::Impl::AxpbyUnificationAttemptTraits::performChecks()"
                      ", onHost, xyRank1Case: InternalTypeY is wrong");
      } else {
        constexpr bool internalTypeA_isOk = (internalTypeA_is_scalar || internalTypeA_is_r1d);
        static_assert(internalTypeA_isOk,
                      "KokkosBlas::Impl::AxpbyUnificationAttemptTraits::performChecks()"
                      ", onHost, xyRank2Case: InternalTypeA is wrong");

        constexpr bool internalTypeX_isOk =
            std::is_same_v<InternalTypeX,
                           Kokkos::View<const AtInputScalarTypeX_nonConst**, InternalLayoutX, typename XMV::device_type,
                                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>>;
        static_assert(internalTypeX_isOk,
                      "KokkosBlas::Impl::AxpbyUnificationAttemptTraits::performChecks()"
                      ", onHost, xyRank2Case: InternalTypeX is wrong");

        constexpr bool internalTypeB_isOk = (internalTypeB_is_scalar || internalTypeB_is_r1d);
        static_assert(internalTypeB_isOk,
                      "KokkosBlas::Impl::AxpbyUnificationAttemptTraits::performChecks()"
                      ", onHost, xyRank2Case: InternalTypeB is wrong");

        constexpr bool internalTypeY_isOk =
            std::is_same_v<InternalTypeY,
                           Kokkos::View<AtInputScalarTypeY_nonConst**, InternalLayoutY, typename YMV::device_type,
                                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>>;
        static_assert(internalTypeY_isOk,
                      "KokkosBlas::Impl::AxpbyUnificationAttemptTraits::performChecks()"
                      ", onHost, xyRank2Case: InternalTypeY is wrong");
      }
    } else {
      if constexpr (xyRank1Case) {
        constexpr bool internalTypeA_isOk =
            internalTypeA_is_r1d || (a_is_scalar && b_is_scalar && internalTypeA_is_scalar);
        static_assert(internalTypeA_isOk,
                      "KokkosBlas::Impl::AxpbyUnificationAttemptTraits::performChecks()"
                      ", onDevice, xyRank1Case: InternalTypeA is wrong");

        constexpr bool internalTypeX_isOk =
            std::is_same_v<InternalTypeX,
                           Kokkos::View<const AtInputScalarTypeX_nonConst*, InternalLayoutX, typename XMV::device_type,
                                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>>;
        static_assert(internalTypeX_isOk,
                      "KokkosBlas::Impl::AxpbyUnificationAttemptTraits::performChecks()"
                      ", onDevice, xyRank1Case: InternalTypeX is wrong");

        constexpr bool internalTypeB_isOk =
            internalTypeB_is_r1d || (a_is_scalar && b_is_scalar && internalTypeA_is_scalar);
        static_assert(internalTypeB_isOk,
                      "KokkosBlas::Impl::AxpbyUnificationAttemptTraits::performChecks()"
                      ", onDevice, xyRank1Case: InternalTypeB is wrong");

        constexpr bool internalTypeY_isOk =
            std::is_same_v<InternalTypeY,
                           Kokkos::View<AtInputScalarTypeY_nonConst*, InternalLayoutY, typename YMV::device_type,
                                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>>;
        static_assert(internalTypeY_isOk,
                      "KokkosBlas::Impl::AxpbyUnificationAttemptTraits::performChecks()"
                      ", onDevice, xyRank1Case: InternalTypeY is wrong");
      } else {
        constexpr bool internalTypeA_isOk =
            internalTypeA_is_r1d || (a_is_scalar && b_is_scalar && internalTypeA_is_scalar);
        static_assert(internalTypeA_isOk,
                      "KokkosBlas::Impl::AxpbyUnificationAttemptTraits::performChecks()"
                      ", onDevice, xyRank2Case: InternalTypeA is wrong");

        constexpr bool internalTypeX_isOk =
            std::is_same_v<InternalTypeX,
                           Kokkos::View<const AtInputScalarTypeX_nonConst**, InternalLayoutX, typename XMV::device_type,
                                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>>;
        static_assert(internalTypeX_isOk,
                      "KokkosBlas::Impl::AxpbyUnificationAttemptTraits::performChecks()"
                      ", onDevice, xyRank2Case: InternalTypeX is wrong");

        constexpr bool internalTypeB_isOk =
            internalTypeB_is_r1d || (a_is_scalar && b_is_scalar && internalTypeB_is_scalar);
        static_assert(internalTypeB_isOk,
                      "KokkosBlas::Impl::AxpbyUnificationAttemptTraits::performChecks()"
                      ", onDevice, xyRank2Case: InternalTypeB is wrong");

        constexpr bool internalTypeY_isOk =
            std::is_same_v<InternalTypeY,
                           Kokkos::View<AtInputScalarTypeY_nonConst**, InternalLayoutY, typename YMV::device_type,
                                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>>;
        static_assert(internalTypeY_isOk,
                      "KokkosBlas::Impl::AxpbyUnificationAttemptTraits::performChecks()"
                      ", onDevice, xyRank2Case: InternalTypeY is wrong");
      }
    }

    if constexpr (onHost) {
      // ****************************************************************
      // We are in the 'onHost' case, with 2 possible subcases::
      //
      // 1) xyRank1Case, with the following possible situations:
      // - [InternalTypeA, B] = [S_a, S_b], or
      // - [InternalTypeA, B] = [view<S_a*,1>, view<S_b*,1>]
      //
      // or
      //
      // 2) xyRank2Case, with the following possible situations:
      // - [InternalTypeA, B] = [S_a, S_b], or
      // - [InternalTypeA, B] = [view<S_a*,1 / m>, view<S_b*,1 / m>]
      // ****************************************************************
      static_assert(internalTypesAB_bothScalars || internalTypesAB_bothViews,
                    "KokkosBlas::Impl::AxpbyUnificationAttemptTraits::performChecks()"
                    ", onHost, invalid combination of types");
    }  // If onHost
    else if constexpr (onDevice) {
      // ****************************************************************
      // We are in the 'onDevice' case, with 2 possible subcases:
      //
      // 1) xyRank1Case, with the following possible situations:
      // - [InternalTypeA, B] = [S_a, S_b], or
      // - [InternalTypeA, B] = [view<S_a*,1>, view<S_b*,1>]
      //
      // or
      //
      // 2) xyRank2Case, with the following possible situations:
      // - [InternalTypeA, B] = [S_a, S_b], or
      // - [InternalTypeA, B] = [view<S_a*,1 / m>, view<S_b*,1 / m>]
      // ****************************************************************
      static_assert(internalTypesAB_bothViews || (a_is_scalar && b_is_scalar && internalTypesAB_bothScalars),
                    "KokkosBlas::Impl::AxpbyUnificationAttemptTraits::performChecks()"
                    ", onDevice, invalid combination of types");
    }

    if constexpr (xyRank2Case && (a_is_r1d || a_is_r1s) && atInputLayoutA_isStride) {
      static_assert(std::is_same_v<typename getLayoutFromView<InternalTypeA, Kokkos::is_view_v<InternalTypeA>>::type,
                                   Kokkos::LayoutStride>,
                    "KokkosBlas::Impl::AxpbyUnificationAttemptTraits::performChecks()"
                    ", xyRank2Case: coeff 'a' is rank-1 and has LayoutStride at input"
                    ", but no LayoutStride internally");
    }

    if constexpr (xyRank2Case && (b_is_r1d || b_is_r1s) && atInputLayoutB_isStride) {
      static_assert(std::is_same_v<typename getLayoutFromView<InternalTypeB, Kokkos::is_view_v<InternalTypeB>>::type,
                                   Kokkos::LayoutStride>,
                    "KokkosBlas::Impl::AxpbyUnificationAttemptTraits::performChecks()"
                    ", xyRank2Case: coeff 'b' is rank-1 and has LayoutStride at input"
                    ", but no LayoutStride internally");
    }
  }  // Constructor

  // ********************************************************************
  // Routine to print information on input variables and internal variables
  // ********************************************************************
#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
  static void printInformation(std::ostream& os, std::string const& headerMsg) {
    os << headerMsg << ": AV = "
       << typeid(AV).name()
       //<< ", AV::const_data_type = "     << typeid(AV::const_data_type).name()
       //<< ", AV::non_const_data_type = " <<
       // typeid(AV::non_const_data_type).name()
       << ", AtInputScalarTypeA = " << typeid(AtInputScalarTypeA).name()
       << ", isConst = " << std::is_const_v<AtInputScalarTypeA> << ", isComplex = "
       << Kokkos::ArithTraits<AtInputScalarTypeA_nonConst>::is_complex
       << ", AtInputScalarTypeA_nonConst = " << typeid(AtInputScalarTypeA_nonConst).name()
       << ", InternalTypeA = " << typeid(InternalTypeA).name() << "\n"
       << ", InternalTypeA_managed = " << typeid(InternalTypeA_managed).name() << "\n"
       << "\n"
       << "XMV = " << typeid(XMV).name() << "\n"
       << "XMV::value_type = " << typeid(typename XMV::value_type).name() << "\n"
       << "XMV::const_data_type = " << typeid(typename XMV::const_data_type).name() << "\n"
       << "XMV::non_const_data_type = " << typeid(typename XMV::non_const_data_type).name() << "\n"
       << "AtInputScalarTypeX = " << typeid(AtInputScalarTypeX).name() << "\n"
       << "isConst = " << std::is_const_v<AtInputScalarTypeX> << "\n"
       << "isComplex = " << Kokkos::ArithTraits<AtInputScalarTypeX_nonConst>::is_complex << "\n"
       << "AtInputScalarTypeX_nonConst = " << typeid(AtInputScalarTypeX_nonConst).name() << "\n"
       << "InternalTypeX = " << typeid(InternalTypeX).name() << "\n"
       << "\n"
       << "BV = "
       << typeid(BV).name()
       //<< ", BV::const_data_type = "     << typeid(BV::const_data_type).name()
       //<< ", BV::non_const_data_type = " <<
       // typeid(BV::non_const_data_type).name()
       << ", AtInputScalarTypeB = " << typeid(AtInputScalarTypeB).name()
       << ", isConst = " << std::is_const_v<AtInputScalarTypeB> << ", isComplex = "
       << Kokkos::ArithTraits<AtInputScalarTypeB_nonConst>::is_complex
       << ", AtInputScalarTypeB_nonConst = " << typeid(AtInputScalarTypeB_nonConst).name()
       << ", InternalTypeB = " << typeid(InternalTypeB).name() << "\n"
       << ", InternalTypeB_managed = " << typeid(InternalTypeB_managed).name() << "\n"
       << "\n"
       << "YMV = " << typeid(YMV).name() << "\n"
       << "YMV::value_type = " << typeid(typename YMV::value_type).name() << "\n"
       << "YMV::const_data_type = " << typeid(typename YMV::const_data_type).name() << "\n"
       << "YMV::non_const_data_type = " << typeid(typename YMV::non_const_data_type).name() << "\n"
       << "AtInputScalarTypeY = " << typeid(AtInputScalarTypeY).name() << "\n"
       << "isConst = " << std::is_const_v<AtInputScalarTypeY> << "\n"
       << "isComplex = " << Kokkos::ArithTraits<AtInputScalarTypeY_nonConst>::is_complex << "\n"
       << "AtInputScalarTypeY_nonConst = " << typeid(AtInputScalarTypeY_nonConst).name() << "\n"
       << "InternalTypeY = " << typeid(InternalTypeY).name() << "\n"
       << std::endl;
  }
#endif

};  // struct AxpbyUnificationAttemptTraits

// --------------------------------

template <typename T, int rankT>
struct getScalarValueFromVariableAtHost {
  getScalarValueFromVariableAtHost() {
    static_assert((rankT == -1) || (rankT == 0) || (rankT == 1), "Generic struct should not have been invoked!");
  }
};

template <typename T>
struct getScalarValueFromVariableAtHost<T, -1> {
  static T getValue(T const& var) { return var; }
};

template <typename T>
struct getScalarValueFromVariableAtHost<T, 0> {
  static typename T::value_type getValue(T const& var) { return var(); }
};

template <class T>
struct getScalarValueFromVariableAtHost<T, 1> {
  static typename T::value_type getValue(T const& var) { return var[0]; }
};

// --------------------------------

template <typename T>
size_t getAmountOfScalarsInCoefficient(T const& coeff) {
  size_t result = 1;
  if constexpr (Kokkos::is_view_v<T>) {
    if constexpr (T::rank == 1) {
      result = coeff.extent(0);
    }
  }
  return result;
}

// --------------------------------

template <typename T>
size_t getStrideInCoefficient(T const& coeff) {
  size_t result = 1;
  if constexpr (Kokkos::is_view_v<T>) {
    if constexpr ((T::rank == 1) && (std::is_same_v<typename T::array_layout, Kokkos::LayoutStride>)) {
      result = coeff.stride_0();
    }
  }
  return result;
}

// --------------------------------

template <class T_in, class T_out>
static void populateRank1Stride1ViewWithScalarOrNonStrideView(T_in const& coeff_in, T_out& coeff_out) {
  // ***********************************************************************
  // 'coeff_out' is assumed to be rank-1, of LayoutLeft or LayoutRight
  //
  // One has to be careful with situations like the following:
  // - a coeff_in that deals with 'double', and
  // - a coeff_out deals with 'complex<double>'
  // ***********************************************************************
  using ScalarOutType = typename std::remove_const<typename T_out::value_type>::type;

  if constexpr (!Kokkos::is_view_v<T_in>) {
    // *********************************************************************
    // 'coeff_in' is scalar
    // *********************************************************************
    ScalarOutType scalarValue(coeff_in);
    Kokkos::deep_copy(coeff_out, scalarValue);
  } else if constexpr (T_in::rank == 0) {
    // *********************************************************************
    // 'coeff_in' is rank-0
    // *********************************************************************
    typename T_in::HostMirror h_coeff_in("h_coeff_in");
    Kokkos::deep_copy(h_coeff_in, coeff_in);
    ScalarOutType scalarValue(h_coeff_in());
    Kokkos::deep_copy(coeff_out, scalarValue);
  } else {
    // *********************************************************************
    // 'coeff_in' is also rank-1
    // *********************************************************************
    if (coeff_out.extent(0) != coeff_in.extent(0)) {
      std::ostringstream msg;
      msg << "In populateRank1Stride1ViewWithScalarOrNonStrideView()"
          << ": 'in' and 'out' should have the same extent(0)"
          << ", T_in = " << typeid(T_in).name() << ", coeff_in.label() = " << coeff_in.label()
          << ", coeff_in.extent(0) = " << coeff_in.extent(0) << ", T_out = " << typeid(T_out).name()
          << ", coeff_out.label() = " << coeff_out.label() << ", coeff_out.extent(0) = " << coeff_out.extent(0);
      KokkosKernels::Impl::throw_runtime_exception(msg.str());
    }

    using ScalarInType = typename std::remove_const<typename T_in::value_type>::type;
    if constexpr (std::is_same_v<ScalarInType, ScalarOutType>) {
      coeff_out = coeff_in;
    } else if (coeff_out.extent(0) == 1) {
      typename T_in::HostMirror h_coeff_in("h_coeff_in");
      Kokkos::deep_copy(h_coeff_in, coeff_in);
      ScalarOutType scalarValue(h_coeff_in[0]);
      Kokkos::deep_copy(coeff_out, scalarValue);
    } else {
      std::ostringstream msg;
      msg << "In populateRank1Stride1ViewWithScalarOrNonStrideView()"
          << ": scalar types 'in' and 'out' should be the same"
          << ", T_in = " << typeid(T_in).name() << ", ScalarInType = " << typeid(ScalarInType).name()
          << ", coeff_in.label() = " << coeff_in.label() << ", coeff_in.extent(0) = " << coeff_in.extent(0)
          << ", T_out = " << typeid(T_out).name() << ", ScalarOutType = " << typeid(ScalarOutType).name()
          << ", coeff_out.label() = " << coeff_out.label() << ", coeff_out.extent(0) = " << coeff_out.extent(0);
      KokkosKernels::Impl::throw_runtime_exception(msg.str());
    }
  }
}  // populateRank1Stride1ViewWithScalarOrNonStrideView()

}  // namespace Impl
}  // namespace KokkosBlas

#endif  // KOKKOS_BLAS1_AXPBY_UNIFICATION_ATTEMPT_TRAITS_HPP_
