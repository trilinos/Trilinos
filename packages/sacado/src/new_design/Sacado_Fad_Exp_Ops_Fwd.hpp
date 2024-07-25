// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_EXP_FAD_OPS_FWD_HPP
#define SACADO_EXP_FAD_OPS_FWD_HPP

#include <type_traits>

namespace Sacado {

  template <typename T> struct IsSimdType;
  template <typename T> struct IsFad;
  template <typename T> struct ValueType;

  namespace Fad {
  namespace Exp {

    template <typename T, typename E> class UnaryPlusOp;
    template <typename T, typename E> class UnaryMinusOp;
    template <typename T, typename E> class ExpOp;
    template <typename T, typename E> class LogOp;
    template <typename T, typename E> class Log10Op;
    template <typename T, typename E> class SqrtOp;
    template <typename T, typename E> class CosOp;
    template <typename T, typename E> class SinOp;
    template <typename T, typename E> class TanOp;
    template <typename T, typename E> class ACosOp;
    template <typename T, typename E> class ASinOp;
    template <typename T, typename E> class ATanOp;
    template <typename T, typename E> class CoshOp;
    template <typename T, typename E> class SinhOp;
    template <typename T, typename E> class TanhOp;
    template <typename T, typename E> class ACoshOp;
    template <typename T, typename E> class ASinhOp;
    template <typename T, typename E> class ATanhOp;
    template <typename T, typename E> class AbsOp;
    template <typename T, typename E> class FAbsOp;
    template <typename T, typename E> class CbrtOp;
    template <typename T, typename E, bool is_simd = IsSimdType<T>::value>
    class SafeSqrtOp;

    template <typename T1, typename T2, bool, bool, typename E>
    class AdditionOp;
    template <typename T1, typename T2, bool, bool, typename E>
    class SubtractionOp;
    template <typename T1, typename T2, bool, bool, typename E>
    class Multiplicationp;
    template <typename T1, typename T2, bool, bool, typename E>
    class DivisionOp;
    template <typename T1, typename T2, bool, bool, typename E> class Atan2Op;
    template <typename T1, typename T2, bool, bool, typename E> class MaxOp;
    template <typename T1, typename T2, bool, bool, typename E> class MinOp;

    namespace PowerImpl {
      struct Simd {};
      struct Nested {};
      struct NestedSimd {};
      struct Scalar {};

      template <typename T1, typename T2>
      struct Selector {
        static constexpr bool is_simd =
          IsSimdType<T1>::value || IsSimdType<T2>::value;
        static constexpr bool is_fad =
          IsFad< typename Sacado::ValueType<T1>::type >::value ||
          IsFad< typename Sacado::ValueType<T2>::type >::value;
        typedef typename std::conditional<
          is_simd && is_fad,
          NestedSimd,
          typename std::conditional<
            is_simd,
            Simd,
            typename std::conditional<
              is_fad,
              Nested,
              Scalar>::type
            >::type
          >::type type;
      };
    }
    template <typename T1, typename T2, bool, bool, typename E,
              typename Impl = typename PowerImpl::Selector<T1,T2>::type >
    class PowerOp;

  } // namespace Exp
  } // namespace Fad

} // namespace Sacado

#endif // SACADO_FAD_EXPRESSION_HPP
