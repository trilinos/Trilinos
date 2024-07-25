// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
//
// ***********************************************************************
//
// The forward-mode AD classes in Sacado are a derivative work of the
// expression template classes in the Fad package by Nicolas Di Cesare.
// The following banner is included in the original Fad source code:
//
// ************ DO NOT REMOVE THIS BANNER ****************
//
//  Nicolas Di Cesare <Nicolas.Dicesare@ann.jussieu.fr>
//  http://www.ann.jussieu.fr/~dicesare
//
//            CEMRACS 98 : C++ courses,
//         templates : new C++ techniques
//            for scientific computing
//
//********************************************************
//
//  A short implementation ( not all operators and
//  functions are overloaded ) of 1st order Automatic
//  Differentiation in forward mode (FAD) using
//  EXPRESSION TEMPLATES.
//
//********************************************************
// @HEADER

#ifndef SACADO_FAD_OPS_FWD_HPP
#define SACADO_FAD_OPS_FWD_HPP

#include <type_traits>

namespace Sacado {

  template <typename T> struct IsSimdType;
  template <typename T> struct IsFad;
  template <typename T> struct ValueType;

  namespace Fad {

    template <typename ExprT> class UnaryPlusOp;
    template <typename ExprT> class UnaryMinusOp;
    template <typename ExprT> class ExpOp;
    template <typename ExprT> class LogOp;
    template <typename ExprT> class Log10Op;
    template <typename ExprT> class SqrtOp;
    template <typename ExprT> class CosOp;
    template <typename ExprT> class SinOp;
    template <typename ExprT> class TanOp;
    template <typename ExprT> class ACosOp;
    template <typename ExprT> class ASinOp;
    template <typename ExprT> class ATanOp;
    template <typename ExprT> class CoshOp;
    template <typename ExprT> class SinhOp;
    template <typename ExprT> class TanhOp;
    template <typename ExprT> class ACoshOp;
    template <typename ExprT> class ASinhOp;
    template <typename ExprT> class ATanhOp;
    template <typename ExprT> class AbsOp;
    template <typename ExprT> class FAbsOp;
    template <typename ExprT> class CbrtOp;
    template <typename ExprT, bool is_simd = IsSimdType<ExprT>::value>
    class SafeSqrtOp;

    template <typename ExprT1, typename ExprT2> class AdditionOp;
    template <typename ExprT1, typename ExprT2> class SubtractionOp;
    template <typename ExprT1, typename ExprT2> class Multiplicationp;
    template <typename ExprT1, typename ExprT2> class DivisionOp;
    template <typename ExprT1, typename ExprT2> class Atan2Op;
    template <typename ExprT1, typename ExprT2> class MaxOp;
    template <typename ExprT1, typename ExprT2> class MinOp;

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
    template <typename ExprT1, typename ExprT2,
              typename Impl = typename PowerImpl::Selector<ExprT1,ExprT2>::type>
    class PowerOp;

  } // namespace Fad

} // namespace Sacado

#endif // SACADO_FAD_EXPRESSION_HPP
