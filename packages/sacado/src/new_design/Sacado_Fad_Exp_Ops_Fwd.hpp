// @HEADER
// ***********************************************************************
//
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef SACADO_EXP_FAD_OPS_FWD_HPP
#define SACADO_EXP_FAD_OPS_FWD_HPP

namespace Sacado {

  template <typename T> struct IsSimdType;

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

    template <typename T1, typename T2, bool, bool, typename E>
    class AdditionOp;
    template <typename T1, typename T2, bool, bool, typename E>
    class SubtractionOp;
    template <typename T1, typename T2, bool, bool, typename E>
    class Multiplicationp;
    template <typename T1, typename T2, bool, bool, typename E>
    class DivisionOp;
    template <typename T1, typename T2, bool, bool, typename E> class Atan2Op;
    template <typename T1, typename T2, bool, bool, typename E,
              bool is_simd = IsSimdType<T1>::value || IsSimdType<T2>::value>
    class PowerOp;
    template <typename T1, typename T2, bool, bool, typename E> class MaxOp;
    template <typename T1, typename T2, bool, bool, typename E> class MinOp;

  } // namespace Exp
  } // namespace Fad

} // namespace Sacado

#endif // SACADO_FAD_EXPRESSION_HPP
