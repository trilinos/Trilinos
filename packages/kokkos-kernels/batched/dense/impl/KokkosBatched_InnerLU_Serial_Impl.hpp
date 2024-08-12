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
#ifndef __KOKKOSBATCHED_INNER_LU_SERIAL_IMPL_HPP__
#define __KOKKOSBATCHED_INNER_LU_SERIAL_IMPL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_InnerLU_Decl.hpp"

namespace KokkosBatched {

///
/// Fixed size LU
/// ================

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerLU<5>::serial_invoke(ValueType *KOKKOS_RESTRICT A) {
  // load
  ValueType a_00 = A[0 * _as0 + 0 * _as1], a_01 = A[0 * _as0 + 1 * _as1], a_02 = A[0 * _as0 + 2 * _as1],
            a_03 = A[0 * _as0 + 3 * _as1], a_04 = A[0 * _as0 + 4 * _as1], a_10 = A[1 * _as0 + 0 * _as1],
            a_11 = A[1 * _as0 + 1 * _as1], a_12 = A[1 * _as0 + 2 * _as1], a_13 = A[1 * _as0 + 3 * _as1],
            a_14 = A[1 * _as0 + 4 * _as1], a_20 = A[2 * _as0 + 0 * _as1], a_21 = A[2 * _as0 + 1 * _as1],
            a_22 = A[2 * _as0 + 2 * _as1], a_23 = A[2 * _as0 + 3 * _as1], a_24 = A[2 * _as0 + 4 * _as1],
            a_30 = A[3 * _as0 + 0 * _as1], a_31 = A[3 * _as0 + 1 * _as1], a_32 = A[3 * _as0 + 2 * _as1],
            a_33 = A[3 * _as0 + 3 * _as1], a_34 = A[3 * _as0 + 4 * _as1], a_40 = A[4 * _as0 + 0 * _as1],
            a_41 = A[4 * _as0 + 1 * _as1], a_42 = A[4 * _as0 + 2 * _as1], a_43 = A[4 * _as0 + 3 * _as1],
            a_44 = A[4 * _as0 + 4 * _as1];

  // 0 iteration
  a_10 /= a_00;
  a_11 -= a_10 * a_01;
  a_12 -= a_10 * a_02;
  a_13 -= a_10 * a_03;
  a_14 -= a_10 * a_04;
  a_20 /= a_00;
  a_21 -= a_20 * a_01;
  a_22 -= a_20 * a_02;
  a_23 -= a_20 * a_03;
  a_24 -= a_20 * a_04;
  a_30 /= a_00;
  a_31 -= a_30 * a_01;
  a_32 -= a_30 * a_02;
  a_33 -= a_30 * a_03;
  a_34 -= a_30 * a_04;
  a_40 /= a_00;
  a_41 -= a_40 * a_01;
  a_42 -= a_40 * a_02;
  a_43 -= a_40 * a_03;
  a_44 -= a_40 * a_04;

  // 1 iteration
  a_21 /= a_11;
  a_22 -= a_21 * a_12;
  a_23 -= a_21 * a_13;
  a_24 -= a_21 * a_14;
  a_31 /= a_11;
  a_32 -= a_31 * a_12;
  a_33 -= a_31 * a_13;
  a_34 -= a_31 * a_14;
  a_41 /= a_11;
  a_42 -= a_41 * a_12;
  a_43 -= a_41 * a_13;
  a_44 -= a_41 * a_14;

  // 2 iteration
  a_32 /= a_22;
  a_33 -= a_32 * a_23;
  a_34 -= a_32 * a_24;
  a_42 /= a_22;
  a_43 -= a_42 * a_23;
  a_44 -= a_42 * a_24;

  // 3 iteration
  a_43 /= a_33;
  a_44 -= a_43 * a_34;

  // store
  A[1 * _as0 + 0 * _as1] = a_10;
  A[1 * _as0 + 1 * _as1] = a_11;
  A[1 * _as0 + 2 * _as1] = a_12;
  A[1 * _as0 + 3 * _as1] = a_13;
  A[1 * _as0 + 4 * _as1] = a_14;
  A[2 * _as0 + 0 * _as1] = a_20;
  A[2 * _as0 + 1 * _as1] = a_21;
  A[2 * _as0 + 2 * _as1] = a_22;
  A[2 * _as0 + 3 * _as1] = a_23;
  A[2 * _as0 + 4 * _as1] = a_24;
  A[3 * _as0 + 0 * _as1] = a_30;
  A[3 * _as0 + 1 * _as1] = a_31;
  A[3 * _as0 + 2 * _as1] = a_32;
  A[3 * _as0 + 3 * _as1] = a_33;
  A[3 * _as0 + 4 * _as1] = a_34;
  A[4 * _as0 + 0 * _as1] = a_40;
  A[4 * _as0 + 1 * _as1] = a_41;
  A[4 * _as0 + 2 * _as1] = a_42;
  A[4 * _as0 + 3 * _as1] = a_43;
  A[4 * _as0 + 4 * _as1] = a_44;

  return 0;
}

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerLU<4>::serial_invoke(ValueType *KOKKOS_RESTRICT A) {
  // load
  ValueType a_00 = A[0 * _as0 + 0 * _as1], a_01 = A[0 * _as0 + 1 * _as1], a_02 = A[0 * _as0 + 2 * _as1],
            a_03 = A[0 * _as0 + 3 * _as1], a_10 = A[1 * _as0 + 0 * _as1], a_11 = A[1 * _as0 + 1 * _as1],
            a_12 = A[1 * _as0 + 2 * _as1], a_13 = A[1 * _as0 + 3 * _as1], a_20 = A[2 * _as0 + 0 * _as1],
            a_21 = A[2 * _as0 + 1 * _as1], a_22 = A[2 * _as0 + 2 * _as1], a_23 = A[2 * _as0 + 3 * _as1],
            a_30 = A[3 * _as0 + 0 * _as1], a_31 = A[3 * _as0 + 1 * _as1], a_32 = A[3 * _as0 + 2 * _as1],
            a_33 = A[3 * _as0 + 3 * _as1];

  // 0 iteration
  a_10 /= a_00;
  a_11 -= a_10 * a_01;
  a_12 -= a_10 * a_02;
  a_13 -= a_10 * a_03;
  a_20 /= a_00;
  a_21 -= a_20 * a_01;
  a_22 -= a_20 * a_02;
  a_23 -= a_20 * a_03;
  a_30 /= a_00;
  a_31 -= a_30 * a_01;
  a_32 -= a_30 * a_02;
  a_33 -= a_30 * a_03;

  // 1 iteration
  a_21 /= a_11;
  a_22 -= a_21 * a_12;
  a_23 -= a_21 * a_13;
  a_31 /= a_11;
  a_32 -= a_31 * a_12;
  a_33 -= a_31 * a_13;

  // 2 iteration
  a_32 /= a_22;
  a_33 -= a_32 * a_23;

  // store
  A[1 * _as0 + 0 * _as1] = a_10;
  A[1 * _as0 + 1 * _as1] = a_11;
  A[1 * _as0 + 2 * _as1] = a_12;
  A[1 * _as0 + 3 * _as1] = a_13;
  A[2 * _as0 + 0 * _as1] = a_20;
  A[2 * _as0 + 1 * _as1] = a_21;
  A[2 * _as0 + 2 * _as1] = a_22;
  A[2 * _as0 + 3 * _as1] = a_23;
  A[3 * _as0 + 0 * _as1] = a_30;
  A[3 * _as0 + 1 * _as1] = a_31;
  A[3 * _as0 + 2 * _as1] = a_32;
  A[3 * _as0 + 3 * _as1] = a_33;

  return 0;
}

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerLU<3>::serial_invoke(ValueType *KOKKOS_RESTRICT A) {
  // load
  ValueType a_00 = A[0 * _as0 + 0 * _as1], a_01 = A[0 * _as0 + 1 * _as1], a_02 = A[0 * _as0 + 2 * _as1],
            a_10 = A[1 * _as0 + 0 * _as1], a_11 = A[1 * _as0 + 1 * _as1], a_12 = A[1 * _as0 + 2 * _as1],
            a_20 = A[2 * _as0 + 0 * _as1], a_21 = A[2 * _as0 + 1 * _as1], a_22 = A[2 * _as0 + 2 * _as1];

  // 0 iteration
  a_10 /= a_00;
  a_11 -= a_10 * a_01;
  a_12 -= a_10 * a_02;
  a_20 /= a_00;
  a_21 -= a_20 * a_01;
  a_22 -= a_20 * a_02;

  // 1 iteration
  a_21 /= a_11;
  a_22 -= a_21 * a_12;

  // store
  A[1 * _as0 + 0 * _as1] = a_10;
  A[1 * _as0 + 1 * _as1] = a_11;
  A[1 * _as0 + 2 * _as1] = a_12;
  A[2 * _as0 + 0 * _as1] = a_20;
  A[2 * _as0 + 1 * _as1] = a_21;
  A[2 * _as0 + 2 * _as1] = a_22;

  return 0;
}

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerLU<2>::serial_invoke(ValueType *KOKKOS_RESTRICT A) {
  // load
  ValueType a_00 = A[0 * _as0 + 0 * _as1], a_01 = A[0 * _as0 + 1 * _as1], a_10 = A[1 * _as0 + 0 * _as1],
            a_11 = A[1 * _as0 + 1 * _as1];

  // 0 iteration
  a_10 /= a_00;
  a_11 -= a_10 * a_01;

  // store
  A[1 * _as0 + 0 * _as1] = a_10;
  A[1 * _as0 + 1 * _as1] = a_11;

  return 0;
}

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerLU<1>::serial_invoke(ValueType *KOKKOS_RESTRICT /* A */) {
  return 0;
}

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerLU<5>::serial_invoke(const int m, ValueType *KOKKOS_RESTRICT A) {
  if (m > 5) Kokkos::abort("InnerLU<5>::serial_invoke, assert failure (m<=5)");
  if (m <= 0) return 0;

  switch (m) {
    case 5: {
      InnerLU<5> inner(_as0, _as1);
      inner.serial_invoke(A);
      break;
    }
    case 4: {
      InnerLU<4> inner(_as0, _as1);
      inner.serial_invoke(A);
      break;
    }
    case 3: {
      InnerLU<3> inner(_as0, _as1);
      inner.serial_invoke(A);
      break;
    }
    case 2: {
      InnerLU<2> inner(_as0, _as1);
      inner.serial_invoke(A);
      break;
    }
    case 1: {
      InnerLU<1> inner(_as0, _as1);
      inner.serial_invoke(A);
      break;
    }
  }
  return 0;
}

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerLU<4>::serial_invoke(const int m, ValueType *KOKKOS_RESTRICT A) {
  if (m > 4) Kokkos::abort("InnerLU<4>::serial_invoke, assert failure (m<=4)");
  if (m <= 0) return 0;

  switch (m) {
    case 4: {
      InnerLU<4> inner(_as0, _as1);
      inner.serial_invoke(A);
      break;
    }
    case 3: {
      InnerLU<3> inner(_as0, _as1);
      inner.serial_invoke(A);
      break;
    }
    case 2: {
      InnerLU<2> inner(_as0, _as1);
      inner.serial_invoke(A);
      break;
    }
    case 1: {
      InnerLU<1> inner(_as0, _as1);
      inner.serial_invoke(A);
      break;
    }
  }
  return 0;
}

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerLU<3>::serial_invoke(const int m, ValueType *KOKKOS_RESTRICT A) {
  if (m > 3) Kokkos::abort("InnerLU<3>::serial_invoke, assert failure (m<=3)");
  if (m <= 0) return 0;

  switch (m) {
    case 3: {
      InnerLU<3> inner(_as0, _as1);
      inner.serial_invoke(A);
      break;
    }
    case 2: {
      InnerLU<2> inner(_as0, _as1);
      inner.serial_invoke(A);
      break;
    }
    case 1: {
      InnerLU<1> inner(_as0, _as1);
      inner.serial_invoke(A);
      break;
    }
  }
  return 0;
}

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerLU<2>::serial_invoke(const int m, ValueType *KOKKOS_RESTRICT A) {
  if (m > 2) Kokkos::abort("InnerLU<2>::serial_invoke, assert failure (m<=2)");
  if (m <= 0) return 0;

  switch (m) {
    case 2: {
      InnerLU<2> inner(_as0, _as1);
      inner.serial_invoke(A);
      break;
    }
    case 1: {
      InnerLU<1> inner(_as0, _as1);
      inner.serial_invoke(A);
      break;
    }
  }
  return 0;
}

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerLU<1>::serial_invoke(const int m, ValueType *KOKKOS_RESTRICT A) {
  if (m > 1) Kokkos::abort("InnerLU<1>::serial_invoke, assert failure (m<=1)");
  if (m <= 0) return 0;

  switch (m) {
    case 1: {
      InnerLU<1> inner(_as0, _as1);
      inner.serial_invoke(A);
      break;
    }
  }
  return 0;
}

// template<int bmn>
// template<typename ValueType>
// KOKKOS_INLINE_FUNCTION
// int
// InnerLU<bmn>::
// serial_invoke(const int m, const int n,
//               ValueType *KOKKOS_RESTRICT A) {
//   if (m <= 0 || n <= 0) return 0;
//   const int k = m < n ? m : n;
//   for (int p=0;p<k;++p) {
//     const ValueType
//       // inv_alpha11 = 1.0/A[p*_as0+p*_as1],
//       alpha11 = A[p*_as0+p*_as1],
//       *KOKKOS_RESTRICT a12t = A + (p  )*_as0 + (p+1)*_as1;

//     ValueType
//       *KOKKOS_RESTRICT a21  = A + (p+1)*_as0 + (p  )*_as1,
//       *KOKKOS_RESTRICT A22  = A + (p+1)*_as0 + (p+1)*_as1;

//     const int
//       iend = m-p-1,
//       jend = n-p-1;

//     for (int i=0;i<iend;++i) {
//       // a21[i*_as0] *= inv_alpha11;
//       a21[i*_as0] /= alpha11;
//       for (int j=0;j<jend;++j)
//         A22[i*_as0+j*_as1] -= a21[i*_as0] * a12t[j*_as1];
//     }
//   }
//   return 0;
// }

}  // namespace KokkosBatched

#endif
