// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBATCHED_NRM_INTERNAL_HPP_
#define KOKKOSBATCHED_NRM_INTERNAL_HPP_

#include <KokkosBatched_Util.hpp>

namespace KokkosBatched {
namespace Impl {

template <typename NrmValueType, typename ValueType>
KOKKOS_INLINE_FUNCTION NrmValueType asum(const ValueType &x) {
  if constexpr (KokkosKernels::ArithTraits<ValueType>::is_complex) {
    return KokkosKernels::ArithTraits<NrmValueType>::abs(KokkosKernels::ArithTraits<ValueType>::real(x)) +
           KokkosKernels::ArithTraits<NrmValueType>::abs(KokkosKernels::ArithTraits<ValueType>::imag(x));
  } else {
    return KokkosKernels::ArithTraits<NrmValueType>::abs(x);
  }
}

template <typename NrmValueType, typename ValueType>
KOKKOS_INLINE_FUNCTION NrmValueType square(const ValueType &x) {
  if constexpr (KokkosKernels::ArithTraits<ValueType>::is_complex) {
    return KokkosKernels::ArithTraits<ValueType>::real(KokkosKernels::ArithTraits<ValueType>::conj(x) * x);
  } else {
    return x * x;
  }
}

template <typename NrmValueType, typename ExecutionSpace>
struct NormAccumulator {
 public:
  using value_type = Kokkos::pair<NrmValueType, NrmValueType>;
  using reducer    = NormAccumulator<NrmValueType, ExecutionSpace>;

 private:
  value_type &m_accum;  // (scale, sumsq)

 public:
  KOKKOS_INLINE_FUNCTION NormAccumulator(value_type &val) : m_accum(val) {}
  KOKKOS_INLINE_FUNCTION
  void join(value_type &dest, const value_type &src) const {
    if (src.first == KokkosKernels::ArithTraits<NrmValueType>::zero()) {
      return;
    }
    if (dest.first == KokkosKernels::ArithTraits<NrmValueType>::zero()) {
      dest = src;
      return;
    }

    if (dest.first >= src.first) {
      auto r = src.first / dest.first;
      dest.second += src.second * square<NrmValueType>(r);
    } else {
      auto r      = dest.first / src.first;
      dest.second = src.second + dest.second * square<NrmValueType>(r);
      dest.first  = src.first;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type &val) const {
    val.first  = KokkosKernels::ArithTraits<NrmValueType>::zero();
    val.second = KokkosKernels::ArithTraits<NrmValueType>::one();
  }

  KOKKOS_INLINE_FUNCTION
  value_type &reference() const { return m_accum; }
};

/// \brief Helper function to compute scaled square for ScaledL2 norm. It updates the scale and sumsq values based on
/// the input x.
/// Reference:
/// dnrm2: https://www.netlib.org/lapack/lapack-3.1.1/html/dnrm2.f.html,
/// dznrm2: https://www.netlib.org/lapack/lapack-3.1.1/html/dznrm2.f.html
/// \tparam NrmValueType: Value type for the norm (e.g., double)
/// \tparam ValueType: Value type for the input (e.g., Kokkos::complex<double>)
/// \param[in] x Input value for which the scaled square is computed
/// \param[in,out] scale Current scale value, which is updated if the absolute value of x
/// is greater than the current scale
/// \param[in,out] sumsq Current sum of squares, which is updated if the absolute
/// value of x is greater than the current scale
template <typename NrmValueType, typename ValueType>
KOKKOS_INLINE_FUNCTION void scaled_square(const ValueType &x, NrmValueType &scale, NrmValueType &sumsq) {
  if constexpr (KokkosKernels::ArithTraits<ValueType>::is_complex) {
    using mag_type = typename KokkosKernels::ArithTraits<ValueType>::mag_type;
    auto x_real    = KokkosKernels::ArithTraits<ValueType>::real(x);
    if (x_real != KokkosKernels::ArithTraits<mag_type>::zero()) {
      const NrmValueType abs_val = KokkosKernels::ArithTraits<mag_type>::abs(x_real);
      if (scale < abs_val) {
        sumsq = KokkosKernels::ArithTraits<mag_type>::one() + sumsq * square<NrmValueType>(scale / abs_val);
        scale = abs_val;
      } else {
        sumsq += square<NrmValueType>(abs_val / scale);
      }
    }
    auto x_imag = KokkosKernels::ArithTraits<ValueType>::imag(x);
    if (x_imag != KokkosKernels::ArithTraits<mag_type>::zero()) {
      const NrmValueType abs_val = KokkosKernels::ArithTraits<mag_type>::abs(x_imag);
      if (scale < abs_val) {
        sumsq = KokkosKernels::ArithTraits<mag_type>::one() + sumsq * square<NrmValueType>(scale / abs_val);
        scale = abs_val;
      } else {
        sumsq += square<NrmValueType>(abs_val / scale);
      }
    }
  } else {
    if (x != KokkosKernels::ArithTraits<ValueType>::zero()) {
      const NrmValueType abs_val = KokkosKernels::ArithTraits<ValueType>::abs(x);
      if (scale < abs_val) {
        sumsq = KokkosKernels::ArithTraits<NrmValueType>::one() + sumsq * square<NrmValueType>(scale / abs_val);
        scale = abs_val;
      } else {
        sumsq += square<NrmValueType>(abs_val / scale);
      }
    }
  }
}

///
/// Serial Internal Impl
/// ====================
template <typename NrmType>
struct SerialNrmInternal {
  template <typename ValueType, typename NrmValueType>
  KOKKOS_INLINE_FUNCTION static void invoke(const int n, const ValueType *KOKKOS_RESTRICT x, const int xs0,
                                            NrmValueType *KOKKOS_RESTRICT norm);
};

template <typename NrmType>
template <typename ValueType, typename NrmValueType>
KOKKOS_INLINE_FUNCTION void SerialNrmInternal<NrmType>::invoke(const int n, const ValueType *KOKKOS_RESTRICT x,
                                                               const int xs0, NrmValueType *KOKKOS_RESTRICT norm) {
  NrmValueType nrm = 0;

  if constexpr (std::is_same_v<NrmType, Norm::L1>) {
    for (int i = 0; i < n; ++i) {
      nrm += asum<NrmValueType>(x[i * xs0]);
    }
  } else if constexpr (std::is_same_v<NrmType, Norm::L2>) {
    for (int i = 0; i < n; ++i) {
      nrm += square<NrmValueType>(x[i * xs0]);
    }
    nrm = KokkosKernels::ArithTraits<NrmValueType>::sqrt(nrm);
  } else if constexpr (std::is_same_v<NrmType, Norm::LInf>) {
    for (int i = 0; i < n; ++i) {
      const NrmValueType abs_val = KokkosKernels::ArithTraits<ValueType>::abs(x[i * xs0]);
      if (abs_val > nrm) nrm = abs_val;
    }
  } else if constexpr (std::is_same_v<NrmType, Norm::ScaledL2>) {
    if (n == 1) {
      nrm = KokkosKernels::ArithTraits<ValueType>::abs(x[0]);
    } else {
      NrmValueType scale = KokkosKernels::ArithTraits<NrmValueType>::zero();
      NrmValueType sumsq = KokkosKernels::ArithTraits<NrmValueType>::one();
      for (int i = 0; i < n; ++i) {
        scaled_square(x[i * xs0], scale, sumsq);
      }
      nrm = scale * KokkosKernels::ArithTraits<NrmValueType>::sqrt(sumsq);
    }
  }

  *norm = nrm;
}

///
/// Team Internal Impl
/// ==================
template <typename MemberType, typename NrmType>
struct TeamNrmInternal {
  template <typename ValueType, typename NrmValueType>
  KOKKOS_INLINE_FUNCTION static void invoke(const MemberType &member, const int n, const ValueType *KOKKOS_RESTRICT x,
                                            const int xs0, NrmValueType *KOKKOS_RESTRICT norm);
};

template <typename MemberType, typename NrmType>
template <typename ValueType, typename NrmValueType>
KOKKOS_INLINE_FUNCTION void TeamNrmInternal<MemberType, NrmType>::invoke(const MemberType &member, const int n,
                                                                         const ValueType *KOKKOS_RESTRICT x,
                                                                         const int xs0,
                                                                         NrmValueType *KOKKOS_RESTRICT norm) {
  NrmValueType nrm = 0;

  if constexpr (std::is_same_v<NrmType, Norm::L1>) {
    Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(member, n),
        [&](const int i, NrmValueType &thread_nrm) { thread_nrm += asum<NrmValueType>(x[i * xs0]); }, nrm);
  } else if constexpr (std::is_same_v<NrmType, Norm::L2>) {
    Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(member, n),
        [&](const int i, NrmValueType &thread_nrm) { thread_nrm += square<NrmValueType>(x[i * xs0]); }, nrm);
    nrm = KokkosKernels::ArithTraits<NrmValueType>::sqrt(nrm);
  } else if constexpr (std::is_same_v<NrmType, Norm::LInf>) {
    Kokkos::Max<NrmValueType, typename MemberType::execution_space> max_nrm(nrm);
    Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(member, n),
        [&](const int i, NrmValueType &thread_nrm) {
          const NrmValueType abs_val = KokkosKernels::ArithTraits<ValueType>::abs(x[i * xs0]);
          if (abs_val > thread_nrm) thread_nrm = abs_val;
        },
        max_nrm);
  } else if constexpr (std::is_same_v<NrmType, Norm::ScaledL2>) {
    if (n == 1) {
      nrm = KokkosKernels::ArithTraits<ValueType>::abs(x[0]);
    } else {
      Kokkos::pair<NrmValueType, NrmValueType> init_val{KokkosKernels::ArithTraits<NrmValueType>::zero(),
                                                        KokkosKernels::ArithTraits<NrmValueType>::one()};
      NormAccumulator<NrmValueType, typename MemberType::execution_space> accumulator(init_val);
      Kokkos::parallel_reduce(
          Kokkos::TeamThreadRange(member, n),
          [&](const int i,
              typename NormAccumulator<NrmValueType, typename MemberType::execution_space>::value_type &thread_accum) {
            scaled_square(x[i * xs0], thread_accum.first, thread_accum.second);
          },
          accumulator);
      const auto [scale, sumsq] = accumulator.reference();
      nrm                       = scale * KokkosKernels::ArithTraits<NrmValueType>::sqrt(sumsq);
    }
  }

  *norm = nrm;
}

///
/// TeamVector Internal Impl
/// ==================
template <typename MemberType, typename NrmType>
struct TeamVectorNrmInternal {
  template <typename ValueType, typename NrmValueType>
  KOKKOS_INLINE_FUNCTION static void invoke(const MemberType &member, const int n, const ValueType *KOKKOS_RESTRICT x,
                                            const int xs0, NrmValueType *KOKKOS_RESTRICT norm);
};

template <typename MemberType, typename NrmType>
template <typename ValueType, typename NrmValueType>
KOKKOS_INLINE_FUNCTION void TeamVectorNrmInternal<MemberType, NrmType>::invoke(const MemberType &member, const int n,
                                                                               const ValueType *KOKKOS_RESTRICT x,
                                                                               const int xs0,
                                                                               NrmValueType *KOKKOS_RESTRICT norm) {
  NrmValueType nrm = 0;

  if constexpr (std::is_same_v<NrmType, Norm::L1>) {
    Kokkos::parallel_reduce(
        Kokkos::TeamVectorRange(member, n),
        [&](const int i, NrmValueType &thread_nrm) { thread_nrm += asum<NrmValueType>(x[i * xs0]); }, nrm);
  } else if constexpr (std::is_same_v<NrmType, Norm::L2>) {
    Kokkos::parallel_reduce(
        Kokkos::TeamVectorRange(member, n),
        [&](const int i, NrmValueType &thread_nrm) { thread_nrm += square<NrmValueType>(x[i * xs0]); }, nrm);
    nrm = KokkosKernels::ArithTraits<NrmValueType>::sqrt(nrm);
  } else if constexpr (std::is_same_v<NrmType, Norm::LInf>) {
    Kokkos::Max<NrmValueType, typename MemberType::execution_space> max_nrm(nrm);
    Kokkos::parallel_reduce(
        Kokkos::TeamVectorRange(member, n),
        [&](const int i, NrmValueType &thread_nrm) {
          const NrmValueType abs_val = KokkosKernels::ArithTraits<ValueType>::abs(x[i * xs0]);
          if (abs_val > thread_nrm) thread_nrm = abs_val;
        },
        max_nrm);
  } else if constexpr (std::is_same_v<NrmType, Norm::ScaledL2>) {
    if (n == 1) {
      nrm = KokkosKernels::ArithTraits<ValueType>::abs(x[0]);
    } else {
      Kokkos::pair<NrmValueType, NrmValueType> init_val{KokkosKernels::ArithTraits<NrmValueType>::zero(),
                                                        KokkosKernels::ArithTraits<NrmValueType>::one()};
      NormAccumulator<NrmValueType, typename MemberType::execution_space> accumulator(init_val);
      Kokkos::parallel_reduce(
          Kokkos::TeamVectorRange(member, n),
          [&](const int i,
              typename NormAccumulator<NrmValueType, typename MemberType::execution_space>::value_type &thread_accum) {
            scaled_square(x[i * xs0], thread_accum.first, thread_accum.second);
          },
          accumulator);
      const auto [scale, sumsq] = accumulator.reference();
      nrm                       = scale * KokkosKernels::ArithTraits<NrmValueType>::sqrt(sumsq);
    }
  }
  *norm = nrm;
}

}  // namespace Impl
}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_NRM_INTERNAL_HPP_
