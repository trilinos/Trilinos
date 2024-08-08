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
#ifndef __KOKKOSBATCHED_VECTOR_HPP__
#define __KOKKOSBATCHED_VECTOR_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

// forward declaration
namespace KokkosBatched {

template <typename T, int l>
class Vector;

template <typename T, int l>
struct is_vector<Vector<SIMD<T>, l>> : public std::true_type {};

template <typename ValueType, typename MemorySpace>
struct DefaultVectorLength {
  enum : int { value = 1 };
};

template <>
struct DefaultVectorLength<float, Kokkos::HostSpace> {
#if defined(__AVX512F__)
  enum : int{value = 16};
#elif defined(__AVX__) || defined(__AVX2__)
  enum : int{value = 8};
#elif defined(__ARM_ARCH)
  enum : int{value = 8};
#else
  enum : int { value = 8 };
#endif
};
template <>
struct DefaultVectorLength<double, Kokkos::HostSpace> {
#if defined(__AVX512F__)
  enum : int{value = 8};
#elif defined(__AVX__) || defined(__AVX2__)
  enum : int{value = 4};
#elif defined(__ARM_ARCH)
  enum : int{value = 4};
#else
  enum : int { value = 4 };
#endif
};
template <>
struct DefaultVectorLength<Kokkos::complex<float>, Kokkos::HostSpace> {
#if defined(__AVX512F__)
  enum : int{value = 8};
#elif defined(__AVX__) || defined(__AVX2__)
  enum : int{value = 4};
#elif defined(__ARM_ARCH)
  enum : int{value = 4};
#else
  enum : int { value = 4 };
#endif
};
template <>
struct DefaultVectorLength<Kokkos::complex<double>, Kokkos::HostSpace> {
#if defined(__AVX512F__)
  enum : int{value = 4};
#elif defined(__AVX__) || defined(__AVX2__)
  enum : int{value = 2};
#elif defined(__ARM_ARCH)
  enum : int{value = 2};
#else
  enum : int { value = 2 };
#endif
};

#if defined(KOKKOS_ENABLE_CUDA)
template <>
struct DefaultVectorLength<float, Kokkos::CudaSpace> {
  enum : int { value = 8 };
};
template <>
struct DefaultVectorLength<double, Kokkos::CudaSpace> {
  enum : int { value = 8 };
};
template <>
struct DefaultVectorLength<Kokkos::complex<float>, Kokkos::CudaSpace> {
  enum : int { value = 8 };
};
template <>
struct DefaultVectorLength<Kokkos::complex<double>, Kokkos::CudaSpace> {
  enum : int { value = 8 };
};
template <>
struct DefaultVectorLength<float, Kokkos::CudaUVMSpace> {
  enum : int { value = 8 };
};
template <>
struct DefaultVectorLength<double, Kokkos::CudaUVMSpace> {
  enum : int { value = 8 };
};
template <>
struct DefaultVectorLength<Kokkos::complex<float>, Kokkos::CudaUVMSpace> {
  enum : int { value = 8 };
};
template <>
struct DefaultVectorLength<Kokkos::complex<double>, Kokkos::CudaUVMSpace> {
  enum : int { value = 8 };
};
#endif

#if defined(KOKKOS_ENABLE_HIP)
template <>
struct DefaultVectorLength<float, Kokkos::HIPSpace> {
  enum : int { value = 16 };
};
template <>
struct DefaultVectorLength<double, Kokkos::HIPSpace> {
  enum : int { value = 16 };
};
template <>
struct DefaultVectorLength<Kokkos::complex<float>, Kokkos::HIPSpace> {
  enum : int { value = 16 };
};
template <>
struct DefaultVectorLength<Kokkos::complex<double>, Kokkos::HIPSpace> {
  enum : int { value = 16 };
};
#endif

template <typename ValueType, typename MemorySpace>
struct DefaultInternalVectorLength {
  enum : int { value = 1 };
};
template <typename ValueType>
struct DefaultInternalVectorLength<ValueType, Kokkos::HostSpace> {
  enum : int { value = DefaultVectorLength<ValueType, Kokkos::HostSpace>::value };
};

#if defined(KOKKOS_ENABLE_CUDA)
template <>
struct DefaultInternalVectorLength<float, Kokkos::CudaSpace> {
  enum : int { value = 4 };
};
template <>
struct DefaultInternalVectorLength<double, Kokkos::CudaSpace> {
  enum : int { value = 2 };
};
template <>
struct DefaultInternalVectorLength<Kokkos::complex<float>, Kokkos::CudaSpace> {
  enum : int { value = 2 };
};
template <>
struct DefaultInternalVectorLength<Kokkos::complex<double>, Kokkos::CudaSpace> {
  enum : int { value = 1 };
};
template <>
struct DefaultInternalVectorLength<float, Kokkos::CudaUVMSpace> {
  enum : int { value = 4 };
};
template <>
struct DefaultInternalVectorLength<double, Kokkos::CudaUVMSpace> {
  enum : int { value = 2 };
};
template <>
struct DefaultInternalVectorLength<Kokkos::complex<float>, Kokkos::CudaUVMSpace> {
  enum : int { value = 2 };
};
template <>
struct DefaultInternalVectorLength<Kokkos::complex<double>, Kokkos::CudaUVMSpace> {
  enum : int { value = 1 };
};
#endif

#if defined(KOKKOS_ENABLE_HIP)
template <>
struct DefaultInternalVectorLength<float, Kokkos::HIPSpace> {
  enum : int { value = 8 };
};
template <>
struct DefaultInternalVectorLength<double, Kokkos::HIPSpace> {
  enum : int { value = 4 };
};
template <>
struct DefaultInternalVectorLength<Kokkos::complex<float>, Kokkos::HIPSpace> {
  enum : int { value = 4 };
};
template <>
struct DefaultInternalVectorLength<Kokkos::complex<double>, Kokkos::HIPSpace> {
  enum : int { value = 2 };
};
#endif

template <typename T>
struct MagnitudeScalarType;

template <>
struct MagnitudeScalarType<float> {
  typedef float type;
};
template <>
struct MagnitudeScalarType<double> {
  typedef double type;
};
template <>
struct MagnitudeScalarType<Kokkos::complex<float>> {
  typedef float type;
};
template <>
struct MagnitudeScalarType<Kokkos::complex<double>> {
  typedef double type;
};

template <int l>
struct MagnitudeScalarType<Vector<SIMD<float>, l>> {
  typedef float type;
};
template <int l>
struct MagnitudeScalarType<Vector<SIMD<double>, l>> {
  typedef double type;
};
template <int l>
struct MagnitudeScalarType<Vector<SIMD<Kokkos::complex<float>>, l>> {
  typedef float type;
};
template <int l>
struct MagnitudeScalarType<Vector<SIMD<Kokkos::complex<double>>, l>> {
  typedef double type;
};

}  // namespace KokkosBatched

#include "KokkosBatched_Vector_SIMD.hpp"

// arith traits overload for vector types
namespace Kokkos {

// do not use Vector alone as other can use the name.

template <typename T, int l>
class ArithTraits<KokkosBatched::Vector<KokkosBatched::SIMD<T>, l>> {
 public:
  typedef typename ArithTraits<T>::val_type val_scalar_type;
  typedef typename ArithTraits<T>::mag_type mag_scalar_type;

  typedef KokkosBatched::Vector<KokkosBatched::SIMD<val_scalar_type>, l> val_type;
  typedef KokkosBatched::Vector<KokkosBatched::SIMD<mag_scalar_type>, l> mag_type;

  static KOKKOS_FORCEINLINE_FUNCTION mag_type real(const val_type &val) { return val; }

  static KOKKOS_FORCEINLINE_FUNCTION val_type conj(const val_type &val) { return val; }

  static KOKKOS_FORCEINLINE_FUNCTION val_type abs(const val_type &val) {
    using KAT = ArithTraits<typename val_type::value_type>;
    val_type v{};
    for (int i = 0; i < l; ++i) {
      v[i] = KAT::abs(val[i]);
    }
    return v;
  }

  static const bool is_specialized = ArithTraits<T>::is_specialized;
  static const bool is_signed      = ArithTraits<T>::is_signed;
  static const bool is_integer     = ArithTraits<T>::is_integer;
  static const bool is_exact       = ArithTraits<T>::is_exact;
  static const bool is_complex     = ArithTraits<T>::is_complex;
};

template <typename T, int l>
class ArithTraits<KokkosBatched::Vector<KokkosBatched::SIMD<Kokkos::complex<T>>, l>> {
 public:
  typedef typename ArithTraits<T>::val_type val_scalar_type;
  typedef typename ArithTraits<T>::mag_type mag_scalar_type;

  typedef KokkosBatched::Vector<KokkosBatched::SIMD<Kokkos::complex<val_scalar_type>>, l> val_type;
  typedef KokkosBatched::Vector<KokkosBatched::SIMD<mag_scalar_type>, l> mag_type;

  static KOKKOS_FORCEINLINE_FUNCTION mag_type real(const val_type &val) {
    mag_type r_val;
    for (int i = 0; i < l; ++i) {
      r_val[i] = val[i].real();
    }
    return r_val;
  }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type imag(const val_type &val) {
    mag_type r_val;
    for (int i = 0; i < l; ++i) {
      r_val[i] = val[i].imag();
    }
    return r_val;
  }

  static KOKKOS_FORCEINLINE_FUNCTION val_type conj(const val_type &val) {
    using KAT = ArithTraits<typename val_type::value_type>;
    val_type v{};
    for (int i = 0; i < l; ++i) {
      v[i] = KAT::conj(val[i]);
    }
    return v;
  }

  static KOKKOS_FORCEINLINE_FUNCTION val_type abs(const val_type &val) {
    using KAT = ArithTraits<typename val_type::value_type>;
    val_type v{};
    for (int i = 0; i < l; ++i) {
      v[i] = KAT::abs(val[i]);
    }
    return v;
  }
};

}  // namespace Kokkos

#endif
