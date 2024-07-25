// clang-format off
// @HEADER
// *****************************************************************************
//                            Tacho package
//
// Copyright 2022 NTESS and the Tacho contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
// clang-format on
#ifndef __TACHO_HPP__
#define __TACHO_HPP__

#include "Tacho_config.h"

#include "Kokkos_Core.hpp"
#include "Kokkos_Timer.hpp"

#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

/// \file Tacho.hpp
/// \brief Header to be included by users
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

///
/// default ordinal and size type
///

#if defined(TACHO_USE_INT_INT)
typedef int ordinal_type;
typedef int size_type;
#elif defined(TACHO_USE_INT_SIZE_T)
typedef int ordinal_type;
typedef size_t size_type;
#else
typedef int ordinal_type;
typedef size_t size_type;
#endif

///
/// default Kokkos types (non-specialized code path is error)
///
template <typename ExecSpace> struct UseThisDevice;

template <typename ExecSpace> struct UseThisScheduler;

template <typename T, typename ExecSpace> struct UseThisFuture;

///
/// dummy objects when kokkos tasking is not used
///
template <typename ExecSpace> struct DummyTaskScheduler {
  static_assert(Kokkos::is_execution_space<ExecSpace>::value, "Error: ExecSpace is not an execution space");
  using execution_space = ExecSpace;
};

template <typename T, typename ExecSpace> struct DummyFuture {
  DummyFuture() = default;
  DummyFuture(const DummyFuture<T, ExecSpace> &b) = default;

  void clear() {}
};

/// until kokkos dual view issue is resolved, we follow the default space in Trilinos (uvm)
#if defined(KOKKOS_ENABLE_CUDA)
template <> struct UseThisDevice<Kokkos::Cuda> {
  using type = Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>;
  using device_type = type;
};
template <> struct UseThisScheduler<Kokkos::Cuda> {
  using type = DummyTaskScheduler<Kokkos::Cuda>;
  using scheduler_type = type;
};
template <typename T> struct UseThisFuture<T, Kokkos::Cuda> {
  using type = DummyFuture<T, Kokkos::Cuda>;
  using future_type = type;
};
#endif
#if defined(KOKKOS_ENABLE_HIP)
template <> struct UseThisDevice<Kokkos::HIP> {
  using type = Kokkos::Device<Kokkos::HIP, Kokkos::HIPSpace>;
  using device_type = type;
};
template <> struct UseThisScheduler<Kokkos::HIP> {
  using type = DummyTaskScheduler<Kokkos::HIP>;
  using scheduler_type = type;
};
template <typename T> struct UseThisFuture<T, Kokkos::HIP> {
  using type = DummyFuture<T, Kokkos::HIP>;
  using future_type = type;
};
#endif
#if defined(KOKKOS_ENABLE_SYCL)
template <> struct UseThisDevice<Kokkos::Experimental::SYCL> {
  using type = Kokkos::Device<Kokkos::Experimental::SYCL, Kokkos::Experimental::SYCLDeviceUSMSpace>;
  using device_type = type;
};
template <> struct UseThisScheduler<Kokkos::Experimental::SYCL> {
  using type = DummyTaskScheduler<Kokkos::Experimental::SYCL>;
  using scheduler_type = type;
};
template <typename T> struct UseThisFuture<T, Kokkos::Experimental::SYCL> {
  using type = DummyFuture<T, Kokkos::Experimental::SYCL>;
  using future_type = type;
};
#endif
#if defined(KOKKOS_ENABLE_OPENMP)
template <> struct UseThisDevice<Kokkos::OpenMP> {
  using type = Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace>;
  using device_type = type;
};
template <> struct UseThisScheduler<Kokkos::OpenMP> {
#if defined(KOKKOS_ENABLE_TASKDAG) && false
  using type = Kokkos::TaskSchedulerMultiple<Kokkos::OpenMP>;
#else
  using type = DummyTaskScheduler<Kokkos::OpenMP>;
#endif
  using scheduler_type = type;
};
template <typename T> struct UseThisFuture<T, Kokkos::OpenMP> {
#if defined(KOKKOS_ENABLE_TASKDAG) && false
  using type = Kokkos::BasicFuture<T, Kokkos::OpenMP>;
#else
  using type = DummyFuture<T, Kokkos::OpenMP>;
#endif
  using future_type = type;
};
#endif
#if defined(KOKKOS_ENABLE_SERIAL)
template <> struct UseThisDevice<Kokkos::Serial> {
  using type = Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>;
  using device_type = type;
};
template <> struct UseThisScheduler<Kokkos::Serial> {
  using type = DummyTaskScheduler<Kokkos::Serial>;
  using scheduler_type = type;
};
template <typename T> struct UseThisFuture<T, Kokkos::Serial> {
  using type = DummyFuture<T, Kokkos::Serial>;
  using future_type = type;
};
#endif

///
/// print execution spaces
///
template <typename SpT> void printExecSpaceConfiguration(std::string name, const bool detail = false) {
  if (!Kokkos::is_space<SpT>::value) {
    std::string msg("SpT is not Kokkos execution space");
    fprintf(stderr, ">> Error in file %s, line %d\n", __FILE__, __LINE__);
    fprintf(stderr, "   %s\n", msg.c_str());
    throw std::logic_error(msg.c_str());
  }
  bool is_printed(false);
#if defined(KOKKOS_ENABLE_SERIAL)
  if (std::is_same<SpT, Kokkos::Serial>::value) {
    is_printed = true;
    std::cout << std::setw(16) << name << ":: Serial \n";
  }
#endif
#if defined(KOKKOS_ENABLE_OPENMP)
  if (std::is_same<SpT, Kokkos::OpenMP>::value) {
    is_printed = true;
    std::cout << std::setw(16) << name << ":: OpenMP \n";
  }
#endif
#if defined(KOKKOS_ENABLE_CUDA)
  if (std::is_same<SpT, Kokkos::Cuda>::value) {
    is_printed = true;
    std::cout << std::setw(16) << name << ":: Cuda \n";
  }
#endif
#if defined(KOKKOS_ENABLE_HIP)
  if (std::is_same<SpT, Kokkos::HIP>::value) {
    is_printed = true;
    std::cout << std::setw(16) << name << ":: HIP \n";
  }
#endif
  if (!is_printed) {
    std::cout << std::setw(16) << name << ":: not supported Kokkos execution space\n";
    SpT().print_configuration(std::cout, detail);
    throw std::logic_error("Error: not supported Kokkos execution space");
  }
  if (detail)
    SpT().print_configuration(std::cout, true);
}

template <typename T> struct ArithTraits;

template <> struct ArithTraits<float> {
  typedef float val_type;
  typedef float mag_type;

  enum : bool { is_complex = false };
  static KOKKOS_FORCEINLINE_FUNCTION mag_type abs(const val_type &x) { return x > 0 ? x : -x; }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type real(const val_type &x) { return x; }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type imag(const val_type &x) { return 0; }
  static KOKKOS_FORCEINLINE_FUNCTION val_type conj(const val_type &x) { return x; }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type epsilon() { return FLT_EPSILON; }
  static KOKKOS_FORCEINLINE_FUNCTION void set_real(val_type &x, const mag_type &val) { x = val; }
  static KOKKOS_FORCEINLINE_FUNCTION void set_imag(val_type &x, const mag_type &val) {}
};

template <> struct ArithTraits<double> {
  typedef double val_type;
  typedef double mag_type;

  enum : bool { is_complex = false };
  static KOKKOS_FORCEINLINE_FUNCTION mag_type abs(const val_type &x) { return x > 0 ? x : -x; }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type real(const val_type &x) { return x; }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type imag(const val_type &x) { return 0; }
  static KOKKOS_FORCEINLINE_FUNCTION val_type conj(const val_type &x) { return x; }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type epsilon() { return DBL_EPSILON; }
  static KOKKOS_FORCEINLINE_FUNCTION void set_real(val_type &x, const mag_type &val) { x = val; }
  static KOKKOS_FORCEINLINE_FUNCTION void set_imag(val_type &x, const mag_type &val) {}
};

template <> struct ArithTraits<std::complex<float>> {
  typedef std::complex<float> val_type;
  typedef float mag_type;

  enum : bool { is_complex = true };
  static inline mag_type abs(const val_type &x) { return std::abs(x); }
  static inline mag_type real(const val_type &x) { return x.real(); }
  static inline mag_type imag(const val_type &x) { return x.imag(); }
  static inline val_type conj(const val_type &x) { return std::conj(x); }
  static inline mag_type epsilon() { return FLT_EPSILON; }
  static inline void set_real(val_type &x, const mag_type &val) { x.real(val); }
  static inline void set_imag(val_type &x, const mag_type &val) { x.imag(val); }
};

template <> struct ArithTraits<std::complex<double>> {
  typedef std::complex<double> val_type;
  typedef double mag_type;

  enum : bool { is_complex = true };
  static inline mag_type abs(const val_type &x) { return std::abs(x); }
  static inline mag_type real(const val_type &x) { return x.real(); }
  static inline mag_type imag(const val_type &x) { return x.imag(); }
  static inline val_type conj(const val_type &x) { return std::conj(x); }
  static inline mag_type epsilon() { return DBL_EPSILON; }
  static inline void set_real(val_type &x, const mag_type &val) { x.real(val); }
  static inline void set_imag(val_type &x, const mag_type &val) { x.imag(val); }
};

template <> struct ArithTraits<Kokkos::complex<float>> {
  typedef Kokkos::complex<float> val_type;
  typedef float mag_type;

  enum : bool { is_complex = true };
  static KOKKOS_FORCEINLINE_FUNCTION mag_type abs(const val_type &x) { return Kokkos::abs(x); }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type real(const val_type &x) { return x.real(); }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type imag(const val_type &x) { return x.imag(); }
  static KOKKOS_FORCEINLINE_FUNCTION val_type conj(const val_type &x) { return Kokkos::conj(x); }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type epsilon() { return FLT_EPSILON; }
  static KOKKOS_FORCEINLINE_FUNCTION void set_real(val_type &x, const mag_type &val) { x.real(val); }
  static KOKKOS_FORCEINLINE_FUNCTION void set_imag(val_type &x, const mag_type &val) { x.imag(val); }
};

template <> struct ArithTraits<Kokkos::complex<double>> {
  typedef Kokkos::complex<double> val_type;
  typedef double mag_type;

  enum : bool { is_complex = true };
  static KOKKOS_FORCEINLINE_FUNCTION mag_type abs(const val_type &x) { return Kokkos::abs(x); }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type real(const val_type &x) { return x.real(); }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type imag(const val_type &x) { return x.imag(); }
  static KOKKOS_FORCEINLINE_FUNCTION val_type conj(const val_type &x) { return Kokkos::conj(x); }
  static KOKKOS_FORCEINLINE_FUNCTION mag_type epsilon() { return DBL_EPSILON; }
  static KOKKOS_FORCEINLINE_FUNCTION void set_real(val_type &x, const mag_type &val) { x.real(val); }
  static KOKKOS_FORCEINLINE_FUNCTION void set_imag(val_type &x, const mag_type &val) { x.imag(val); }
};

} // namespace Tacho

#endif
