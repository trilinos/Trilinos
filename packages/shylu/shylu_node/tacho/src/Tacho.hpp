#ifndef __TACHO_HPP__
#define __TACHO_HPP__

#include "Tacho_config.h" 

#include "Kokkos_Core.hpp"
#include "Kokkos_Timer.hpp"

#include <cstddef>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <memory>
#include <string>
#include <stdexcept>
#include <vector>

/// \file Tacho.hpp
/// \brief Header to be included by users
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

  ///
  /// default ordinal and size type
  ///

#if defined( TACHO_USE_INT_INT )
  typedef int ordinal_type;
  typedef int size_type;
#elif defined( TACHO_USE_INT_SIZE_T )
  typedef int ordinal_type;
  typedef size_t size_type;
#else
  typedef int ordinal_type;
  typedef size_t size_type;
#endif

  ///
  /// default device type used in tacho
  ///
  template<typename ExecSpace>
  struct UseThisDevice {
    using exec_space = ExecSpace;
    using memory_space = typename exec_space::memory_space;
    using type = Kokkos::Device<exec_space,memory_space>;
    using device_type = type;
  };

  template<typename ExecSpace>
  struct UseThisScheduler {
    using type = Kokkos::TaskSchedulerMultiple<ExecSpace>;
    using scheduler_type = type;
  };

  /// until kokkos dual view issue is resolved, we follow the default space in Trilinos (uvm)
#if defined(KOKKOS_ENABLE_CUDA)
  template<>
  struct UseThisDevice<Kokkos::Cuda> { 
    using type = Kokkos::Device<Kokkos::Cuda,Kokkos::CudaSpace>; 
    using device_type = type;
  };
#endif
#if defined(KOKKOS_ENABLE_OPENMP)
  template<>
  struct UseThisDevice<Kokkos::OpenMP> { 
    using type = Kokkos::Device<Kokkos::OpenMP,Kokkos::HostSpace>; 
    using device_type = type; 
  };
#endif
#if defined(KOKKOS_ENABLE_SERIAL)
  template<>
  struct UseThisDevice<Kokkos::Serial> { 
    using type = Kokkos::Device<Kokkos::Serial,Kokkos::HostSpace>;
    using device_type = type;
  };
#endif

  ///
  /// print execution spaces
  ///
  template<typename SpT>
  void printExecSpaceConfiguration(std::string name, const bool detail = false) {
    if (!Kokkos::is_space<SpT>::value) {
      std::string msg("SpT is not Kokkos execution space");
      fprintf(stderr, ">> Error in file %s, line %d\n",__FILE__,__LINE__);
      fprintf(stderr, "   %s\n", msg.c_str());
      throw std::logic_error(msg.c_str());
    }
    std::cout << std::setw(16) << name << "::  ";
    SpT::print_configuration(std::cout, detail);
  }

  template<typename T>
  struct ArithTraits;

  template<>
  struct ArithTraits<float> {
    typedef float val_type;
    typedef float mag_type;

    enum : bool { is_complex = false };
    static KOKKOS_FORCEINLINE_FUNCTION mag_type abs (const val_type& x) { return x > 0 ? x : -x; }
    static KOKKOS_FORCEINLINE_FUNCTION mag_type real(const val_type& x) { return x; }
    static KOKKOS_FORCEINLINE_FUNCTION mag_type imag(const val_type& x) { return x; }
    static KOKKOS_FORCEINLINE_FUNCTION val_type conj(const val_type& x) { return x; }
  };

  template<>
  struct ArithTraits<double> {
    typedef double val_type;
    typedef double mag_type;

    enum : bool { is_complex = false };
    static KOKKOS_FORCEINLINE_FUNCTION mag_type abs (const val_type& x) { return x > 0 ? x : -x; }
    static KOKKOS_FORCEINLINE_FUNCTION mag_type real(const val_type& x) { return x; }
    static KOKKOS_FORCEINLINE_FUNCTION mag_type imag(const val_type& x) { return x; }
    static KOKKOS_FORCEINLINE_FUNCTION val_type conj(const val_type& x) { return x; }
  };

  template<>
  struct ArithTraits<std::complex<float> > {
    typedef std::complex<float> val_type;
    typedef float mag_type;

    enum : bool { is_complex = true };
    static inline mag_type abs (const val_type& x) { return std::abs(x); }
    static inline mag_type real(const val_type& x) { return x.real(); }
    static inline mag_type imag(const val_type& x) { return x.imag(); }
    static inline val_type conj(const val_type& x) { return std::conj(x); }
  };

  template<>
  struct ArithTraits<std::complex<double> > {
    typedef std::complex<double> val_type;
    typedef double mag_type;

    enum : bool { is_complex = true };
    static inline mag_type abs (const val_type& x) { return std::abs(x); }
    static inline mag_type real(const val_type& x) { return x.real(); }
    static inline mag_type imag(const val_type& x) { return x.imag(); }
    static inline val_type conj(const val_type& x) { return std::conj(x); }
  };
    
  template<>
  struct ArithTraits<Kokkos::complex<float> > {
    typedef Kokkos::complex<float> val_type;
    typedef float mag_type;

    enum : bool { is_complex = true };
    static KOKKOS_FORCEINLINE_FUNCTION mag_type abs (const val_type& x) { return Kokkos::abs(x); }
    static KOKKOS_FORCEINLINE_FUNCTION mag_type real(const val_type& x) { return x.real(); }
    static KOKKOS_FORCEINLINE_FUNCTION mag_type imag(const val_type& x) { return x.imag(); }
    static KOKKOS_FORCEINLINE_FUNCTION val_type conj(const val_type& x) { return Kokkos::conj(x); }
  };

  template<>
  struct ArithTraits<Kokkos::complex<double> > {
    typedef Kokkos::complex<double> val_type;
    typedef double mag_type;

    enum : bool { is_complex = true };
    static KOKKOS_FORCEINLINE_FUNCTION mag_type abs (const val_type& x) { return Kokkos::abs(x); }
    static KOKKOS_FORCEINLINE_FUNCTION mag_type real(const val_type& x) { return x.real(); }
    static KOKKOS_FORCEINLINE_FUNCTION mag_type imag(const val_type& x) { return x.imag(); }
    static KOKKOS_FORCEINLINE_FUNCTION val_type conj(const val_type& x) { return Kokkos::conj(x); }
  };

}

#endif
