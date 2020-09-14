#ifndef __TACHO_HPP__
#define __TACHO_HPP__

#include "Tacho_config.h" 

#include "Kokkos_Core.hpp"
#include "impl/Kokkos_Timer.hpp"

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
    using default_exec_space = Kokkos::DefaultExecutionSpace;
    using default_memory_space = typename default_exec_space::memory_space;
    using device_type = Kokkos::Device<default_exec_space,default_memory_space>;
  };

  /// until kokkos dual view issue is resolved, we follow the default space in Trilinos (uvm)
#if defined(KOKKOS_ENABLE_CUDA)
  template<>
  struct UseThisDevice<Kokkos::Cuda> { using device_type = Kokkos::Device<Kokkos::Cuda,Kokkos::CudaSpace>; };
#endif
#if defined(KOKKOS_ENABLE_OPENMP)
  template<>
  struct UseThisDevice<Kokkos::OpenMP> { using device_type = Kokkos::Device<Kokkos::OpenMP,Kokkos::HostSpace>; };
#endif
#if defined(KOKKOS_ENABLE_SERIAL)
  template<>
  struct UseThisDevice<Kokkos::Serial> { using device_type = Kokkos::Device<Kokkos::Serial,Kokkos::HostSpace>; };
#endif

  ///
  /// print execution spaces
  ///
  template<typename SpT>
  void printExecSpaceConfiguration(std::string name, const bool detail = false) {
    if (!Kokkos::Impl::is_space<SpT>::value) {
      std::string msg("SpT is not Kokkos execution space");
      fprintf(stderr, ">> Error in file %s, line %d\n",__FILE__,__LINE__);
      fprintf(stderr, "   %s\n", msg.c_str());
      throw std::logic_error(msg.c_str());
    }
    std::cout << std::setw(16) << name << "::  ";
    SpT::print_configuration(std::cout, detail);
  }

}

#endif
