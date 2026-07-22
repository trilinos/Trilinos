// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHALANX_VIRTUAL_FUNCTION_ON_DEVICE_HPP
#define PHALANX_VIRTUAL_FUNCTION_ON_DEVICE_HPP

#include "Phalanx_KokkosDeviceTypes.hpp"
#include <memory>

namespace PHX {

  /// Struct for deleting device instantiation
  template<typename Device>
  struct DeviceDeleter {
    template<typename T>
    void operator()(T* ptr) {
      Kokkos::parallel_for(Kokkos::RangePolicy<typename Device::execution_space>(0,1),
                           KOKKOS_LAMBDA (const int i) { ptr->~T(); });
      typename Device::execution_space().fence();
      Kokkos::kokkos_free<typename Device::memory_space>(ptr);
    }
  };

  /// Function for creating a vtable on device (requires copy ctor for
  /// derived object). Allocates device memory and must be called from
  /// host.
  template<typename Device,typename Derived>
  std::unique_ptr<Derived,DeviceDeleter<Device>>
  copy_virtual_class_to_device(const Derived& host_source)
  {
    auto* p = static_cast<Derived*>(Kokkos::kokkos_malloc<typename Device::memory_space>(sizeof(Derived)));
    Kokkos::parallel_for(Kokkos::RangePolicy<typename Device::execution_space>(0,1),
                         KOKKOS_LAMBDA (const int i) {new (p) Derived(host_source); });
    typename Device::execution_space().fence();
    return std::unique_ptr<Derived,DeviceDeleter<Device>>(p);
  }

  /// Struct for holding pointers to objects in a Kokkos::View. Used
  /// for putting virtual functions on device. We can't create a
  /// pointer as the Scalar type since the "*" is used to show
  /// rank. Need to wrap pointers in a struct.
  template<typename T>
  struct DevicePtrWrapper {
    T* ptr;
  };

}

#endif
