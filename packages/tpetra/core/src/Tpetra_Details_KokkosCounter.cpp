// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// clang-format off
#include "Tpetra_Details_KokkosCounter.hpp"
#include "TpetraCore_config.h"
#include "Kokkos_Core.hpp"
#include "Teuchos_TestForException.hpp"
#include <cstring>
#include <string>

namespace Tpetra {
namespace Details {


  /***************************** Deep Copy *****************************/
  namespace DeepCopyCounterDetails {
    // Static variables
    bool is_initialized=true;
    size_t count_same=0;
    size_t count_different=0;
    bool count_active=false;
       
    void kokkosp_begin_deep_copy(Kokkos::Tools::SpaceHandle dst_handle, const char* dst_name, const void* dst_ptr,                                 
                                 Kokkos::Tools::SpaceHandle src_handle, const char* src_name, const void* src_ptr,
                                 uint64_t size) {

      if(count_active) {
        if(strcmp(dst_handle.name,src_handle.name))
          count_different++;
        else
          count_same++;        
      }
    }

  }// end DeepCopyCounterDetails


  void DeepCopyCounter::start() {
    DeepCopyCounterDetails::count_active=true;
    Kokkos::Tools::Experimental::set_begin_deep_copy_callback(DeepCopyCounterDetails::kokkosp_begin_deep_copy);
  }

  void DeepCopyCounter::reset() {
    DeepCopyCounterDetails::count_same=0;
    DeepCopyCounterDetails::count_different=0;
  }

  void DeepCopyCounter::stop() {
    DeepCopyCounterDetails::count_active=false;
  }

  size_t DeepCopyCounter::get_count_same_space() {
    return DeepCopyCounterDetails::count_same;
  }

  size_t DeepCopyCounter::get_count_different_space() {
    return DeepCopyCounterDetails::count_different;
  }



  /***************************** Fence *****************************/


  namespace FenceCounterDetails {

    // Static variables
    bool is_initialized=false;
    bool count_active=false;
    std::vector<size_t> count_instance;
    std::vector<size_t> count_global;
    int num_devices=0;


    void kokkosp_begin_fence(const char* name, const uint32_t deviceId,
                             uint64_t* handle) {

      if(count_active) {
        using namespace Kokkos::Tools::Experimental;
        ExecutionSpaceIdentifier eid = identifier_from_devid(deviceId);
        
        // Figure out what count bin to stick this in
        int idx = (int) eid.type;
#if KOKKOS_VERSION >= 40499
        if(eid.instance_id == int_for_synchronization_reason(SpecialSynchronizationCases::GlobalDeviceSynchronization))
#else
        if(eid.instance_id == Impl::int_for_synchronization_reason(SpecialSynchronizationCases::GlobalDeviceSynchronization))
#endif
          count_global[idx]++;
        else
          count_instance[idx]++;
      }
    }


    std::string get_label(int i) {
      using namespace Kokkos::Tools::Experimental;
      DeviceType i_type = devicetype_from_uint32t(i);
      std::string device_label;
      if      (i_type == DeviceType::Serial)       device_label="Serial";
      else if (i_type == DeviceType::OpenMP)       device_label="OpenMP";
      else if (i_type == DeviceType::Cuda)         device_label="Cuda";
      else if (i_type == DeviceType::HIP)          device_label="HIP";
      else if (i_type == DeviceType::OpenMPTarget) device_label="OpenMPTarget";
      else if (i_type == DeviceType::HPX)          device_label="HPX";
      else if (i_type == DeviceType::Threads)      device_label="Threats";
      else if (i_type == DeviceType::SYCL)         device_label="SYCL";
      else if (i_type == DeviceType::OpenACC)      device_label="OpenACC";
      else if (i_type == DeviceType::Unknown)      device_label="Unknown";
      
      return device_label;
    }

    void initialize() {
      using namespace Kokkos::Tools::Experimental;
      num_devices = (int) DeviceType::Unknown;
      count_instance.resize(num_devices);
      count_instance.assign(num_devices,0);
      count_global.resize(num_devices);
      count_global.assign(num_devices,0);
      is_initialized=true;
    }

  }// end FenceCounterDetails


  

  void FenceCounter::start() {
    if(!FenceCounterDetails::is_initialized) 
      FenceCounterDetails::initialize();
    FenceCounterDetails::count_active=true;
    Kokkos::Tools::Experimental::set_begin_fence_callback(FenceCounterDetails::kokkosp_begin_fence);
  }

  void FenceCounter::reset() {
    FenceCounterDetails::count_instance.assign(FenceCounterDetails::num_devices,0);
    FenceCounterDetails::count_global.assign(FenceCounterDetails::num_devices,0);
  }
  
  void FenceCounter::stop() {
    FenceCounterDetails::count_active=false;
  }

  size_t FenceCounter::get_count_global(const std::string & device) {
    using namespace Kokkos::Tools::Experimental;
    for(int i=0;i<FenceCounterDetails::num_devices; i++) {
      std::string device_label = FenceCounterDetails::get_label(i);      

      if(device == device_label)
        return FenceCounterDetails::count_global[i];
    }

    // Haven't found a device by this name
    TEUCHOS_TEST_FOR_EXCEPTION(1,std::runtime_error,std::string("Error: ") + device + std::string(" is not a device known to Tpetra"));
  }


  size_t FenceCounter::get_count_instance(const std::string & device) {
    using namespace Kokkos::Tools::Experimental;
    for(int i=0;i<FenceCounterDetails::num_devices; i++) {
      std::string device_label = FenceCounterDetails::get_label(i);      

      if(device == device_label)
        return FenceCounterDetails::count_instance[i];
    }

    // Haven't found a device by this name
    TEUCHOS_TEST_FOR_EXCEPTION(1,std::runtime_error,std::string("Error: ") + device + std::string(" is not a device known to Tpetra"));
  }

// clang-format on
namespace KokkosRegionCounterDetails {
std::vector<std::string> regions;

void push_region_callback(const char *label) { regions.push_back(label); }
static_assert(std::is_same_v<decltype(&push_region_callback),
                             Kokkos_Profiling_pushFunction>,
              "Unexpected Kokkos profiling interface API. This is an internal "
              "Tpetra developer error, please report this.");

} // namespace KokkosRegionCounterDetails

void KokkosRegionCounter::start() {
  Kokkos::Tools::Experimental::set_push_region_callback(
      KokkosRegionCounterDetails::push_region_callback);
}

void KokkosRegionCounter::reset() {
  KokkosRegionCounterDetails::regions.clear();
}

void KokkosRegionCounter::stop() {
  Kokkos::Tools::Experimental::set_push_region_callback(nullptr);
}

size_t
KokkosRegionCounter::get_count_region_contains(const std::string &needle) {
  size_t count = 0;
  for (const auto &region : KokkosRegionCounterDetails::regions) {
    count += (region.find(needle) != std::string::npos);
  }
  return count;
}

void KokkosRegionCounter::dump_regions(Teuchos::FancyOStream &os) {
  for (const auto &region : KokkosRegionCounterDetails::regions) {
    os << region << "\n";
  }
}

void KokkosRegionCounter::dump_regions(std::ostream &os) {
  for (const auto &region : KokkosRegionCounterDetails::regions) {
    os << region << "\n";
  }
}


// clang-format off


} // namespace Details
} // namespace Tpetra

