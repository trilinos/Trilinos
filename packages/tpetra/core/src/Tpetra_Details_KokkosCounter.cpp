/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// ************************************************************************
// @HEADER
*/
#include "Tpetra_Details_KokkosCounter.hpp"
#include "TpetraCore_config.h"
#include "Kokkos_Core.hpp"
#include "Teuchos_TestForException.hpp"
#include <cstring>

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
        if(eid.instance_id == Impl::int_for_synchronization_reason(SpecialSynchronizationCases::GlobalDeviceSynchronization))
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



} // namespace Details
} // namespace Tpetra

