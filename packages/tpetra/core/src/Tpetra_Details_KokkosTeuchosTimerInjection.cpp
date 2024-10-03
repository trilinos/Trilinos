// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_Details_KokkosTeuchosTimerInjection.hpp"
#include "TpetraCore_config.h"
#include "Tpetra_Details_Behavior.hpp"
#include "Kokkos_Core.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_RCP.hpp"
#ifdef HAVE_TEUCHOS_ADD_TIME_MONITOR_TO_STACKED_TIMER
#include "Teuchos_StackedTimer.hpp"
#include <sstream>
#endif
#include <string>

#include "Tpetra_Details_Behavior.hpp"



namespace {
  // Get a useful label from the deviceId
  // NOTE: Relevant code is in: kokkos/core/src/impl/Kokkos_Profiling_Interface.hpp
  std::string deviceIdToString(const uint32_t deviceId) {
    using namespace Kokkos::Tools::Experimental;
    std::string device_label("(");
    ExecutionSpaceIdentifier eid = identifier_from_devid(deviceId);
    if      (eid.type == DeviceType::Serial)       device_label+="Serial";
    else if (eid.type == DeviceType::OpenMP)       device_label+="OpenMP";
    else if (eid.type == DeviceType::Cuda)         device_label+="Cuda";
    else if (eid.type == DeviceType::HIP)          device_label+="HIP";
    else if (eid.type == DeviceType::OpenMPTarget) device_label+="OpenMPTarget";
    else if (eid.type == DeviceType::HPX)          device_label+="HPX";
    else if (eid.type == DeviceType::Threads)      device_label+="Threads";
    else if (eid.type == DeviceType::SYCL)         device_label+="SYCL";
    else if (eid.type == DeviceType::OpenACC)      device_label+="OpenACC";
    else if (eid.type == DeviceType::Unknown)      device_label+="Unknown";
    else                                           device_label+="Unknown to Tpetra";
#if KOKKOS_VERSION >= 40499
    if(eid.instance_id == int_for_synchronization_reason(SpecialSynchronizationCases::GlobalDeviceSynchronization))
      device_label += " All Instances)";
    else if(eid.instance_id == int_for_synchronization_reason(SpecialSynchronizationCases::DeepCopyResourceSynchronization))
      device_label += " DeepCopyResource)";
#else
    if(eid.instance_id == Impl::int_for_synchronization_reason(SpecialSynchronizationCases::GlobalDeviceSynchronization))
      device_label += " All Instances)";
    else if(eid.instance_id == Impl::int_for_synchronization_reason(SpecialSynchronizationCases::DeepCopyResourceSynchronization))
      device_label += " DeepCopyResource)";
#endif
    else
      device_label += " Instance " + std::to_string(eid.instance_id) + ")";

    return device_label;
  }

  void overlappingWarning() {
    std::ostringstream warning;
    warning <<
      "\n*********************************************************************\n"
      "WARNING: Overlapping timers detected!\n"
      "A TimeMonitor timer was stopped before a nested subtimer was\n"
      "stopped. This is not allowed by the StackedTimer. This corner case\n"
      "typically occurs if the TimeMonitor is stored in an RCP and the RCP is\n"
      "assigned to a new timer. To disable this warning, either fix the\n"
            "ordering of timer creation and destuction or disable the StackedTimer\n";
    std::cout << warning.str() << std::endl;
  }
    
}// anonymous space


namespace Tpetra {
namespace Details {

  namespace DeepCopyTimerInjection {
    Teuchos::RCP<Teuchos::Time> timer_;
    bool initialized_ = false;

    void kokkosp_begin_deep_copy(Kokkos::Tools::SpaceHandle dst_handle, const char* dst_name, const void* dst_ptr,                                 
                                 Kokkos::Tools::SpaceHandle src_handle, const char* src_name, const void* src_ptr,
                                 uint64_t size) {      
      // In verbose mode, we add the src/dst names as well
      std::string extra_label;
      if(Tpetra::Details::Behavior::timeKokkosDeepCopyVerbose1()) {
        extra_label = std::string(" {") + src_name + "=>" + dst_name + "}";
      } else if(Tpetra::Details::Behavior::timeKokkosDeepCopyVerbose2()) {
        extra_label = std::string(" {") + src_name + "=>" + dst_name + "," + std::to_string(size)+"}";
      }    

      if(timer_ != Teuchos::null)
        std::cout << "WARNING: Kokkos::deep_copy() started within another Kokkos::deep_copy().  Timers will be in error"<<std::endl;

      // If the src_name is "Scalar" or "(none)" then we're doing a "Fill" style copy from host to devices, which we want to record separately.  
      if(!strcmp(src_name,"Scalar") || !strcmp(src_name,"(none)")) 
        timer_ = Teuchos::TimeMonitor::getNewTimer(std::string("Kokkos::deep_copy_scalar [")+src_handle.name+"=>"+dst_handle.name+"]" + extra_label);
      // If the size is under 65 bytes, we're going to flag this as "small" to make it easier to watch the big stuff
      else if(size <= 64)
        timer_ = Teuchos::TimeMonitor::getNewTimer(std::string("Kokkos::deep_copy_small [")+src_handle.name+"=>"+dst_handle.name+"]" + extra_label);
      else
        timer_ = Teuchos::TimeMonitor::getNewTimer(std::string("Kokkos::deep_copy [")+src_handle.name+"=>"+dst_handle.name+"]" + extra_label);
      timer_->start();
      timer_->incrementNumCalls();
#ifdef HAVE_TEUCHOS_ADD_TIME_MONITOR_TO_STACKED_TIMER
      const auto stackedTimer = Teuchos::TimeMonitor::getStackedTimer();
      if (nonnull(stackedTimer))
        stackedTimer->start(timer_->name());
#endif
    }

    void kokkosp_end_deep_copy() {
      if (timer_ != Teuchos::null) {
        timer_->stop();
#ifdef HAVE_TEUCHOS_ADD_TIME_MONITOR_TO_STACKED_TIMER
        try {
          const auto stackedTimer = Teuchos::TimeMonitor::getStackedTimer();
          if (nonnull(stackedTimer))
            stackedTimer->stop(timer_->name());
        }
        catch (std::runtime_error&) {
          overlappingWarning();
          Teuchos::TimeMonitor::setStackedTimer(Teuchos::null);
        }
#endif
      }      
      timer_ = Teuchos::null;
    }

  }// end DeepCopyTimerInjection

  void AddKokkosDeepCopyToTimeMonitor(bool force) {
    if (!DeepCopyTimerInjection::initialized_) {
      if (force || Tpetra::Details::Behavior::timeKokkosDeepCopy() || Tpetra::Details::Behavior::timeKokkosDeepCopyVerbose1()
                                                                   || Tpetra::Details::Behavior::timeKokkosDeepCopyVerbose2()) {
        Kokkos::Tools::Experimental::set_begin_deep_copy_callback(DeepCopyTimerInjection::kokkosp_begin_deep_copy);
        Kokkos::Tools::Experimental::set_end_deep_copy_callback(DeepCopyTimerInjection::kokkosp_end_deep_copy);
        DeepCopyTimerInjection::initialized_=true;
      }
    }
  }

  
  namespace FenceTimerInjection {
    Teuchos::RCP<Teuchos::Time> timer_;
    bool initialized_ = false;
    uint64_t active_handle;

    void kokkosp_begin_fence(const char* name, const uint32_t deviceId,
                             uint64_t* handle) {

      // Nested fences are not allowed
      if(timer_ != Teuchos::null)
        return;        
      active_handle = (active_handle+1) % 1024;
      *handle = active_handle;

      std::string device_label = deviceIdToString(deviceId);

      timer_ = Teuchos::TimeMonitor::getNewTimer(std::string("Kokkos::fence ")+name + " " + device_label);
      timer_->start();
      timer_->incrementNumCalls();
#ifdef HAVE_TEUCHOS_ADD_TIME_MONITOR_TO_STACKED_TIMER
      const auto stackedTimer = Teuchos::TimeMonitor::getStackedTimer();
      if (nonnull(stackedTimer))
        stackedTimer->start(timer_->name());
#endif
      
    }


    void kokkosp_end_fence(const uint64_t handle) {
      if(handle == active_handle) {
        if (timer_ != Teuchos::null) {
          timer_->stop();
#ifdef HAVE_TEUCHOS_ADD_TIME_MONITOR_TO_STACKED_TIMER
          try {
            const auto stackedTimer = Teuchos::TimeMonitor::getStackedTimer();
            if (nonnull(stackedTimer))
              stackedTimer->stop(timer_->name());
          }
          catch (std::runtime_error&) {
            overlappingWarning();
            Teuchos::TimeMonitor::setStackedTimer(Teuchos::null);
          }
#endif
        }        
        timer_ = Teuchos::null;        
      }
      // Else: We've nested our fences, and we need to ignore the inner fences
    }


  }//end FenceTimerInjection

  void AddKokkosFenceToTimeMonitor(bool force) {
    if (!FenceTimerInjection::initialized_) {
      if (force || Tpetra::Details::Behavior::timeKokkosFence()) {
        Kokkos::Tools::Experimental::set_begin_fence_callback(FenceTimerInjection::kokkosp_begin_fence);
        Kokkos::Tools::Experimental::set_end_fence_callback(FenceTimerInjection::kokkosp_end_fence);
        FenceTimerInjection::initialized_=true;
      }
    }
  }


  namespace FunctionsTimerInjection {
    Teuchos::RCP<Teuchos::Time> timer_;
    bool initialized_ = false;

    void kokkosp_begin_kernel(const char* kernelName, const char* kernelPrefix, const uint32_t devID,
                              uint64_t* kernelID) {
      // Nested fences are not allowed
      if(timer_ != Teuchos::null)
        return;        
      std::string device_label = deviceIdToString(devID);

      timer_ = Teuchos::TimeMonitor::getNewTimer(std::string("Kokkos::")+ kernelName + " " +kernelPrefix + " " + device_label);
      timer_->start();
      timer_->incrementNumCalls();
#ifdef HAVE_TEUCHOS_ADD_TIME_MONITOR_TO_STACKED_TIMER
      const auto stackedTimer = Teuchos::TimeMonitor::getStackedTimer();
      if (nonnull(stackedTimer))
        stackedTimer->start(timer_->name());
#endif
      
    }

    void kokkosp_begin_for(const char* kernelPrefix, const uint32_t devID, uint64_t* kernelID) {
      kokkosp_begin_kernel("parallel_for",kernelPrefix,devID,kernelID);
    }

    void kokkosp_begin_scan(const char* kernelPrefix, const uint32_t devID, uint64_t* kernelID) {
      kokkosp_begin_kernel("parallel_scan",kernelPrefix,devID,kernelID);
    }

    void kokkosp_begin_reduce(const char* kernelPrefix, const uint32_t devID, uint64_t* kernelID) {
      kokkosp_begin_kernel("parallel_reduce",kernelPrefix,devID,kernelID);
    }
                              
    void kokkosp_end_kernel(const uint64_t handle) {
      if (timer_ != Teuchos::null) {
        timer_->stop();
#ifdef HAVE_TEUCHOS_ADD_TIME_MONITOR_TO_STACKED_TIMER
        try {
          const auto stackedTimer = Teuchos::TimeMonitor::getStackedTimer();
          if (nonnull(stackedTimer))
            stackedTimer->stop(timer_->name());
        }
        catch (std::runtime_error&) {
          overlappingWarning();
          Teuchos::TimeMonitor::setStackedTimer(Teuchos::null);
        }
#endif
      }
      
      timer_ = Teuchos::null;      
    }
  }//end FunctionsInjection

  void AddKokkosFunctionsToTimeMonitor(bool force) {
    if (!FunctionsTimerInjection::initialized_) {
      if (force || Tpetra::Details::Behavior::timeKokkosFunctions()) {
        Kokkos::Tools::Experimental::set_begin_parallel_for_callback(FunctionsTimerInjection::kokkosp_begin_for);
        Kokkos::Tools::Experimental::set_begin_parallel_reduce_callback(FunctionsTimerInjection::kokkosp_begin_reduce);
        Kokkos::Tools::Experimental::set_begin_parallel_scan_callback(FunctionsTimerInjection::kokkosp_begin_scan);

        // The end-call is generic, even though the start-call is not.
        Kokkos::Tools::Experimental::set_end_parallel_for_callback(FunctionsTimerInjection::kokkosp_end_kernel);
        Kokkos::Tools::Experimental::set_end_parallel_reduce_callback(FunctionsTimerInjection::kokkosp_end_kernel);
        Kokkos::Tools::Experimental::set_end_parallel_scan_callback(FunctionsTimerInjection::kokkosp_end_kernel);
        FunctionsTimerInjection::initialized_=true;
      }
    }
  }



} // namespace Details
} // namespace Tpetra

