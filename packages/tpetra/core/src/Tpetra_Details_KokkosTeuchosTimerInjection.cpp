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

      // If the src_name is "Scalar" then we're doing a "Fill" style copy from host to devices, which we want to record separately.  
      if(!strcmp(src_name,"Scalar")) 
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

      // Get a useful label from the deviceId
      // NOTE: Relevant code is in: kokkos/core/src/impl/Kokkos_Profiling_Interface.hpp
      std::string device_label("(");
      {
        using namespace Kokkos::Tools::Experimental;

        ExecutionSpaceIdentifier eid = identifier_from_devid(deviceId);
        if      (eid.type == DeviceType::Serial)       device_label+="Serial";
        else if (eid.type == DeviceType::OpenMP)       device_label+="OpenMP";
        else if (eid.type == DeviceType::Cuda)         device_label+="Cuda";
        else if (eid.type == DeviceType::HIP)          device_label+="HIP";
        else if (eid.type == DeviceType::OpenMPTarget) device_label+="OpenMPTarget";
        else if (eid.type == DeviceType::HPX)          device_label+="HPX";
        else if (eid.type == DeviceType::Threads)      device_label+="Threats";
        else if (eid.type == DeviceType::SYCL)         device_label+="SYCL";
        else if (eid.type == DeviceType::OpenACC)      device_label+="OpenACC";
        else if (eid.type == DeviceType::Unknown)      device_label+="Unknown";
        else                                           device_label+="Unknown to Tpetra";

        if(eid.instance_id == Impl::int_for_synchronization_reason(SpecialSynchronizationCases::GlobalDeviceSynchronization))
           device_label += " All Instances)";
        else if(eid.instance_id == Impl::int_for_synchronization_reason(SpecialSynchronizationCases::DeepCopyResourceSynchronization))
          device_label += " DeepCopyResource)";
        else
          device_label += " Instance " + std::to_string(eid.instance_id) + ")";
      }

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



} // namespace Details
} // namespace Tpetra

