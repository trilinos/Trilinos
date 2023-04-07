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
#include "Tpetra_Details_DeepCopyTeuchosTimerInjection.hpp"
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

namespace Tpetra {
namespace Details {

  namespace DeepCopyTimerInjection {
    Teuchos::RCP<Teuchos::Time> timer_;
    bool initialized_ = false;

    void kokkosp_begin_deep_copy(Kokkos::Tools::SpaceHandle dst_handle, const char* dst_name, const void* dst_ptr,                                 
                                 Kokkos::Tools::SpaceHandle src_handle, const char* src_name, const void* src_ptr,
                                 uint64_t size) {      
      if(timer_ != Teuchos::null)
        std::cout << "WARNING: Kokkos::deep_copy() started within another Kokkos::deep_copy().  Timers will be in error"<<std::endl;

      timer_ = Teuchos::TimeMonitor::getNewTimer(std::string("Kokkos::deep_copy [")+src_handle.name+"=>"+dst_handle.name+"]");
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
      if (force || Tpetra::Details::Behavior::timeKokkosDeepCopy()) {
        Kokkos::Tools::Experimental::set_begin_deep_copy_callback(DeepCopyTimerInjection::kokkosp_begin_deep_copy);
        Kokkos::Tools::Experimental::set_end_deep_copy_callback(DeepCopyTimerInjection::kokkosp_end_deep_copy);
        DeepCopyTimerInjection::initialized_=true;
      }
    }
  }


} // namespace Details
} // namespace Tpetra

