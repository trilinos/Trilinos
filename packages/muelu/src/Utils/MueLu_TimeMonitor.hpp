// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_TIMEMONITOR_HPP
#define MUELU_TIMEMONITOR_HPP

#include <string>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_VerboseObject.hpp"
#include "MueLu_MutuallyExclusiveTime.hpp"

namespace MueLu {

  // Helper function. Similar to Teuchos::TimeMonitor::summarize().
  ArrayRCP<double> ReduceMaxMinAvg(double localValue, Teuchos::Comm<int> const &comm, int rootNode = 0);

  /*! @class TimeMonitor
  
      @brief Integrates Teuchos::TimeMonitor with MueLu verbosity system.
  */
  class TimeMonitor : public BaseClass {

  public:

    TimeMonitor(const BaseClass& object, const std::string& msg, MsgType timerLevel = Timings0)
    {

      // Inherit verbosity from 'object'
      SetVerbLevel(object.GetVerbLevel());
      setOStream(object.getOStream());

      if (IsPrint(timerLevel) &&
          /* disable timer if never printed: */ (IsPrint(RuntimeTimings) || (!IsPrint(NoTimeReport)))) {

        if (!IsPrint(NoTimeReport)) {
          // TODO: there is no function to register a timer in Teuchos::TimeMonitor after the creation of the timer. But would be useful...
          timer_ = Teuchos::TimeMonitor::getNewTimer("MueLu: " + msg);
        } else {
          timer_ = rcp(new Teuchos::Time("MueLu: " + msg));
        }

        // Start the timer (this is what is done by Teuchos::TimeMonitor)
        timer_->start();
        timer_->incrementNumCalls();
      }
    }

    ~TimeMonitor() {
      if (timer_ != Teuchos::null) {

        // Stop the timer
        timer_->stop();

        if (IsPrint(RuntimeTimings)) {
          //FIXME: creates lot of barriers. An option to report time of proc0 only instead would be nice
          //FIXME: MPI_COMM_WORLD only... BTW, it is also the case in Teuchos::TimeMonitor...
          //
          // mfh 11 Nov 2012: Actually, Teuchos::TimeMonitor::summarize() has multiple overloads that take a Teuchos::Comm.
          ArrayRCP<double> stats = ReduceMaxMinAvg(timer_->totalElapsedTime(), *Teuchos::DefaultComm<int>::getComm ());

          //FIXME: Not very important for now, but timer will be printed even if verboseLevel of Monitor/Object changed
          //       between Monitor constructor and destructor.
          GetOStream(RuntimeTimings, 0) << "Timer: " << " max=" << stats[0] << " min=" << stats[1] << " avg=" << stats[2] << std::endl;
        }
      }
    }

  protected:
    TimeMonitor() { }

  private:
    RCP<Teuchos::Time> timer_;
  }; //class TimeMonitor

  //TODO: code duplication MutuallyExclusiveTimeMonitor / TimeMonitor

  /*! @class MutuallyExclusiveTimeMonitor

      @brief Similar to TimeMonitor, but uses MutuallyExclusiveTime objects.

  */
  template <class TagName>
  class MutuallyExclusiveTimeMonitor : public BaseClass {

  public:

    MutuallyExclusiveTimeMonitor(const BaseClass& object, const std::string& msg, MsgType timerLevel = Timings0)
    {
      // Inherit verbosity from 'object'
      SetVerbLevel(object.GetVerbLevel());
      setOStream(object.getOStream());

      if (IsPrint(timerLevel) &&
          /* disable timer if never printed: */ (IsPrint(RuntimeTimings) || (!IsPrint(NoTimeReport)))) {

        if (!IsPrint(NoTimeReport)) {
          timer_ = MutuallyExclusiveTime<TagName>::getNewTimer("MueLu: " + msg /*+ " (MutuallyExclusive)" */);
        } else {
          timer_ = rcp(new MutuallyExclusiveTime<TagName>     ("MueLu: " + msg /*+ " (MutuallyExclusive)" */));
        }

        timer_->start();
        timer_->incrementNumCalls();
      }
    }

    ~MutuallyExclusiveTimeMonitor() {
      if (timer_ != Teuchos::null) {

        // Stop the timer
        timer_->stop();

        if (IsPrint(RuntimeTimings)) {
          //FIXME: creates lot of barriers. An option to report time of proc0 only instead would be nice
          //FIXME: MPI_COMM_WORLD only... BTW, it is also the case in Teuchos::TimeMonitor...
          //TODO          ArrayRCP<double> stats = ReduceMaxMinAvg(timer_->totalElapsedTime(), *Teuchos::DefaultComm<int>::getComm ());
          //
          // mfh 11 Nov 2012: Actually, Teuchos::TimeMonitor::summarize() has multiple overloads that take a Teuchos::Comm.

          //FIXME: Not very important for now, but timer will be printed even if verboseLevel of Monitor/Object changed
          //       between Monitor constructor and destructor.
          //TODO GetOStream(RuntimeTimings, 0) << "Timer: " << " max=" << stats[0] << " min=" << stats[1] << " avg=" << stats[2] << std::endl;
        }
      }
    }

  protected:
    MutuallyExclusiveTimeMonitor() { }

  private:
    RCP<MutuallyExclusiveTime<TagName> > timer_; // keep a reference on the timer to print stats if RuntimeTimings=ON //TODO:use base class instead
  }; //class MutuallyExclusiveTimeMonitor

} // namespace MueLu

#endif // MUELU_TIMEMONITOR_HPP


