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
#ifndef MUELU_MONITOR_HPP
#define MUELU_MONITOR_HPP

#include <string>
#include <Teuchos_Time.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_VerboseObject.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_MutuallyExclusiveTime.hpp"
#include "MueLu_TimeMonitor.hpp"

namespace MueLu {

  //! Manage indentation of output using OSTab and verbosity level
  class PrintMonitor : public BaseClass {

  public:

    PrintMonitor(const BaseClass& object, const std::string& msg, MsgType msgLevel = Runtime0) {

      // Inherit verbosity from 'object'
      SetVerbLevel(object.GetVerbLevel());
      setOStream(object.getOStream());

      // Print description and new indent
      if (IsPrint(msgLevel)) {
        GetOStream(msgLevel, 0) << msg << std::endl;
        tab_ = rcp(new Teuchos::OSTab(getOStream()));
      }

    }

    ~PrintMonitor() { }

  protected:
    PrintMonitor() { }

  private:
    RCP<Teuchos::OSTab> tab_;
  };

  // Main Monitor
  // A timer is created only if 'timerLevel' (Timings0 by default) is enabled
  class Monitor: public BaseClass {
  public:
    Monitor(const BaseClass& object, const std::string & msg, MsgType msgLevel = Runtime0, MsgType timerLevel = Timings0)
      : printMonitor_(object, msg + " (" + object.description() + ")", msgLevel),
        timerMonitor_(object, object.ShortClassName() + ": " + msg + " (total)",    timerLevel)
    { }

  private:
    PrintMonitor printMonitor_;
    TimeMonitor timerMonitor_;
  };

  // Another version of the Monitor, that does not repeat the object description
  // A timer is created only if 'timerLevel' (Timings1 by default) is enabled
  class SubMonitor: public BaseClass {
  public:
    SubMonitor(const BaseClass& object, const std::string & msg, MsgType msgLevel = Runtime1, MsgType timerLevel = Timings1)
      : printMonitor_(object, msg, msgLevel),
        timerMonitor_(object, object.ShortClassName() + ": " + msg + " (sub, total)", timerLevel)
    { }

  private:
    PrintMonitor printMonitor_;
    TimeMonitor timerMonitor_;
  };

  // Factory monitor
  // Similar to Monitor but add a timer level by level
  class FactoryMonitor: public Monitor {
  public:
    FactoryMonitor(const BaseClass& object, const std::string & msg, int levelID, MsgType msgLevel = Runtime0, MsgType timerLevel = Timings0)
      : Monitor(object, msg, msgLevel, timerLevel),
        timerMonitorExclusive_(object, object.ShortClassName() + ": " + msg, timerLevel)
    {
      if (IsPrint(TimingsByLevel)) {
        levelTimeMonitor_ = rcp(new TimeMonitor(object, object.ShortClassName() + ": " + msg + " (total, level=" + Teuchos::Utils::toString(levelID) + ")", timerLevel));
	levelTimeMonitorExclusive_ = rcp(new MutuallyExclusiveTimeMonitor<Level>(object, object.ShortClassName() + ": " + msg + " (level=" + Teuchos::Utils::toString(levelID) + ")", timerLevel));
      }
    }

    //TODO: code factorization
    FactoryMonitor(const BaseClass& object, const std::string & msg, const Level & level, MsgType msgLevel = Runtime0, MsgType timerLevel = Timings0)
      : Monitor(object, msg, msgLevel, timerLevel),
        timerMonitorExclusive_(object, object.ShortClassName() + ": " + msg, timerLevel)
    {
      if (IsPrint(TimingsByLevel)) {
        levelTimeMonitor_ = rcp(new TimeMonitor(object, object.ShortClassName() + ": " +  msg + " (total, level=" + Teuchos::Utils::toString(level.GetLevelID()) + ")", timerLevel));
	levelTimeMonitorExclusive_ = rcp(new MutuallyExclusiveTimeMonitor<Level>(object, object.ShortClassName() + ": " + msg + " (level=" + Teuchos::Utils::toString(level.GetLevelID()) + ")", timerLevel));
      }
    }

  private:
    RCP<TimeMonitor> levelTimeMonitor_;

    MutuallyExclusiveTimeMonitor<FactoryBase>  timerMonitorExclusive_;
    RCP<MutuallyExclusiveTimeMonitor<Level> >  levelTimeMonitorExclusive_;
  };

  // Factory monitor
  // Similar to SubMonitor but add a timer level by level
  class SubFactoryMonitor: public SubMonitor {
  public:
    SubFactoryMonitor(const BaseClass& object, const std::string & msg, int levelID, MsgType msgLevel = Runtime1, MsgType timerLevel = Timings1)
      : SubMonitor(object, msg, msgLevel, timerLevel)
    {
      if (IsPrint(TimingsByLevel)) {
        levelTimeMonitor_ = rcp(new TimeMonitor(object, object.ShortClassName() + ": " + msg + " (sub, total, level=" + Teuchos::Utils::toString(levelID) + ")", timerLevel));
      }
    }

    SubFactoryMonitor(const BaseClass& object, const std::string & msg, const Level & level, MsgType msgLevel = Runtime1, MsgType timerLevel = Timings1)
      : SubMonitor(object, msg, msgLevel, timerLevel)
    {
      if (IsPrint(TimingsByLevel)) {
        levelTimeMonitor_ = rcp(new TimeMonitor(object, object.ShortClassName() + ": " +  msg + " (sub, total, level=" + Teuchos::Utils::toString(level.GetLevelID()) + ")", timerLevel));
      }
    }
  private:
    RCP<TimeMonitor> levelTimeMonitor_;
  };

} // namespace MueLu

#endif // MUELU_MONITOR_HPP
