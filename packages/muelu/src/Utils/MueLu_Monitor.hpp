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

  /*! @class PrintMonitor
    Manages indentation of output using Teuchos::OSTab and verbosity level
  */
  class PrintMonitor : public BaseClass {

  public:

    //! Constructor
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

  /*! @class Monitor
     @brief Main MueLu time monitor.

     A timer is created only if 'timerLevel' (Timings0 by default) is enabled.
     This class manages verbosity, as well as times this object and all its children over all levels.

     The timer yields output such as
      \verbatim
       MueLu: SaPFactory: Prolongator smoothing (total)             6.508 (6)        6.551 (6)        6.557 (6)        1.092 (6)
      \endverbatim
     Note that the keyword \a total denotes timing of the object and its children.
  */
  class Monitor: public BaseClass {
  public:
    Monitor(const BaseClass& object, const std::string & msg, MsgType msgLevel = Runtime0, MsgType timerLevel = Timings0)
      : printMonitor_(object, msg + " (" + object.description() + ")", msgLevel),
        timerMonitor_(object, object.ShortClassName() + ": " + msg + " (total)",    timerLevel)
    { }

  private:
    //! Manages printing.
    PrintMonitor printMonitor_;
    //! Records total time spent  in this object and all its children, over all levels.
    TimeMonitor timerMonitor_;
  };

  //---------------------------------------------------------------------------------------------------

  /*! @class SubMonitor
    @brief Another version of Monitor that doesn't repeat the object description.

    Times object and all its children over all levels.
    A timer is created only if 'timerLevel' (Timings1 by default) is enabled.

    The timer yields output such as
     \verbatim
     MueLu: SaPFactory: Eigenvalue estimate (sub, total)        0.149 (6)        0.1615 (6)       0.1643 (6)       0.02692 (6)
     \endverbatim
     Note that the keyword \a sub denotes that this is output from a SubFactoryMonitor.
     Note that the keyword \a total denotes timing of the object and its children.
  */
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

  //---------------------------------------------------------------------------------------------------

  /*! @class FactoryMonitor
      @brief Similar to Monitor but with additional timers.

     This class provides the following three timers for an object:
     - Timer for this object over all levels, excluding calls to children.
     - Timer for this object per level, excluding calls to children.
     - Timer for this object over all levels, including calls to children.

     The three timers above yield output like the following:
     \verbatim
       MueLu: SaPFactory: Prolongator smoothing                     2.599 (6)        2.602 (6)        2.607 (6)        0.4336 (6)
       MueLu: SaPFactory: Prolongator smoothing (level=5)           0.4731 (1)       0.4735 (1)       0.4736 (1)       0.4735 (1)
       MueLu: SaPFactory: Prolongator smoothing (total)             6.508 (6)        6.551 (6)        6.557 (6)        1.092 (6)
     \endverbatim

     Note that the keyword \a total denotes timing of the object and its children.
  */
  class FactoryMonitor: public Monitor {
  public:

    //! Constructor
    FactoryMonitor(const BaseClass& object, const std::string & msg, int levelID, MsgType msgLevel = Runtime0, MsgType timerLevel = Timings0)
      : Monitor(object, msg, msgLevel, timerLevel),
        timerMonitorExclusive_(object, object.ShortClassName() + ": " + msg, timerLevel)
    {
      if (IsPrint(TimingsByLevel)) {
        levelTimeMonitor_ = rcp(new TimeMonitor(object, object.ShortClassName() + ": " + msg + " (total, level=" + Teuchos::Utils::toString(levelID) + ")", timerLevel));
        levelTimeMonitorExclusive_ = rcp(new MutuallyExclusiveTimeMonitor<Level>(object, object.ShortClassName() + ": " + msg + " (level=" + Teuchos::Utils::toString(levelID) + ")", timerLevel));
      }
    }

    /*! Constructor

      TODO: code factorization
    */
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
    //! Total time spent on this level in this object and all its children.
    RCP<TimeMonitor>                           levelTimeMonitor_;
    //! Total time spent on all levels in this object only, excluding all children.
    MutuallyExclusiveTimeMonitor<FactoryBase>  timerMonitorExclusive_;
    //! Total time spent on this level in this object only, excluding all children.
    RCP<MutuallyExclusiveTimeMonitor<Level> >  levelTimeMonitorExclusive_;
  };

  //---------------------------------------------------------------------------------------------------

  /*! @class SubFactoryMonitor
    @brief Similar to SubMonitor but adds a timer level by level.

    Times an object and all its children on a level-by-level basis.

    The timer yields output such as
     \verbatim
       MueLu: SaPFactory: Eigenvalue estimate (sub, total, level=4)    0.01999 (1)      0.02103 (1)      0.02121 (1)      0.02103 (1)
     \endverbatim
     Note that the keyword \a sub denotes that this is output from a SubFactoryMonitor.
     Note that the keyword \a total denotes timing of the object and its children.
  */
  class SubFactoryMonitor: public SubMonitor {
  public:

    //!Constructor
    SubFactoryMonitor(const BaseClass& object, const std::string & msg, int levelID, MsgType msgLevel = Runtime1, MsgType timerLevel = Timings1)
      : SubMonitor(object, msg, msgLevel, timerLevel)
    {
      if (IsPrint(TimingsByLevel)) {
        levelTimeMonitor_ = rcp(new TimeMonitor(object, object.ShortClassName() + ": " + msg + " (sub, total, level=" + Teuchos::Utils::toString(levelID) + ")", timerLevel));
      }
    }

    //!Constructor
    SubFactoryMonitor(const BaseClass& object, const std::string & msg, const Level & level, MsgType msgLevel = Runtime1, MsgType timerLevel = Timings1)
      : SubMonitor(object, msg, msgLevel, timerLevel)
    {
      if (IsPrint(TimingsByLevel)) {
        levelTimeMonitor_ = rcp(new TimeMonitor(object, object.ShortClassName() + ": " +  msg + " (sub, total, level=" + Teuchos::Utils::toString(level.GetLevelID()) + ")", timerLevel));
      }
    }
  private:
    //! Total time spent on this level in this object and all children.
    RCP<TimeMonitor> levelTimeMonitor_;
  };

} // namespace MueLu

#endif // MUELU_MONITOR_HPP
