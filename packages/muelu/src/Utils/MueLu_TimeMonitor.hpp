// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
#ifdef HAVE_TEUCHOS_ADD_TIME_MONITOR_TO_STACKED_TIMER
#include "Teuchos_StackedTimer.hpp"
#include <sstream>
#endif

namespace MueLu {

/*! @class TimeMonitor

    @brief Integrates Teuchos::TimeMonitor with MueLu verbosity system.
*/
class TimeMonitor : public BaseClass {
 public:
  TimeMonitor(const BaseClass& object, const std::string& msg, MsgType timerLevel = Timings0);

  ~TimeMonitor();

 protected:
  TimeMonitor();

 private:
  RCP<Teuchos::Time> timer_;
  bool useStackedTimer_;
  std::string label_;
};  // class TimeMonitor

}  // namespace MueLu

#endif  // MUELU_TIMEMONITOR_HPP
