// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file Zoltan2_TimerManager.hpp
 *  \brief Declarations for TimerManager.
 */

#ifndef ZOLTAN2_TIMERMANAGER_HPP
#define ZOLTAN2_TIMERMANAGER_HPP

#include <Zoltan2_Standards.hpp>
#include <Zoltan2_Parameters.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <map>
#include <iostream>
#include <fstream>
#include <bitset>

namespace Zoltan2{

class TimerManager {

private:
  RCP<const Comm<int> > comm_;
  std::ostream *myOS_;
  std::ofstream *fileOS_;
  TimerType ttype_;
  std::bitset<NUM_TIMING_OPTIONS> typeSelector_;

  Array<RCP<Teuchos::Time> > timers_;
  std::map<std::string, int> timerMap_;
  int stopHint_;

public:

  /*! \brief Constructor for output to a file.
   */
  TimerManager(const RCP<const Comm<int> > & comm, std::ofstream *of,
    TimerType tt);

  /*! \brief Constructor for output to an std::ostream.
   */
  TimerManager(const RCP<const Comm<int> > & comm, std::ostream *os,
    TimerType tt);

  /*! \brief Destructor.
   */
  ~TimerManager();

  /*! \brief Start the named timer.
   */
  void start(TimerType tt, const std::string &name);

  /*! \brief Stop the named timer.
   */
  void stop(TimerType tt, const std::string &name);

  /*! \brief Print out global summary of timers and reset timers to zero.
   */
  void printAndResetToZero();

  /*! \brief Print out global summary, do not reset timers. 
   */
  void print() const;
};

}  // namespace Zoltan2

#endif
