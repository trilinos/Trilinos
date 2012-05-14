// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_TimerManager.hpp
 *  \brief Declarations for TimerManager.
 */

#ifndef ZOLTAN2_TIMERMANAGER_HPP
#define ZOLTAN2_TIMERMANAGER_HPP

#include <Zoltan2_Standards.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <map>
#include <iostream>
#include <fstream>

namespace Zoltan2{

class TimerManager {

private:
  RCP<const Comm<int> > comm_;
  std::ostream *myOS_;
  std::ofstream *fileOS_;

  Array<RCP<Teuchos::Time> > timers;
  std::map<string, int> timerMap;
  int stopHint;

public:

  /*! \brief Constructor for output to a file.
   */
  TimerManager(const RCP<const Comm<int> > & comm, std::ofstream *of);

  /*! \brief Constructor for output to an std::ostream.
   */
  TimerManager(const RCP<const Comm<int> > & comm, std::ostream *os);

  /*! \brief Destructor.
   */
  ~TimerManager();

  /*! \brief Start the named timer.
   */
  void start(const string &name);

  /*! \brief Stop the named timer.
   */
  void stop(const string &name);

  /*! \brief Print out global summary of timers and reset timers to zero.
   */
  void printAndResetToZero();

  /*! \brief Print out global summary, do not reset timers. 
   */
  void print() const;
};

}  // namespace Zoltan2

#endif
