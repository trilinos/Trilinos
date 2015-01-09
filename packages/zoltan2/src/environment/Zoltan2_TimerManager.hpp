// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
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
