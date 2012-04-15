/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/


#ifndef _fei_Logger_hpp_
#define _fei_Logger_hpp_

#include <fei_fwd.hpp>
#include <fei_iosfwd.hpp>
#include <set>

namespace fei {
/** Class to be inherited by fei classes that wish to write
    to the fei debug-log file. */
class Logger {
 public:
  /** constructor */
  Logger();
  /** destructor */
  virtual ~Logger();

  /** set specified output-level. */
  void setOutputLevel(OutputLevel olevel);

  void addLogID(int ID);
  void addLogEqn(int eqn);

  bool isLogID(int ID);
  bool isLogEqn(int eqn);

  std::set<int>& getLogIDs();
  std::set<int>& getLogEqns();

 protected:
  /** output level
    Note that the OutputLevel enum is defined in fei_fwd.hpp.
  */
  OutputLevel output_level_;
  /** output stream */
  FEI_OSTREAM* output_stream_;

  std::set<int> logIDs_;
  std::set<int> logEqns_;

 private:
  Logger(const Logger& src);
  Logger& operator=(const Logger& src);
};//class Logger
}//namespace fei
#endif

