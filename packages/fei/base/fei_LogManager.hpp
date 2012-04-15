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


#ifndef _fei_LogManager_hpp_
#define _fei_LogManager_hpp_

#include <fei_fwd.hpp>

#include <string>
#include <vector>

namespace fei {

/** Singleton class to manage attributes controlling the type and
   amount of data that should be written to the fei log file.
*/
class LogManager {
 public:
  /** destructor */
  virtual ~LogManager();

  /** Accessor for the one-and-only instance of LogManager.
      Constructs a LogManager instance on the first call, returns
      that same instance on the first and all subsequent calls.
  */
  static LogManager& getLogManager();

  /** Query output-level. Result is an enumeration. The enumeration is
   defined in fei_fwd.hpp. */
  OutputLevel getOutputLevel();

  /** Set output-level, using an enumeration. The enumeration is
   defined in fei_fwd.hpp. */
  void setOutputLevel(OutputLevel olevel);

  /** Set output-level, using a string. Valid values are strings that
   match the names of the enumeration values. e.g., "MATRIX_FILES", etc.
   */
  void setOutputLevel(const char* olevel);

  /** Specify path where debug-log files should be written. */
  void setOutputPath(const std::string& opath);

  /** Query for string specifying path to where debug-log files should
      be written. */
  const std::string& getOutputPath();

  /** Set numProcs and localProc (which will be used in the log-file-name).
  */
  void setNumProcs(int nprocs, int localproc);

 private:
  /** constructor */
  LogManager();

  OutputLevel output_level_;
  std::string output_path_;
  int numProcs_;
  int localProc_;
}; //class LogManager
}//namespace fei
#endif

