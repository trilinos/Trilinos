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


#ifndef _fei_LogFile_hpp_
#define _fei_LogFile_hpp_

#include "fei_iosfwd.hpp"

namespace fei {

/** Singleton class to manage (open, close, etc.) the one-and-only
   fei log file.
*/
class LogFile {
 public:
  /** destructor */
  virtual ~LogFile();

  /** Open a log-file ostream. If one is already open, it is closed
    before the new one is opened.

    The file name is 'fei_log.counter.nprocs.localproc', where
    counter is the number of times this function has been called,
    and nprocs and localproc are specified in the arguments.

    @param path Path, not including file-name, to log-file.
    @param nprocs Number of processors.
    @param localproc Rank of local processor.
  */
  void openOutputStream(const char* path=NULL,
                        int nprocs=1,
                        int localproc=0);

  /** Query for the log-file ostream. */
  FEI_OSTREAM* getOutputStream();

  /** Destroy the log-file ostream (closes the file).
  */
  void closeOutputStream();

  /** Accessor for the one-and-only instance of LogFile.
      Constructs a LogFile instance on the first call, returns
      that same instance on the first and all subsequent calls.
  */
  static LogFile& getLogFile();

 private:
  /** constructor */
  LogFile();

  FEI_OSTREAM* output_stream_;
  unsigned counter_;
}; //class LogFile
}//namespace fei
#endif
