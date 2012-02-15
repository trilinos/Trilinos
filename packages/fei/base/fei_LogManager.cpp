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


#include <fei_utils.hpp>
#include <fei_LogManager.hpp>
#include <fei_Logger.hpp>
#include <fei_LogFile.hpp>

fei::LogManager::LogManager()
 : output_level_(NONE),
   output_path_("./")
{
}

fei::LogManager::~LogManager()
{
}

fei::LogManager& fei::LogManager::getLogManager()
{
  static fei::LogManager log_manager;
  return(log_manager);
}

fei::OutputLevel fei::LogManager::getOutputLevel()
{
  return(output_level_);
}

void fei::LogManager::setOutputLevel(fei::OutputLevel olevel)
{
  if (output_level_ == olevel) {
    return;
  }

  bool no_existing_output_stream = output_level_ < fei::BRIEF_LOGS;

  output_level_ = olevel;

  bool need_output_stream = output_level_ >= fei::BRIEF_LOGS;

  if (need_output_stream && no_existing_output_stream) {
    fei::LogFile::getLogFile().openOutputStream(output_path_.c_str(),
                                                  numProcs_, localProc_);
  }
}

void fei::LogManager::setOutputLevel(const char* olevel)
{
  setOutputLevel(fei::utils::string_to_output_level(olevel));
}

void fei::LogManager::setOutputPath(const std::string& opath)
{
  output_path_ = opath;
}

const std::string& fei::LogManager::getOutputPath()
{
  return(output_path_);
}

void fei::LogManager::setNumProcs(int nprocs, int localproc)
{
  numProcs_ = nprocs;
  localProc_ = localproc;
}

