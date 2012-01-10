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


#include <fei_Logger.hpp>
#include <fei_LogManager.hpp>
#include <fei_LogFile.hpp>

fei::Logger::Logger()
 : output_level_(NONE),
   output_stream_(0),
   logIDs_(),
   logEqns_()
{
  fei::LogFile& log_file = fei::LogFile::getLogFile();
  output_stream_ = log_file.getOutputStream();
}

fei::Logger::~Logger()
{
}

void fei::Logger::setOutputLevel(OutputLevel olevel)
{
  output_level_ = olevel;
  fei::LogFile& log_file = fei::LogFile::getLogFile();
  output_stream_ = log_file.getOutputStream();
}

void fei::Logger::addLogID(int ID)
{
  logIDs_.insert(ID);
}

void fei::Logger::addLogEqn(int eqn)
{
  logEqns_.insert(eqn);
}

bool fei::Logger::isLogID(int ID)
{
  return(logIDs_.find(ID) != logIDs_.end());
}

bool fei::Logger::isLogEqn(int eqn)
{
  return(logEqns_.find(eqn) != logEqns_.end());
}

std::set<int>& fei::Logger::getLogIDs()
{
  return(logIDs_);
}

std::set<int>& fei::Logger::getLogEqns()
{
  return(logEqns_);
}

