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


#include "fei_LogFile.hpp"
#include "fei_iostream.hpp"
#include "fei_fstream.hpp"
#include "fei_sstream.hpp"
#include <string>

fei::LogFile::LogFile()
 : output_stream_(0),
   counter_(0)
{
}

fei::LogFile::~LogFile()
{
  counter_ = 0;
  closeOutputStream();
}

void fei::LogFile::openOutputStream(const char* path,
                                    int nprocs,
                                    int localproc)
{
  closeOutputStream();

  std::string pathstr("./");
  if (path != NULL) {
    pathstr = path;
  }

  if (pathstr[pathstr.size()] != '/') {
    pathstr = pathstr+"/";
  }

  FEI_OSTRINGSTREAM osstr;
  osstr << pathstr << "fei_log."<<counter_<<"."<<nprocs<<"."<<localproc;
  std::string filename = osstr.str();

  ++counter_;

  output_stream_ = new FEI_OFSTREAM(filename.c_str(), IOS_OUT);

  if (output_stream_ == NULL || output_stream_->bad()) {
    fei::console_out() << "couldn't open debug output file: " << filename << FEI_ENDL;
    delete output_stream_;
    output_stream_ = 0;
  }
}

FEI_OSTREAM* fei::LogFile::getOutputStream()
{
  return( output_stream_ );
}

void fei::LogFile::closeOutputStream()
{
  delete output_stream_;
  output_stream_ = 0;
}

fei::LogFile& fei::LogFile::getLogFile()
{
  static fei::LogFile log_file;
  return(log_file);
}

