// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <stk_util/stk_config.h>
#include <stk_util/environment/EnvData.hpp>
#if defined( STK_HAS_MPI)
#  include <mpi.h>                        // for MPI_COMM_NULL, MPI_Comm, etc
#endif
#include <time.h>                       // for time, NULL
#include <iostream>                     // for cout, cerr
#include <stk_util/environment/OutputLog.hpp>  // for register_ostream, etc
#include <stk_util/environment/ProgramOptions.hpp>
#include <stk_util/util/IndentStreambuf.hpp>  // for indent_streambuf
#include "stk_util/util/Null_Streambuf.hpp"  // for null_streambuf


namespace stk {

  EnvData &EnvData::instance() {
    static EnvData s_env;
    
    return s_env;
  }

  EnvData::EnvData()
    : m_productName("not specified"),
      m_vm(stk::get_variables_map()),
      m_nullBuf(),
      m_outputNull(&m_nullBuf),
      m_outputP0(&std::cout),
      m_output(),
      m_startTime((double) ::time(NULL)),
      m_executablePath(),
      m_shutdownRequested(false),
      m_inputFileRequired(true),
      m_checkSubCycle(false),
      m_checkSmRegion(false),
      m_isZapotec(false),
      m_usingDiffingTool(false),
      m_worldComm(MPI_COMM_NULL),
      m_parallelComm(MPI_COMM_NULL),
      m_parallelSize(-1),
      m_parallelRank(-1),
      m_emptyString(),
      m_onString("on"),
      m_inputFile("")
  {
    m_execMap[sierra::Env::EXEC_TYPE_LAG].m_master      = -1;
    m_execMap[sierra::Env::EXEC_TYPE_LAG].m_groupComm   = MPI_COMM_NULL;
    m_execMap[sierra::Env::EXEC_TYPE_LAG].m_interComm   = MPI_COMM_NULL;
    m_execMap[sierra::Env::EXEC_TYPE_FLUID].m_master    = -1;
    m_execMap[sierra::Env::EXEC_TYPE_FLUID].m_groupComm = MPI_COMM_NULL;
    m_execMap[sierra::Env::EXEC_TYPE_FLUID].m_interComm = MPI_COMM_NULL;
    stk::register_log_ostream(std::cout, "cout");
    stk::register_log_ostream(std::cerr, "cerr");

    stk::register_ostream(sierra::out(), "out");
    stk::register_ostream(sierra::pout(), "pout");
    stk::register_ostream(sierra::dout(), "dout");
    stk::register_ostream(sierra::tout(), "tout");

    static_cast<stk::indent_streambuf *>(sierra::dwout().rdbuf())->redirect(sierra::dout().rdbuf());
  }

  EnvData::~EnvData()
  {
    static_cast<stk::indent_streambuf *>(sierra::dwout().rdbuf())->redirect(std::cout.rdbuf());

    stk::unregister_ostream(sierra::tout());
    stk::unregister_ostream(sierra::dout());
    stk::unregister_ostream(sierra::pout());
    stk::unregister_ostream(sierra::out());

    stk::unregister_log_ostream(std::cerr);
    stk::unregister_log_ostream(std::cout);
  }

  void EnvData::setInputFileName(std::string name) {
    instance().m_inputFile = name;
  }

  std::string EnvData::getInputFileName() {
    return instance().m_inputFile;
  }

  int EnvData::parallel_size() {
    return instance().m_parallelSize;
  }

  int EnvData::parallel_rank() {
    return instance().m_parallelRank;
  }

  MPI_Comm EnvData::parallel_comm()
  {
    return instance().m_parallelComm;
  }
}

