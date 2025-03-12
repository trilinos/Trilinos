// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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

#include "stk_util/environment/EnvData.hpp"
#include "stk_util/environment/OutputLog.hpp"       // for register_ostream, unregister_ostream
#include "stk_util/environment/ProgramOptions.hpp"  // for get_parsed_options
#include "stk_util/parallel/Parallel.hpp"           // for MPI_COMM_NULL, MPI_Comm, ompi_communi...
#include "stk_util/util/IndentStreambuf.hpp"        // for indent_streambuf
#include "stk_util/util/Null_Streambuf.hpp"         // for null_streambuf
#include <time.h>                                   // for time
#include <iostream>                                 // for cout, cerr


namespace stk {

  EnvData &EnvData::instance() {
    static EnvData s_env;
    
    return s_env;
  }

  EnvData::EnvData()
    : m_productName("not specified"),
      m_parsedOptions(stk::get_parsed_options()),
      m_nullBuf(),
      m_outputNull(&m_nullBuf),
      m_outputP0(&std::cout),
      m_output(),
      m_startTime((double) ::time(nullptr)),
      m_executablePath(),
      m_shutdownRequested(false),
      m_inputFileRequired(true),
      m_checkSubCycle(false),
      m_checkSmRegion(false),
      m_worldComm(MPI_COMM_NULL),
      m_parallelComm(MPI_COMM_NULL),
      m_parallelSize(-1),
      m_parallelRank(-1),
      m_emptyString(),
      m_onString("on"),
      m_inputFile(""),
      m_outPath("")
  {
    m_execMap[sierra::Env::EXEC_TYPE_LAG].m_rootProcessor      = -1;
    m_execMap[sierra::Env::EXEC_TYPE_LAG].m_groupComm   = MPI_COMM_NULL;
    m_execMap[sierra::Env::EXEC_TYPE_LAG].m_interComm   = MPI_COMM_NULL;
    m_execMap[sierra::Env::EXEC_TYPE_FLUID].m_rootProcessor    = -1;
    m_execMap[sierra::Env::EXEC_TYPE_FLUID].m_groupComm = MPI_COMM_NULL;
    m_execMap[sierra::Env::EXEC_TYPE_FLUID].m_interComm = MPI_COMM_NULL;
    stk::register_log_ostream(std::cout, "cout");
    stk::register_log_ostream(std::cerr, "cerr");

    stk::register_ostream(sierra::out(), "out");
    stk::register_ostream(sierra::pout(), "pout");
    stk::register_ostream(sierra::dout(), "dout");
    stk::register_ostream(sierra::tout(), "tout");

    static_cast<stk::indent_streambuf *>(sierra::dwout().rdbuf())->redirect(sierra::dout().rdbuf());
    stk::bind_output_streams("dout>null");
  }

  void EnvData::initialize(MPI_Comm worldComm)
  {
    m_worldComm = worldComm;
    m_parallelComm = worldComm;
    m_parallelSize = stk::parallel_machine_size(m_parallelComm);
    m_parallelRank = stk::parallel_machine_rank(m_parallelComm);
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

