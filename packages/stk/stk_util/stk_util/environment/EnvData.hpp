/*--------------------------------------------------------------------*/
/*    Copyright (c) 2013, Sandia Corporation.
/*    Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
/*    the U.S. Governement retains certain rights in this software.
/*    
/*    Redistribution and use in source and binary forms, with or without
/*    modification, are permitted provided that the following conditions are
/*    met:
/*    
/*        * Redistributions of source code must retain the above copyright
/*          notice, this list of conditions and the following disclaimer.
/*    
/*        * Redistributions in binary form must reproduce the above
/*          copyright notice, this list of conditions and the following
/*          disclaimer in the documentation and/or other materials provided
/*          with the distribution.
/*    
/*        * Neither the name of Sandia Corporation nor the names of its
/*          contributors may be used to endorse or promote products derived
/*          from this software without specific prior written permission.
/*    
/*    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
/*    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
/*    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
/*    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
/*    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
/*    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
/*    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
/*    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
/*    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
/*    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
/*    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
/*    
/*--------------------------------------------------------------------*/

#ifndef STK_UTIL_ENVIRONMENT_EnvData_h
#define STK_UTIL_ENVIRONMENT_EnvData_h

#include <stk_util/stk_config.h>
#if defined( STK_HAS_MPI)
#  include <mpi.h>
#endif
#include <iosfwd>
#include <map>
#include <string>
#include <boost/program_options.hpp>
#include <stk_util/util/Null_Streambuf.hpp>

namespace sierra {
  namespace Env {
    /**
     * @brief Enumeration ExecutableType defines the known types of coordinated executables that operate
     * with a sierra application.  Unfortunately, this scheme for coordination is currently defined by
     * Gemini whose implementation forces a limit of two executables, namely it and a fluid code.  The
     * startup_multi_exec() function handles the creation of groups which are contiguous processor
     * groups, each with lead processor being the least ranked processor in the group.
     *
     * Modification of the startup_multi_exec() function would need to be made to enable more than the
     * two executable types.
     */
    enum ExecType {
      EXEC_TYPE_WORLD = 0,            ///< Generic application using entire communicator (MPI_COMM_WORLD)
      EXEC_TYPE_FLUID = 1,            ///< Gemini Euler application
      EXEC_TYPE_LAG   = 2,            ///< Sierra Lagrangian application
      EXEC_TYPE_PEER  = 3             ///< Split communicator application; non-Gemini
    };

    struct ExecInfo
    {
      MPI_Comm              m_groupComm;
      int                   m_master;
      MPI_Comm              m_worldComm;
    };
  }
}

namespace stk {

struct EnvData
{
  typedef std::map<sierra::Env::ExecType, sierra::Env::ExecInfo>    ExecMap;

  static EnvData &instance();
  EnvData();

  ~EnvData();

  static void setInputFileName(std::string value);
  static std::string getInputFileName();
  static MPI_Comm parallel_comm();
  static int parallel_size();
  static int parallel_rank();

  std::string           m_productName;

  boost::program_options::variables_map & m_vm;

  null_streambuf	m_nullBuf;
  std::ostream		m_outputNull;
  std::ostream *	m_outputP0;
  std::ostringstream	m_output;

  double		m_startTime;
  std::string		m_executablePath;

  bool			m_shutdownRequested;
  bool                  m_inputFileRequired;
  bool                  m_checkSubCycle;
  bool                  m_isZapotec;

  MPI_Comm		m_worldComm;

  MPI_Comm		m_parallelComm;
  int			m_parallelSize;
  int			m_parallelRank;

  ExecMap               m_execMap;

  const std::string	m_emptyString;
  const std::string	m_onString;

  std::string           m_inputFile;
};
}
#endif
