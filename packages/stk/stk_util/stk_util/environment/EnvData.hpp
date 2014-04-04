/*--------------------------------------------------------------------*/
/*    Copyright 2014 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef STK_UTIL_ENVIRONMENT_EnvData_h
#define STK_UTIL_ENVIRONMENT_EnvData_h

#include <mpi.h>
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
