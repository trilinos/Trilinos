#include <stk_util/environment/EnvData.hpp>
#include <stk_util/environment/OutputLog.hpp>
#include <stk_util/environment/ProgramOptions.hpp>

#include <stk_util/util/IndentStreambuf.hpp>

#include <mpi.h>
#include <ostream>

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
      m_isZapotec(false),
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
    m_execMap[sierra::Env::EXEC_TYPE_FLUID].m_master    = -1;
    m_execMap[sierra::Env::EXEC_TYPE_FLUID].m_groupComm = MPI_COMM_NULL;
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

  MPI_Comm parallel_comm();
  int EnvData::parallel_size();
  int EnvData::parallel_rank();
}

