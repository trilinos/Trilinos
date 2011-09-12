// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER
/*! \file Zoltan2_Environment.cpp

    \brief The definition of the Environment object.
*/

#ifndef _ZOLTAN2_ENVIRONMENT_HPP_
#define _ZOLTAN2_ENVIRONMENT_HPP_

/*! \file Zoltan2_Environment.hpp
  
  \brief Defines the Zoltan2::Environment class.

*/

#include <ostream>
#include <Zoltan2_Environment.hpp>
#include <Zoltan2_Exceptions.hpp>
#include <Teuchos_RCP.hpp>

namespace Zoltan2 {

Environment::Environment( Teuchos::ParameterList &problemParams,
  Teuchos::RCP<Teuchos::Comm<int> > &comm):
  _params(problemParams), _validParams(), 
  _myRank(0), _numProcs(0),
  _printDebugMessages(false), _printProfilingMessages(false),
  _errorCheckLevel(), _debugDepthLevel(), _profilingIndicator(),
  _committed(false), _comm(comm)
{
  _myRank = comm->getRank();
  _numProcs = comm->getSize();
}

Environment::Environment():
  _params(), _validParams(), 
  _myRank(0), _numProcs(0),
  _printDebugMessages(false), _printProfilingMessages(false),
  _errorCheckLevel(), _debugDepthLevel(), _profilingIndicator(),
  _committed(false)
{
  _comm = Teuchos::DefaultComm<int>::getComm();
  _myRank = _comm->getRank();
  _numProcs = _comm->getSize();
}

Environment::~Environment()
{
}

Environment::Environment(const Environment &env)
{
  _params = env._params;
  _validParams = env._validParams;
  _printDebugMessages = env._printDebugMessages;
  _printProfilingMessages = env._printProfilingMessages;
  _errorCheckLevel = env._errorCheckLevel;
  _debugDepthLevel = env._debugDepthLevel;
  _profilingIndicator = env._profilingIndicator;
  _comm = env._comm;
  _dbg = env._dbg;
}

Environment &Environment::operator=(const Environment &env)
{
  if (this == &env) return *this;
  this->_params = env._params;
  this->_validParams = env._validParams;
  this->_printDebugMessages = env._printDebugMessages;
  this->_printProfilingMessages = env._printProfilingMessages;
  this->_errorCheckLevel = env._errorCheckLevel;
  this->_debugDepthLevel = env._debugDepthLevel;
  this->_profilingIndicator = env._profilingIndicator;
  this->_dbg = env._dbg;
  this->_comm = env._comm;
  return *this;
}

void Environment::setCommunicator(Teuchos::RCP<Teuchos::Comm<int> > &comm)
{
  Z2_LOCAL_INPUT_ASSERTION(*comm, *this, 
    "parameters are already committed",
    _committed, BASIC_ASSERTION);

  _comm = comm;
  _myRank = _comm->getRank();
  _numProcs = _comm->getSize();
}

void Environment::setParameters(Teuchos::ParameterList &params)
{
  Z2_LOCAL_INPUT_ASSERTION(*_comm, *this, 
    "parameters are already committed",
    _committed, BASIC_ASSERTION);

  _params = params;
}

void Environment::addParameters(Teuchos::ParameterList &params)
{
  Z2_LOCAL_INPUT_ASSERTION(*_comm, *this, 
    "parameters are already committed",
    _committed, BASIC_ASSERTION);

  _params.setParameters(params);
}

void Environment::commitParameters()
{
  if (_numProcs == 0){
    std::string msg("Can not commit parameters because ");
    msg += std::string("setCommunicator has not been called.");
    throw std::runtime_error(msg);
  }

  createValidParameterList(_validParams);
  _params.validateParametersAndSetDefaults(_validParams);

  Teuchos::Array<int> reporters = 
    _params.get<Teuchos::Array<int> >(std::string("debug_procs"));

  _printDebugMessages = IsInRangeList(_myRank, reporters);

  reporters = 
    _params.get<Teuchos::Array<int> >(std::string("profiling_procs")); 

  _printProfilingMessages = IsInRangeList(_myRank, reporters);

  _errorCheckLevel = _params.get<int>(std::string("error_check_level"));
  _debugDepthLevel = _params.get<int>(std::string("debug_level"));
  _profilingIndicator = _params.get<int>(std::string("timing_level"));

  std::string &fname = _params.get<std::string>(std::string("debug_output_file"));
  std::ofstream dbgFile;
  if (fname.size() > 0){
    try{
      dbgFile.open(fname.c_str(), std::ios::out|std::ios::trunc);
    }
    catch(std::exception &e){
      // TODO
    }
    _dbg = Teuchos::rcp(new DebugManager(
      _myRank, _printDebugMessages,  dbgFile, _debugDepthLevel));
  }
  else{
    std::string &osname = _params.get<std::string>(std::string("debug_output_stream"));
    if (osname == std::string("std::cout"))
      _dbg = Teuchos::rcp(new DebugManager(
        _myRank, _printDebugMessages,  std::cout, _debugDepthLevel));
    else if (osname == std::string("std::cerr"))
      _dbg = Teuchos::rcp(new DebugManager(
        _myRank, _printDebugMessages,  std::cerr, _debugDepthLevel));
  }

  // TODO - do the same thing for profiling message data.

  _committed = true;
}
  
}  //namespace Zoltan2

#endif
