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

#ifndef _ZOLTAN2_ENVIRONMENT_DEF_HPP_
#define _ZOLTAN2_ENVIRONMENT_DEF_HPP_

/*! \file Zoltan2_Environment_def.hpp
  
  \brief Defines the Zoltan2::Environment class.

*/

#include <ostream>
#include <Zoltan2_Util.hpp>
#include <Zoltan2_Environment.hpp>
#include <Teuchos_RCP.hpp>

namespace Zoltan2 {

Environment::Environment(
  Teuchos::ParameterList &problemParams, Teuchos::ParameterList &libraryConfig):
    _errorCheckLevel(0), _errorOStream(&std::cerr)
{
  setProblemParameters(problemParams);
  setLibraryConfiguration(libraryConfig);
}

Environment::Environment()
{
  Teuchos::ParameterList emptyList;
  Environment(emptyList, emptyList);
}

Environment::~Environment()
{
}

Environment::Environment(const Environment &env)
{
  _params = env._params;
  _config = env._config;
  _dbg.setOStream(env._dbg.getOStream());
  _dbg.setDebugLevel(env._dbg.getDebugLevel());
  _errorCheckLevel = env._errorCheckLevel;
  _errorOStream = env._errorOStream;
}

Environment &Environment::operator=(const Environment &env)
{
  if (this == &env) return *this;
  this->_params = env._params;
  this->_config = env._config;
  this->_dbg.setOStream(env._dbg.getOStream());
  this->_dbg.setDebugLevel(env._dbg.getDebugLevel());
  this->_errorCheckLevel = env._errorCheckLevel;
  this->_errorOStream = env._errorOStream;
  return *this;
}

void Environment::setProblemParameters(Teuchos::ParameterList &problemParams)
{
  _params = problemParams;

  const std::string error_check_level("ERROR_CHECK_LEVEL");

  int *checkLevel = _config.getPtr<int>(error_check_level);
  if (checkLevel)
    _errorCheckLevel = *checkLevel;

  if (_errorCheckLevel< 0)
    _errorCheckLevel = 0;
  else if (_errorCheckLevel > Z2_MAX_CHECK_LEVEL)
    _errorCheckLevel = Z2_MAX_CHECK_LEVEL;
}

void Environment::setLibraryConfiguration(Teuchos::ParameterList &libraryConfig)
{
  _config = libraryConfig;

#if 0
  //Zoltan2::getOutputStreamFromParameterList(libraryConfig, "ERROR_OSTREAM", _errorOStream, std::cerr);

  std::ostream *os;
  Zoltan2::getOutputStreamFromParameterList(libraryConfig, "DEBUG_OSTREAM", os, std::cout);
  _dbg.setOStream(os);

  const std::string debug_level("DEBUG_LEVEL");
  int *debugLevel = _config.getPtr<int>(debug_level);
  if (debugLevel)
    _dbg.setDebugLevel(*debugLevel);
#endif

#if 0
  // Same for timing when we add profiling abilities
  std::ostream *timingStream;
  Z2::getOutputStreamFromParameterList(_config, "TIMING_OSTREAM", timingStream);
  Teuchos::RCP<std::ostream> timingOut = Teuchos::rcp(timingStream);
  const std::string timing_level("TIMING_LEVEL");
  int *timingLevel = _config.getPtr<int>(timing_level);
  if (timingLevel)
    tmg.setTimingLevel(*timingLevel);
#endif
}
  
}  //namespace Zoltan2

#endif
