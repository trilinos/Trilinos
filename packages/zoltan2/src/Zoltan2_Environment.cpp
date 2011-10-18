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

#include <Zoltan2_Environment.hpp>

namespace Zoltan2 {

Environment::Environment( ParameterList &problemParams,
  const RCP<const Comm<int> > &comm):
  params_(problemParams), validParams_(), 
  myRank_(0), numProcs_(1),
  printDebugMessages_(false), printProfilingMessages_(false),
  errorCheckLevel_(), debugDepthLevel_(), profilingIndicator_(),
  committed_(false), comm_(comm), dbg_()
{
  myRank_ = comm->getRank();
  numProcs_ = comm->getSize();
}

Environment::Environment():
  params_(), validParams_(), 
  myRank_(0), numProcs_(1),
  printDebugMessages_(false), printProfilingMessages_(false),
  errorCheckLevel_(), debugDepthLevel_(), profilingIndicator_(),
  committed_(false), dbg_()
{
  comm_ = DefaultComm<int>::getComm();
  myRank_ = comm_->getRank();
  numProcs_ = comm_->getSize();
}

Environment::~Environment()
{
}

Environment::Environment(const Environment &env)
{
  params_ = env.params_;
  validParams_ = env.validParams_;
  printDebugMessages_ = env.printDebugMessages_;
  printProfilingMessages_ = env.printProfilingMessages_;
  errorCheckLevel_ = env.errorCheckLevel_;
  debugDepthLevel_ = env.debugDepthLevel_;
  profilingIndicator_ = env.profilingIndicator_;
  comm_ = env.comm_;
  dbg_ = env.dbg_;
}

Environment &Environment::operator=(const Environment &env)
{
  if (this == &env) return *this;
  this->params_ = env.params_;
  this->validParams_ = env.validParams_;
  this->printDebugMessages_ = env.printDebugMessages_;
  this->printProfilingMessages_ = env.printProfilingMessages_;
  this->errorCheckLevel_ = env.errorCheckLevel_;
  this->debugDepthLevel_ = env.debugDepthLevel_;
  this->profilingIndicator_ = env.profilingIndicator_;
  this->dbg_ = env.dbg_;
  this->comm_ = env.comm_;
  return *this;
}

void Environment::setCommunicator(const RCP<const Comm<int> > &comm)
{
  Z2_LOCAL_INPUT_ASSERTION(*comm, *this, 
    "parameters are already committed",
    !committed_, BASIC_ASSERTION);

  comm_ = comm;
  myRank_ = comm_->getRank();
  numProcs_ = comm_->getSize();
}

void Environment::setParameters(ParameterList &params)
{
  Z2_LOCAL_INPUT_ASSERTION(*comm_, *this, 
    "parameters are already committed",
    !committed_, BASIC_ASSERTION);

  params_ = params;
}

void Environment::addParameters(ParameterList &params)
{
  Z2_LOCAL_INPUT_ASSERTION(*comm_, *this, 
    "parameters are already committed",
    !committed_, BASIC_ASSERTION);

  params_.setParameters(params);
}

void Environment::commitParameters()
{
  Z2_LOCAL_INPUT_ASSERTION(*comm_, *this, 
    "Can not commit parameters because setCommunicator has not been called.",
    numProcs_ > 0, BASIC_ASSERTION);

  createValidParameterList(validParams_);
  params_.validateParametersAndSetDefaults(validParams_);

  Array<int> reporters = 
    params_.get<Array<int> >(std::string("debug_procs"));

  printDebugMessages_ = IsInRangeList(myRank_, reporters);

  reporters = 
    params_.get<Array<int> >(std::string("profiling_procs")); 

  printProfilingMessages_ = IsInRangeList(myRank_, reporters);

  errorCheckLevel_ = params_.get<int>(std::string("error_check_level"));
  debugDepthLevel_ = params_.get<int>(std::string("debug_level"));
  profilingIndicator_ = params_.get<int>(std::string("profiling_level"));

  std::string &fname = params_.get<std::string>(std::string("debug_output_file"));
  std::ofstream dbgFile;
  if (fname.size() > 0){
    try{
      dbgFile.open(fname.c_str(), std::ios::out|std::ios::trunc);
    }
    catch(std::exception &e){
      // TODO
    }
    // TODO DebugManager destructor should close dbgFile if necessary
    dbg_ = rcp(new DebugManager(
      myRank_, printDebugMessages_,  dbgFile, debugDepthLevel_));
  }
  else{
    std::string &osname = params_.get<std::string>(std::string("debug_output_stream"));
    if (osname == std::string("std::cout"))
      dbg_ = rcp(new DebugManager(
        myRank_, printDebugMessages_,  std::cout, debugDepthLevel_));
    else if (osname == std::string("std::cerr"))
      dbg_ = rcp(new DebugManager(
        myRank_, printDebugMessages_,  std::cerr, debugDepthLevel_));
  }

  // TODO - do the same thing for profiling message data.

  committed_ = true;
}
  
}  //namespace Zoltan2

#endif
