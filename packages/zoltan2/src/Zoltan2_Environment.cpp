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

#ifndef _ZOLTAN2_ENVIRONMENT_CPP_
#define _ZOLTAN2_ENVIRONMENT_CPP_

/*! \file Zoltan2_Environment.cpp
  
  \brief Defines the Zoltan2::Environment class.

*/

#include <Zoltan2_Environment.hpp>

#include <sstream>
#include <ostream>

namespace Zoltan2 {


Environment::Environment( Teuchos::ParameterList &problemParams,
  const Teuchos::RCP<const Teuchos::Comm<int> > &comm):
  myRank_(0), numProcs_(1), comm_(comm),
  debugOut_(), timerOut_(), memoryOut_(),
  params_(problemParams), validParams_(), committed_(false)
{
  myRank_ = comm->getRank();
  numProcs_ = comm->getSize();
}

Environment::Environment():
  myRank_(0), numProcs_(1), comm_(),
  debugOut_(), timerOut_(), memoryOut_(),
  params_(), validParams_(), committed_(false)
{
  comm_ = Teuchos::DefaultComm<int>::getComm();
  myRank_ = comm_->getRank();
  numProcs_ = comm_->getSize();
}

Environment::~Environment()
{
}

void Environment::setCommunicator(
  const Teuchos::RCP<const Teuchos::Comm<int> > &comm)
{
  if (committed_){
    std::ostringstream msg;
    msg<<myRank_<<" "<<__FILE__<<", "<<__LINE__; 
    msg<<", error: parameters are already committed";
    throw std::runtime_error(msg.str()); 
  }

  comm_ = comm;
  myRank_ = comm_->getRank();
  numProcs_ = comm_->getSize();
}

void Environment::setParameters(Teuchos::ParameterList &params)
{
  if (committed_){
    std::ostringstream msg;
    msg<<myRank_<<" "<<__FILE__<<", "<<__LINE__; 
    msg<<", error: parameters are already committed";
    throw std::runtime_error(msg.str()); 
  }

  params_ = params;
}

void Environment::addParameters(Teuchos::ParameterList &params)
{
  if (committed_){
    std::ostringstream msg;
    msg<<myRank_<<" "<<__FILE__<<", "<<__LINE__; 
    msg<<", error: parameters are already committed";
    throw std::runtime_error(msg.str()); 
  }

  params_.setParameters(params);
}

static void makeDebugManager(int rank, bool iPrint, int level,
  std::string fname, std::string osname, Teuchos::RCP<DebugManager> &mgr);

void Environment::commitParameters()
{
  using std::string;
  using Teuchos::Array;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  string noOutStream("/dev/null");
  string noFile();

  createValidParameterList(validParams_, *comm_);
  params_.validateParametersAndSetDefaults(validParams_);

  int level = params_.get<int>(string("debug_level"));

  if (level > NO_STATUS){
    Array<int> reporters = params_.get<Array<int> >(string("debug_procs"));
    bool iPrint = IsInRangeList(myRank_, reporters);
    string &fname = params_.get<string>(string("debug_output_file"));
    string &osname = params_.get<string>(string("debug_output_stream"));

    try{
      makeDebugManager(myRank_, iPrint, level, fname, osname, debugOut_);
    }
    catch (std::exception &e){
      std::ostringstream oss;
      oss << myRank_ << ": unable to create debug output manager";
      oss << " (" << e.what() << ")";
      throw std::runtime_error(oss.str());
    }
  }
  else{
    debugOut_ = rcp(new DebugManager(myRank_, false, std::cout, NO_STATUS));
  }

  level = params_.get<int>(string("timing_level"));

  if (level > NO_STATUS){
    Array<int> reporters = params_.get<Array<int> >(string("timing_procs"));
    bool iPrint = IsInRangeList(myRank_, reporters);
    string &fname = params_.get<string>(string("timing_output_file"));
    string &osname = params_.get<string>(string("timing_output_stream"));

    try{
      makeDebugManager(myRank_, iPrint, level, fname, osname, timerOut_);
    }
    catch (std::exception &e){
      std::ostringstream oss;
      oss << myRank_ << ": unable to create timing output manager";
      oss << " (" << e.what() << ")";
      throw std::runtime_error(oss.str());
    }
  }
  else{
    timerOut_ = rcp(new DebugManager(myRank_, false, std::cout, NO_STATUS));
  }

  level = params_.get<int>(string("memory_profiling_level"));

  if (level > NO_STATUS){
    Array<int> reporters = 
      params_.get<Array<int> >(string("memory_profiling_procs"));
    bool iPrint = IsInRangeList(myRank_, reporters);
    string &fname = params_.get<string>(string("memory_profiling_output_file"));
    string &osname = 
      params_.get<string>(string("memory_profiling_output_stream"));

    try{
      makeDebugManager(myRank_, iPrint, level, fname, osname, memoryOut_);
    }
    catch (std::exception &e){
      std::ostringstream oss;
      oss << myRank_ << ": unable to create memory profiling output manager";
      oss << " (" << e.what() << ")";
      throw std::runtime_error(oss.str());
    }
  }
  else{
    memoryOut_ = rcp(new DebugManager(myRank_, false, std::cout, NO_STATUS));
  }

  errorCheckLevel_ = static_cast<AssertionLevel>( 
    params_.get<int>(string("error_check_level")));
  
  hasAnyParams_ = true;    // due to defaults, we always have parameters
  hasPartitioningParams_ = true;
  hasOrderingParams_ = true;
  hasColoringParams_ = true;
  hasMatchingParams_ = true;

  try{
    params_.sublist(string("partitioning"));
  }
  catch(std::exception &e){
    hasPartitioningParams_ = false;
  }

  try{
    params_.sublist(string("ordering"));
  }
  catch(std::exception &e){
    hasOrderingParams_ = false;
  }

  try{
    params_.sublist(string("coloring"));
  }
  catch(std::exception &e){
    hasColoringParams_ = false;
  }

  try{
    params_.sublist(string("matching"));
  }
  catch(std::exception &e){
    hasMatchingParams_ = false;
  }

  committed_ = true;
}


static void makeDebugManager(int rank, bool iPrint,
  int level, std::string fname, std::string osname,
  Teuchos::RCP<DebugManager> &mgr)
{
  std::ofstream dbgFile;
  if (fname.size() > 0){
    try{
      dbgFile.open(fname.c_str(), std::ios::out|std::ios::trunc);
    }
    catch(std::exception &e){
      throw std::runtime_error(e.what());
    }
  }

  MessageOutputLevel lvl = static_cast<MessageOutputLevel>(level);

  if (osname == std::string("std::cout"))
    mgr = Teuchos::rcp(new DebugManager(rank, iPrint, std::cout, lvl));
  else if (osname == std::string("std::cerr"))
    mgr = Teuchos::rcp(new DebugManager(rank, iPrint, std::cerr, lvl));
  else if (osname != std::string("/dev/null"))
    mgr = Teuchos::rcp(new DebugManager(rank, iPrint, dbgFile, lvl));
  else
    mgr = Teuchos::rcp(new DebugManager(rank, false, std::cout, lvl));
}

// A helper function to return a default environment.

Teuchos::RCP<const Environment> getDefaultEnvironment()
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = 
    Teuchos::DefaultComm<int>::getComm();

  Teuchos::ParameterList params;   // an empty parameter list;

  Environment *defaultEnv = NULL;

  try{
    defaultEnv = new Environment(params, comm);
  }
  catch (std::exception &e){
    throw std::runtime_error("getDefaultEnvironment");
  }

  defaultEnv->commitParameters();   // will set parameters to defaults

  Teuchos::RCP<const Environment> env = Teuchos::rcp(defaultEnv);
  return env;
}
  
}  //namespace Zoltan2

#endif
