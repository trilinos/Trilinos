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

#include <Zoltan2_Environment.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_RCP.hpp>

#include <sstream>
#include <ostream>

namespace Zoltan2 {

template <typename metric_t>
  static void makeMetricOutputManager(int rank, bool iPrint, std::string fname,
    std::string osname, Teuchos::RCP<MetricOutputManager<metric_t> > &mgr);

static void makeDebugManager(int rank, bool iPrint,
  int level, std::string fname, std::string osname,
  Teuchos::RCP<DebugManager> &mgr);

#define UNSET std::string("notSet")

Environment::Environment( Teuchos::ParameterList &problemParams,
  const Teuchos::RCP<const Teuchos::Comm<int> > &comm):
  myRank_(0), numProcs_(1), comm_(comm),
  debugOut_(), timerOut_(), memoryOut_(),
  errorCheckLevel_(BASIC_ASSERTION),
  params_(problemParams), committed_(false)

{
  myRank_ = comm->getRank();
  numProcs_ = comm->getSize();
}

Environment::Environment():
  myRank_(0), numProcs_(1), comm_(),
  debugOut_(), timerOut_(), memoryOut_(),
  errorCheckLevel_(BASIC_ASSERTION),
  params_(), committed_(false)
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

void Environment::convertStringToInt(Teuchos::ParameterList &params)
{
  using Teuchos::ParameterList;
  using Teuchos::ParameterEntry;
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using std::string;
  ParameterList::ConstIterator next = params.begin();

  // Data type of these parameters will now change from string to int

  string validatorName("StringIntegralValidator(int)");
  typedef Teuchos::StringToIntegralParameterEntryValidator<int> s2i_t;

  while (next != params.end()){

    const std::string &name = next->first;
    ParameterEntry &entry = params.getEntry(name);

    if (entry.isList()){
      ParameterList &pl = entry.getValue<ParameterList>(&pl);
      convertStringToInt(pl);
    }
    else{
      if ((entry.validator()).get()){
        if (entry.validator()->getXMLTypeName() == validatorName){
          string &entryValue = entry.getValue<string>(&entryValue);
          RCP<const s2i_t> s2i =
            Teuchos::rcp_dynamic_cast<const s2i_t>(entry.validator(), true);
          int val = s2i->getIntegralValue(entryValue);
          entry.setValue<int>(val);
        }
      }
    }
    ++next;
  }
}

void Environment::commitParameters()
{
  using std::string;
  using Teuchos::Array;
  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  Array<int> nodeZeroOnly(2);
  nodeZeroOnly[0] = 0;
  nodeZeroOnly[1] = RANGE_IS_LISTED;

  ParameterList validParams;

  createValidatorList(params_, validParams);

  // Note that since validParams only contains parameters that
  // appear in params_, there are no defaults being set.  We
  // call ParameterList::validateParametersAndSetDefaults() instead of
  // ParameterList::validateParameters() because we want the
  // the validators' validateAndModify() to be called instead
  // of validate().  validateAndModify() "fixes" some of the
  // parameters for us.

  params_.validateParametersAndSetDefaults(validParams);

  // For all of the true/false, yes/no, 0/1 string parameters, 
  // convert them to integers zero or one.  I would have
  // expected validateAndModify() to do this.

  convertStringToInt(params_);

  /////////////////////////////////////////////////////////////////////

  int &level = params_.get<int>("debug_level", BASIC_STATUS);

  if (level > NO_STATUS){
    Array<int> &reporters = 
      params_.get<Array<int> >("debug_procs", nodeZeroOnly);
    bool iPrint = IsInRangeList(myRank_, reporters);
    string &fname = params_.get<string>("debug_output_file", UNSET);
    string &osname = params_.get<string>("debug_output_stream", UNSET);

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

  // TODO instead of level, either do timing or don't do timing
  level = params_.get<int>("timing_level", BASIC_STATUS);

  if (level > NO_STATUS){
    Array<int> &reporters = 
      params_.get<Array<int> >("timing_procs", nodeZeroOnly);
    bool iPrint = IsInRangeList(myRank_, reporters);
    string &fname = params_.get<string>("timing_output_file", UNSET);
    string &osname = params_.get<string>("timing_output_stream", UNSET);

    try{
      makeMetricOutputManager<double>(myRank_, iPrint, fname, osname, timerOut_);
    }
    catch (std::exception &e){
      std::ostringstream oss;
      oss << myRank_ << ": unable to create timing output manager";
      oss << " (" << e.what() << ")";
      throw std::runtime_error(oss.str());
    }
  }
  else{
    timerOut_ = 
      rcp(new MetricOutputManager<double>(myRank_, false, std::cout, false));
  }

  // TODO instead of level, either output memory used info or don't
  level = params_.get<int>("memory_profiling_level", BASIC_STATUS);

  if (level > NO_STATUS){
    Array<int> &reporters = 
      params_.get<Array<int> >("memory_profiling_procs", nodeZeroOnly);
    bool iPrint = IsInRangeList(myRank_, reporters);
    string &fname = 
      params_.get<string>("memory_profiling_output_file", UNSET);
    string &osname = 
      params_.get<string>("memory_profiling_output_stream", UNSET);

    try{
      makeMetricOutputManager<long>(myRank_, iPrint, fname, osname, memoryOut_);
    }
    catch (std::exception &e){
      std::ostringstream oss;
      oss << myRank_ << ": unable to create memory profiling output manager";
      oss << " (" << e.what() << ")";
      throw std::runtime_error(oss.str());
    }
  }
  else{
    memoryOut_ = 
      rcp(new MetricOutputManager<long>(myRank_, false, std::cout, false));
  }

  errorCheckLevel_ = static_cast<AssertionLevel>( 
    params_.get<int>("error_check_level", BASIC_ASSERTION));

  committed_ = true;
}
  
bool Environment::hasSublist(const Teuchos::ParameterList &pl, 
                             const std::string &listName) const
{
  const Teuchos::ParameterEntry *entry = pl.getEntryPtr(listName);
  if (entry && entry->isList())
    return true;

  return false;
}

template<typename metric_t>
  static void makeMetricOutputManager(int rank, bool iPrint, std::string fname, 
    std::string osname, Teuchos::RCP<MetricOutputManager<metric_t> > &mgr)
{
  std::ofstream oFile;
  if (fname != UNSET){
    try{
      oFile.open(fname.c_str(), std::ios::out|std::ios::trunc);
    }
    catch(std::exception &e){
      throw std::runtime_error(e.what());
    }
  }

  typedef MetricOutputManager<metric_t> manager_t;

  if (osname == std::string("std::cout"))
    mgr = Teuchos::rcp(new manager_t(rank, iPrint, std::cout, true));
  else if (osname == std::string("std::cerr"))
    mgr = Teuchos::rcp(new manager_t(rank, iPrint, std::cerr, true));
  else if (osname != std::string("/dev/null"))
    mgr = Teuchos::rcp(new manager_t(rank, iPrint, oFile, true));
  else
    mgr = Teuchos::rcp(new manager_t(rank, false, std::cout, true));
}

static void makeDebugManager(int rank, bool iPrint,
  int level, std::string fname, std::string osname,
  Teuchos::RCP<DebugManager> &mgr)
{
  std::ofstream dbgFile;
  if (fname != UNSET){
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
