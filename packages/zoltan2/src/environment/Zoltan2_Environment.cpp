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

//////////////////////////////////////////////////////////////////////
// Prototypes for static namespace methods used by the Environment

/*! \brief A value to indicate a string parameter that was not set by the user.
 */
#define Z2_UNSET std::string("notSet")

/*! \brief Create an output manager for a metric value.
 *
 *  \param rank  the MPI rank of the calling process in the application
 *  \param iPrint   true if this process should print metric information
 *  \param fname    name of file to which output is to be appended, or
 *                      or Z2_UNSET
 *  \param osname   "std::cout", "std::cerr", "/dev/null", or Z2_UNSET
 *  \param mgr     on return, a pointer to the created output manager
 *
 * The template parameter is the data type of the entity being measured.
 */

template <typename metric_t>
  static void makeMetricOutputManager(int rank, bool iPrint, std::string fname,
    std::string osname, Teuchos::RCP<MetricOutputManager<metric_t> > &mgr);

/*! \brief Create an output manager for debugging or status information
 *
 *  \param rank  the MPI rank of the calling process in the application
 *  \param iPrint   true if this process should output information
 *  \param fname    name of file to which output is to be appended, or
 *                      or Z2_UNSET
 *  \param osname   "std::cout", "std::cerr", "/dev/null", or Z2_UNSET
 *  \param mgr     on return, a pointer to the created output manager
 */

static void makeDebugManager(int rank, bool iPrint,
  int level, std::string fname, std::string osname,
  Teuchos::RCP<DebugManager> &mgr);

/*! \brief Helper method to add process rank to a file name.
 *    \param rank   the rank of the calling process
 *    \param fname  the file name to modify
 *    \param newf   on return newf is fname with the rank added to the name.
 *
 *   If fname has no dot in it, then the rank is added to the end of the name.
 *   Otherwise the rank is added before the first dot in fname.
 */
static void addRankToFileName(int rank, std::string fname, std::string &newf);

//////////////////////////////////////////////////////////////////////
// Environment definitions

Environment::Environment( Teuchos::ParameterList &problemParams,
  const Teuchos::RCP<const Teuchos::Comm<int> > &comm):
  myRank_(0), numProcs_(1), comm_(comm), errorCheckLevel_(BASIC_ASSERTION),
  unvalidatedParams_(problemParams), params_(problemParams),
  debugOut_(), timerOut_(), memoryOut_()
{
  myRank_ = comm->getRank();
  numProcs_ = comm->getSize();
  commitParameters();
}

Environment::Environment():
  myRank_(0), numProcs_(1), comm_(), errorCheckLevel_(BASIC_ASSERTION),
  unvalidatedParams_(), params_(),
  debugOut_(), timerOut_(), memoryOut_()
{
  comm_ = Teuchos::DefaultComm<int>::getComm();
  myRank_ = comm_->getRank();
  numProcs_ = comm_->getSize();

  commitParameters();
}

Environment::~Environment()
{
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

  bool emptyList = (params_.begin() == params_.end());

  if (!emptyList){

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
  }

  /////////////////////////////////////////////////////////////////////

  // Set up for debugging/status output.

  int &level = params_.get<int>("debug_level", BASIC_STATUS);

  if (level > NO_STATUS){
    Array<int> &reporters = 
      params_.get<Array<int> >("debug_procs", nodeZeroOnly);
    bool iPrint = IsInRangeList(myRank_, reporters);
    string &fname = params_.get<string>("debug_output_file", Z2_UNSET);
    string &osname = params_.get<string>("debug_output_stream", Z2_UNSET);

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

  // Set up for timing output.

  string &f1 = params_.get<string>("timing_output_file", Z2_UNSET);
  string &os1 = params_.get<string>("timing_output_stream", Z2_UNSET);

  if (f1 != Z2_UNSET || os1 != Z2_UNSET){
    Array<int> &reporters = 
      params_.get<Array<int> >("timing_procs", nodeZeroOnly);
    bool iPrint = IsInRangeList(myRank_, reporters);

    try{
      makeMetricOutputManager<double>(myRank_, iPrint, f1, os1, timerOut_);
    }
    catch (std::exception &e){
      std::ostringstream oss;
      oss << myRank_ << ": unable to create timing output manager";
      oss << " (" << e.what() << ")";
      throw std::runtime_error(oss.str());
    }
  }
  else{
    // Says no processes output timing metrics.
    timerOut_ = 
      rcp(new MetricOutputManager<double>(myRank_, false, std::cout, false));
  }

  // Set up for memory usage output.
  
  string &f2 = params_.get<string>("memory_profiling_output_file", Z2_UNSET);
  string &os2 = params_.get<string>("memory_profiling_output_stream", Z2_UNSET);

  if (f2 != Z2_UNSET || os2 != Z2_UNSET){
    Array<int> &reporters = 
      params_.get<Array<int> >("memory_profiling_procs", nodeZeroOnly);
    bool iPrint = IsInRangeList(myRank_, reporters);

    try{
      makeMetricOutputManager<long>(myRank_, iPrint, f2, os2, memoryOut_);
    }
    catch (std::exception &e){
      std::ostringstream oss;
      oss << myRank_ << ": unable to create memory profiling output manager";
      oss << " (" << e.what() << ")";
      throw std::runtime_error(oss.str());
    }
  }
  else{
    // Says no processes output memory usage metrics.
    memoryOut_ = 
      rcp(new MetricOutputManager<long>(myRank_, false, std::cout, false));
  }

  errorCheckLevel_ = static_cast<AssertionLevel>( 
    params_.get<int>("error_check_level", BASIC_ASSERTION));
}
  
bool Environment::hasSublist(const Teuchos::ParameterList &pl, 
                             const std::string &listName) const
{
  const Teuchos::ParameterEntry *entry = pl.getEntryPtr(listName);
  if (entry && entry->isList())
    return true;

  return false;
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


//////////////////////////////////////////////////////////////////////
// Definitions static methods used by the Environment

static void addRankToFileName(int rank, std::string fname, std::string &newf)
{
  std::ostringstream id;
  id << "_";
  id.width(6);
  id.fill('0');
  id << rank;

  std::ostringstream localFileName;
  std::string::size_type loc = fname.find('.');

  if (loc == std::string::npos)
    localFileName << fname << id.str();
  else
    localFileName << fname.substr(0, loc) << id.str() << fname.substr(loc);

  newf = localFileName.str();
}

template<typename metric_t>
  static void makeMetricOutputManager(int rank, bool iPrint, std::string fname, 
    std::string osname, Teuchos::RCP<MetricOutputManager<metric_t> > &mgr)
{
  std::ofstream oFile;

  if (fname != Z2_UNSET){
    std::string newFname;
    addRankToFileName(rank, fname, newFname);

    try{
      oFile.open(newFname.c_str(), std::ios::out|std::ios::trunc);
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
  if (fname != Z2_UNSET){
    std::string newFname;
    addRankToFileName(rank, fname, newFname);
    try{
      dbgFile.open(newFname.c_str(), std::ios::out|std::ios::trunc);
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


}  //namespace Zoltan2

#endif
