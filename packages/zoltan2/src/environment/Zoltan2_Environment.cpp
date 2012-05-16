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
#include <Zoltan2_Util.hpp>

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_RCP.hpp>

#include <sstream>
#include <ostream>

namespace Zoltan2 {

//////////////////////////////////////////////////////////////////////
// Namespace definitions used by this class.

/*! \brief Create an output manager for debugging or status information
 *
 *  \param rank  the MPI rank of the calling process in the application
 *  \param iPrint   true if this process should output information
 *  \param fname    name of file to which output is to be appended, or
 *                      or Z2_UNSET
 *  \param osname   "std::cout", "std::cerr", "/dev/null", or Z2_UNSET
 *  \param mgr     on return, a pointer to the created output manager
 */

void makeDebugManager(int rank, bool iPrint,
  int level, std::string fname, std::string osname,
  Teuchos::RCP<DebugManager> &mgr)
{
  MessageOutputLevel lvl = static_cast<MessageOutputLevel>(level);

  bool haveFname = (fname != Z2_UNSET);
  bool haveStreamName = (!haveFname && (osname != Z2_UNSET));

  if (!haveFname && !haveStreamName){
    mgr = Teuchos::rcp(new DebugManager(rank, iPrint, std::cout, lvl));
    return;
  }

  if (haveFname){
    std::ofstream *dbgFile = new std::ofstream;
    if (iPrint){
      std::string newFname;
      addNumberToFileName(rank, fname, newFname);
      try{
        dbgFile->open(newFname.c_str(), std::ios::out|std::ios::trunc);
      }
      catch(std::exception &e){
        throw std::runtime_error(e.what());
      }
    }
    mgr = Teuchos::rcp(new DebugManager(rank, iPrint, *dbgFile, lvl));
    return;
  }

  if (osname == std::string("std::cout"))
    mgr = Teuchos::rcp(new DebugManager(rank, iPrint, std::cout, lvl));
  else if (osname == std::string("std::cerr"))
    mgr = Teuchos::rcp(new DebugManager(rank, iPrint, std::cerr, lvl));
  else if (osname == std::string("/dev/null"))
    mgr = Teuchos::rcp(new DebugManager(rank, false, std::cout, lvl));
  else   // should never happen
    throw std::logic_error("invalid debug output stream was not caught");
}

//////////////////////////////////////////////////////////////////////
// Environment definitions

Environment::Environment( Teuchos::ParameterList &problemParams,
  const Teuchos::RCP<const Teuchos::Comm<int> > &comm):
  myRank_(comm->getRank()), numProcs_(comm->getSize()), comm_(comm), 
  errorCheckLevel_(BASIC_ASSERTION),
  unvalidatedParams_(problemParams), params_(problemParams),
  debugOut_(), timerOut_(), timingOn(false), memoryOut_()
{
  commitParameters();
}

Environment::Environment():
  myRank_(0), numProcs_(1), comm_(), errorCheckLevel_(BASIC_ASSERTION),
  unvalidatedParams_("emptyList"), params_("emptyList"), 
  debugOut_(), timerOut_(), timingOn(false), memoryOut_()
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
  //   By default: if no output stream is specified, then node zero
  //       outputs BASIC_STATUS.

  int &level = params_.get<int>("debug_level", BASIC_STATUS);

  if (level > NO_STATUS){
    bool iPrint = (myRank_ == 0);   // default
    const Array<int> *reporters = 
      params_.getPtr<Array<int> >("debug_procs");
    if (reporters)
      iPrint = IsInRangeList(myRank_, *reporters);
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

  // Set up for memory usage output. 
  
  string &f2 = params_.get<string>("memory_profiling_output_file", Z2_UNSET);
  string &os2 = params_.get<string>("memory_profiling_output_stream", Z2_UNSET);


  if (f2 != Z2_UNSET || os2 != Z2_UNSET){
    long numKbytes = 0;
    if (myRank_ == 0)
      numKbytes = getProcessKilobytes();

    Teuchos::broadcast<int, long>(*comm_, 0, 1, &numKbytes);

    if (numKbytes == 0){
      // This is not a Linux system with proc/pid/statm.
      f2 = Z2_UNSET;
      os2 = Z2_UNSET;
      debugOut_->print(BASIC_ASSERTION, 
        "Warning: memory profiling requests but not available.");
    }
  }
  
  if (f2 != Z2_UNSET || os2 != Z2_UNSET){
    bool iPrint = (myRank_ == 0);   // default
    const Array<int> *reporters = 
      params_.getPtr<Array<int> >("memory_profiling_procs");
    if (reporters)
      iPrint = IsInRangeList(myRank_, *reporters);

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

}  //namespace Zoltan2

#endif
