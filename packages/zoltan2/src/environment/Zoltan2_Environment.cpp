// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
/*! \file Zoltan2_Environment.cpp
    \brief The definition of the Environment object.
*/

#ifndef _ZOLTAN2_ENVIRONMENT_CPP_
#define _ZOLTAN2_ENVIRONMENT_CPP_

#include <Zoltan2_Environment.hpp>
#include <Zoltan2_IntegerRangeList.hpp>
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
 *                      or Z2_UNSET_STRING
 *  \param ost     output stream type
 *  \param mgr     on return, a pointer to the created output manager
 */

void makeDebugManager(int rank, bool iPrint,
  int level, std::string fname, int ost,
  Teuchos::RCP<DebugManager> &mgr)
{
  MessageOutputLevel lvl = static_cast<MessageOutputLevel>(level);

  if (fname != Z2_UNSET_STRING){
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

  OSType os = static_cast<OSType>(ost);

  if (os == COUT_STREAM)
    mgr = Teuchos::rcp(new DebugManager(rank, iPrint, std::cout, lvl));
  else if (os == CERR_STREAM)
    mgr = Teuchos::rcp(new DebugManager(rank, iPrint, std::cerr, lvl));
  else if (os == NULL_STREAM)
    mgr = Teuchos::rcp(new DebugManager(rank, false, std::cout, lvl));
}

//////////////////////////////////////////////////////////////////////
// Environment definitions

Environment::Environment( const Teuchos::ParameterList &problemParams,
  const Teuchos::ParameterList &sourceParams,
  const Teuchos::RCP<const Teuchos::Comm<int> > &comm):
  myRank_(comm->getRank()), numProcs_(comm->getSize()), comm_(comm), 
  errorCheckLevel_(BASIC_ASSERTION),
  unvalidatedParams_(problemParams), params_(problemParams),
  sourceParams_(sourceParams),
  debugOut_(), timerOut_(), timingOn_(false), memoryOut_(), memoryOn_(false),
  memoryOutputFile_()
{
  try{
    commitParameters();
  }
  Z2_FORWARD_EXCEPTIONS
}

Environment::Environment():
  myRank_(0), numProcs_(1), comm_(), errorCheckLevel_(BASIC_ASSERTION),
  unvalidatedParams_("emptyList"), params_("emptyList"), 
  sourceParams_("emptyList"),
  debugOut_(), timerOut_(), timingOn_(false), memoryOut_(), memoryOn_(false),
  memoryOutputFile_()
{
  comm_ = Teuchos::DefaultComm<int>::getComm();
  myRank_ = comm_->getRank();
  numProcs_ = comm_->getSize();

  try{
    commitParameters();
  }
  Z2_FORWARD_EXCEPTIONS
}

Environment::~Environment()
{
  if (!memoryOutputFile_.is_null())
    memoryOutputFile_->close();
}

void Environment::getBaseParameters(ParameterList & pl)
{
  // these are parameters which are generic to all environments - timers, debugging, etc

  // we set the name here because this always happens
  pl.setName("zoltan2ValidatingParameters");

  // error_check_level
  RCP<Teuchos::StringToIntegralParameterEntryValidator<int> > error_check_level_Validator = Teuchos::rcp( new Teuchos::StringToIntegralParameterEntryValidator<int>(
    Teuchos::tuple<std::string>( "no_assertions", "basic_assertions", "complex_assertions", "debug_mode_assertions" ),
    Teuchos::tuple<std::string>( "no assertions will be performed", "typical checks of argument validity (fast, default)", "additional checks, i.e. is input graph a valid graph)", "check for everything including logic errors (slowest)" ),
    Teuchos::tuple<int>( 0, 1, 2, 3 ),
      "error_check_level") );
  pl.set("error_check_level", "basic_assertions", "  the amount of error checking performed     (If the compile flag Z2_OMIT_ALL_ERROR_CHECKING was set,     then error checking code is not executed at runtime.)");
  pl.getEntryRCP("error_check_level")->setValidator(error_check_level_Validator);

  // basic_status
  RCP<Teuchos::StringToIntegralParameterEntryValidator<int> > debug_level_Validator = Teuchos::rcp( new Teuchos::StringToIntegralParameterEntryValidator<int>(
    Teuchos::tuple<std::string>( "no_status", "basic_status", "detailed_status", "verbose_detailed_status" ),
    Teuchos::tuple<std::string>( "library outputs no status information", "library outputs basic status information (default)", "library outputs detailed information", "library outputs very detailed information" ),
    Teuchos::tuple<int>( 0, 1, 2, 3 ),
    "basic_status") );
  pl.set("debug_level", "basic_status", "  the amount of status/debugging/warning information to print");
  pl.getEntryRCP("debug_level")->setValidator(debug_level_Validator);

  // timer_type
  RCP<Teuchos::StringToIntegralParameterEntryValidator<int> > timer_type_Validator = Teuchos::rcp( new Teuchos::StringToIntegralParameterEntryValidator<int>(
    Teuchos::tuple<std::string>( "no_timers", "macro_timers", "micro_timers", "both_timers", "test_timers" ),
    Teuchos::tuple<std::string>( "No timing data will be collected (the default).", "Time an algorithm (or other entity) as a whole.", "Time the substeps of an entity.", "Run both MACRO and MICRO timers.", "Run timers added to code for testing, removed later" ),
    Teuchos::tuple<int>( 0, 1, 2, 3, 4 ),
      "no_timers") );
  pl.set("timer_type", "no_timers", "  the type of timing information to collect     (If the compile flag Z2_OMIT_ALL_PROFILING was set,     then the timing code is not executed at runtime.)");
  pl.getEntryRCP("timer_type")->setValidator(timer_type_Validator);

  // debug_output_stream
  RCP<Teuchos::StringToIntegralParameterEntryValidator<int> > debug_output_stream_Validator = Teuchos::rcp( new Teuchos::StringToIntegralParameterEntryValidator<int>(
    Teuchos::tuple<std::string>( "std::cout", "cout", "stdout", "std::cerr", "cerr", "stderr", "/dev/null", "null" ),
    Teuchos::tuple<std::string>( "", "", "", "", "", "", "", "" ), // original did not have this documented - left this here for building later
    Teuchos::tuple<int>( 0, 0, 0, 1, 1, 1, 2, 2 ),
    "cout") );
  pl.set("debug_output_stream", "cout", "  output stream for debug/status/warning messages (default cout)");
  pl.getEntryRCP("debug_output_stream")->setValidator(debug_output_stream_Validator);

  // timer_output_stream
  RCP<Teuchos::StringToIntegralParameterEntryValidator<int> > timer_output_stream_Validator = Teuchos::rcp( new Teuchos::StringToIntegralParameterEntryValidator<int>(
    Teuchos::tuple<std::string>( "std::cout", "cout", "stdout", "std::cerr", "cerr", "stderr", "/dev/null", "null" ),
    Teuchos::tuple<std::string>( "", "", "", "", "", "", "", "" ), // original did not have this documented - left this here for building later
    Teuchos::tuple<int>( 0, 0, 0, 1, 1, 1, 2, 2 ),
    "cout") );
  pl.set("timer_output_stream", "cout", "  output stream for timing report (default cout)");
  pl.getEntryRCP("timer_output_stream")->setValidator(timer_output_stream_Validator);

  // memory_output_stream
  RCP<Teuchos::StringToIntegralParameterEntryValidator<int> > memory_output_stream_Validator = Teuchos::rcp( new Teuchos::StringToIntegralParameterEntryValidator<int>(
    Teuchos::tuple<std::string>( "std::cout", "cout", "stdout", "std::cerr", "cerr", "stderr", "/dev/null", "null" ),
    Teuchos::tuple<std::string>( "", "", "", "", "", "", "", "" ), // original did not have this documented - left this here for building later
    Teuchos::tuple<int>( 0, 0, 0, 1, 1, 1, 2, 2 ),
    "cout") );
  pl.set("memory_output_stream", "cout", "  output stream for memory usage messages (default cout)");
  pl.getEntryRCP("memory_output_stream")->setValidator(memory_output_stream_Validator);

  // validator for file does not have to exist
  RCP<Teuchos::FileNameValidator> file__not_required_validator = Teuchos::rcp( new Teuchos::FileNameValidator(false) ); // file does not have to exist

  // debug_output_file
  pl.set("memory_output_file", "/dev/null", "  name of file to which memory profiling information should be written     (process rank will be included in file name)");
  pl.getEntryRCP("memory_output_file")->setValidator(file__not_required_validator);

  // timer_output_file
  pl.set("timer_output_file", "/dev/null", "  name of file to which timing information should be written     (process rank will be included in file name)");
  pl.getEntryRCP("timer_output_file")->setValidator(file__not_required_validator);

  // debug_output_file
  pl.set("debug_output_file", "/dev/null", "  name of file to which debug/status messages should be written     (process rank will be included in file name)");
  pl.getEntryRCP("debug_output_file")->setValidator(file__not_required_validator);

  // debug_procs
  pl.set("debug_procs", "0", "  list of ranks that output debugging/status messages (default \"0\")");
  RCP<Zoltan2::IntegerRangeListValidator<int>> debug_procs_Validator = Teuchos::rcp( new Zoltan2::IntegerRangeListValidator<int>(false) ); // unsorted false
  pl.getEntryRCP("debug_procs")->setValidator(debug_procs_Validator);

  // memory_procs
  RCP<Zoltan2::IntegerRangeListValidator<int>> memory_procs_Validator = Teuchos::rcp( new Zoltan2::IntegerRangeListValidator<int>(false) ); // unsorted false
  pl.set("memory_procs", "0", "  list of ranks that do memory profiling information (default \"0\")");
  pl.getEntryRCP("memory_procs")->setValidator(memory_procs_Validator);

  // random_seed
  RCP<Teuchos::AnyNumberParameterEntryValidator> random_seed_Validator = Teuchos::rcp( new Teuchos::AnyNumberParameterEntryValidator() );  // default is DOUBLE and accept Double, Int, String
  pl.set("random_seed", "0.5", "  random seed");
  pl.getEntryRCP("random_seed")->setValidator(random_seed_Validator);
}

void Environment::commitParameters()
{
  using Teuchos::Array;
  using Teuchos::ParameterList;

  // first get the base parameters - this is the core set that environment creates and always exists with any environment regardless of whether a Problem was involved or not
  ParameterList sourceParameters;
  getBaseParameters(sourceParameters); // these are the base set params which come from

  // Add in all the source parameters - these were generated by a Problem + Model + Algorithm - all the derived classes, etc
  // Note that this could be empty
  sourceParameters.setParameters(sourceParams_);

  bool emptyList = (params_.begin() == params_.end());

  if (!emptyList){

    ParameterList validParams;

    try{
      setValidatorsInList(params_, validParams, sourceParameters);
    }
    Z2_FORWARD_EXCEPTIONS
  
    // Note that since validParams only contains parameters that
    // appear in params_, there are no defaults being set.  We
    // call ParameterList::validateParametersAndSetDefaults() instead of
    // ParameterList::validateParameters() because we want the
    // the validators' validateAndModify() to be called instead
    // of validate().  validateAndModify() "fixes" some of the
    // parameters for us.
    // Note:  depth==0 --> do not validate sublists, 
    //                     since they are for TPL parameters
  
    params_.validateParametersAndSetDefaults(validParams, 0);

    // For all of the string to integer parameters, convert
    // them to the integer.  I would have
    // expected validateAndModify() to do this.
  
    convertStringToInt(params_);
  }

  /////////////////////////////////////////////////////////////////////

  // Set up for debugging/status output.
  //   By default: if no output stream is specified, then node zero
  //       outputs BASIC_STATUS to std::cout.

#ifndef Z2_OMIT_ALL_STATUS_MESSAGES
  int &level = params_.get<int>("debug_level", NUM_STATUS_OUTPUT_LEVELS);
  std::string &fname = params_.get<std::string>("debug_output_file", Z2_UNSET_STRING);
  int &os = params_.get<int>("debug_output_stream", NUM_OUTPUT_STREAMS);

  if (fname==Z2_UNSET_STRING && os==NUM_OUTPUT_STREAMS)
    os = COUT_STREAM;        // default output target
  if (level == NUM_STATUS_OUTPUT_LEVELS)
    level = BASIC_STATUS;    // default level of verbosity

  bool iPrint = (myRank_ == 0);   // default reporter

  const Array<int> *reporters = 
    params_.getPtr<Array<int> >("debug_procs");
  if (reporters)
    iPrint = IsInRangeList(myRank_, *reporters);

  try{
    makeDebugManager(myRank_, iPrint, level, fname, os, debugOut_);
  }
  catch (std::exception &e){
    std::ostringstream oss;
    oss << myRank_ << ": unable to create debug output manager";
    oss << " (" << e.what() << ")";
    throw std::runtime_error(oss.str());
  }
#endif

  // Set up for memory usage output. 
  
#ifndef Z2_OMIT_ALL_PROFILING
  std::string &f2 = 
    params_.get<std::string>("memory_output_file", Z2_UNSET_STRING);
  int &os2 = 
    params_.get<int>("memory_output_stream", NUM_OUTPUT_STREAMS);

  const Array<int> *reporters2 = 
    params_.getPtr<Array<int> >("memory_procs");

  bool doMemory = true;

  if (f2 != Z2_UNSET_STRING || os2 != NUM_OUTPUT_STREAMS || reporters2 != NULL){
    // user indicated they want memory usage information
    long numKbytes = 0;
    if (myRank_ == 0)
      numKbytes = getProcessKilobytes();

    Teuchos::broadcast<int, long>(*comm_, 0, 1, &numKbytes);

    if (numKbytes == 0){
      // This is not a Linux system with proc/pid/statm.
      f2 = Z2_UNSET_STRING;
      os2 = NUM_OUTPUT_STREAMS;
      reporters2 = NULL;
      this->debug(BASIC_STATUS, 
        std::string("Warning: memory profiling requested but not available."));
      doMemory = false;   // can't do it
    }
  }
  else{
    doMemory = false;   // not requested
  }

  if (doMemory){
    iPrint = (myRank_ == 0);   // default
    if (reporters2)
      iPrint = IsInRangeList(myRank_, *reporters2);

    try{
      makeMetricOutputManager<long>(myRank_, iPrint, f2, os2, memoryOut_,
        std::string("KB"), 10, memoryOutputFile_);
    }
    catch (std::exception &e){
      std::ostringstream oss;
      oss << myRank_ << ": unable to create memory profiling output manager";
      oss << " (" << e.what() << ")";
      throw std::runtime_error(oss.str());
    }

    memoryOn_ = true;
  }
#endif

#ifdef Z2_OMIT_ALL_ERROR_CHECKING
  errorCheckLevel_ = NO_ASSERTIONS;
#else
  errorCheckLevel_ = static_cast<AssertionLevel>( 
    params_.get<int>("error_check_level", BASIC_ASSERTION));
#endif
}
  
void Environment::convertStringToInt(Teuchos::ParameterList &params)
{
  using Teuchos::ParameterList;
  using Teuchos::ParameterEntry;
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  ParameterList::ConstIterator next = params.begin();

  // Data type of these parameters will now change from string to int

  std::string validatorName("StringIntegralValidator(int)");
  typedef Teuchos::StringToIntegralParameterEntryValidator<int> s2i_t;

  while (next != params.end()){

    const std::string &name = next->first;
    ParameterEntry &entry = params.getEntry(name);

    if (entry.isList()){
      ParameterList *dummy = NULL;
      ParameterList &pl = entry.getValue<ParameterList>(dummy);
      convertStringToInt(pl);
    }
    else{
      if ((entry.validator()).get()){
        if (entry.validator()->getXMLTypeName() == validatorName){
          std::string dummy("");
          std::string &entryValue = entry.getValue<std::string>(&dummy);
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
