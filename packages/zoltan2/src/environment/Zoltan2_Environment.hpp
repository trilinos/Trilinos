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

/*! \file Zoltan2_Environment.hpp
    \brief Defines the Environment class.
*/

#ifndef _ZOLTAN2_ENVIRONMENT_HPP_
#define _ZOLTAN2_ENVIRONMENT_HPP_

#include <Zoltan2_config.h>
#include <Zoltan2_Util.hpp>
#include <Zoltan2_IO.hpp>
#include <Zoltan2_Parameters.hpp>
#include <Zoltan2_DebugManager.hpp>
#include <Zoltan2_TimerManager.hpp>
#include <Zoltan2_MetricOutputManager.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_CommHelpers.hpp>

namespace Zoltan2 {

/*!  \brief The user parameters, debug, timing and memory profiling output objects, 
            and error checking methods.

  This is object is passed to almost every method in the library. It
  has the problem parameters and the configuration information that governs
  how the library should behave when reporting status information,
  testing for errors, and so on.

  The environment holds the application's default communicator.  Note that 
  this communicator may differ from the problem communicator for any
  given problem.
*/

class Environment{

public:

  typedef long   memory_t;     /*!< \brief data type for Kilobytes */

  typedef Teuchos::RCP<const Teuchos::Comm<int> > Comm_t;
  typedef Teuchos::RCP<DebugManager>     DebugManager_t;
  typedef Teuchos::RCP<MetricOutputManager<memory_t> > MemoryProfilerManager_t;
  typedef Teuchos::RCP<TimerManager> Timer_t;

  int  myRank_;   /*!< \brief mpi rank (relative to comm_) */

  int  numProcs_; /*!< \brief number of processes (relative to comm_) */

  Comm_t comm_;   /*!< \brief communicator for environment*/

  AssertionLevel errorCheckLevel_; /*!< \brief level of error checking to do */

  /*! \brief Constructor
   *
   *  \param problemParams  the parameters supplied by the user, and
   *                          not yet validated by the Environment
   *  \param comm           the default communicator for the application
   *
   *   Note that the communicator is for the application, not the problem.
   */
  Environment(Teuchos::ParameterList &problemParams,
    const Teuchos::RCP<const Teuchos::Comm<int> > &comm );

  /*! \brief Default Constructor
   *
   *    The default constructor uses the Teuchos default communicator,
   *    BASIC_STATUS for debug_level, and does not timing or memory profiling.
   *    It has error_check_level BASIC_ASSERTION. It has no other parameters.
   */
  Environment();

  /*! \brief Destructor
   */
  ~Environment();

  /*! \brief Provide the Timer object to the Environment.
   *
   *  Having a Timer implies that the user asked for for timing.
   *  The Problem creates and holds the Timer.
   */

  void setTimer(RCP<TimerManager> &timer) { timerOut_=timer; timingOn_=true;}

#ifdef Z2_OMIT_ALL_ERROR_CHECKING

  void localInputAssertion(const char *file, int lineNum,
    const char *msg, bool ok, AssertionLevel level) const {}

  void globalInputAssertion(const char *file, int lineNum,
    const char *msg, bool ok, AssertionLevel level, 
    const Comm_t &comm=comm_) const {}

  void localBugAssertion(const char *file, int lineNum,
    const char *msg, bool ok, AssertionLevel level) const {}

  void globalBugAssertion(const char *file, int lineNum,
    const char *msg, bool ok, AssertionLevel level, 
    const Comm_t &comm=comm_) const {}

  void localMemoryAssertion(const char *file, int lineNum,
    size_t nobj, bool ok) const {}

  void globalMemoryAssertion(const char *file, int lineNum,
    size_t nobj, bool ok, const Comm_t &comm=comm_) const {}

#else

  /*! \brief %Test for valid user input on local process only.
   *
   *   \param file     the __FILE__ value of the caller.
   *   \param lineNum  the __LINE__ value of the caller.
   *   \param msg      an optional descriptive message
   *   \param ok       a boolean which if false indicates an error
   *   \param level    a AssertionLevel value
   *
   *  If the \c level does not exceed the \c error_check_level parameter
   *  set by the user, then the assertion is tested and a std::runtime  error
   *  is thrown if it is false.
   */

  void localInputAssertion(const char *file, int lineNum,
    const char *msg, bool ok, AssertionLevel level) const {

    if (level <= errorCheckLevel_ && !ok){
      std::ostringstream emsg;
      emsg << myRank_<< ": " << file << "," << lineNum<< std::endl;
      if (msg)
        emsg << myRank_ << ": error: " << msg << std::endl;
      throw std::runtime_error(emsg.str()); 
    }
  }
  /*! \brief %Test globally for valid user input.
   *
   *   \param file     the __FILE__ value of the caller.
   *   \param lineNum  the __LINE__ value of the caller.
   *   \param msg      an optional descriptive message
   *   \param ok       a boolean which if false indicates an error
   *   \param level    a AssertionLevel value
   *   \param comm     a RCP<const Comm<int> > for the global check,
   *        if not specified we use the Environment's communicator.
   *
   *  If the \c level does not exceed the \c error_check_level parameter
   *  set by the user, then the assertion is tested on all processes in
   *  the \c comm.  If it fails on any process a std::runtime  error
   *  is thrown.
   */

  void globalInputAssertion(const char *file, int lineNum,
    const char *msg, bool ok, AssertionLevel level, 
    const Comm_t &comm) const {

    if (level <= errorCheckLevel_){
      int anyFail=0, fail = (!ok ? 1 : 0); 
      Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_MAX, 1, &fail, 
        &anyFail); 
      if (anyFail > 0){
        std::ostringstream emsg;
        emsg << myRank_<< ": " << file << "," << lineNum<< std::endl;
        if (msg && !ok)
          emsg << myRank_ << ": error: " << msg << std::endl;
        else
          emsg << myRank_ << ": exiting" << std::endl;
  
        throw std::runtime_error(emsg.str()); 
      }
    }
  }

  /*! \brief %Test for valid library behavior on local process only.
   *
   *   \param file     the __FILE__ value of the caller.
   *   \param lineNum  the __LINE__ value of the caller.
   *   \param msg      an optional descriptive message
   *   \param ok       a boolean which if false indicates an error
   *   \param level    a AssertionLevel value
   *
   *  If the \c level does not exceed the \c error_check_level parameter
   *  set by the user, then the assertion is tested and a std::logic_error
   *  error is thrown if it is false.
   *
   *  A failure of this test indicates a bug in Zoltan2.  Because of this
   *  consider doing these tests at the level of COMPLEX_ASSERTION, 
   *  so they they only get checked when we specifically ask for this higher
   *  level of checking.  An exception would be a test that is unlikely 
   *  to be executed (the default in a switch for example).
   */

  void localBugAssertion(const char *file, int lineNum,
    const char *msg, bool ok, AssertionLevel level) const {

    if (level <= errorCheckLevel_ && !ok){
      std::ostringstream emsg;
      emsg << myRank_<< ": " << file << "," << lineNum<< std::endl;
      if (msg)
        emsg << myRank_ << ": bug: " << msg << std::endl;
      throw std::logic_error(emsg.str()); 
    }
  }

  /*! \brief %Test for valid library behavior on every process.
   *
   *   \param file     the __FILE__ value of the caller.
   *   \param lineNum  the __LINE__ value of the caller.
   *   \param msg      an optional descriptive message
   *   \param ok       a boolean which if false indicates an error
   *   \param level    a AssertionLevel value
   *   \param comm     a RCP<const Comm<int> > for the global check,
   *        if not specified we use the Environment's communicator.
   *
   *  If the \c level does not exceed the \c error_check_level parameter
   *  set by the user, then the assertion is tested and a std::logic_error
   *  error is thrown if it is false on any process.
   *
   *  A failure of this test indicates a bug in Zoltan2.  Because of this
   *  consider doing these tests at the level of COMPLEX_ASSERTION, 
   *  so they they only get checked when we specifically ask for this higher
   *  level of checking.  An exception would be a test that is unlikely 
   *  to be executed (the default in a switch for example).
   */

  void globalBugAssertion(const char *file, int lineNum,
    const char *msg, bool ok, AssertionLevel level, 
   const Comm_t &comm) const {

    if (level <= errorCheckLevel_){
      int anyFail=0, fail = (!ok ? 1 : 0); 
      Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_MAX, 1, &fail, 
        &anyFail); 
      if (anyFail > 0){

        std::ostringstream emsg;
        emsg << myRank_<< ": " << file << "," << lineNum<< std::endl;
        if (msg && !ok)
          emsg << myRank_ << ": bug: " << msg << std::endl;
        else
          emsg << myRank_ << ": exiting" << std::endl;
  
        throw std::logic_error(emsg.str()); 
      }
    }
  }

  /*! \brief %Test for successful memory allocation on local process only.
   *
   *   \param file     the __FILE__ value of the caller.
   *   \param lineNum  the __LINE__ value of the caller.
   *   \param nobj     a value indicating the amount of memory wanted
   *   \param ok       a boolean which if false indicates failure
   *
   *  If the assertion fails, we throw std::bad_alloc.  There is no
   *  level to this macro because memory assertions are BASIC_ASSERTIONs.
   *  We always check for successful memory allocation unless compiled
   *  with -DZ2_OMIT_ALL_ERROR_CHECKING.
   */

  void localMemoryAssertion(const char *file, int lineNum, size_t nobj, 
    bool ok) const {

   if (!ok){ 
     std::cerr << myRank_ << ": " << file << ", " << lineNum<< std::endl;
     std::cerr << myRank_ << ": " << nobj << " objects" << std::endl;
     throw std::bad_alloc();
    } 
  }

  /*! \brief %Test for successful memory allocation on every process.
   *
   *   \param file     the __FILE__ value of the caller.
   *   \param lineNum  the __LINE__ value of the caller.
   *   \param nobj     a value indicating the amount of memory wanted
   *   \param ok       a boolean which if false indicates failure
   *   \param comm     a RCP<const Comm<int> > for the global check,
   *        if not specified we use the Environment's communicator.
   *
   *  If the assertion fails anywhere, we throw std::bad_alloc.  There is no
   *  level to this macro because memory assertions are BASIC_ASSERTIONs.
   *  We always check for successful memory allocation unless compiled
   *  with -DZ2_OMIT_ALL_ERROR_CHECKING.
   */

  void globalMemoryAssertion(const char *file, int lineNum,
    size_t nobj, bool ok, const Comm_t &comm) const {

    int anyFail=0, fail = (!ok ? 1 : 0); 
    Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_MAX, 1, &fail, &anyFail); 
    if (anyFail > 0){
      std::cerr << myRank_ << ": " << file << ", " << lineNum<< std::endl;
      if (!ok)
        std::cerr << myRank_ << ": " << nobj << " objects" << std::endl;
      else
        std::cerr << myRank_ << ": exiting" << std::endl;

      throw std::bad_alloc();
    }
  }
#endif

  // For debugging and profiling output, we define "char *" versions
  // as well as "string" versions to avoid runtime conversion of "char *"
  // to "string".

  /*! \brief  Send a message to the debug output manager.
   *
   *   \param level  If \c level does not exceed the \c debug_level
   *          parameter set by the user, then if this process is one that 
   *          prints debug messages (as indicated by the parameter 
   *          \c debug_profiling_procs) then the \c msg
   *          will be output to either \c debug_output_stream
   *          or \c debug_output_file.
   *   \param msg   The message to be output.
   */

#ifdef Z2_OMIT_ALL_STATUS_MESSAGES
  void debug(MessageOutputLevel level, const char *msg) const{ return;}
  void debug(MessageOutputLevel level, const std::string msg) const{ return;}
#else
  void debug(MessageOutputLevel level, const char *msg) const{ 
    debugOut_->print(level, msg);}
  void debug(MessageOutputLevel level, const std::string &msg) const{ 
    debugOut_->print(level, msg);}
#endif

#ifdef Z2_OMIT_ALL_PROFILING

  void timerStart(TimerType tt, const char * timerName) const  {return;}
  void timerStart(TimerType tt, const string &timerName) const  {return;}
  void timerStart(TimerType tt, const char * timerName, int, 
    int fieldWidth=0) const  {return;}
  void timerStart(TimerType tt, const string &timerName, int, 
    int fieldWidth=0) const  {return;}

  void timerStop(TimerType tt, const char * timerName) const {return;}
  void timerStop(TimerType tt, const string &timerName) const {return;}
  void timerStop(TimerType tt, const char * timerName, int, 
    int fieldWidth=0) const {return;}
  void timerStop(TimerType tt, const string &timerName, int, 
    int fieldWidth=0) const {return;}

#else
  /*! \brief  Start a named timer.
   */

  void timerStart(TimerType tt, const char *timerName) const  {
    if (timingOn_) timerOut_->start(tt, timerName); }

  void timerStart(TimerType tt, const string &timerName) const  {
    if (timingOn_) timerOut_->start(tt, timerName); }

  /*! \brief  Start a named timer, with a number as part of the name.
   */
  void timerStart(TimerType tt, const char *timerName, int num, 
    int fieldWidth=0) const  {
    if (timingOn_){
      ostringstream oss;
      oss << timerName << " ";
      if (fieldWidth > 0){
        oss.width(fieldWidth);
        oss.fill('0');
      }
      oss << num;
      timerOut_->start(tt, oss.str()); 
    }
  }

  void timerStart(TimerType tt, const string &timerName, int num, 
    int fieldWidth=0) const  {
    if (timingOn_){
      ostringstream oss;
      oss << timerName << " ";
      if (fieldWidth > 0){
        oss.width(fieldWidth);
        oss.fill('0');
      }
      oss << num;
      timerOut_->start(tt, oss.str()); 
    }
  }

  /*! \brief  Stop a named timer.
   */

  void timerStop(TimerType tt, const char *timerName) const {
    if (timingOn_) timerOut_->stop(tt, timerName); }

  void timerStop(TimerType tt, const string &timerName) const {
    if (timingOn_) timerOut_->stop(tt, timerName); }

  /*! \brief  Stop a named timer, with a number as part of the name.
   */

  void timerStop(TimerType tt, const char *timerName, int num, 
    int fieldWidth=0) const {
    if (timingOn_){
      ostringstream oss;
      oss << timerName << " ";
      if (fieldWidth > 0){
        oss.width(fieldWidth);
        oss.fill('0');
      }
      oss << num;
      timerOut_->stop(tt, oss.str()); 
    }
  }

  void timerStop(TimerType tt, const string &timerName, int num, 
    int fieldWidth=0) const {
    if (timingOn_){
      ostringstream oss;
      oss << timerName << " ";
      if (fieldWidth > 0){
        oss.width(fieldWidth);
        oss.fill('0');
      }
      oss << num;
      timerOut_->stop(tt, oss.str()); 
    }
  }
  
#endif

  /*! \brief Print a message and the kilobytes in use by this process.
   *
   *   \param msg   The message to be output. If this process
   *          is one that prints memory profiling messages
   *          (as indicated by the parameter \c memory_procs), 
   *          the \c msg (along with kilobytes currently allocated to
   *          this process) will issued.  The output target is either the
   *       \c memory_output_stream or \c memory_output_file.
   *          If neither was set, it goes to std::cout.
   *
   * Memory profiling is only supported on Linux nodes that 
   * have /proc/pid/statm.  If this is an unsupported node, the call 
   * does nothing.
   */

#ifdef Z2_OMIT_ALL_PROFILING
  void memory(const char *msg) const {return;}

  void memory(const std::string &msg) const {return; }
#else
  void memory(const char *msg) const
    {if (memoryOn_)
       memoryOut_->print(msg, getProcessKilobytes());}

  void memory(const std::string &msg) const
    {if (memoryOn_)
      memoryOut_->print(msg, getProcessKilobytes());}
#endif

  /*! \brief Returns a reference to the user's parameter list.
   *
   *  This is the parameter list after validation and modification.
   */
  const Teuchos::ParameterList &getParameters() const { return params_; }

  /*! \brief Returns a reference to a non-const copy of the parameters.
   *
   *  This is the parameter list after validation and modification.
   */
  Teuchos::ParameterList &getParametersNonConst() { return params_; }

  /*! \brief Return true if timing was requested, even if this
   *    process is not printing out timing messages.
   */
  bool doTiming() const { return timingOn_; }

  /*! \brief Return true if debug output was requested, even if
   *     this process is not printing out debug messages.
   */
#ifdef Z2_OMIT_ALL_STATUS_MESSAGES
  bool doStatus() const { return false;}
  MessageOutputLevel getDebugLevel() const {return NO_STATUS;}
#else
  bool doStatus() const { return debugOut_->getDebugLevel() > NO_STATUS;}

  MessageOutputLevel getDebugLevel() const {return debugOut_->getDebugLevel(); }
#endif

  /*! \brief Return true if memory usage output was requested, even if
   *     this process is not printing out memory used messages.
   */
  bool doMemoryProfiling() const { return memoryOn_;}

  /*! \brief Returns a const reference to the user's original list.
   *
   *  This is the parameter list before it was validated.  It is
   *  not the version supplied to algorithms.
   *  
   *  It is made available in case a Problem wants to create a new
   *  Environment by augmenting the user's original parameters.
   *  See PartitioningProblem::createPartitioningProblem() for an 
   *  example of doing this.
   */

  const Teuchos::ParameterList &getUnvalidatedParameters() const { 
    return unvalidatedParams_; }

  /*! \brief Convert parameters of type 
   *  Teuchos::StringToIntegralParameterEntryValidator<int> to integer.
   *
   * \param params  on input, a list of parameter, on return, all of the
   *                     StringToIntegral parameters have been converted
   *                     to integer values.
   *
   *    Given a parameter list, this function converts all of the entries that
   *    have valiator of type StringToIntegralParameterEntryValidator<int>
   *    from their string value to their int value.
   *
   */

  static void convertStringToInt(Teuchos::ParameterList &params);

private:

  /*! \brief Set up the Environment for the constructor.
   */
  void commitParameters();

  /*! \brief The Zoltan2 parameters supplied by the user.
   *
   *   Parameters lists are relatively small, so we keep a copy.
   */
  Teuchos::ParameterList unvalidatedParams_;

  /*! \brief The validated user parameters.
   *
   *  When the constructor calls commitParameters(), some of the
   *  input parameters are changed. For example, "yes"/"no" types
   *  of parameters are changed to the integers 0 or 1.  So
   *  this list of parameters provided to Problems and algorithms
   *  is this validated/converted list.
   */
  Teuchos::ParameterList params_;

  DebugManager_t debugOut_;    /*!< \brief output for status messages */

  Timer_t timerOut_;             /*!< \brief timer output */
  bool timingOn_;

  MemoryProfilerManager_t memoryOut_;  /*!< \brief memory profiling output */
  bool memoryOn_;
  RCP<std::ofstream> memoryOutputFile_;
};

/*! \brief A value to indicate a string parameter that was 
              not set by the user.
 */
#define Z2_UNSET_STRING std::string("notSet")

//////////////////////////////////////////////////////////////////////
// Templated namespace definitions used by the class

/*! \brief Create an output manager for a metric value.
 *
 *  \param rank  the MPI rank of the calling process in the application
 *  \param iPrint   true if this process should print metric information
 *  \param fname    name of file to which output is to be appended, or
 *                      or Z2_UNSET_STRING
 *  \param ost     output stream type
 *  \param mgr     on return, a pointer to the created output manager
 *  \param fptr    on return, an RCP to an ofstream if one was created.
 *
 * The template parameter is the data type of the entity being measured.
 */
template<typename metric_t>
  void makeMetricOutputManager(int rank, bool iPrint, 
    std::string fname, int ost,
    Teuchos::RCP<MetricOutputManager<metric_t> > &mgr,
    std::string units, int fieldWidth,
    RCP<std::ofstream> &fptr)
{
  typedef MetricOutputManager<metric_t> manager_t;

  OSType os = static_cast<OSType>(ost);

  bool haveFname = (fname != Z2_UNSET_STRING);
  bool haveStreamName = (os != NUM_OUTPUT_STREAMS);

  if (!haveFname && !haveStreamName){
    mgr = Teuchos::rcp(new manager_t(rank, iPrint, std::cout, true,
      units, fieldWidth));
    return;
  }

  if (haveFname){
    std::ofstream *oFile = NULL;
    if (iPrint){
      oFile = new std::ofstream;
      std::string newFname;
      addNumberToFileName(rank, fname, newFname);
      try{
        oFile->open(newFname.c_str(), std::ios::out|std::ios::trunc);
      }
      catch(std::exception &e){
        throw std::runtime_error(e.what());
      }
      fptr = rcp(oFile);
    }
    mgr = Teuchos::rcp(new manager_t(rank, iPrint, *oFile, true,
      units, fieldWidth));
    return;
  }

  if (os == COUT_STREAM)
    mgr = Teuchos::rcp(new manager_t(rank, iPrint, std::cout, true,
      units, fieldWidth));
  else if (os == CERR_STREAM)
    mgr = Teuchos::rcp(new manager_t(rank, iPrint, std::cerr, true,
      units, fieldWidth));
  else if (os == NULL_STREAM)
    mgr = Teuchos::rcp(new manager_t(rank, false, std::cout, true,
      units, fieldWidth));
  else
    throw std::logic_error("invalid metric output stream was not caught");
}

}  // namespace Zoltan2

#endif
