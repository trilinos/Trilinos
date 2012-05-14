// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_Environment.hpp
    \brief Defines the Environment class.
*/

#ifndef _ZOLTAN2_ENVIRONMENT_HPP_
#define _ZOLTAN2_ENVIRONMENT_HPP_

#include <Zoltan2_config.h>
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

/*!  \brief The parameters, debug and profiling output objects, and error checking methods.

  This is object is passed to almost every method in the library. It
  has the problem parameters and the configuration information that governs
  how the library should behave when reporting status information,
  testing for errors, and so on.

  The environment holds the application's default communicator.  Note that 
  this communicator may differ from the problem communicator for any
  given problem.

  \todo The Environment could make it easier for algorithms to get parameter
    values without having to navigate Teuchos::ParameterLists and sublists.
*/

class Environment{

public:

  typedef long   memory_t;

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
  Environment(Teuchos::ParameterList &problemParams, const Comm_t &comm);

  /*! \brief Default Constructor
   *
   *    The default constructor uses the Teuchos default communicator,
   *    BASIC_STATUS for debug_level, and NO_STATUS for timing_level and
   *    memory_profiling_level.  It has error_check_level BASIC_ASSERTION.
   *    It has no other parameters.
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

  void setTimer(RCP<TimerManager> &timer) { timerOut_=timer; timingOn=true;}

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

  /*! \brief Test for valid user input on local process only.
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
  /*! \brief Test globally for valid user input.
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

  /*! \brief Test for valid library behavior on local process only.
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

  /*! \brief Test for valid library behavior on every process.
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

  /*! \brief Test for successful memory allocation on local process only.
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

  /*! \brief Test for successful memory allocation on every process.
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

  void debug(MessageOutputLevel level, const char *msg) const{ 
    debugOut_->print(level, msg);}

#ifdef Z2_OMIT_ALL_ERROR_CHECKING

  // These should inline to a no-op.  
  void timerStart(string timerName) const  {return;}
  void timerStop(string timerName) const {return;}

#else
  /*! \brief  Start a named timer.
   */

  void timerStart(string timerName) const  {
    if (timingOn) timerOut_->start(timerName); }

  /*! \brief  Stop a named timer.
   */

  void timerStop(string timerName) const {
    if (timingOn) timerOut_->stop(timerName); }
#endif

  /*! \brief  Send a message to the memory profiling output manager.
   *
   *   \param msg   The message to be output. If this process
   *          is one that prints memory profiling messages
   *          (as indicated by the parameter \c memory_profiling_procs), 
   *          the \c msg will be output to either the
   *       \c memory_profiling_output_stream or \c memory_profiling_output_file.
   *   \param units A string indicating the memory_profiling units.
   *   \param val   The memory_profiling value, a long.
   */

  void memory(const char *msg, const char *units, const memory_t val) const{ 
    memoryOut_->print(msg, units, val); }

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

  void debug(MessageOutputLevel level, const std::string &msg) const{ 
    debugOut_->print(level, msg);}

  /*! \brief  Send a message to the memory profiling output manager.
   *
   *   \param msg   The message to be output. If this process
   *          is one that prints memory profiling messages
   *          (as indicated by the parameter \c memory_profiling_procs), 
   *          the \c msg will be output to either the
   *       \c memory_profiling_output_stream or \c memory_profiling_output_file.
   *   \param units A string indicating the memory_profiling units.
   *   \param val   The memory_profiling value, a long.
   */

  void memory(const std::string &msg, const std::string &units, 
    const memory_t val) const{ memoryOut_->print(msg, units, val); }

  /*! \brief Returns true if the parameter list has the named sublist.
   */
  bool hasSublist(const Teuchos::ParameterList &pl, 
    const std::string &plname) const;

  /*! \brief Returns true if the parameter list has a "partitioning" sublist.
   */
  bool hasPartitioningParameters() const { 
   return hasSublist(params_, "partitioning"); }

  /*! \brief Returns true if there is a "ordering" parameter sublist.
   */
  bool hasOrderingParameters() const {
   return hasSublist(params_, "ordering"); }

  /*! \brief Returns true if there is a "matching" parameter sublist.
   */
  bool hasMatchingParameters() const {
   return hasSublist(params_, "matching"); }

  /*! \brief Returns true if there is a "coloring" parameter sublist.
   */
  bool hasColoringParameters() const {
   return hasSublist(params_, "coloring"); }

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
  bool doTiming() const { return timingOn; }

  /*! \brief Return true if debug output was requested, even if
   *     this process is not printing out debug messages.
   */
  bool doStatus() const { return debugOut_->getDebugLevel() > NO_STATUS;}

  /*! \brief Return true if memory usage output was requested, even if
   *     this process is not printing out memory used messages.
   */
  bool doMemoryProfiling() const { return memoryOut_->getMetricsOn(); }

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
  bool timingOn;

  MemoryProfilerManager_t memoryOut_;  /*!< \brief memory profiling output */
};

/*! \brief A value to indicate a string parameter that was 
              not set by the user.
 */
#define Z2_UNSET std::string("notSet")

//////////////////////////////////////////////////////////////////////
// Templated namespace definitions used by the class

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
template<typename metric_t>
  void makeMetricOutputManager(int rank, bool iPrint, 
    std::string fname, std::string osname, 
    Teuchos::RCP<MetricOutputManager<metric_t> > &mgr)
{
  typedef MetricOutputManager<metric_t> manager_t;
  bool haveFname = (fname != Z2_UNSET);
  bool haveStreamName = (!haveFname && (osname != Z2_UNSET));

  if (!haveFname && !haveStreamName){
    mgr = Teuchos::rcp(new manager_t(rank, iPrint, std::cout, true));
    return;
  }

  if (haveFname){
    std::ofstream *oFile = new std::ofstream;
    if (iPrint){
      std::string newFname;
      addNumberToFileName(rank, fname, newFname);
      try{
        oFile->open(newFname.c_str(), std::ios::out|std::ios::trunc);
      }
      catch(std::exception &e){
        throw std::runtime_error(e.what());
      }
    }
    mgr = Teuchos::rcp(new manager_t(rank, iPrint, *oFile, true));
    return;
  }

  if (osname == std::string("std::cout"))
    mgr = Teuchos::rcp(new manager_t(rank, iPrint, std::cout, true));
  else if (osname == std::string("std::cerr"))
    mgr = Teuchos::rcp(new manager_t(rank, iPrint, std::cerr, true));
  else if (osname == std::string("/dev/null"))
    mgr = Teuchos::rcp(new manager_t(rank, false, std::cout, true));
  else
    throw std::logic_error("invalid metric output stream was not caught");
}

}  // namespace Zoltan2

#endif
