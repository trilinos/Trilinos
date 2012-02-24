// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_Environment.hpp

    \brief The Environment object, which governs library behavior.
*/

#ifndef _ZOLTAN2_ENVIRONMENT_HPP_
#define _ZOLTAN2_ENVIRONMENT_HPP_

/*! \file Zoltan2_Environment.hpp
  
  \brief Declares the Zoltan2::Environment class.

*/

#include <Zoltan2_config.h>
#include <Zoltan2_Parameters.hpp>
#include <Zoltan2_DebugManager.hpp>
#include <Zoltan2_MetricOutputManager.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultComm.hpp>

namespace Zoltan2 {

/*!  \brief The parameters and other information needed at runtime.

  This is object is passed to almost every method in the library. It
  has the problem parameters and the configuration information that governs
  how the library should behave when reporting status information,
  testing for errors, and so on.

  The environment has the application's default communicator, unless another
  communicator is specified in the constructor or with setCommunicator().
  Note that this communicator may differ from the problem communicator for any
  given problem.
*/

class Environment{

public:

  typedef Teuchos::RCP<const Teuchos::Comm<int> > Comm_t;
  typedef Teuchos::RCP<DebugManager>     DebugManager_t;
  typedef Teuchos::RCP<MetricOutputManager<double> > TimerManager_t;
  typedef Teuchos::RCP<MetricOutputManager<long> > MemoryProfilerManager_t;

  // These are accessed directly rather than through get methods
  // due to the frequency with which we'll be referring to them
  // all over the code.

  int  myRank_;                        /*!< mpi rank (relative to comm_) */

  int  numProcs_;           /*!< number of processes (relative to comm_) */

  Comm_t comm_;                         /*!< communicator for environment*/

  DebugManager_t debugOut_;              /*!< output for status messages */

  TimerManager_t timerOut_;                            /*!< timer output */

  MemoryProfilerManager_t memoryOut_;       /*!< memory profiling output */

  AssertionLevel errorCheckLevel_;      /*!< level of error checking to do */

  /*! \brief Constructor
   */
  Environment(Teuchos::ParameterList &problemParams, const Comm_t &comm);

  /*! \brief Default Constructor
   */
  Environment();

  /*! \brief Destructor
   */
  ~Environment();

  /*! \brief Set or reset the environment's communicator.
   */
  void setCommunicator(const Comm_t &comm);

  /*! \brief Set or reset the user parameters.
   */
  void setParameters(Teuchos::ParameterList &params);

  /*! \brief Add to the user parameters.
   *  User parameters can only be added before calling commitParameters().
   */
  void addParameters(Teuchos::ParameterList &params);

  /*! \brief Finalize the user's parameters.
   *  
   * When commitParameters() is called, the user's parameters are
   * validated and some are changed from parameter strings to the internal
   * values that Zoltan2 will use.  For this reason, commitParameters()
   * should only be called once.
   */
  void commitParameters();

  /*! brief Returns true if the parameter list has the named sublist
   */
  bool hasSublist(const Teuchos::ParameterList &pl, 
    const std::string &plname) const;

  /*! brief Returns true if the parameter list has the named sublist
   */
  bool hasPartitioningParameters() const { 
   return hasSublist(params_, "partitioning"); }

  /*! brief Returns true if there is a "ordering" parameter sublist.
   */
  bool hasOrderingParameters() const {
   return hasSublist(params_, "ordering"); }

  /*! brief Returns true if there is a "matching" parameter sublist.
   */
  bool hasMatchingParameters() const {
   return hasSublist(params_, "matching"); }

  /*! brief Returns true if there is a "coloring" parameter sublist.
   */
  bool hasColoringParameters() const {
   return hasSublist(params_, "coloring"); }

  /*! \brief Returns a reference to the user's parameter list.
   */
  const Teuchos::ParameterList &getParameters() const { return params_; }

  /*! \brief Returns a reference to a non-const copy of the parameters.
   */
  Teuchos::ParameterList &getParametersNonConst() { return params_; }

  /*! \brief Returns true if the user parameters have been processed.
   *
   * The user's parameters are processed in commitParameters(). During
   * this call they are checked for validity, and some are transformed
   * to an internal format.  
   *
   * If parameters are already committed, you can not call
   * addParameters() or setParameters().
   */
  bool parametersAreCommitted() const { return committed_; }

  /*! \brief Return true if timing was requested, even if this
   *    process is not printing out timing messages.
   */
  bool doTiming() const { return timerOut_->getMetricsOn(); }

  /*! \brief Return true if debug output was requested, even if
   *     this process is not printing out debug messages.
   */
  bool doStatus() const { return debugOut_->getDebugLevel() > NO_STATUS;}

  /*! \brief Return true if memory usage output was requested, even if
   *     this process is not printing out memory used messages.
   */
  bool doMemoryProfiling() const { return memoryOut_->getMetricsOn(); }

  /*! \brief Given a parameter list, convert all of the entries that
   *     have valiator of type StringToIntegralParameterEntryValidator<int>
   *     from their string value to their int value.
   *
   *    I expected this to happen in validateAndModify, but
   *    it doesn't.
   */
  static void convertStringToInt(Teuchos::ParameterList &params);

private:

  /*! \brief The Zoltan2 parameters supplied by the user.
   *
   *   Parameters lists are relatively small, so we keep a copy.
   */
  Teuchos::ParameterList params_;

  /*! \brief The user parameters can be validated only once.  True
   *            if this validation has occured.
   */
  bool committed_;
};

/*! A helper function to return a default environment.
 *
 *    A default environment has the default communicator and the
 *    the default parameters.  This is helpful when writing tests
 *    of classes or methods that need an environment.
 */

Teuchos::RCP<const Environment> getDefaultEnvironment();


}  // namespace Zoltan2

#endif
