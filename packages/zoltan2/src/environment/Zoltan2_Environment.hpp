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
*/


class Environment{

public:

  // These are accessed directly rather than through get methods
  // due to the frequency with which we'll be referring to them
  // all over the code.

  int  myRank_;             /*!< mpi rank in original problem */

  int  numProcs_;           /*!< number of processes in original problem */

  Teuchos::RCP<const Teuchos::Comm<int> > comm_; /*!< for application */

  Teuchos::RCP<DebugManager> debugOut_;  /*!< output for status messages */

  Teuchos::RCP<MetricOutputManager<double> > timerOut_;  /*!< output for timing messages */

  Teuchos::RCP<MetricOutputManager<long> > memoryOut_; /*!< output for memory usage messages*/

  AssertionLevel errorCheckLevel_; /*!< level of error checking to do */

  /*! \brief Constructor
   */
  Environment( Teuchos::ParameterList &problemParams, 
    const Teuchos::RCP<const Teuchos::Comm<int> > &comm);

  /*! \brief Default Constructor
   */
  Environment();

  /*! \brief Destructor
   */
  ~Environment();

  /*! \brief Set or reset the application communicator.
   */
  void setCommunicator(const Teuchos::RCP<const Teuchos::Comm<int> > &comm);

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

  bool hasParameters() const { return hasAnyParams_; }

  bool hasPartitioningParameters() const { return hasPartitioningParams_; }

  bool hasOrderingParameters() const {return hasOrderingParams_; }

  bool hasMatchingParameters() const { return hasMatchingParams_; }

  bool hasColoringParameters() const { return hasColoringParams_; }

  /*! \brief Returns a reference to the user's parameter list.
   *
   *   If there are no parameters, this call throws a std::exception.
   */
  const Teuchos::ParameterList &getParams() const { return params_; }

  /*! \brief Returns a reference to a non-const copy of the parameters.
   *   If there are no parameters, this call throws a std::exception.
   */
  Teuchos::ParameterList &getParamsNonConst() { return params_; }

  /*! \brief Returns a reference to the user's partitioning parameters.
   *   If there are no parameters, this call throws a std::exception.
   */
  const Teuchos::ParameterList &getPartitioningParams() const 
  { 
    return params_.sublist("partitioning"); 
  }

  /*! \brief Returns a reference to the user's partitioning parameters.
   *   If there are no parameters, this call throws a std::exception.
   */
  Teuchos::ParameterList &getPartitioningParamsNonConst() 
  { 
    return params_.sublist("partitioning"); 
  }

  /*! \brief Returns a reference to the user's ordering parameters.
   *   If there are no parameters, this call throws a std::exception.
   */
  const Teuchos::ParameterList &getOrderingParams() const 
  { 
    return params_.sublist("ordering"); 
  }

  /*! \brief Returns a reference to the user's ordering parameters.
   *   If there are no parameters, this call throws a std::exception.
   */
  Teuchos::ParameterList &getOrderingParamsNonConst()
  { 
    return params_.sublist("ordering"); 
  }

  /*! \brief Returns a reference to the user's coloring parameters.
   *   If there are no parameters, this call throws a std::exception.
   */
  const Teuchos::ParameterList &getColoringParams() const 
  { 
    return params_.sublist("coloring"); 
  }

  /*! \brief Returns true if the user parameters have been processed.
   *
   * The user's parameters are processed in commitParameters(). During
   * this call they are checked for validity, and some are transformed
   * to an internal format.  (For example "std::cout" will be transformed
   * to a pointer to that ostream.)
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

private:
  /*! \brief The Zoltan2 parameters supplied by the user.
   */
  Teuchos::ParameterList params_;

  /*! \brief The list of valid parameters against which to
   *             check the user parameters.
   */
  Teuchos::ParameterList validParams_;

  /*! \brief The user parameters can be validated only once.  True
   *            if this validation has occured.
   */
  bool committed_;

  bool hasAnyParams_;
  bool hasPartitioningParams_;
  bool hasOrderingParams_;
  bool hasColoringParams_;
  bool hasMatchingParams_;
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
