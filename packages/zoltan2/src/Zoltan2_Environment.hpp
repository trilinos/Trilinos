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

#ifndef _ZOLTAN2_ENVIRONMENT_DECL_HPP_
#define _ZOLTAN2_ENVIRONMENT_DECL_HPP_

/*! \file Zoltan2_Environment_decl.hpp
  
  \brief Declares the Zoltan2::Environment class.

*/

#include <Zoltan2_config.h>
#include <Zoltan2_Parameters.hpp>
#include <Zoltan2_DebugManager.hpp>

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

  Teuchos::RCP<const Teuchos::Comm<int> > comm_; /*!< from original problem */

  Teuchos::RCP<DebugManager> debugOut_;  /*!< output for status messages */

  Teuchos::RCP<DebugManager> timerOut_;  /*!< output for timing messages */

  Teuchos::RCP<DebugManager> memoryOut_; /*!< output for memory usage messages*/

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

  /*! \brief Returns a reference to the user's parameter list.
   */
  const Teuchos::ParameterList &getParams() const { return params_; }

  /*! \brief Returns a reference to a non-const copy of the parameters.
   */
  Teuchos::ParameterList &getParamsNonConst() { return params_; }

  /*! \brief Returns a reference to the user's partitioning parameters.
   */
  const Teuchos::ParameterList &getPartitioningParams() const 
  { 
    return params_.sublist("partitioning"); 
  }

  /*! \brief Returns a reference to the user's partitioning parameters.
   */
  Teuchos::ParameterList &getPartitioningParamsNonConst() 
  { 
    return params_.sublist("partitioning"); 
  }

  /*! \brief Returns a reference to the user's ordering parameters.
   */
  const Teuchos::ParameterList &getOrderingParams() const 
  { 
    return params_.sublist("ordering"); 
  }

  /*! \brief Returns a reference to the user's ordering parameters.
   */
  Teuchos::ParameterList &getOrderingParamsNonConst()
  { 
    return params_.sublist("ordering"); 
  }

  /*! \brief Returns a reference to the user's coloring parameters.
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

  /*! \brief Return true if timing was requested.
   */
  bool doTiming() const { return timerOut_->getDebugLevel() > NO_STATUS;}

  /*! \brief Return true if debug output was requested.
   */
  bool doStatus() const { return timerOut_->getDebugLevel() > NO_STATUS;}

  /*! \brief Return true if memory usage output was requested.
   */
  bool doMemoryProfiling() const { 
    return memoryOut_->getDebugLevel() > NO_STATUS;}

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
};

}  // namespace Zoltan2

#endif
