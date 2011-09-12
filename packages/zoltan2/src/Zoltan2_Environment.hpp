// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_Environment.hpp

    \brief The declarations for the Environment object.
*/


#ifndef _ZOLTAN2_ENVIRONMENT_DECL_HPP_
#define _ZOLTAN2_ENVIRONMENT_DECL_HPP_

/*! \file Zoltan2_Environment_decl.hpp
  
  \brief Declares the Zoltan2::Environment class.

*/

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Kokkos_DefaultNode.hpp>
#include <Zoltan2_Parameters.hpp>
#include <Zoltan2_DebugManager.hpp>

namespace Zoltan2 {

/*! Zoltan2::Environment
    \brief The parameters and other information needed at runtime.

  This is object is passed to almost every method in the library.  We may want
  to have a memory manager here as well.
*/


class Environment{

private:
  /*! The user's parameters
   */
  Teuchos::ParameterList _params;

  /*! The parameters, their defaults, their validators,
   *  and their documentation
   */
  Teuchos::ParameterList _validParams;

public:
  /*! Values needed for quick access
   */
  int  _myRank;
  int  _numProcs;
  bool _printDebugMessages; // true if this proc prints
  bool _printProfilingMessages;// true if this proc prints
  int  _errorCheckLevel;    // how vigilant are we with assertions
  int  _debugDepthLevel;    // how much info do we write out
  int  _profilingIndicator;    // how much profiling (should really be
                           // "what are we profiling", not a level)
  bool _committed;
  Teuchos::RCP<const Teuchos::Comm<int> > _comm;

#if 0
  /*! The node description is not yet implemented.
   */
  Kokkos::CUDANodeMemoryModel     _gpuNode;
  Kokkos::StandardNodeMemoryModel _standardNode;
  
  /*! The machine model is not yet implemented.  It will
      not necessarily be implemented as a ParameterList.
   */
  Teuchos::ParameterList _machine;
#endif

  /*! Constructor 
      Because parameter lists are small, we save a
      copy of them instead of requiring an RCP.
   */
  Environment(Teuchos::ParameterList &prob, Teuchos::RCP<Teuchos::Comm<int> > &comm);

  /*! Constructor
   */
  Environment();

  /*! Destructor */
  virtual ~Environment();

  /*! Copy Constructor */
  Environment(const Environment &env);

  /*! Assignment operator */
  Environment &operator=(const Environment &env);

  /*! Set communicator */
  void setCommunicator(Teuchos::RCP<Teuchos::Comm<int> > &comm);

  /*! Set or reset the problem parameters*/
  void setParameters(Teuchos::ParameterList &Params);

  /*! Add and to or replace the existing parameters with these */
  void addParameters(Teuchos::ParameterList &Params);

  /*! Parameter setting is done, process parameters for use.
      This should be done in the Problem.
   */
  void commitParameters();

  /*! Get a reference to a read-only copy of the parameters. */
  const Teuchos::ParameterList &getParameters() const;

  /*! The debug manager, used by debug statements */
  Teuchos::RCP<Zoltan2::DebugManager> _dbg;

  // TODO: some kind of a manager for profiling, with
  //   labels, hierarchy, and so on.

};
  
  
}  //namespace Zoltan2
  
#endif
