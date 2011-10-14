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

#include <Zoltan2_Parameters.hpp>
#include <Zoltan2_DebugManager.hpp>
#include <Zoltan2_Standards.hpp>

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
  ParameterList params_;

  /*! The parameters, their defaults, their validators,
   *  and their documentation
   */
  ParameterList validParams_;

public:
  /*! Values needed for quick access
   */
  int  myRank_;
  int  numProcs_;
  bool printDebugMessages_; // true if this proc prints
  bool printProfilingMessages_;// true if this proc prints
  int  errorCheckLevel_;    // how vigilant are we with assertions
  int  debugDepthLevel_;    // how much info do we write out
  int  profilingIndicator_;    // how much profiling (should really be
                           // "what are we profiling", not a level)
  bool committed_;
  RCP<const Comm<int> > comm_;

#if 0
  /*! The node description is not yet implemented.
   */
  Kokkos::CUDANodeMemoryModel     gpuNode_;
  Kokkos::StandardNodeMemoryModel standardNode_;
  
  /*! The machine model is not yet implemented.  It will
      not necessarily be implemented as a ParameterList.
   */
  ParameterList machine_;
#endif

  /*! Constructor 
      Because parameter lists are small, we save a
      copy of them instead of requiring an RCP.
   */
  Environment(ParameterList &prob, 
    const RCP<const Comm<int> > &comm);

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
  void setCommunicator(const RCP<const Comm<int> > &comm);

  /*! Set or reset the problem parameters*/
  void setParameters(ParameterList &Params);

  /*! Add and to or replace the existing parameters with these */
  void addParameters(ParameterList &Params);

  /*! Parameter setting is done, process parameters for use.
      This should be done in the Problem.
   */
  void commitParameters();

  /*! Get a reference to a read-only copy of the parameters. */
  const ParameterList &getParameters() const;

  /*! The debug manager, used by debug statements */
  RCP<Zoltan2::DebugManager> dbg_;

  // TODO: some kind of a manager for profiling, with
  //   labels, hierarchy, and so on.

};
  
  
}  //namespace Zoltan2
  
#endif
