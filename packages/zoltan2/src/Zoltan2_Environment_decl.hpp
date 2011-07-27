// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

#ifndef _ZOLTAN2_ENVIRONMENT_DECL_HPP_
#define _ZOLTAN2_ENVIRONMENT_DECL_HPP_

/*! \file Zoltan2_Environment_decl.hpp
  
  \brief Declares the ZOLTAN2::Environment class.

*/

#include <Teuchos_ParameterList.hpp>
#include <Kokkos_DefaultNode.hpp>
#include <Zoltan2_DebugManager.hpp>

/*! \namespace ZOLTAN2

  \brief The ZOLTAN2 namespace contains symbols that are part of the Zoltan2 API.
*/

namespace ZOLTAN2 {

/*! ZOLTAN2::Environment
    \brief The problem parameters, library configuration, and other information.

  This is object is passed to almost every method in the library.  We may want
  to have a memory manager here as well.

  TODO: add teuchos validators

  TODO: Not sure whether we'll need AppGID (fixed vertices?) and the
    others.  If not we can remove them.
*/

template <typename AppLID, typename AppGID, typename LNO=int, 
  typename GNO=AppGID, typename NODE=Kokkos::DefaultNode>
    class Environment{

public:
  /*! Constructor 
   */
  Environment(Teuchos::ParameterList &prob, Teuchos::ParameterList &config,
    NODE &node, Teuchos::ParameterList &machine);

  /*! Constructor 
   */
  Environment(Teuchos::ParameterList &prob, Teuchos::ParameterList &config);


  /*! Constructor
   */
  Environment();

  /*! Destructor */
  ~Environment();

  /*! Copy Constructor */
  Environment(const Environment &env);

  /*! Assignment operator */
  Environment &operator=(const Environment &env);

  /*! Set or reset the list of problem parameters
   */
  void setProblemParameters(Teuchos::ParameterList &p);

  /*! Set or reset the library configuration settings.
   */
  void setLibraryConfiguration(Teuchos::ParameterList &c);
  
  /*! The problem parameters
   */
  Teuchos::ParameterList params;
  
  /*! The library configuration 
   */
  Teuchos::ParameterList config;
  
  /*! The node description
  
    Not yet implemented.
   */
  NODE node;
  
  /*! The machine model
  
    Not yet implemented.  Not necessarily implemented as a ParameterList.
   */
  Teuchos::ParameterList machine;

  /*! The debug manager, used by debug statements
   */

  Z2::DebugManager dbg;

  /*! The level of checking to do at runtime.
   */
  int errorCheckLevel;

private:

static void getOutputStream(Teuchos::ParameterList &pl, char *key, std::ostream *&os);

};
  
  
}  //namespace ZOLTAN2
  
#endif
