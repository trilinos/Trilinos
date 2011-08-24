// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_InputAdapter_def.hpp

    \brief The abstract interface for a Zoltan2 input adapter.
*/

#ifndef _ZOLTAN2_INPUTADAPTER_DECL_HPP_
#define _ZOLTAN2_INPUTADAPTER_DECL_HPP_

#include <Teuchos_Comm.hpp>
#include <Teuchos_RCP.hpp>
#include <Kokkos_DefaultNode.hpp>
#include <Kokkos_StandardNodeMemoryModel.hpp>
#include <Kokkos_CUDANodeMemoryModel.hpp>

namespace Zoltan2 {

/*! Zoltan2::InputAdapter
    \brief The InputAdapter presents a uniform interface to Zoltan2 for many types of user input.

   Zoltan2 defines input adapters for many common types of user input.  In addition, users
   can write their own input adapter.

*/

class InputAdapter{

private:

public:

  /*! Destructor */
  virtual ~InputAdapter(){};

  /*! Copy Constructor */
  virtual InputAdapter(const InputAdapter &env) = 0;

  /*! Assignment operator */
  virtual InputAdapter &operator=(const InputAdapter &env) = 0;

  /*! The caller's communicator */
  const Teuchos::Comm<int> &getComm() const {return comm;}

  /*! TODO perhaps we don't need these here.  If the only callers that
   *   use Kokkos are Tpetra users, then these can be pushed down to
   *   the Tpetra input adapters.

  MOVE THESE

  bool haveMultiCoreNode() const { return !multicoreNode.is_null();}

  bool haveGpuNode() const { return !gpuNode.is_null();}

  Teuchos::RCP<const Kokkos_StandardNodeMemoryModel> &getMultiCoreNode() const { return multicoreNode;}

  Teuchos::RCP<const Kokkos_CUDANodeMemoryModel> &getGpuNode() const { return gpuNode;}
   */

};
  
  
}  //namespace Zoltan2
  
#endif
