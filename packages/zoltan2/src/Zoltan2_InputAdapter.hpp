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

   TODO: add methods that return a view of the vertices and edges, not a copy.

*/

class InputAdapter{

private:
  
protected:

  Teuchos::RCP<Zoltan2::Environment> _env;

  /*! True if the application input to the adapter is finished.
      Some input adapters may be initialized with the application data
      in several steps.  They should indicate with this boolean when
      input is complete.
   */
  bool _inputComplete;

  /*! True if adapter processing of application input is complete.
      Some input adapters may need to process the application input
      before they can respond to queries. They should indicate with 
      this boolean when such process has been done.
   */
  bool _adapterProcessingComplete;

public:

  bool adapterComplete() {return _adapterProcessingComplete;}

  /*! Direct the input adapter to complete processing input data.
      Some input adapters may choose to delay processing the input
      until they get a query. completeProcessing() directs that the
      processing must be done now.
   */
  virtual void completeProcessing() = 0;

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

  MOVE THESE   will probably only be used with Tpetra objects

  bool haveMultiCoreNode() const { return !multicoreNode.is_null();}

  bool haveGpuNode() const { return !gpuNode.is_null();}

  Teuchos::RCP<const Kokkos_StandardNodeMemoryModel> &getMultiCoreNode() const { return multicoreNode;}

  Teuchos::RCP<const Kokkos_CUDANodeMemoryModel> &getGpuNode() const { return gpuNode;}
   */

};
  
  
}  //namespace Zoltan2
  
#endif
