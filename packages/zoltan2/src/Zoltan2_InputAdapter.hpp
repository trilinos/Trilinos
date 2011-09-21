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

#define KDDKDD
#ifndef KDDKDD
  KDDKDD I think the InputAdapter class should NOT depend on Teuchos, Kokkos
  KDDKDD or the Zoltan2_Environment.  It is only an interface to the 
  KDDKDD application data.   
#include <Teuchos_Comm.hpp>
#include <Teuchos_RCP.hpp>
#include <Kokkos_DefaultNode.hpp>
#include <Zoltan2_Environment.hpp>
//#include <Kokkos_StandardNodeMemoryModel.hpp>
//#include <Kokkos_CUDANodeMemoryModel.hpp>
#endif // KDDKDD

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

#ifndef KDDKDD
  KDDKDD Since application developers are implementing the InputAdapter,
  KDDKDD the InputAdapter should not have storage; it should have only
  KDDKDD functions (as is currently in the GraphInput.hpp file.

  Teuchos::RCP<Teuchos::MpiComm<int> > _comm;

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
#endif //KDDKDD

public:

#ifndef KDDKDD
  // SR These virtual functions need to be defined by every implementation.
  // SR Leaving them out now until we have good use case. Reusing Karen's
  // macro here.

  virtual bool adapterComplete() = 0;

  /*! Direct the input adapter to complete processing input data.
      Some input adapters may choose to delay processing the input
      until they get a query. completeProcessing() directs that the
      processing must be done now.
   */
  virtual void completeProcessing() = 0;

  /*! Destructor */
//  virtual ~InputAdapter(){};

  /*! Copy Constructor */
//  virtual InputAdapter(const InputAdapter &env) = 0;

  /*! Assignment operator */
  virtual InputAdapter &operator=(const InputAdapter &rhs) = 0;
#endif

  /*! The caller's communicator */
#ifdef KDDKDD
  // KDDKDD Do we want to force Zoltan2 users to use Teuchos communicators?
  // KDDKDD Maybe we need two options -- getTeuchosComm and getMPIComm??
  // KDDKDD Either way, these functions must be implemented by the
  // KDDKDD application developer, not by us. 
  // KDDKDD The way to get a communicator from, say, a Tpetra matrix will
  // KDDKDD differ from the way to get it from a mesh.
  virtual const Teuchos::Comm<int> &getTeuchosComm() = 0;
  virtual const MPI_Comm &getMpiComm() = 0;
#else //KDDKDD
  const Teuchos::Comm<int> &getComm() const {return *_comm;}
#endif //KDDKDD

  /*! TODO perhaps we don't need these here.  If the only callers that
   *   use Kokkos are Tpetra users, then these can be pushed down to
   *   the Tpetra input adapters.

  MOVE THESE   will probably only be used with Tpetra objects

  bool haveMultiCoreNode() const { return !multicoreNode.is_null();}

  bool haveGpuNode() const { return !gpuNode.is_null();}

  Teuchos::RCP<const Kokkos_StandardNodeMemoryModel> &getMultiCoreNode() const { return multicoreNode;}

  Teuchos::RCP<const Kokkos_CUDANodeMemoryModel> &getGpuNode() const { return gpuNode;}
   */

#ifndef KDDKDD
  KDDKDD I think these functions do not belong in the InputAdapter.

  void setComm(MPI_Comm comm){
    MPI_Comm commCopy;
    MPI_Comm_dup(comm, &commCopy);
    MPI_Comm_set_errhandler(commCopy, MPI_ERRORS_RETURN);
    Teuchos::RCP<Teuchos::MPIComm<int> > tcomm = getTeuchosMPIComm(commCopy);
  }

  void freeComm(){
    MPI_Comm_free(static_cast<MPI_Comm *>(_comm->getRawMpiComm()->getRawPtr()));
  }
#endif

};
  
  
}  //namespace Zoltan2
  
#endif
