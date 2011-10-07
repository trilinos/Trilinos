// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Zoltan2_Environment.hpp>
#include <Zoltan2_Standards.hpp>

namespace Zoltan2{

// Using C-language MPI rather than C++ because it seems to be more portable.

/*! Convert an MPI communicator to a MpiComm object.
 */

#ifdef HAVE_MPI
template <typename Ordinal>
  RCP<MpiComm<Ordinal> >
    getTeuchosMpiComm(const MPI_Comm &comm)
{
  RCP<Teuchos::OpaqueWrapper<MPI_Comm> >handle = Teuchos::opaqueWrapper<MPI_Comm>(comm); 
  RCP<MpiComm<Ordinal> > tcommPtr(new MpiComm<Ordinal>(handle));

  return tcommPtr;
}

#endif

}  //namespace Zoltan2
