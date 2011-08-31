// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Zoltan2_Environment.hpp>
#include <Zoltan2_Exceptions.hpp>
#include <mpi.h>
#include <sstream>
#include <fstream>
#include <ostream>

namespace Zoltan2{

// Using C-language MPI rather than C++ because it seems to be more portable.

// TODO move ostreams to environment methods rather than parameterList entries
void getOutputStreamFromParameterList(
  Teuchos::ParameterList &pl, std::string key, std::ostream *&os,
  std::ostream &defaultValue);

/*! Convert an MPI communicator to a Teuchos::MpiComm object.
 */

template <typename Ordinal>
  Teuchos::RCP<Teuchos::MpiComm<Ordinal> >
    getTeuchosMpiComm(const MPI_Comm comm)
{
  Teuchos::RCP<Teuchos::OpaqueWrapper<MPI_Comm> >handle = Teuchos::opaqueWrapper<MPI_Comm>(comm);
  Teuchos::MpiComm<Ordinal> tcomm(handle);
  Teuchos::RCP<Teuchos::MpiComm<Ordinal> > tcommPtr(&tcomm);
  return tcommPtr;
}


/*! Create a sub communicator. 
    \param comm  The original communicator
    \param members  The ranks in the original communicator of the members in the sub communicator
    \return A communicator containing only the specified members of the original communicator

   All processes in the original communicator must call this.
 */

template <typename Ordinal>
  Teuchos::RCP<Teuchos::MpiComm<Ordinal> >
    getTeuchosMpiSubComm(const Teuchos::RCP<Teuchos::MpiComm<Ordinal> > comm, 
                          const std::vector<Ordinal> members,
                          const Environment &env)
{
  // TODO test this only once during the library execution
  Z2_GLOBAL_BUG_ASSERTION(*comm, env, "size of Ordinal assumption", 
    sizeof(Ordinal) != sizeof(int), Z2_BASIC_ASSERTION);

  Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > superComm = comm->getRawMpiComm();
  MPI_Group mainGroup, subGroup;
  MPI_Comm  mainComm, subComm;

  int rc = MPI_Comm_group(*superComm, &mainGroup);
  Z2_LOCAL_INPUT_ASSERTION(*comm, env, "obtaining group", 
    rc != MPI_SUCCESS, Z2_BASIC_ASSERTION);

  rc = MPI_Group_incl(mainGroup, members.size(), (int *)(&members[0]), &subGroup);
  Z2_LOCAL_INPUT_ASSERTION(*comm, env, "creating subgroup", 
    rc != MPI_SUCCESS, Z2_BASIC_ASSERTION);

  rc = MPI_Comm_create(mainComm, subGroup, &subComm);
  Z2_LOCAL_INPUT_ASSERTION(*comm, env, "creating sub communicator", 
    rc != MPI_SUCCESS, Z2_BASIC_ASSERTION);

  MPI_Comm_set_errhandler(subComm, MPI_ERRORS_RETURN);

  return getTeuchosMpiComm<Ordinal>(subComm);
}


/*! Create a sub communicator.
    \param comm  The original communicator
    \param color All processes supplying the same value for color will be in the same subcommunicator.
    \return A communicator containing only the members of the original communicator that supplied the same color in the call.
 */
template <typename Ordinal>
  Teuchos::RCP<Teuchos::MpiComm<Ordinal> >
    getTeuchosMpiSubComm(const Teuchos::RCP<Teuchos::MpiComm<Ordinal> > comm, 
                          const Ordinal color, const Environment &env)
{
  MPI_Comm subComm;
  Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > superComm = comm->getRawMpiComm();

  int rc = MPI_Comm_split(*superComm, color, 0, &subComm);
  Z2_LOCAL_INPUT_ASSERTION(*comm, env, "creating sub communicator", 
    rc != MPI_SUCCESS, Z2_BASIC_ASSERTION);

  MPI_Comm_set_errhandler(subComm, MPI_ERRORS_RETURN);

  return getTeuchosMpiComm<Ordinal>(subComm);
}


}  //namespace Zoltan2
