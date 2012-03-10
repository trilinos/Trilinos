// @HEADER
// ***********************************************************************
//                Copyright message goes here. 
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_Util.hpp
 *  \brief A gathering of useful namespace methods.
 *  \todo Should each class of utility functions be in a separate source file
 *         instead of having a source file with the unhelpful name of Util?
 */

#ifndef ZOLTAN2_UTIL_HPP
#define ZOLTAN2_UTIL_HPP

#include <Zoltan2_Standards.hpp>
#include <Zoltan2_AlltoAll.hpp>
#include <Zoltan2_PartitioningSolution.hpp>

#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_DefaultComm.hpp>

namespace Zoltan2{

/*! \brief Convert an MPI communicator to a MpiComm object.
 */

#ifdef HAVE_ZOLTAN2_MPI
template <typename Ordinal>
  RCP<MpiComm<Ordinal> >
    getTeuchosMpiComm(const MPI_Comm &comm)
{
  RCP<Teuchos::OpaqueWrapper<MPI_Comm> >handle = Teuchos::opaqueWrapper<MPI_Comm>(comm);
  RCP<MpiComm<Ordinal> > tcommPtr(new MpiComm<Ordinal>(handle));

  return tcommPtr;
}
#endif

/*! \brief Convert an export part list to an import list.
 *
 * Given a PartitioningSolution return a list of all the global IDs 
 * that are mine.  
 *
 * If there are sizes associated with the IDs (like number of non-zeros)
 * include that in xtraInfo array.  Get back new sizes in newXtraInfo.
 * Otherwise xtraInfo.size() must be zero.
 *
 * \todo document params and return 
 */

template <typename User, typename Extra>
  size_t convertSolutionToImportList(
    const PartitioningSolution<User> &solution,
    ArrayRCP<Extra> &xtraInfo,
    ArrayRCP<typename InputTraits<User>::gid_t> &imports,
    ArrayRCP<Extra> &newXtraInfo)
{
  typedef typename InputTraits<User>::lno_t lno_t;
  typedef typename InputTraits<User>::gid_t gid_t;
  typedef Teuchos::Comm<int> comm_t;

  const RCP<const comm_t> &comm     = solution.getCommunicator();
  int numProcs                = comm->getSize();
  const gid_t *myGids         = solution.getIdList();
  size_t localNumIds          = solution.getLocalNumberOfIds();
  const partId_t *partList    = solution.getPartList();
  size_t numParts             = solution.getGlobalNumberOfParts();

  //
  // If procList is NULL, then part P goes to process P for all P.
  //
  const int *procList         = solution.getProcList();

  int localSend = ((xtraInfo.size() == localNumIds) ? 1 : 0);
  int globalSend =0;
  reduceAll<int, int>(*comm, Teuchos::REDUCE_SUM, 1, &localSend, &globalSend);

  bool sendExtra = (globalSend == numProcs);

  Array<lno_t> counts(numProcs, 0);
  if (procList)
    for (size_t i=0; i < localNumIds; i++)
      counts[procList[i]]++;
  else
    for (size_t i=0; i < localNumIds; i++)
      counts[partList[i]]++;

  Array<lno_t> offsets(numProcs+1, 0);
  for (size_t i=1; i <= numProcs; i++){
    offsets[i] = offsets[i-1] + counts[i-1];
  }

  Array<gid_t> gidList(localNumIds);
  Array<Extra> numericInfo;

  if (sendExtra)
    numericInfo.resize(localNumIds);

  if (procList){
    for (size_t i=0; i < localNumIds; i++){
      lno_t idx = offsets[procList[i]];
      gidList[idx] = myGids[i];
      if (sendExtra)
        numericInfo[idx] = xtraInfo[i];
      offsets[procList[i]] = idx + 1;
    }
  }
  else{
    for (size_t i=0; i < localNumIds; i++){
      lno_t idx = offsets[partList[i]];
      gidList[idx] = myGids[i];
      if (sendExtra)
        numericInfo[idx] = xtraInfo[i];
      offsets[partList[i]] = idx + 1;
    }
  }

  ArrayRCP<lno_t> recvCounts;
  RCP<const Environment> env = rcp(new Environment);

  try{
    AlltoAllv<gid_t, lno_t>(*comm, *env, gidList.view(0,localNumIds), 
      counts.view(0, numProcs), imports, recvCounts);
  }
  catch (std::exception &e){
    throw std::runtime_error("alltoallv 1");
  }

  if (sendExtra){
    try{
      AlltoAllv<Extra, lno_t>(*comm, *env, xtraInfo.view(0, localNumIds), 
        counts.view(0, numProcs), newXtraInfo, recvCounts);
    }
    catch (std::exception &e){
      throw std::runtime_error("alltoallv 2");
    }
  }

  return imports.size();
}

} // namespace Zoltan2

#endif
