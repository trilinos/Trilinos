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

template <typename scalar_t>
  inline bool outsideRegion(scalar_t val, scalar_t mark, double epsilon){
    return ((val < mark-epsilon) || (val > mark+epsilon));
}

#ifdef HAVE_ZOLTAN2_MPI

/*! \brief Convert an MPI communicator to a MpiComm object.
 */

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
 * If there are values of interest associated with the IDs (like number of non-zeros
 * or a global permutation) include that in xtraInfo array.  Get back new 
 * values in newXtraInfo. Otherwise xtraInfo.size() must be zero.
 *
 * Global Ids appear in the import list in process rank order.  So all ids to be
 * sent from rank 0 are listed first, then all ids to be sent by rank 1, and so on.
 *
 *  \param solution a partitioning solution
 *  \param xtraInfo an array of data cooresponding to the objects for which parts are
 *                         supplied the \c solution, or an empty array.
 *  \param imports on return lists the global Ids assigned to this process by the solution.
 *  \param newXtraInfo if \c xtraInfo was supplied, this array on return 
 *         has the extra information associated with the global Ids in the
 *         \c imports list.  Otherwise it is an empty list.
 *
 *  \todo  In the case where processes are not synomomous with parts, do we want to
 *              return part Id for each global Id as well?
 */

template <typename Adapter, typename Extra>
  size_t convertSolutionToImportList(
    const PartitioningSolution<Adapter> &solution,
    ArrayRCP<Extra> &xtraInfo,
    ArrayRCP<typename Adapter::gid_t> &imports,
    ArrayRCP<Extra> &newXtraInfo)
{
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gid_t gid_t;
  typedef Teuchos::Comm<int> comm_t;

  const RCP<const comm_t> &comm     = solution.getCommunicator();
  int numProcs                = comm->getSize();
  const gid_t *myGids         = solution.getIdList();
  size_t localNumIds          = solution.getLocalNumberOfIds();
  const partId_t *partList    = solution.getPartList();

  //
  // If procList is NULL, then part P goes to process P for all P.
  //
  const int *procList         = solution.getProcList();

  int localSend = ((size_t(xtraInfo.size()) == localNumIds) ? 1 : 0);
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
  for (int i=1; i <= numProcs; i++){
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
