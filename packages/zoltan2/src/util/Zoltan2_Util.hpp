// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

#ifndef ZOLTAN2_UTIL_HPP
#define ZOLTAN2_UTIL_HPP

#include <Zoltan2_Standards.hpp>
#include <Zoltan2_AlltoAll.hpp>

#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_DefaultComm.hpp>

namespace Zoltan2{

/*! Convert an MPI communicator to a MpiComm object.
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
 * Given a list of global IDs and their assigned parts, return
 * a list of all the global IDs that are mine.  Assumption is that
 * parts are 0 through nprocs-1, and process p gets part p.
 *
 * If there are sizes associated with the IDs (like number of non-zeros)
 * include that in xtraInfo array.  Get back new sizes in newXtraInfo.
 * Otherwise xtraInfo.size() must be zero.
 *
 * return the size of the import list
 */

template <typename GID, typename LNO, typename EXTRA>
  size_t convertPartListToImportList(
    const Teuchos::Comm<int> &comm,
    ArrayRCP<size_t> &part,
    ArrayRCP<GID> &gid,
    ArrayRCP<EXTRA> &xtraInfo,
    ArrayRCP<GID> &imports,
    ArrayRCP<EXTRA> &newXtraInfo)
{
  size_t numParts = comm.getSize();
  size_t localNumIds = gid.size();

  int localSend = (xtraInfo.size() == gid.size() ? 1 : 0);
  int globalSend =0;
  reduceAll<int, int>(comm, Teuchos::REDUCE_SUM, 1,
    &localSend, &globalSend);

  bool sendSizes = (globalSend == comm.getSize());

  Array<LNO> counts(numParts, 0);
  for (size_t i=0; i < localNumIds; i++){
    counts[part[i]]++;
  }

  Array<LNO> offsets(numParts+1, 0);
  for (size_t i=1; i <= numParts; i++){
    offsets[i] = offsets[i-1] + counts[i-1];
  }

  Array<GID> gidList(localNumIds);
  Array<EXTRA> numericInfo(localNumIds);

  for (size_t i=0; i < localNumIds; i++){
    LNO idx = offsets[part[i]];
    gidList[idx] = gid[i];
    if (sendSizes)
      numericInfo[idx] = xtraInfo[i];
    offsets[part[i]] = idx + 1;
  }

  ArrayRCP<LNO> recvCounts;
  RCP<const Environment> env = rcp(new Environment);

  try{
    AlltoAllv<GID, LNO>(comm, *env, gidList(), counts(), imports, recvCounts);
  }
  catch (std::exception &e){
    throw std::runtime_error("alltoallv 1");
  }

  if (sendSizes){
    try{
      AlltoAllv<EXTRA, LNO>(comm, *env, xtraInfo(), counts(), newXtraInfo, recvCounts);
    }
    catch (std::exception &e){
      throw std::runtime_error("alltoallv 2");
    }
  }

  return imports.size();
}

} // namespace Zoltan2

#endif
