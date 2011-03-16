// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

#ifndef _ZOLTAN2_IDENTIFIERMAP_DEF_HPP_
#define _ZOLTAN2_IDENTIFIERMAP_DEF_HPP_

/*! \file Zoltan2_IdentifierMap.hpp
*/

#include <stdexcept>
#include <vector>
#include <set>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Map.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_TestForException.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Zoltan2_IdentifierTraits.hpp>
#include <Zoltan2_Hash.hpp>
#include <Zoltan2_AlltoAll.hpp>

namespace Z2
{

template<typename AppGID, typename AppLID, typename GNO, typename LNO>
  IdentifierMap<AppGID,AppLID,GNO,LNO>::IdentifierMap(
    Teuchos::RCP<const Teuchos::Comm<int> > &in_comm, 
    Teuchos::RCP<std::vector<AppGID> > &gids, 
    Teuchos::RCP<std::vector<AppLID> > &lids) 
         : comm(in_comm), myGids(gids), myLids(lids),
           consecutive(false), base(0), 
           totalGlobalIds(0), totalLocalIds(0), localNumberOfIds(0){

  int nprocs = comm->getSize(); 
  GNO numIds[2] = {static_cast<GNO>(gids->size()), static_cast<GNO>(lids->size())};
  size_t nvals = 2 * (nprocs+1);

  Teuchos::ArrayRCP<GNO> totals = Teuchos::arcp(new GNO [nvals], 0, nvals);

  totals[0] = totals[1] = GNO(0);

  Teuchos::scan(*comm, Teuchos::REDUCE_SUM, int(2), numIds, &(totals[2]));

  GNO totalGlobalIds = totals[nvals-2];
  GNO totalLocalIds = totals[nvals-1];
  LNO localNumberOfIds = static_cast<LNO>(gids->size());

  TEST_FOR_EXCEPTION((totalLocalIds > 0) && (totalLocalIds != totalGlobalIds),
                      std::runtime_error,
                     "local IDs are provided but number of global IDs"
                     " does not equal number of local IDs");

  if (lids->size() > 0){
    lidHash = Teuchos::rcp(new Teuchos::Hashtable<std::string, LNO>(localNumberOfIds));

    for (LNO i=0; i < localNumberOfIds; i++){
      lidHash->put(Z2::IdentifierTraits<AppLID>::key((*lids)[i]), i);
    }
  }

#ifdef APPGID_IS_NOT_GNO

  consecutive = true;
  base = 0;

  gnoDist.reserve(nprocs+1);
  gnoDist[0] = 0;
  for (LNO i=1; i < nprocs+1; i++){
    gnoDist[i] = totals[i*2];
  }

  gidHash = Teuchos::rcp(new Teuchos::Hashtable<std::string, LNO>(localNumberOfIds));

  for (LNO i=0; i < localNumberOfIds; i++){
    gidHash->put(Z2::IdentifierTraits<AppGID>::key((*gids)[i]), i);
  }

#else

  GNO min(0), max(0), globalMin(0), globalMax(0);
  min = max = static_cast<GNO>((*gids)[0]);

  for (LNO i=1; i < localNumberOfIds; i++){
    if ((*gids)[i] < min)
      min = static_cast<GNO>((*gids)[i]);
    else if ((*gids)[i] > max)
      max = static_cast<GNO>((*gids)[i]);
  }

  Teuchos::reduceAll(comm, Teuchos::REDUCE_MIN, int(1), &min, &globalMin);
  Teuchos::reduceAll(comm, Teuchos::REDUCE_MAX, int(1), &max, &globalMax);

  if (globalMax - globalMin + 1 == globalNumGids){
    consecutive=true;
    base = globalMin;
  }

  if (totalLocalIds > 0){
    gidHash = Teuchos::rcp(new Teuchos::Hashtable<AppGID, LNO>(localNumberOfIds));

    for (LNO i=0; i < localNumberOfIds; i++){
      gidHash->put((*gids)[i], i);
    }
  }
#endif
}

  /*! Constructor */
template<typename AppGID, typename AppLID, typename GNO, typename LNO>
  IdentifierMap<AppGID,AppLID,GNO,LNO>::IdentifierMap()  
         : comm(), myGids(), myLids(),
           consecutive(false), base(0), 
           totalGlobalIds(0), totalLocalIds(0), localNumberOfIds(0)
{
}

  /*! Destructor */
template<typename AppGID, typename AppLID, typename GNO, typename LNO>
  IdentifierMap<AppGID,AppLID,GNO,LNO>::~IdentifierMap() 
  {
  }

  /*! Copy Constructor */
template<typename AppGID, typename AppLID, typename GNO, typename LNO>
  IdentifierMap<AppGID,AppLID,GNO,LNO>::IdentifierMap(const IdentifierMap &id)
{
    //TODO
}

  /*! Assignment operator */
template<typename AppGID, typename AppLID, typename GNO, typename LNO>
  IdentifierMap<AppGID,AppLID,GNO,LNO> &
    IdentifierMap<AppGID,AppLID,GNO,LNO>::operator=(const IdentifierMap<AppGID,AppLID,GNO,LNO> &id)
{
    //TODO
}

template<typename AppGID, typename AppLID, typename GNO, typename LNO>
  void IdentifierMap<AppGID,AppLID,GNO,LNO>::initialize(
    Teuchos::RCP<const Teuchos::Comm<int> > &in_comm, 
    Teuchos::RCP<std::vector<AppGID> > &gids, 
    Teuchos::RCP<std::vector<AppLID> > &lids) 
{
    reset();
    this->IdentifierMap(in_comm, gids, lids);
}

template<typename AppGID, typename AppLID, typename GNO, typename LNO>
  void IdentifierMap<AppGID, AppLID, GNO, LNO>::reset()
{
  if (!comm.is_null())
    comm.release();

  if (!myGids.is_null())
    myGids.release();

  if (!myLids.is_null())
    myLids.release();

  if (!gidHash.is_null())
    gidHash.release();

  if (!lidHash.is_null())
    lidHash.release();

  consecutive =false;
  base = totalGlobalIds = totalLocalIds = GNO(0);
  localNumberOfIds = LNO(0);
}

template<typename AppGID, typename AppLID, typename GNO, typename LNO>
  bool IdentifierMap<AppGID, AppLID, GNO, LNO>::gnosAreGids()
{
#ifdef APPGID_IS_NOT_GNO
  return false;
#else
  return true;
#endif
}

template<typename AppGID, typename AppLID, typename GNO, typename LNO>
  void IdentifierMap<AppGID, AppLID, GNO, LNO>::gidTranslate(
    std::vector<AppGID> &gid, std::vector<GNO> &gno)
{
  if ((gid.size() == 0) && (gno.size() == 0)){
    return;
  }

  bool getGids = (gno.size() > 0);

  if (getGids)
    gid.reserve(gno.size());
  else
    gno.reserve(gid.size());

#ifdef APPGID_IS_NOT_GNO

  GNO gnoStart = gnoDist[comm->getRank()];

  if (getGids){
    for (LNO i=0; i < gno.size(); i++){
      gid[i] = myGids[gno[i] - gnoStart];
    }
  }
  else{
    for (LNO i=0; i < gid.size(); i++){
      const LNO &idx = gidHash->get(Z2::IdentifierTraits<AppGID>::key(gid[i]));
      gno[i] = gnoStart + idx;
    }
  }

#else

  if (getGids)
    gid = gno;
  else
    gno = gid;

#endif
}

template<typename AppGID, typename AppLID, typename GNO, typename LNO>
  void IdentifierMap<AppGID, AppLID, GNO, LNO>::lidTranslate(
    std::vector<AppLID> &lid, std::vector<GNO> &gno)
{
  if ((lid.size() == 0) && (gno.size() == 0)){
    return;
  }

  bool getLids = (gno.size() > 0);

  if (getLids)
    lid.reserve(gno.size());
  else
    gno.reserve(lid.size());

  TEST_FOR_EXCEPTION(myLids.size() == 0,
                      std::runtime_error,
                     "local ID translation is requested but none were provided");

#ifdef APPGID_IS_NOT_GNO

  GNO gnoStart = gnoDist[comm->getRank()];

  if (getLids){
    for (LNO i=0; i < gno.size(); i++){
      lid[i] = myLids[gno[i] - gnoStart];
    }
  }
  else{
    for (LNO i=0; i < lid.size(); i++){
      const LNO &idx = lidHash->get(Z2::IdentifierTraits<AppLID>::key(lid[i]));
      gno[i] = gnoStart + idx;
    }
  }

#else

  if (getLids){
    for (LNO i=0; i < gno.size(); i++)
      const LNO &idx = gidHash->get(gno[i]);
      lid[i] = myLids[idx];
  }
  else{
    for (LNO i=0; i < lid.size(); i++){
      const LNO &idx = lidHash->get(Z2::IdentifierTraits<AppLID>::key(lid[i]));
      gno[i] = myGids[idx];
    }
  }
#endif
}

template<typename AppGID, typename AppLID, typename GNO, typename LNO>
  void IdentifierMap<AppGID, AppLID, GNO, LNO>::gidGlobalTranslate(
    std::vector<AppGID> &gid,                                 /* input  */
    std::vector<GNO> &gno, std::vector<int> &proc)            /* output */
{
  // Problem: the following does not compile:
  //    std::map<std::string, std::vector<LNO> >::iterator
  // but this does:
  //    std::map<std::string, std::vector<int> >::iterator

  TEST_FOR_EXCEPTION ( sizeof(LNO) > sizeof(int),
                      std::runtime_error,
                     "int is not sufficient for vector type");

#ifdef APPGID_IS_NOT_GNO
  // AppGID doesn't implement "<" and ">", so we use std::string for key
  std::map<std::string, std::vector<int> > remoteAppGIDs; 
  std::map<std::string, std::vector<int> >::iterator next, last;
  int rank = comm->getRank();
  int nprocs = comm->getSize();
#else
  std::map<AppGID, std::vector<int> > remoteAppGIDs;
  std::map<AppGID, std::vector<int> >::iterator next, last;
#endif

  gno.reserve(gid.size());
  proc.reserve(gid.size());

  for (int i=0; i < gid.size(); i++){

#ifdef APPGID_IS_NOT_GNO
    std::string key(Z2::IdentifierTraits<AppGID>::key(gid[i]));
#else
    AppGID key(gid[i]);
#endif

    if (!gidHash.containsKey(key)){

      proc[i] = -1;       // not ours

      next = remoteAppGIDs.find(key);
      last = remoteAppGIDs.end();

      if (next == last){
        std::vector<int> indices(10);
        indices.push_back(i);
        remoteAppGIDs[key] = indices;
      }
      else{
        std::vector<int> &indices = next->second;
        indices.push_back(i);
      }
    }
    else{

#ifdef APPGID_IS_NOT_GNO
      gno[i] = gnoDist[rank] + gidHash.get(key);
#else
      gno[i] = GNO(gid[i]);
#endif
      proc[i] = rank;

    }
  }

  LNO queryCount = remoteAppGIDs.size();

  GNO localCount(queryCount), globalCount(0);

  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, int(1), &localCount, &globalCount);

  if (globalCount == 0)
    return;


#ifdef APPGID_IS_NOT_GNO

  // Hash all my GID strings to processes that will answer for them

  std::vector<std::string> gidStrings(localNumberOfIds);
  std::vector<int> answerProcs(localNumberOfIds);

  for (LNO i=0; i < (*myGids).size(); i++){
    gidStrings.push_back(Z2::IdentifierTraits<AppGID>::key((*myGids)[i]));
    answerProcs.push_back(static_cast<int>(Zoltan2_Hash(gidStrings[i], nprocs)));
  }

  Teuchos::ArrayView<std::string> gidArray(gidStrings);

  // Export my GID strings to the processes that will answer for them.  

  // Create a search structure for the GIDs I answer for.

  // Request information about my remote GIDs from the processes that answer for them.

#else

  Teuchos::ArrayView<AppGID> gidArray((*myGids));
  Tpetra::Map<LNO, AppGID> map(totalGlobalIds, gidArray, base, comm);

  std::vector<AppGID> queryGids(queryCount);
  std::vector<int> procIds(queryCount);

  next = remoteAppGIDs.begin();
  last = remoteAppGIDs.end();

  while (next != last){
    queryGids.push_back(next->first);
    next++;
  }

  Teuchos::ArrayView<AppGID> queryGidsArray(queryGids);
  Teuchos::ArrayView<int> procIdsArray(procIds);

  Tpetra::LookupStatus status = map.getRemoteIndexList(queryGidsArray, procIdsArray);

  TEST_FOR_EXCEPTION( status != Tpetra::AllIDsPresent,
                      std::runtime_error,
                     "invalid remote application global ID in lookup call");

  next = remoteAppGIDs.begin();
  last = remoteAppGIDs.end();

  std::vector<int>::iterator nextProc = procIds.begin();

  while (next != last){
    int gidOwner = *nextProc++;
    AppGID &gno = next->first;
    std::vector<int> &indices = next->second;
    for (int i=0; i < indices.size(); i++){
      gno[indices[i]] = gno;
      proc[indices[i]] = gidOwner;
    }
    next++;
  }
#endif


  
}





}   // end namespace Z2

#endif /* _ZOLTAN2_IDENTIFIERMAP_DEF_HPP_ */
