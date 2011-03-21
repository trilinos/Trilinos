//@HEADER
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
           globalNumberOfIds(0), localNumberOfIds(0), haveLocalIds(false)
{

  int nprocs = comm->getSize(); 

  localNumberOfIds = gids.get()->size();

  GNO *val = new GNO [4];
  val[0] = lids.get()->size();
  val[1] = localNumberOfIds;
  val[2] = 0;
  val[3] = 0;

  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 2, val+0, val+2);

  globalNumberOfIds = val[3];
  haveLocalIds = (val[2] > 0);

 TEST_FOR_EXCEPTION(haveLocalIds && 
             (static_cast<LNO>(lids.get()->size()) != localNumberOfIds),
             std::runtime_error,
          "local IDs are provided but number of global IDs"
          " does not equal number of local IDs");

  if (haveLocalIds){
    lidHash = Teuchos::rcp(new Teuchos::Hashtable<std::string, LNO>
                (localNumberOfIds));

    for (LNO i=0; i < localNumberOfIds; i++){
      lidHash->put(Z2::IdentifierTraits<AppLID>::key((*lids)[i]), i);
    }
  }

#ifdef APPGID_IS_NOT_GNO

  consecutive = true;
  base = 0;

  gnoDist = Teuchos::ArrayRCP<GNO>(nprocs + 1, 0);
  GNO *distArray = gnoDist.get();

  distArray[nprocs] = localNumberOfIds;

  Teuchos::gatherAll(*comm, 1, distArray + nprocs, 1, distArray);

  gnoDist.get()[nprocs] = globalNumberOfIds;

  for (int p=nprocs-1; p >=0; p--){
    gnoDist.get()[p] = gnoDist.get()[p+1] - gnoDist.get()[p];
  }

  gidHash = Teuchos::rcp(new Teuchos::Hashtable<std::string, LNO>
              (localNumberOfIds));

  for (LNO i=0; i < localNumberOfIds; i++){
    gidHash->put(Z2::IdentifierTraits<AppGID>::key((*gids)[i]), i);
  }

#else

  GNO min(0), max(0), globalMin(0), globalMax(0);
  std::vector<AppGID> &ids = gids.get();

  min = max = static_cast<GNO>(ids[0]);

  for (LNO i=1; i < localNumberOfIds; i++){
    if (ids[i] < min)
      min = static_cast<GNO>(ids[i]);
    else if (ids[i] > max)
      max = static_cast<GNO>(ids[i]);
  }

  Teuchos::reduceAll(comm, Teuchos::REDUCE_MIN, int(1), &min, &globalMin);
  Teuchos::reduceAll(comm, Teuchos::REDUCE_MAX, int(1), &max, &globalMax);

  if (globalMax - globalMin + 1 == globalNumberOfIds){
    consecutive=true;
  }
  else{
    consecutive=false;
  }

  base = globalMin;

  gidHash = Teuchos::rcp(new Teuchos::Hashtable<AppGID, LNO>(localNumberOfIds));

  for (LNO i=0; i < localNumberOfIds; i++){
    gidHash->put((*gids)[i], i);
  }
#endif
}

  /*! Constructor */
template<typename AppGID, typename AppLID, typename GNO, typename LNO>
  IdentifierMap<AppGID,AppLID,GNO,LNO>::IdentifierMap()  
         : comm(), myGids(), myLids(),
           consecutive(false), base(0), 
           globalNumberOfIds(0), localNumberOfIds(0), haveLocalIds(false)
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
    IdentifierMap<AppGID,AppLID,GNO,LNO>::operator=(const IdentifierMap<AppGID,
                  AppLID,GNO,LNO> &id)
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
  base = globalNumberOfIds = GNO(0);
  localNumberOfIds = LNO(0);
  haveLocalIds = false;
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
    std::vector<AppGID> &gid,          // input or output, take your pick
    std::vector<GNO> &gno)             // output or input
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

  GNO firstGNO = gnoDist[comm->getRank()];

  if (getGids){
    for (LNO i=0; i < gno.size(); i++){
      gid.push_back(myGids.get()[gno[i] - firstGNO]);
    }
  }
  else{
    for (LNO i=0; i < gid.size(); i++){
      const LNO &idx = gidHash->get(Z2::IdentifierTraits<AppGID>::key(gid[i]));
      gno.push_back(firstGNO + idx);
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
    std::vector<AppLID> &lid,          // input or output, take your pick
    std::vector<GNO> &gno)             // output or input
{
  if ((lid.size() == 0) && (gno.size() == 0)){
    return;
  }

  TEST_FOR_EXCEPTION(!haveLocalIds,
                  std::runtime_error,
                 "local ID translation is requested but none were provided");

  bool getLids = (gno.size() > 0);

  if (getLids)
    lid.reserve(gno.size());
  else
    gno.reserve(lid.size());


#ifdef APPGID_IS_NOT_GNO

  GNO firstGNO = gnoDist[comm->getRank()];

  if (getLids){
    for (LNO i=0; i < gno.size(); i++){
      lid.push_back(myLids.get()[gno[i] - firstGNO]);
    }
  }
  else{
    for (LNO i=0; i < lid.size(); i++){
      const LNO &idx = lidHash->get(Z2::IdentifierTraits<AppLID>::key(lid[i]));
      gno.push_back(firstGNO + idx);
    }
  }

#else

  if (getLids){
    for (LNO i=0; i < gno.size(); i++)
      const LNO &idx = gidHash->get(gno[i]);
      lid.push_back(myLids.get()[idx]);
  }
  else{
    for (LNO i=0; i < lid.size(); i++){
      const LNO &idx = lidHash->get(Z2::IdentifierTraits<AppLID>::key(lid[i]));
      gno.push_back(myGids.get()[idx]);
    }
  }
#endif
}

template<typename AppGID, typename AppLID, typename GNO, typename LNO>
  void IdentifierMap<AppGID, AppLID, GNO, LNO>::gidGlobalTranslate(
    std::vector<AppGID> &gid,    /* input, gids may be repeated in list  */
    std::vector<GNO> &gno,       /* output */ 
    std::vector<int> &proc)      /* output */
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
  GNO firstGNO = gnoDist[rank];
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
      gno[i] = firstGNO + gidHash.get(key);
#else
      gno[i] = GNO(gid[i]);
#endif
      proc[i] = rank;
    }
  }

  LNO queryCount = remoteAppGIDs.size();

  GNO localCount(queryCount), globalCount(0);

  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &localCount, &globalCount);

  if (globalCount == 0)
    return;

#ifdef APPGID_IS_NOT_GNO

  // Hash all my GID strings to processes that will answer for them

  std::vector<AppGID> &ids = myGids.get();
  std::vector<std::string> gidStrings(localNumberOfIds);
  std::vector<int> answerProcs(localNumberOfIds);
  std::vector<int> msgCount(nprocs,0);

  for (LNO i=0; i < ids.size(); i++){
    std::string key(Z2::IdentifierTraits<AppGID>::key(ids[i]));
    int p = Z2::Hash<AppGID, int>(key, nprocs);
    answerProcs.push_back(p);
    msgCount[p]++;
  }

  // Send my GIDs, plus their local index, 
  // to the processes that will answer for them

  AppGID *sendGidBuf = new AppGID [localNumberOfIds];
  AppGID *sendIndexBuf = new LNO [localNumberOfIds];
  std::vector<int> sendGidOffsets(nprocs+1,0);
  for (int i=0; i < nprocs; i++){
    sendGidOffsets[i+1] = sendGidOffsets[i] + msgCount[i];
  }

  for (LNO i=0; i < ids.size(); i++){
    int p = answerProcs[i];
    sendGidBuf[sendGidOffsets[p]] = ids[i];
    sendIndexBuf[sendGidOffsets[p]] = i;
    sendGidOffsets[p]++;
  }

#else

  Teuchos::ArrayView<AppGID> gidArray((*myGids));
  Tpetra::Map<LNO, AppGID> map(globalNumberOfIds, gidArray, base, comm);

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

  Tpetra::LookupStatus status = map.getRemoteIndexList(queryGidsArray, 
   procIdsArray);

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
