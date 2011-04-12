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
#include <map>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Map.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_TestForException.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Zoltan2_IdentifierTraits.hpp>
#include <Zoltan2_Util.hpp>
#include <Zoltan2_AlltoAll.hpp>

namespace Z2
{

template<typename AppLID, typename AppGID, typename LNO, typename GNO> 
  IdentifierMap<AppLID,AppGID,LNO,GNO>::IdentifierMap(
    Teuchos::RCP<const Teuchos::Comm<int> > &in_comm, 
    Teuchos::RCP<std::vector<AppGID> > &gids, 
    Teuchos::RCP<std::vector<AppLID> > &lids) 
         : _comm(in_comm), _myGids(gids), _myLids(lids),
           _globalNumberOfIds(0), _localNumberOfIds(0), _haveLocalIds(false),
           _myRank(0), _numProcs(0)
{
  _numProcs = _comm->getSize(); 
  _myRank = _comm->getRank(); 

  std::vector<AppGID> &ids = *_myGids;

  _localNumberOfIds = ids.size();

  GNO *val = new GNO [4];
  val[0] = (*lids).size();
  val[1] = _localNumberOfIds;
  val[2] = val[3] = 0;

  Teuchos::reduceAll(*_comm, Teuchos::REDUCE_SUM, 2, val+0, val+2);

  _haveLocalIds = (val[2] > 0);
  _globalNumberOfIds = val[3];

  TEST_FOR_EXCEPTION(_haveLocalIds && (val[0] != _localNumberOfIds),
             std::runtime_error,
          "local IDs are provided but number of global IDs"
          " does not equal number of local IDs");

  _gnoDist=Teuchos::null;
  _gidHash=Teuchos::null;
  _lidHash=Teuchos::null;

  if (_haveLocalIds){   // hash LID to index in LID vector
    _lidHash = Teuchos::rcp(
      new Teuchos::Hashtable<std::string, LNO> (_localNumberOfIds));

    for (LNO i=0; i < _localNumberOfIds; i++){
      _lidHash->put(Z2::IdentifierTraits<AppLID>::key((*lids)[i]), i);
    }
  }

  // If the application's global ID data type (AppGID) is a Teuchos Ordinal,
  // we will be able to use it as our internal global numbers (GNO).  
  // Otherwise we will have to map their global IDs to valid global numbers.

  if (IdentifierTraits<AppGID>::isGlobalOrdinalType()){

    // Are the AppGIDs consecutive and increasing with process rank? 
    // If so GID/proc lookups are optimized.

    GNO *startGID = NULL;
    GNO min(0), max(0), globalMin(0), globalMax(0);
  
    min = max = static_cast<GNO>(ids[0]);
    GNO checkVal = min;
    bool consecutive = true;
  
    for (LNO i=1; i < _localNumberOfIds; i++){
      if (consecutive && (ids[i] != ++checkVal)){
        consecutive=false;
        break;
      }
      if (ids[i] < min)
        min = static_cast<GNO>(ids[i]);
      else if (ids[i] > max)
        max = static_cast<GNO>(ids[i]);
    }
  
    val[0] = (consecutive ? 1 : 0);
    val[1] = min;
    val[2] = val[3] = 0;

    Teuchos::reduceAll<int, GNO>(*_comm, Teuchos::REDUCE_MIN, 2, val+0, val+2);

    if (val[2] != 1)       // min of consecutive flags
      consecutive=false;

    if (consecutive){
      globalMin = val[3];
      Teuchos::reduceAll<int, GNO>(*_comm, Teuchos::REDUCE_MAX, 1, &max, 
        &globalMax);
      if (globalMax - globalMin + 1 != _globalNumberOfIds)
        consecutive = false;   // there are gaps in the gids
  
      if (consecutive){
        GNO myStart = ids[0];
        startGID = new GNO [_numProcs+1];
      
        Teuchos::gatherAll<int, GNO>(*_comm, 1, &myStart, _numProcs, startGID);
      
        for (int p=1; p < _numProcs; p++){
          if (startGID[p] < startGID[p-1]){
            consecutive = false;  // gids do not increase with process rank
            break;
          }
        }
        if (consecutive){
          startGID[_numProcs] = globalMax + 1;
          _gnoDist = Teuchos::ArrayRCP<GNO>(startGID, 0, _numProcs+1, true);
        }
        else{
          delete [] startGID;
          startGID = NULL;
        }
      }
    }
  }
  else{

    // AppGIDs are not Ordinals.  We map them to consecutive 
    // global numbers starting with 0. 

    GNO *numScan = new GNO [_numProcs+1];
    GNO myNum = _localNumberOfIds;

    Teuchos::scan<int, GNO>(*_comm, Teuchos::REDUCE_SUM, 1, &myNum, 
      numScan + 1);

    numScan[0] = 0;
    _gnoDist = Teuchos::ArrayRCP<GNO>(numScan, 0, _numProcs + 1, true);
  }

  if (_gnoDist.is_null()){

    // We need a hash table mapping the application global ID
    // to its index in _myGids.

    _gidHash = Teuchos::rcp(new Teuchos::Hashtable<std::string, LNO>
              (_localNumberOfIds));

    for (LNO i=0; i < _localNumberOfIds; i++){
      _gidHash->put(IdentifierTraits<AppGID>::key(ids[i]), i);
    }
  }
}

  /*! Constructor */
template<typename AppLID, typename AppGID, typename LNO, typename GNO>
  IdentifierMap<AppLID,AppGID,LNO,GNO>::IdentifierMap()  
         : _comm(), _myGids(), _myLids(),
           _globalNumberOfIds(0), _localNumberOfIds(0), _haveLocalIds(false),
           _myRank(0), _numProcs(0)
{
}

  /*! Destructor */
template<typename AppLID, typename AppGID, typename LNO, typename GNO>
  IdentifierMap<AppLID,AppGID,LNO,GNO>::~IdentifierMap() 
  {
  }

  /*! Copy Constructor */
template<typename AppLID, typename AppGID, typename LNO, typename GNO>
  IdentifierMap<AppLID,AppGID,LNO,GNO>::IdentifierMap(const IdentifierMap &id)
{
    //TODO
}

  /*! Assignment operator */
template<typename AppLID, typename AppGID, typename LNO, typename GNO>
  IdentifierMap<AppLID,AppGID,LNO,GNO> &
    IdentifierMap<AppLID,AppGID,LNO,GNO>::operator=(const IdentifierMap<AppLID,
                  AppGID,LNO,GNO> &id)
{
    //TODO
}

template<typename AppLID, typename AppGID, typename LNO, typename GNO>
  void IdentifierMap<AppLID,AppGID,LNO,GNO>::initialize(
    Teuchos::RCP<const Teuchos::Comm<int> > &in_comm, 
    Teuchos::RCP<std::vector<AppGID> > &gids, 
    Teuchos::RCP<std::vector<AppLID> > &lids) 
{
    this->IdentifierMap(in_comm, gids, lids);
}

template<typename AppLID, typename AppGID, typename LNO, typename GNO>
  bool IdentifierMap<AppLID, AppGID, LNO, GNO>::gnosAreGids()
{
  return IdentifierTraits<AppGID>::isGlobalOrdinalType();
}

template<typename AppLID, typename AppGID, typename LNO, typename GNO>
  void IdentifierMap<AppLID, AppGID, LNO, GNO>::gidTranslate(
    std::vector<AppGID> &gid,          // input or output, take your pick
    std::vector<GNO> &gno)             // output or input
{
  if ((gid.size() == 0) && (gno.size() == 0)){
    return;
  }

  bool getGids;
  int len=0;

  if (gno.size() > 0){
    getGids = true;
    len = gno.size();
    gid.erase();
    gid.reserve(len);
  }
  else{
    getGids = false;
    len = gid.size();
    gno.erase();
    gno.reserve(len);
  }

  if (IdentifierTraits<AppGID>::isGlobalOrdinalType()){   
    // our gnos are the app gids
    if (getGids)
      for (LNO i=0; i < len; i++)
        gid.push_back(static_cast<AppGID>(gno[i]));
    else
      for (LNO i=0; i < len; i++)
        gno.push_back(static_cast<GNO>(gid[i]));
  }
  else{              // we mapped gids to consecutive gnos
  
    GNO firstGNO = _gnoDist[_myRank];
    if (getGids)
      for (LNO i=0; i < len; i++)
        gid.push_back((*_myGids)[gno[i] - firstGNO]);
    else
      for (LNO i=0; i < len; i++){
        const LNO &idx = _gidHash->get(
          Z2::IdentifierTraits<AppGID>::key(gid[i]));
        gno.push_back(firstGNO + idx);
      }
  }
  return;
}

template<typename AppLID, typename AppGID, typename LNO, typename GNO>
  void IdentifierMap<AppLID, AppGID, LNO, GNO>::lidTranslate(
    std::vector<AppLID> &lid,          // input or output, take your pick
    std::vector<GNO> &gno)             // output or input
{
  TEST_FOR_EXCEPTION(!_haveLocalIds,
                  std::runtime_error,
                 "local ID translation is requested but none were provided");

  if ((lid.size() == 0) && (gno.size() == 0)){
    return;
  }

  bool getLids;
  int len=0;

  if (gno.size() > 0){
    getLids = true;
    len = gno.size();
    lid.erase();
    lid.reserve(len);
  }
  else{
    getLids = false;
    len = lid.size();
    gno.erase();
    gno.reserve(len);
  }

  GNO firstGNO = - 1;
  if (!_gnoDist.is_null())
    firstGNO = _gnoDist[_myRank];
  

  if (getLids){
    for (LNO i=0; i < len; i++){
      int idx = 0;
      if (firstGNO >= 0)       // gnos are consecutive
        idx = gno[i] - firstGNO;
      else                     // gnos must be the app gids
        const LNO &idx = _gidHash->get(
          Z2::IdentifierTraits<AppGID>::key(static_cast<AppGID>(gno[i])));
      
      lid.push_back((*_myLids)[idx]);
    }
  }
  else{
    for (LNO i=0; i < len; i++){
      const LNO &idx = _lidHash->get(
        Z2::IdentifierTraits<AppLID>::key(lid[i]));

      if (firstGNO >= 0)       // gnos are consecutive
        gno.push_back(static_cast<GNO>(idx + firstGNO));
      else                     // gnos must be the app gids
        gno.push_back(static_cast<GNO>((*_myGids)[idx]));
    }
  }
}

template<typename AppLID, typename AppGID, typename LNO, typename GNO>
  void IdentifierMap<AppLID, AppGID, LNO, GNO>::gidGlobalTranslate(
    const std::vector<AppGID> &in_gid,
    std::vector<GNO> &out_gno,
    std::vector<int> &out_proc)
{
  LNO numGids = in_gid.size();

  out_gno.clear();
  out_gno.reserve(numGids);
  out_proc.clear();
  out_proc.reserve(numGids);

  if (IdentifierTraits<AppGID>::isGlobalOrdinalType() && !_gnoDist.is_null()){

    // Easy case - communication is not needed.
    // Global numbers are the application global IDs and
    // they are increasing consecutively with rank.
 
    typename std::map<GNO, int> firstGnoToProc;
    typename std::map<GNO, int>::iterator pos;

    for (int p=0; p <= _numProcs; p++){
      firstGnoToProc[_gnoDist[p]] = p;
    }

    for (LNO i=0; i < numGids; i++){
      GNO globalNumber = static_cast<GNO>(in_gid[i]);
      out_gno[i] = globalNumber;
      pos = firstGnoToProc.upper_bound(globalNumber);
      out_proc[i] = pos->first - 1;
    }

    return;
  }

  AppGID *gidOutBuf = NULL;
  GNO *gnoOutBuf = NULL;
  int *procOutBuf = NULL;
  LNO *countOutBuf = NULL;

  Teuchos::ArrayRCP<AppGID> gidInBuf();
  Teuchos::ArrayRCP<GNO> gnoInBuf();
  Teuchos::ArrayRCP<int> procInBuf();
  Teuchos::ArrayRCP<LNO> countInBuf();

  LNO *offsetBuf = NULL;

  bool needGnoInfo = !IdentifierTraits<AppGID>::isGlobalOrdinalType();

  ///////////////////////////////////////////////////////////////////////
  // First: Hash each of my AppGIDs to a process that will answer
  // for it.  Send my AppGIDs (and the Gnos if they are different)
  // to their assigned processes.  Build a search structure for
  // the AppGIDs that were assigned to me, so I can reply with
  // with the process owning them (and their Gnos if they are different).
  ///////////////////////////////////////////////////////////////////////

  int *hashProc = NULL;
  countOutBuf = new LNO [_numProcs];
  offsetBuf = new int [_numProcs+1];

  if (_localNumberOfIds > 0){

    hashProc = new int [_localNumberOfIds];
  
    for (int p=0; p < _numProcs; p++)
      countOutBuf[p] = 0;
  
    for (LNO i=0; i < _localNumberOfIds; i++){
      hashProc[i] = 
        IdentifierTraits<AppGID>::hashCode((*_myGids)[i]) % _numProcs;
      countOutBuf[hashProc[i]]++;
    }
  
    offsetBuf[0] = 0;
  
    for (int p=1; p <= _numProcs; p++){
      offsetBuf[p] = offsetBuf[p-1] + countOutBuf[p-1];
    }
  
    gidOutBuf = new AppGID [_localNumberOfIds];
  
    if (needGnoInfo){            // the gnos are not the gids
      gnoOutBuf = new GNO [_localNumberOfIds];
    }
  
    for (LNO i=0; i < _localNumberOfIds; i++){
      LNO offset = offsetBuf[hashProc[i]];
      gidOutBuf[offset] = (*_myGids)[i];
      if (needGnoInfo)
        gnoOutBuf[offset] = _gnoDist[_myRank] + i;
      offsetBuf[hashProc[i]] = offset + 1;
    }
  
    delete [] hashProc;
    hashProc = NULL;
  }
  else{
    for (int p=0; p < _numProcs; p++)
      countOutBuf[p] = 0;
  }

  Teuchos::ArrayView<LNO> countOutView(countOutBuf, _numProcs);
  Teuchos::ArrayView<AppGID> gidOutView(Teuchos::null);
  Teuchos::ArrayView<GNO> gnoOutView(Teuchos::null);

  if (_localNumberOfIds){
    gidOutView = Teuchos::ArrayView<AppGID>(gidOutBuf, _localNumberOfIds);
    if (needGnoInfo)
      gnoOutView = Teuchos::ArrayView<GNO>(gnoOutBuf, _localNumberOfIds);
  }

  AlltoAllv(*_comm, gidOutView, countOutView, gidInBuf, countInBuf);

  if (gidOutBuf != NULL){
    delete [] gidOutBuf;
    gidOutBuf = NULL;
  }

  if (needGnoInfo){
    countInBuf.release();
    AlltoAllv(*_comm, gnoOutView, countOutView, gnoInBuf, countInBuf);
    if (gnoOutBuf){
      delete [] gnoOutBuf;
      gnoOutBuf = NULL;
    }
  }

  //
  // Save the information that was hashed to me so I can do lookups.
  //

  std::map<LNO, int> firstIndexToProc;
  LNO total = 0;

  for (int p=0; p < _numProcs; p++){
    firstIndexToProc[total] = p;
    total += countInBuf[p];
  }

  firstIndexToProc[total] = _numProcs;

  Teuchos::Hashtable<std::string, LNO> *gidToAnswerIndex = 
    new Teuchos::Hashtable<std::string, LNO>(total);

  total = 0;
  for (int p=0; p < _numProcs; p++){
    for (LNO i=countInBuf[p]; i < countInBuf[p+1]; i++, total++){
      gidToAnswerIndex->put( 
        IdentifierTraits<AppGID>::key(gidInBuf[total]), total);
    }
  }

  // Keep gnoInBuf.  We're done with the others.
  gidInBuf.release();
  countInBuf.release();

  ///////////////////////////////////////////////////////////////////////
  // Send a request for information to the "answer process" for each 
  // of the GIDs in in_gid.  First remove duplicate GIDs from the list.
  ///////////////////////////////////////////////////////////////////////
  
  Teuchos::Array<std::string> uniqueGidQueries;
  Teuchos::Array<std::vector<int> > uniqueGidQueryIndices;
  LNO numberOfUniqueGids = 0;
  LNO *gidLocation = NULL;

  if (numGids > 0){
    LNO tableSize = numGids / 2;
    if (!tableSize) tableSize = 1;
  
    // In an input adapter, how many objects will have the same neighbor?
    LNO sizeChunk = 4;
  
    Teuchos::Hashtable<std::string, std::vector<LNO> > *gidIndices = 
      new Teuchos::Hashtable<std::string, std::vector<LNO> >(tableSize);
  
    for (LNO i=0; i < numGids; i++){
      std::string uniqueKey(IdentifierTraits<AppGID>::key(in_gid[i]));
      if (gidIndices->containsKey(uniqueKey)){
        std::vector<LNO> &v = gidIndices->get(uniqueKey);
        if (v.size() % sizeChunk == 0){
          v.reserve(v.size() + sizeChunk);
        }
        v.push_back(i);
      }
      else{
        std::vector<LNO> v;
        v.reserve(sizeChunk);
        v.push_back(i);
        gidIndices->put(uniqueKey, v);
      }
    }
  
    numberOfUniqueGids = gidIndices->size();
  
    gidIndices->arrayify(uniqueGidQueries, uniqueGidQueryIndices);
  
    delete gidIndices;
  
    gidOutBuf = new AppGID [numberOfUniqueGids];
    for (int p=0; p < _numProcs; p++){
      countOutBuf[p] = 0;
    }
  
    hashProc = new int [numberOfUniqueGids];
  
    for (LNO i=0; i < numberOfUniqueGids; i++){
      hashProc[i] = Teuchos::hashCode(uniqueGidQueries[i]) % _numProcs;
      countOutBuf[hashProc[i]]++;
    }
  
    offsetBuf[0] = 0;
  
    for (int p=0; p < _numProcs; p++){
      offsetBuf[p+1] = offsetBuf[p] + countOutBuf[p];
    }
  
    gidLocation = new LNO [numberOfUniqueGids];
  
    for (LNO i=0; i < numberOfUniqueGids; i++){
      AppGID gid = IdentifierTraits<AppGID>::keyToGid(uniqueGidQueries[i]);
      LNO offset = offsetBuf[hashProc[i]];
      gidLocation[i] = offset;
      gidOutBuf[offset] = gid;
      offsetBuf[hashProc[i]] = offset + 1;
    }

    delete [] hashProc;
    hashProc = NULL;

    gidOutView = Teuchos::ArrayView<AppGID>(gidOutBuf, numberOfUniqueGids);
  }
  else{
    for (int p=0; p < _numProcs; p++)
      countOutBuf[p] = 0;

    gidOutView = Teuchos::ArrayView<AppGID>(Teuchos::null);
  }

  if (offsetBuf){
    delete [] offsetBuf;
    offsetBuf = NULL;
  }

  AlltoAllv(*_comm, gidOutView, countOutView, gidInBuf, countInBuf);

  if (gidOutBuf){
    delete [] gidOutBuf;
    gidOutBuf = NULL;
  }

  ///////////////////////////////////////////////////////////////////////
  // Create and send answers to the processes that made requests of me.
  ///////////////////////////////////////////////////////////////////////

  total = 0;
  for (int p=0; p < _numProcs; p++){
    countOutBuf[p] = (*countInBuf)[p];
    total += countOutBuf[p];
  }
  gnoOutBuf = NULL;
  procOutBuf = NULL;
  Teuchos::ArrayView<int> procOutView;

  if (total > 0){
    procOutBuf = new int [total];
  
    if (needGnoInfo){
      gnoOutBuf = new GNO [total];
    }
  
    total=0;
  
    for (int p=0; p < _numProcs; p++){
      for (LNO i=0; i < countInBuf[p]; i++, total++){
        std::string s(IdentifierTraits<AppGID>::key(gidInBuf[total]));
        LNO index = gidToAnswerIndex->get(s);
        int proc = firstIndexToProc.upper_bound(index);
        procOutBuf[total] = proc-1;
  
        if (needGnoInfo){
          gnoOutBuf[total] = gnoInBuf[index];
        }
      }
    }
  
    gidInBuf.release();
    delete gidToAnswerIndex;
  
    if (needGnoInfo){
      gnoInBuf.release();
    }

    procOutView = Teuchos::ArrayView<int>(procOutBuf, total);
  }
  else{
    procOutView = Teuchos::ArrayView<int>(Teuchos::null);
  }

  AlltoAllv(*_comm, procOutView, countOutView, procInBuf, countInBuf);

  if (procOutBuf){
    delete [] procOutBuf;
    procOutBuf = NULL;
  }

  if (needGnoInfo){
    if (total > 0)
      gnoOutView = Teuchos::ArrayView<GNO>(gnoOutBuf, total);
    else
      gnoOutView = Teuchos::ArrayView<GNO>(Teuchos::null);

    AlltoAllv(*_comm, gnoOutView, countOutView, gnoInBuf, countInBuf);

    delete [] gnoOutBuf;
  }

  delete [] countOutBuf;

  ///////////////////////////////////////////////////////////////////////
  // Done.  Process the replies to my queries
  ///////////////////////////////////////////////////////////////////////

  for (LNO i=0; i < numberOfUniqueGids; i++){

    std::string s(uniqueGidQueries[i]);
    std::vector<LNO> v(uniqueGidQueryIndices[i]);

    int gidProc = (*procInBuf)[gidLocation[i]];

    LNO gno;
    if (needGnoInfo){
      gno = (*gnoInBuf)[gidLocation[i]];
    }
    else{
      gno = static_cast<GNO>(IdentifierTraits<AppGID>::keyToGid(s));
    }

    for (LNO j=0; j < v.size(); j++){
      out_gno[v[j]] = gno;
      out_proc[v[j]] = gidProc;
    }
  }
}

}   // end namespace Z2

#endif /* _ZOLTAN2_IDENTIFIERMAP_DEF_HPP_ */
