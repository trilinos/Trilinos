//@HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_IdentifierMap_def.hpp

    \brief Definitions for the IdentifierMap class.
*/


#ifndef _ZOLTAN2_IDENTIFIERMAP_DEF_HPP_
#define _ZOLTAN2_IDENTIFIERMAP_DEF_HPP_

#include <vector>
#include <map>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_TestForException.hpp>
#include <Zoltan2_IdentifierTraits.hpp>
#include <Zoltan2_Environment.hpp>
#include <Zoltan2_AlltoAll.hpp>

namespace Z2
{

template<typename AppLID, typename AppGID, typename LNO, typename GNO> 
  IdentifierMap<AppLID,AppGID,LNO,GNO>::IdentifierMap(
    Teuchos::RCP<const Teuchos::Comm<int> > &in_comm, 
    Teuchos::RCP<Zoltan2::Environment> &env,
    typename Teuchos::ArrayRCP<AppGID> &gids, 
    typename Teuchos::ArrayRCP<AppLID> &lids) 
         : _comm(in_comm),  _env(env), _myGids(gids), _myLids(lids),
           _globalNumberOfIds(0), _localNumberOfIds(0), _haveLocalIds(false),
           _myRank(0), _numProcs(0)
{
  setupMap();
}

  /*! Constructor */
template<typename AppLID, typename AppGID, typename LNO, typename GNO>
  IdentifierMap<AppLID,AppGID,LNO,GNO>::IdentifierMap()  
         : _comm(), _env(), _myGids(), _myLids(),
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
    //TODO    default should work, shouldn't it?
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
    Teuchos::RCP<Zoltan2::Environment> &env,
    Teuchos::ArrayRCP<AppGID> &gids, 
    Teuchos::ArrayRCP<AppLID> &lids)
{
  _gnoDist.release();
  _gidHash.release();
  _lidHash.release();
  _comm.release();
  _env.release();
  _myGids.release();
  _myLids.release();

  _comm = in_comm;
  _env = env;
  _myGids = gids;
  _myLids = lids;
  _globalNumberOfIds = 0; 
  _localNumberOfIds = 0; 
  _haveLocalIds=false;;
  _myRank = 0; 
  _numProcs = 0;

  setupMap();
}

template<typename AppLID, typename AppGID, typename LNO, typename GNO>
  bool IdentifierMap<AppLID, AppGID, LNO, GNO>::gnosAreGids() const
{
  return IdentifierTraits<AppGID>::isGlobalOrdinalType();
}

template<typename AppLID, typename AppGID, typename LNO, typename GNO>
  void IdentifierMap<AppLID, AppGID, LNO, GNO>::gidTranslate(
    Teuchos::ArrayView<AppGID> gid, 
    Teuchos::ArrayView<GNO> gno,
    TranslationType tt) const
{
  typedef typename Teuchos::Array<GNO>::size_type teuchos_size_t;
  teuchos_size_t len=gid.size();

  if (len == 0){
    return;
  }

  Z2_LOCAL_INPUT_ASSERTION(*_comm, *_env, 
    "invalid TranslationType", 
    (tt==TRANSLATE_GNO_TO_GID) || (tt==TRANSLATE_GID_TO_GNO), 
    Z2_BASIC_ASSERTION);

  Z2_LOCAL_INPUT_ASSERTION(*_comm, *_env, 
    "Destination array is too small",
    ((tt==TRANSLATE_GNO_TO_GID) && (gid.size() >= gno.size())) || 
     ((tt==TRANSLATE_GID_TO_GNO) && (gno.size() >= gid.size())),
    Z2_BASIC_ASSERTION);

  if (IdentifierTraits<AppGID>::isGlobalOrdinalType()){   
                              // our gnos are the app gids
    if (tt == TRANSLATE_GNO_TO_GID)
      for (teuchos_size_t i=0; i < len; i++)
        gid[i] = static_cast<AppGID>(gno[i]);
    else
      for (teuchos_size_t i=0; i < len; i++)
        gno[i] = static_cast<GNO>(gid[i]);
  }
  else{              // we mapped gids to consecutive gnos
  
    GNO firstGno = _gnoDist[_myRank];
    GNO endGno = _gnoDist[_myRank + 1];

    if (tt == TRANSLATE_GNO_TO_GID){
      for (teuchos_size_t i=0; i < len; i++){

        Z2_LOCAL_INPUT_ASSERTION(*_comm, *_env, "invalid global number", 
        (gno[i] < firstGno) || (gno[i] >= endGno), Z2_BASIC_ASSERTION);

        gid[i] = _myGids[gno[i] - firstGno];
      }
    }
    else{
      LNO idx;
      for (teuchos_size_t i=0; i < len; i++){
        try{
          idx = _gidHash->get(Z2::IdentifierTraits<AppGID>::key(gid[i]));
        }
        catch (const std::exception &e) {
          Z2_THROW_OUTSIDE_ERROR(*_env, e);
        }
        
        gno[i] = firstGno + idx;
      }
    }
  }
  return;
}

template<typename AppLID, typename AppGID, typename LNO, typename GNO>
  void IdentifierMap<AppLID, AppGID, LNO, GNO>::lidTranslate(
    Teuchos::ArrayView<AppLID> lid, 
    Teuchos::ArrayView<GNO> gno, 
    TranslationType tt) const
{
  typedef typename Teuchos::Array<GNO>::size_type teuchos_size_t;
  teuchos_size_t len=lid.size();

  if (len == 0){
    return;
  }

  Z2_LOCAL_INPUT_ASSERTION(*_comm, *_env, 
    "invalid TranslationType", 
    (tt==TRANSLATE_GNO_TO_LID) || (tt==TRANSLATE_LID_TO_GNO), 
    Z2_BASIC_ASSERTION);

  Z2_LOCAL_INPUT_ASSERTION(*_comm, *_env, 
    "Destination array is too small",
    ((tt==TRANSLATE_GNO_TO_LID) && (lid.size() >= gno.size())) || ((tt==TRANSLATE_LID_TO_GNO) && (gno.size() >= lid.size())),
    Z2_BASIC_ASSERTION);

  Z2_LOCAL_INPUT_ASSERTION(*_comm, *_env, 
    "local ID translation is requested but none were provided",
     _haveLocalIds,
    Z2_BASIC_ASSERTION);

  GNO firstGno(0), endGno(0);
  if (_gnoDist.size() > 0){
    firstGno = _gnoDist[_myRank];
    endGno = _gnoDist[_myRank+1];
  }
  
  if (tt == TRANSLATE_GNO_TO_LID){
    for (teuchos_size_t i=0; i < len; i++){
      LNO idx = 0;
      if (_gnoDist.size() > 0) {// gnos are consecutive

        Z2_LOCAL_INPUT_ASSERTION(*_comm, *_env, 
          "invalid global number", 
          (gno[i] >= firstGno) && (gno[i] < endGno),
          Z2_BASIC_ASSERTION);

        idx = gno[i] - firstGno;
      }
      else {                    // gnos must be the app gids
        try{
          idx = _gidHash->get(
            IdentifierTraits<AppGID>::key(static_cast<AppGID>(gno[i])));
        }
        catch (const std::exception &e) {
          Z2_THROW_OUTSIDE_ERROR(*_env, e);
        }
      }
      
      lid[i] = _myLids[idx];
    }
  }
  else{
    for (teuchos_size_t i=0; i < len; i++){
      LNO idx(0);
      try{
        idx = _lidHash->get(IdentifierTraits<AppLID>::key(lid[i]));
      }
      catch (const std::exception &e) {
        Z2_THROW_OUTSIDE_ERROR(*_env, e);
      }

      if (_gnoDist.size() > 0)  // gnos are consecutive
        gno[i] = firstGno + idx;
      else                     // gnos must be the app gids
        gno[i] = static_cast<GNO>(_myGids[idx]);
    }
  }
}

template<typename AppLID, typename AppGID, typename LNO, typename GNO>
  void IdentifierMap<AppLID, AppGID, LNO, GNO>::gidGlobalTranslate(
    Teuchos::ArrayView<const AppGID> in_gid,
    Teuchos::ArrayView<GNO> out_gno,
    Teuchos::ArrayView<int> out_proc) const
{
  typedef typename Teuchos::Array<GNO>::size_type teuchos_size_t;
  typedef typename Teuchos::Hashtable<double, LNO> id2index_hash_t;
  typedef typename Teuchos::Hashtable<double, Teuchos::Array<LNO> > 
    id2index_array_hash_t;

  teuchos_size_t len=in_gid.size();

  if (len == 0){
    return;
  }

  Z2_LOCAL_INPUT_ASSERTION(*_comm, *_env, 
    "Destination array is too small", 
    (out_gno.size() >= len) && (out_proc.size() >= len),
    Z2_BASIC_ASSERTION);

  if (IdentifierTraits<AppGID>::isGlobalOrdinalType() && (_gnoDist.size() > 0)){

    // Easy case - communication is not needed.
    // Global numbers are the application global IDs and
    // they are increasing consecutively with rank.
 
    typename std::map<GNO, int> firstGnoToProc;
    typename std::map<GNO, int>::iterator pos;

    for (int p=0; p <= _numProcs; p++){
      firstGnoToProc[_gnoDist[p]] = p;
    }

    for (teuchos_size_t i=0; i < len; i++){
      GNO globalNumber = static_cast<GNO>(in_gid[i]);
      out_gno[i] = globalNumber;
      pos = firstGnoToProc.upper_bound(globalNumber);
      out_proc[i] = pos->first - 1;
    }

    return;
  }

  bool needGnoInfo = !IdentifierTraits<AppGID>::isGlobalOrdinalType();

  ///////////////////////////////////////////////////////////////////////
  // First: Hash each of my AppGIDs to a process that will answer
  // for it.  Send my AppGIDs (and the Gnos if they are different)
  // to their assigned processes.  Build a search structure for
  // the AppGIDs that were assigned to me, so I can reply with
  // with the process owning them (and their Gnos if they are different).
  ///////////////////////////////////////////////////////////////////////

  Teuchos::Array<int> hashProc(0);
  Teuchos::Array<AppGID> gidOutBuf(0);
  Teuchos::Array<GNO> gnoOutBuf(0);
  Teuchos::Array<LNO> countOutBuf(_numProcs, 0);
  Teuchos::Array<LNO> offsetBuf(_numProcs + 1, 0);

  Teuchos::ArrayRCP<AppGID> gidInBuf();
  Teuchos::ArrayRCP<GNO> gnoInBuf();
  Teuchos::ArrayRCP<LNO> countInBuf();

  if (_localNumberOfIds > 0){

    try{ 
      hashProc.reserve(_localNumberOfIds); 
    }
    catch(...)
      Z2_LOCAL_MEMORY_ASSERTION(*_comm, *_env, _localNumberOfIds, false); 

    try{ 
      gidOutBuf.reserve(_localNumberOfIds); 
    }
    catch(...)
      Z2_LOCAL_MEMORY_ASSERTION(*_comm, *_env, _localNumberOfIds, false); 

    for (teuchos_size_t i=0; i < _localNumberOfIds; i++){
      hashProc[i] = IdentifierTraits<AppGID>::hashCode(_myGids[i]) % _numProcs;
      countOutBuf[hashProc[i]]++;
    }
  
    for (int p=1; p <= _numProcs; p++){
      offsetBuf[p] = offsetBuf[p-1] + countOutBuf[p-1];
    }
  
    if (needGnoInfo){   
      // The gnos are not the gids, which also implies that
      // gnos are consecutive numbers given by _gnoDist.
      gnoOutBuf.resize(_localNumberOfIds);
    }
  
    for (teuchos_size_t i=0; i < _localNumberOfIds; i++){
      LNO offset = offsetBuf[hashProc[i]];
      gidOutBuf[offset] = _myGids[i];
      if (needGnoInfo)
        gnoOutBuf[offset] = _gnoDist[_myRank] + i;
      offsetBuf[hashProc[i]] = offset + 1;
    }
    hashProc.clear();
  }

  // Teuchos comment #1: An Array can be passed for an ArrayView parameter.
  // Teuchos comment #2: AppGID need not be a Teuchos Packet type,
  //                     so we wrote our own AlltoAllv.
  // Z2::AlltoAllv comment: Buffers are in process rank order.

  try{
    AlltoAllv(*_comm, *_env, gidOutBuf, countOutBuf, gidInBuf, countInBuf);
  }
  catch (const std::exception &e)
    Z2_THROW_ZOLTAN2_ERROR(*_env, e);

  gidOutBuf.clear();
  
  if (needGnoInfo){
    countInBuf.release();
    try{
      AlltoAllv(*_comm, *_env, gnoOutBuf, countOutBuf, gnoInBuf, countInBuf);
    }
    catch (const std::exception &e)
      Z2_THROW_ZOLTAN2_ERROR(*_env, e);
  }

  gnoOutBuf.clear();
  countOutBuf.clear();

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

  id2index_hash_t gidToIndex(total);

  total = 0;
  for (int p=0; p < _numProcs; p++){
    for (LNO i=countInBuf[p]; i < countInBuf[p+1]; i++, total++){
      try{
        gidToIndex.put(IdentifierTraits<AppGID>::key(gidInBuf[total]), total);
      }
      catch (const std::exception &e) 
        Z2_THROW_OUTSIDE_ERROR(*_env, e);
    }
  }

  // Keep gnoInBuf.  We're done with the others.

  gidInBuf.release();
  countInBuf.release();

  ///////////////////////////////////////////////////////////////////////
  // Send a request for information to the "answer process" for each 
  // of the GIDs in in_gid.  First remove duplicate GIDs from the list.
  ///////////////////////////////////////////////////////////////////////

  
  Teuchos::Array<double> uniqueGidQueries(0);
  Teuchos::Array<Teuchos::Array<LNO> > uniqueGidQueryIndices(0);
  teuchos_size_t numberOfUniqueGids = 0;
  Teuchos::Array<LNO> gidLocation(0);

  countOutBuf.resize(_numProcs, 0);

  if (len > 0){
  
    // For efficiency, guess a reasonable size for the array.
    // (In an input adapter, how many objects will have the same neighbor?)

    teuchos_size_t sizeChunk = 4;

    teuchos_size_t tableSize = len / sizeChunk;

    tableSize =  (tableSize < 1) ? 1 : tableSize;

    id2index_array_hash_t *gidIndices = NULL;
    Z2_ASYNC_MEMORY_ALLOC(*_comm, *_env, id2index_array_hash_t, gidIndices, tableSize);

    for (LNO i=0; i < len; i++){

      double uniqueKey(IdentifierTraits<AppGID>::key(in_gid[i]));

      if (gidIndices->containsKey(uniqueKey)){
        Teuchos::Array<LNO> v;
        try{
          v = gidIndices->get(uniqueKey);
        }
        catch (const std::exception &e) 
          Z2_THROW_OUTSIDE_ERROR(*_env, e);
        
        teuchos_size_t n = v.size();
        if (n % sizeChunk == 0){
          v.reserve(n + sizeChunk);
        }
        v.push_back(i);
      }
      else{
        Teuchos::Array<LNO> v(sizeChunk);
        v[0] = i;
        try{
          gidIndices->put(uniqueKey, v);
        }
        catch (const std::exception &e) 
          Z2_THROW_OUTSIDE_ERROR(*_env, e);
      }
    }
  
    numberOfUniqueGids = gidIndices->size();

    gidIndices->arrayify(uniqueGidQueries, uniqueGidQueryIndices);
  
    delete gidIndices;
  
    try{ 
      gidOutBuf.reserve(numberOfUniqueGids); 
    }
    catch(...)
      Z2_LOCAL_MEMORY_ASSERTION(*_comm, *_env, numberOfUniqueGids, false); 

    try{ 
      hashProc.reserve(numberOfUniqueGids);
    }
    catch(...)
      Z2_LOCAL_MEMORY_ASSERTION(*_comm, *_env, numberOfUniqueGids, false); 
  
    for (teuchos_size_t i=0; i < numberOfUniqueGids; i++){
      hashProc[i] = Teuchos::hashCode(uniqueGidQueries[i]) % _numProcs;
      countOutBuf[hashProc[i]]++;
    }
  
    offsetBuf[0] = 0;
  
    for (int p=0; p < _numProcs; p++){
      offsetBuf[p+1] = offsetBuf[p] + countOutBuf[p];
    }
  
    try{ 
      gidLocation.reserve(numberOfUniqueGids);
    }
    catch(...)
      Z2_LOCAL_MEMORY_ASSERTION(*_comm, *_env, numberOfUniqueGids, false); 
  
    for (teuchos_size_t i=0; i < numberOfUniqueGids; i++){
      AppGID gid = IdentifierTraits<AppGID>::keyToGid(uniqueGidQueries[i]);
      gidLocation[i] = offsetBuf[hashProc[i]];
      gidOutBuf[gidLocation[i]] = gid;
      offsetBuf[hashProc[i]] = gidLocation[i] + 1;
    }

    hashProc.clear();
  }

  try{
    AlltoAllv(*_comm, *_env, gidOutBuf, countOutBuf, gidInBuf, countInBuf);
  }
  catch (const std::exception &e)
    Z2_THROW_ZOLTAN2_ERROR(*_env, e)

  gidOutBuf.clear();

  ///////////////////////////////////////////////////////////////////////
  // Create and send answers to the processes that made requests of me.
  ///////////////////////////////////////////////////////////////////////

  total = 0;

  for (int p=0; p < _numProcs; p++){
    countOutBuf[p] = countInBuf[p];
    total += countOutBuf[p];
  }

  Teuchos::Array<int> procOutBuf(total);
  Teuchos::ArrayRCP<int> procInBuf();

  if (needGnoInfo){
    try{ 
      gnoOutBuf.reserve(total);
    }
    catch(...)
      Z2_LOCAL_MEMORY_ASSERTION(*_comm, *_env, total, false); 
  }

  if (total > 0){
  
    total=0;
  
    for (int p=0; p < _numProcs; p++){
      for (LNO i=0; i < countInBuf[p]; i++, total++){
        double k(IdentifierTraits<AppGID>::key(gidInBuf[total]));
        LNO index(0);
        try{
          index = gidToIndex.get(k);
        }
        catch (const std::exception &e) 
          Z2_THROW_OUTSIDE_ERROR(*_env, e);
        
        int proc = firstIndexToProc.upper_bound(index);
        procOutBuf[total] = proc-1;
  
        if (needGnoInfo){
          gnoOutBuf[total] = gnoInBuf[index];
        }
      }
    }
  }

  gidInBuf.release();
  if (needGnoInfo){
    gnoInBuf.release();
  }

  try{
    AlltoAllv(*_comm, *_env, procOutBuf, countOutBuf, procInBuf, countInBuf);
  }
  catch (const std::exception &e)
    Z2_THROW_ZOLTAN2_ERROR(*_env, e);

  procOutBuf.clear();

  if (needGnoInfo){
    try{
      AlltoAllv(*_comm, *_env, gnoOutBuf, countOutBuf, gnoInBuf, countInBuf);
    }
    catch (const std::exception &e)
      Z2_THROW_ZOLTAN2_ERROR(*_env, e);

    gnoOutBuf.clear();
  }

  countOutBuf.clear();

  ///////////////////////////////////////////////////////////////////////
  // Done.  Process the replies to my queries
  ///////////////////////////////////////////////////////////////////////

  for (teuchos_size_t i=0; i < numberOfUniqueGids; i++){

    double key(uniqueGidQueries[i]);
    Teuchos::Array<LNO> v(uniqueGidQueryIndices[i]);

    int gidProc = procInBuf[gidLocation[i]];

    LNO gno;
    if (needGnoInfo){
      gno = gnoInBuf[gidLocation[i]];
    }
    else{
      gno = static_cast<GNO>(IdentifierTraits<AppGID>::keyToGid(key));
    }

    for (teuchos_size_t j=0; j < v.size(); j++){
      out_gno[v[j]] = gno;
      out_proc[v[j]] = gidProc;
    }
  }
}


template<typename AppLID, typename AppGID, typename LNO, typename GNO> 
  void IdentifierMap<AppLID,AppGID,LNO,GNO>::setupMap(void)
{
  _numProcs = _comm->getSize(); 
  _myRank = _comm->getRank(); 

  Z2_GLOBAL_INPUT_ASSERTION( *_comm, *_env, 
           "application global ID type is not supported",
           IdentifierTraits<AppGID>::is_valid_id_type() == true,
           Z2_BASIC_ASSERTION);

  _localNumberOfIds = _myGids.size();

  typedef typename Teuchos::Array<GNO>::size_type teuchos_size_t;
  typedef typename Teuchos::Hashtable<double, LNO> id2index_hash_t;

  Teuchos::Tuple<teuchos_size_t, 4> counts;

  counts[0] = _myLids.size();
  counts[1] = _localNumberOfIds;
  counts[2] = counts[3] = 0;

  try{
    Teuchos::reduceAll<int, teuchos_size_t>(*_comm, Teuchos::REDUCE_SUM, 
      2, counts.getRawPtr(), counts.getRawPtr()+2);
  } 
  catch (const std::exception &e) {
    Z2_THROW_OUTSIDE_ERROR(*_env, e);
  }

  _haveLocalIds = (counts[2] > 0);
  _globalNumberOfIds = counts[3];

  Z2_GLOBAL_INPUT_ASSERTION( *_comm, *_env, 
       "number of global IDs does not equal number of local IDs",
      !_haveLocalIds || (counts[0] == _localNumberOfIds),
       Z2_BASIC_ASSERTION);

  if (_haveLocalIds){   // hash LID to index in LID vector
    id2index_hash_t *p = NULL;
    if (_localNumberOfIds){
      try{
        p = new id2index_hash_t(_localNumberOfIds);
      }
      catch (const std::exception &e) 
        Z2_LOCAL_MEMORY_ASSERTION(*_comm, *_env, _localNumberOfIds, false); 
    }

    AppLID *lidPtr = _myLids.get();  // for performance

    for (teuchos_size_t i=0; i < _localNumberOfIds; i++){
      try{
        p->put(IdentifierTraits<AppLID>::key(lidPtr[i]), i);
      }
      catch (const std::exception &e) 
        Z2_THROW_OUTSIDE_ERROR(*_env, e);
    }

    _lidHash = Teuchos::RCP<id2index_hash_t>(p);
  }

  // If the application's global ID data type (AppGID) is a Teuchos Ordinal,
  // we will be able to use it as our internal global numbers (GNO).  
  // Otherwise we will have to map their global IDs to valid global numbers.

  if (IdentifierTraits<AppGID>::isGlobalOrdinalType()){

    // Are the AppGIDs consecutive and increasing with process rank? 
    // If so GID/proc lookups can be optimized.

    AppGID min(0), max(0), globalMin(0), globalMax(0);
    AppGID *gidPtr = _myGids.get();  // for performance
    min = max = gidPtr[0];
    AppGID checkVal = min;
    bool consecutive = true;
  
    for (teuchos_size_t i=1; i < _localNumberOfIds; i++){
      if (consecutive && (gidPtr[i] != ++checkVal)){
        consecutive=false;
        break;
      }
      if (gidPtr[i] < min)
        min = gidPtr[i];
      else if (gidPtr[i] > max)
        max = gidPtr[i];
    }

    Teuchos::Tuple<AppGID, 4> results;

    results[0] = static_cast<AppGID>(consecutive ? 1 : 0);
    results[1] = min;
    results[2] = results[3] = 0;

    try{
      Teuchos::reduceAll<int, AppGID>(*_comm, Teuchos::REDUCE_MIN, 2, 
        results.getRawPtr(), results.getRawPtr()+2);
    }
    catch (const std::exception &e) {
      Z2_THROW_OUTSIDE_ERROR(*_env, e);
    }

    if (results[2] != 1)       // min of consecutive flags
      consecutive=false;

    if (consecutive){
      globalMin = results[3];
      try{
        Teuchos::reduceAll<int, AppGID>(*_comm, Teuchos::REDUCE_MAX, 1, &max, 
          &globalMax);
      }
      catch (const std::exception &e) {
        Z2_THROW_OUTSIDE_ERROR(*_env, e);
      }

      if (globalMax - globalMin + 1 != static_cast<AppGID>(_globalNumberOfIds))
        consecutive = false;   // there are gaps in the gids
  
      if (consecutive){
        AppGID myStart = _myGids[0];

        _gnoDist = Teuchos::ArrayRCP<GNO>(_numProcs + 1);

        AppGID *startGID = static_cast<AppGID *>(_gnoDist.getRawPtr());
      
        try{
          Teuchos::gatherAll<int, AppGID>(*_comm, 1, &myStart, _numProcs, 
            startGID);
        }
        catch (const std::exception &e) {
          Z2_THROW_OUTSIDE_ERROR(*_env, e);
        }
      
        for (int p=1; p < _numProcs; p++){
          if (startGID[p] < startGID[p-1]){
            consecutive = false;  // gids do not increase with process rank
            break;
          }
        }
        if (consecutive){
          startGID[_numProcs] = globalMax + 1;
        }
      }
    }
  }
  else{

    // AppGIDs are not Ordinals.  We map them to consecutive 
    // global numbers starting with 0. 

    try{
      _gnoDist = Teuchos::ArrayRCP<GNO>(_numProcs + 1, 0);
    }
    catch (const std::exception &e) {
      Z2_THROW_OUTSIDE_ERROR(*_env, e);
    }

    GNO myNum = static_cast<GNO>(_localNumberOfIds);

    try{
      Teuchos::scan<int, GNO>(*_comm, Teuchos::REDUCE_SUM, 1, &myNum, 
        _gnoDist.getRawPtr() + 1);
    }
    catch (const std::exception &e) {
      Z2_THROW_OUTSIDE_ERROR(*_env, e);
    }
  }

  if (_gnoDist.size() == 0){

    // We need a hash table mapping the application global ID
    // to its index in _myGids.

    id2index_hash_t *p = NULL;
    if (_localNumberOfIds){
      try{
        p = new id2index_hash_t(_localNumberOfIds);
      }
      catch (const std::exception &e) 
        Z2_LOCAL_MEMORY_ASSERTION(*_comm, *_env, _localNumberOfIds, false); 
    }

    AppGID *gidPtr = _myGids.get();  // for performance

    for (teuchos_size_t i=0; i < _localNumberOfIds; i++){
      try{
        p->put(IdentifierTraits<AppGID>::key(gidPtr[i]), i);
      }
      catch (const std::exception &e) {
        Z2_THROW_OUTSIDE_ERROR(*_env, e);
      }
    }
    _gidHash = Teuchos::RCP<id2index_hash_t>(p);
  }
}

}   // end namespace Z2

#endif /* _ZOLTAN2_IDENTIFIERMAP_DEF_HPP_ */
