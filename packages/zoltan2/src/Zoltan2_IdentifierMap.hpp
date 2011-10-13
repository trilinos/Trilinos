// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

#ifndef _ZOLTAN2_IDENTIFIERMAP_HPP_
#define _ZOLTAN2_IDENTIFIERMAP_HPP_

/*! \file Zoltan2_IdentifierMap.hpp

    \brief IdentifierMap class.
*/

#include <vector>
#include <map>
#include <Teuchos_Hashtable.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_TestForException.hpp>
#include <Zoltan2_Standards.hpp>
#include <Zoltan2_IdentifierTraits.hpp>

namespace Zoltan2
{

/*! \brief Identifier translations available from IdentifierMap.

  GNO and LNO are the internal Zoltan2 global and local IDs.
  GID and LID are the global and local IDs used by the caller.
  The internal GNOs may be the same as the GIDs if the caller
  used a data type that can be used as a Teuchos Ordinal.
*/

enum TranslationType {
  TRANSLATE_GNO_TO_GID,
  TRANSLATE_GID_TO_GNO,
  TRANSLATE_GNO_TO_LID,
  TRANSLATE_LID_TO_GNO
};

/*! Z2::IdentifierMap
    \brief An IdentifierMap manages a global space of unique object identifiers.

    TODO - trim down comments and code once we get it all straight
           replace new/delete with memory wrappers
           exception handling
           use Kokkos node 

    TODO - methods called by Problem to redistribute in place or to
          create a redistributed copy.
*/

////////////////////////////////////////////////////////////////////
// Declarations
////////////////////////////////////////////////////////////////////

template<typename AppLID, typename AppGID, typename LNO=int, 
  typename GNO=AppGID> 
    class IdentifierMap{

private:

  // Input communicator

  RCP<const Comm<int> > comm_;

  // Problem parameters, library configuration.

  RCP<Environment> env_;

  // Application global and local IDs

  ArrayRCP<AppGID> myGids_; 
  ArrayRCP<AppLID> myLids_;

  // In the case of consecutive ordinal application global IDs,
  // gnoDist[p] is the first global number on process p, and
  // we don't need the gidHash_.

  ArrayRCP<GNO> gnoDist_;

  // A hash table from application global ID key to our local index.

  RCP<Teuchos::Hashtable<double, LNO> >  gidHash_;

  // A hash table from application local ID key to our local index.

  RCP<Teuchos::Hashtable<double, LNO> >  lidHash_;

  typename Array<GNO>::size_type globalNumberOfIds_;
  typename Array<GNO>::size_type localNumberOfIds_;
  bool haveLocalIds_;
  int myRank_;
  int numProcs_;
  bool consecutiveIdsRequired_;   // TODO

  void setupMap();

public:

  /*! Constructor - Must be called by all processes 
   *
   *  We are not taking a Kokkos node here.  In the future we should.
   *
   *  Problem: we have to template the communicator on int, because Tpetra_Map
   *  only takes Comm<int>. 
   */

  explicit IdentifierMap( RCP<const Comm<int> > &incomm_, 
                          RCP<Environment > &env, 
                          ArrayRCP<AppGID> &gids, 
                          ArrayRCP<AppLID> &lids,
                          bool idsMustBeConsecutive=false);

  /*! Constructor 
      This constructor does not need to be called by all processes.
   */
  IdentifierMap();

  /*! Destructor */
  ~IdentifierMap();

  /*! Copy Constructor */
  IdentifierMap(const IdentifierMap &id);

  /*! Assignment operator */
  IdentifierMap &operator=(const IdentifierMap &id);

  /*! Initialize object if not done in the constructor */
  void initialize(RCP<const Comm<int> > &incomm_,
                  RCP<Environment > &env,
                  ArrayRCP<AppGID> &gids,
                  ArrayRCP<AppLID> &lids,
                  bool idsMustBeConsecutive=false);

  /*! Return true if we are using the application global IDs 
   *  for our internal global numbers 
   */
  bool gnosAreGids() const;

  /*! Map the application global IDs to internal global numbers or vice versa.

      \param gid an array of caller's global IDs
      \param gno an array of Zoltan2 global numbers
      \param tt should be TRANSLATE_GNO_TO_GID or TRANSLATE_GID_TO_GNO 

      This is a local call.  If gid is a vector of application global IDs, then
      gno will be set to the corresponding internal global numbers.  If the
      gno vector contains internal global numbers, they will be translated
      to application global IDs.  The application global IDs must be from
      those supplied by this process.

      TODO: Why does code fail to compile if I pass gid and gno by reference?
   */
  void gidTranslate(ArrayView<AppGID> gid, 
                    ArrayView<GNO> gno,
                    TranslationType tt) const;

  /*! Map application local IDs to internal global numbers or vice versa.

      \param lid an array of caller's local IDs
      \param gno an array of Zoltan2 global numbers
      \param tt should be TRANSLATE_GNO_TO_LID or TRANSLATE_LID_TO_GNO 

      This is a local call.  If lid is a vector of application local IDs, then
      gno will be set to the corresponding internal global numbers.  If the
      gno vector contains internal global numbers, they will be translated
      to application local IDs.  The application local IDs must be from
      those supplied by this process.
   */
  void lidTranslate(ArrayView<AppLID> lid, 
                    ArrayView<GNO> gno,
                    TranslationType tt) const;

  /*! Map application global IDs to internal global numbers.

      \param in_gid input, an array of the Caller's global IDs
      \param out_gno output, an optional array of the corresponding 
          global numbers used by Zoltan2.  If out_gno.size() is zero,
          we assume global numbers are not needed.
      \param out_proc output, an array of the process ranks corresponding with
                         the in_gid and out_gno, out_proc[i] is the process
                         that supplied in_gid[i] to Zoltan2.

      All processes must call this.  The global IDs 
      supplied may belong to another process.  
   */
  void gidGlobalTranslate( ArrayView<const AppGID> in_gid,
                           ArrayView<GNO> out_gno,
                           ArrayView<int> out_proc) const;
};

////////////////////////////////////////////////////////////////////
// Definitions
////////////////////////////////////////////////////////////////////

template<typename AppLID, typename AppGID, typename LNO, typename GNO> 
  IdentifierMap<AppLID,AppGID,LNO,GNO>::IdentifierMap(
    RCP<const Comm<int> > &incomm_, RCP<Zoltan2::Environment> &env,
    ArrayRCP<AppGID> &gids, ArrayRCP<AppLID> &lids,
    bool idsMustBeConsecutive) 
         : comm_(incomm_),  env_(env), myGids_(gids), myLids_(lids),
           globalNumberOfIds_(0), localNumberOfIds_(0), haveLocalIds_(false),
           myRank_(0), numProcs_(0), 
           consecutiveIdsRequired_(idsMustBeConsecutive)
{
  setupMap();
}

  /*! Constructor */
template<typename AppLID, typename AppGID, typename LNO, typename GNO>
  IdentifierMap<AppLID,AppGID,LNO,GNO>::IdentifierMap()  
         : comm_(), env_(), myGids_(), myLids_(),
           globalNumberOfIds_(0), localNumberOfIds_(0), haveLocalIds_(false),
           myRank_(0), numProcs_(0), consecutiveIdsRequired_(false)
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
    RCP<const Comm<int> > &incomm_, 
    RCP<Zoltan2::Environment> &env,
    ArrayRCP<AppGID> &gids, ArrayRCP<AppLID> &lids,
    bool idsMustBeConsecutive) 
{
  gnoDist_.release();
  gidHash_.release();
  lidHash_.release();
  comm_.release();
  env_.release();
  myGids_.release();
  myLids_.release();

  comm_ = incomm_;
  env_ = env;
  myGids_ = gids;
  myLids_ = lids;
  globalNumberOfIds_ = 0; 
  localNumberOfIds_ = 0; 
  haveLocalIds_=false;;
  myRank_ = 0; 
  numProcs_ = 0;
  consecutiveIdsRequired_ = idsMustBeConsecutive;

  setupMap();
}

template<typename AppLID, typename AppGID, typename LNO, typename GNO>
  bool IdentifierMap<AppLID, AppGID, LNO, GNO>::gnosAreGids() const
{
  return IdentifierTraits<AppGID>::isGlobalOrdinalType();
}

template<typename AppLID, typename AppGID, typename LNO, typename GNO>
  void IdentifierMap<AppLID, AppGID, LNO, GNO>::gidTranslate(
    ArrayView<AppGID> gid, 
    ArrayView<GNO> gno,
    TranslationType tt) const
{
  typedef typename Array<GNO>::size_type teuchos_size_t;
  teuchos_size_t len=gid.size();

  if (len == 0){
    return;
  }

  Z2_LOCAL_INPUT_ASSERTION(*comm_, *env_, 
    "invalid TranslationType", 
    (tt==TRANSLATE_GNO_TO_GID) || (tt==TRANSLATE_GID_TO_GNO), 
    BASIC_ASSERTION);

  Z2_LOCAL_INPUT_ASSERTION(*comm_, *env_, 
    "Destination array is too small",
    ((tt==TRANSLATE_GNO_TO_GID) && (gid.size() >= gno.size())) || 
     ((tt==TRANSLATE_GID_TO_GNO) && (gno.size() >= gid.size())),
    BASIC_ASSERTION);

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
  
    GNO firstGno = gnoDist_[myRank_];
    GNO endGno = gnoDist_[myRank_ + 1];

    if (tt == TRANSLATE_GNO_TO_GID){
      for (teuchos_size_t i=0; i < len; i++){

        Z2_LOCAL_INPUT_ASSERTION(*comm_, *env_, "invalid global number", 
        (gno[i] < firstGno) || (gno[i] >= endGno), BASIC_ASSERTION);

        gid[i] = myGids_[gno[i] - firstGno];
      }
    }
    else{
      LNO idx;
      for (teuchos_size_t i=0; i < len; i++){
        try{
          idx = gidHash_->get(Zoltan2::IdentifierTraits<AppGID>::key(gid[i]));
        }
        catch (const std::exception &e) {
          Z2_THROW_OUTSIDE_ERROR(*env_, e);
        }
        
        gno[i] = firstGno + idx;
      }
    }
  }
  return;
}

template<typename AppLID, typename AppGID, typename LNO, typename GNO>
  void IdentifierMap<AppLID, AppGID, LNO, GNO>::lidTranslate(
    ArrayView<AppLID> lid, 
    ArrayView<GNO> gno, 
    TranslationType tt) const
{
  typedef typename Array<GNO>::size_type teuchos_size_t;
  teuchos_size_t len=lid.size();

  if (len == 0){
    return;
  }

  Z2_LOCAL_INPUT_ASSERTION(*comm_, *env_, 
    "invalid TranslationType", 
    (tt==TRANSLATE_GNO_TO_LID) || (tt==TRANSLATE_LID_TO_GNO), 
    BASIC_ASSERTION);

  Z2_LOCAL_INPUT_ASSERTION(*comm_, *env_, 
    "Destination array is too small",
    ((tt==TRANSLATE_GNO_TO_LID) && (lid.size() >= gno.size())) || ((tt==TRANSLATE_LID_TO_GNO) && (gno.size() >= lid.size())),
    BASIC_ASSERTION);

  Z2_LOCAL_INPUT_ASSERTION(*comm_, *env_, 
    "local ID translation is requested but none were provided",
     haveLocalIds_,
    BASIC_ASSERTION);

  GNO firstGno(0), endGno(0);
  if (gnoDist_.size() > 0){
    firstGno = gnoDist_[myRank_];
    endGno = gnoDist_[myRank_+1];
  }
  
  if (tt == TRANSLATE_GNO_TO_LID){
    for (teuchos_size_t i=0; i < len; i++){
      LNO idx = 0;
      if (gnoDist_.size() > 0) {// gnos are consecutive

        Z2_LOCAL_INPUT_ASSERTION(*comm_, *env_, 
          "invalid global number", 
          (gno[i] >= firstGno) && (gno[i] < endGno),
          BASIC_ASSERTION);

        idx = gno[i] - firstGno;
      }
      else {                    // gnos must be the app gids
        try{
          idx = gidHash_->get(
            IdentifierTraits<AppGID>::key(static_cast<AppGID>(gno[i])));
        }
        catch (const std::exception &e) {
          Z2_THROW_OUTSIDE_ERROR(*env_, e);
        }
      }
      
      lid[i] = myLids_[idx];
    }
  }
  else{
    for (teuchos_size_t i=0; i < len; i++){
      LNO idx(0);
      try{
        idx = lidHash_->get(IdentifierTraits<AppLID>::key(lid[i]));
      }
      catch (const std::exception &e) {
        Z2_THROW_OUTSIDE_ERROR(*env_, e);
      }

      if (gnoDist_.size() > 0)  // gnos are consecutive
        gno[i] = firstGno + idx;
      else                     // gnos must be the app gids
        gno[i] = static_cast<GNO>(myGids_[idx]);
    }
  }
}

template<typename AppLID, typename AppGID, typename LNO, typename GNO>
  void IdentifierMap<AppLID, AppGID, LNO, GNO>::gidGlobalTranslate(
    ArrayView<const AppGID> in_gid,
    ArrayView<GNO> out_gno,
    ArrayView<int> out_proc) const
{
  typedef typename Array<GNO>::size_type teuchos_size_t;
  typedef typename Teuchos::Hashtable<double, LNO> id2index_hash_t;
  typedef typename Teuchos::Hashtable<double, Array<LNO> > 
    id2index_array_hash_t;

  teuchos_size_t len=in_gid.size();

  if (len == 0){
    return;
  }

  bool skipGno = (out_gno.size() == 0);

  Z2_LOCAL_INPUT_ASSERTION(*comm_, *env_, 
    "Destination array is too small", 
    (out_proc.size() >= len) && (skipGno || (out_gno.size() >= len)),
    BASIC_ASSERTION);

  if (IdentifierTraits<AppGID>::isGlobalOrdinalType() && (gnoDist_.size() > 0)){

    // Easy case - communication is not needed.
    // Global numbers are the application global IDs and
    // they are increasing consecutively with rank.
 
    typename std::map<GNO, int> firstGnoToProc;
    typename std::map<GNO, int>::iterator pos;

    for (int p=0; p <= numProcs_; p++){
      firstGnoToProc[gnoDist_[p]] = p;
    }

    for (teuchos_size_t i=0; i < len; i++){
      GNO globalNumber = static_cast<GNO>(in_gid[i]);
      if (!skipGno)
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

  Array<int> hashProc(0);
  Array<AppGID> gidOutBuf(0);
  Array<GNO> gnoOutBuf(0);
  Array<LNO> countOutBuf(numProcs_, 0);
  Array<LNO> offsetBuf(numProcs_ + 1, 0);

  ArrayRCP<AppGID> gidInBuf();
  ArrayRCP<GNO> gnoInBuf();
  ArrayRCP<LNO> countInBuf();

  if (localNumberOfIds_ > 0){

    try{ 
      hashProc.reserve(localNumberOfIds_); 
    }
    catch(...)
      Z2_LOCAL_MEMORY_ASSERTION(*comm_, *env_, localNumberOfIds_, false); 

    try{ 
      gidOutBuf.reserve(localNumberOfIds_); 
    }
    catch(...)
      Z2_LOCAL_MEMORY_ASSERTION(*comm_, *env_, localNumberOfIds_, false); 

    for (teuchos_size_t i=0; i < localNumberOfIds_; i++){
      hashProc[i] = IdentifierTraits<AppGID>::hashCode(myGids_[i]) % numProcs_;
      countOutBuf[hashProc[i]]++;
    }
  
    for (int p=1; p <= numProcs_; p++){
      offsetBuf[p] = offsetBuf[p-1] + countOutBuf[p-1];
    }
  
    if (needGnoInfo){   
      // The gnos are not the gids, which also implies that
      // gnos are consecutive numbers given by gnoDist_.
      gnoOutBuf.resize(localNumberOfIds_);
    }
  
    for (teuchos_size_t i=0; i < localNumberOfIds_; i++){
      LNO offset = offsetBuf[hashProc[i]];
      gidOutBuf[offset] = myGids_[i];
      if (needGnoInfo)
        gnoOutBuf[offset] = gnoDist_[myRank_] + i;
      offsetBuf[hashProc[i]] = offset + 1;
    }
    hashProc.clear();
  }

  // Teuchos comment #1: An Array can be passed for an ArrayView parameter.
  // Teuchos comment #2: AppGID need not be a Teuchos Packet type,
  //                     so we wrote our own AlltoAllv.
  // Z2::AlltoAllv comment: Buffers are in process rank order.

  try{
    AlltoAllv(*comm_, *env_, gidOutBuf, countOutBuf, gidInBuf, countInBuf);
  }
  catch (const std::exception &e)
    Z2_THROW_ZOLTAN2_ERROR(*env_, e);

  gidOutBuf.clear();
  
  if (needGnoInfo){
    countInBuf.release();
    try{
      AlltoAllv(*comm_, *env_, gnoOutBuf, countOutBuf, gnoInBuf, countInBuf);
    }
    catch (const std::exception &e)
      Z2_THROW_ZOLTAN2_ERROR(*env_, e);
  }

  gnoOutBuf.clear();
  countOutBuf.clear();

  //
  // Save the information that was hashed to me so I can do lookups.
  //

  std::map<LNO, int> firstIndexToProc;
  LNO total = 0;

  for (int p=0; p < numProcs_; p++){
    firstIndexToProc[total] = p;
    total += countInBuf[p];
  }

  firstIndexToProc[total] = numProcs_;

  id2index_hash_t gidToIndex(total);

  total = 0;
  for (int p=0; p < numProcs_; p++){
    for (LNO i=countInBuf[p]; i < countInBuf[p+1]; i++, total++){
      try{
        gidToIndex.put(IdentifierTraits<AppGID>::key(gidInBuf[total]), total);
      }
      catch (const std::exception &e) 
        Z2_THROW_OUTSIDE_ERROR(*env_, e);
    }
  }

  // Keep gnoInBuf.  We're done with the others.

  gidInBuf.release();
  countInBuf.release();

  ///////////////////////////////////////////////////////////////////////
  // Send a request for information to the "answer process" for each 
  // of the GIDs in in_gid.  First remove duplicate GIDs from the list.
  ///////////////////////////////////////////////////////////////////////

  
  Array<double> uniqueGidQueries(0);
  Array<Array<LNO> > uniqueGidQueryIndices(0);
  teuchos_size_t numberOfUniqueGids = 0;
  Array<LNO> gidLocation(0);

  countOutBuf.resize(numProcs_, 0);

  if (len > 0){
  
    // For efficiency, guess a reasonable size for the array.
    // (In an input adapter, how many objects will have the same neighbor?)

    teuchos_size_t sizeChunk = 4;

    teuchos_size_t tableSize = len / sizeChunk;

    tableSize =  (tableSize < 1) ? 1 : tableSize;

    id2index_array_hash_t *gidIndices = NULL;
    Z2_ASYNC_MEMORY_ALLOC(*comm_, *env_, id2index_array_hash_t, 
      gidIndices, tableSize);

    for (LNO i=0; i < len; i++){

      double uniqueKey(IdentifierTraits<AppGID>::key(in_gid[i]));

      if (gidIndices->containsKey(uniqueKey)){
        Array<LNO> v;
        try{
          v = gidIndices->get(uniqueKey);
        }
        catch (const std::exception &e) 
          Z2_THROW_OUTSIDE_ERROR(*env_, e);
        
        teuchos_size_t n = v.size();
        if (n % sizeChunk == 0){
          v.reserve(n + sizeChunk);
        }
        v.push_back(i);
      }
      else{
        Array<LNO> v(sizeChunk);
        v[0] = i;
        try{
          gidIndices->put(uniqueKey, v);
        }
        catch (const std::exception &e) 
          Z2_THROW_OUTSIDE_ERROR(*env_, e);
      }
    }
  
    numberOfUniqueGids = gidIndices->size();

    gidIndices->arrayify(uniqueGidQueries, uniqueGidQueryIndices);
  
    delete gidIndices;
  
    try{ 
      gidOutBuf.reserve(numberOfUniqueGids); 
    }
    catch(...)
      Z2_LOCAL_MEMORY_ASSERTION(*comm_, *env_, numberOfUniqueGids, false); 

    try{ 
      hashProc.reserve(numberOfUniqueGids);
    }
    catch(...)
      Z2_LOCAL_MEMORY_ASSERTION(*comm_, *env_, numberOfUniqueGids, false); 
  
    for (teuchos_size_t i=0; i < numberOfUniqueGids; i++){
      hashProc[i] = Teuchos::hashCode(uniqueGidQueries[i]) % numProcs_;
      countOutBuf[hashProc[i]]++;
    }
  
    offsetBuf[0] = 0;
  
    for (int p=0; p < numProcs_; p++){
      offsetBuf[p+1] = offsetBuf[p] + countOutBuf[p];
    }
  
    try{ 
      gidLocation.reserve(numberOfUniqueGids);
    }
    catch(...)
      Z2_LOCAL_MEMORY_ASSERTION(*comm_, *env_, numberOfUniqueGids, false); 
  
    for (teuchos_size_t i=0; i < numberOfUniqueGids; i++){
      AppGID gid = IdentifierTraits<AppGID>::keyToGid(uniqueGidQueries[i]);
      gidLocation[i] = offsetBuf[hashProc[i]];
      gidOutBuf[gidLocation[i]] = gid;
      offsetBuf[hashProc[i]] = gidLocation[i] + 1;
    }

    hashProc.clear();
  }

  try{
    AlltoAllv(*comm_, *env_, gidOutBuf, countOutBuf, gidInBuf, countInBuf);
  }
  catch (const std::exception &e)
    Z2_THROW_ZOLTAN2_ERROR(*env_, e)

  gidOutBuf.clear();

  ///////////////////////////////////////////////////////////////////////
  // Create and send answers to the processes that made requests of me.
  ///////////////////////////////////////////////////////////////////////

  total = 0;

  for (int p=0; p < numProcs_; p++){
    countOutBuf[p] = countInBuf[p];
    total += countOutBuf[p];
  }

  Array<int> procOutBuf(total);
  ArrayRCP<int> procInBuf();

  if (needGnoInfo){
    try{ 
      gnoOutBuf.reserve(total);
    }
    catch(...)
      Z2_LOCAL_MEMORY_ASSERTION(*comm_, *env_, total, false); 
  }

  if (total > 0){
  
    total=0;
  
    for (int p=0; p < numProcs_; p++){
      for (LNO i=0; i < countInBuf[p]; i++, total++){
        double k(IdentifierTraits<AppGID>::key(gidInBuf[total]));
        LNO index(0);
        try{
          index = gidToIndex.get(k);
        }
        catch (const std::exception &e) 
          Z2_THROW_OUTSIDE_ERROR(*env_, e);
        
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
    AlltoAllv(*comm_, *env_, procOutBuf, countOutBuf, procInBuf, countInBuf);
  }
  catch (const std::exception &e)
    Z2_THROW_ZOLTAN2_ERROR(*env_, e);

  procOutBuf.clear();

  if (needGnoInfo){
    try{
      AlltoAllv(*comm_, *env_, gnoOutBuf, countOutBuf, gnoInBuf, countInBuf);
    }
    catch (const std::exception &e)
      Z2_THROW_ZOLTAN2_ERROR(*env_, e);

    gnoOutBuf.clear();
  }

  countOutBuf.clear();

  ///////////////////////////////////////////////////////////////////////
  // Done.  Process the replies to my queries
  ///////////////////////////////////////////////////////////////////////

  for (teuchos_size_t i=0; i < numberOfUniqueGids; i++){

    double key(uniqueGidQueries[i]);
    Array<LNO> v(uniqueGidQueryIndices[i]);

    int gidProc = procInBuf[gidLocation[i]];

    LNO gno;
    if (needGnoInfo){
      gno = gnoInBuf[gidLocation[i]];
    }
    else{
      gno = static_cast<GNO>(IdentifierTraits<AppGID>::keyToGid(key));
    }

    for (teuchos_size_t j=0; j < v.size(); j++){
      out_proc[v[j]] = gidProc;
      if (!skipGno)
        out_gno[v[j]] = gno;
    }
  }
}


template<typename AppLID, typename AppGID, typename LNO, typename GNO> 
  void IdentifierMap<AppLID,AppGID,LNO,GNO>::setupMap(void)
{
  numProcs_ = comm_->getSize(); 
  myRank_ = comm_->getRank(); 

  Z2_GLOBAL_INPUT_ASSERTION( *comm_, *env_, 
           "application global ID type is not supported",
           IdentifierTraits<AppGID>::is_valid_id_type() == true,
           BASIC_ASSERTION);

  localNumberOfIds_ = myGids_.size();

  typedef typename Array<GNO>::size_type teuchos_size_t;
  typedef typename Teuchos::Hashtable<double, LNO> id2index_hash_t;

  Teuchos::Tuple<teuchos_size_t, 4> counts;

  counts[0] = myLids_.size();
  counts[1] = localNumberOfIds_;
  counts[2] = counts[3] = 0;

  try{
    Teuchos::reduceAll<int, teuchos_size_t>(*comm_, Teuchos::REDUCE_SUM, 
      2, counts.getRawPtr(), counts.getRawPtr()+2);
  } 
  catch (const std::exception &e) {
    Z2_THROW_OUTSIDE_ERROR(*env_, e);
  }

  haveLocalIds_ = (counts[2] > 0);
  globalNumberOfIds_ = counts[3];

  Z2_GLOBAL_INPUT_ASSERTION( *comm_, *env_, 
       "number of global IDs does not equal number of local IDs",
      !haveLocalIds_ || (counts[0] == localNumberOfIds_),
       BASIC_ASSERTION);

  if (haveLocalIds_){   // hash LID to index in LID vector
    id2index_hash_t *p = NULL;
    if (localNumberOfIds_){
      try{
        p = new id2index_hash_t(localNumberOfIds_);
      }
      catch (const std::exception &e) 
        Z2_LOCAL_MEMORY_ASSERTION(*comm_, *env_, localNumberOfIds_, false); 
    }

    AppLID *lidPtr = myLids_.get();  // for performance

    for (teuchos_size_t i=0; i < localNumberOfIds_; i++){
      try{
        p->put(IdentifierTraits<AppLID>::key(lidPtr[i]), i);
      }
      catch (const std::exception &e) 
        Z2_THROW_OUTSIDE_ERROR(*env_, e);
    }

    lidHash_ = RCP<id2index_hash_t>(p);
  }

  // If the application's global ID data type (AppGID) is a Teuchos Ordinal,
  // we will be able to use it as our internal global numbers (GNO).  
  // Otherwise we will have to map their global IDs to valid global numbers.

  if (IdentifierTraits<AppGID>::isGlobalOrdinalType()){

    // Are the AppGIDs consecutive and increasing with process rank? 
    // If so GID/proc lookups can be optimized.

    AppGID min(0), max(0), globalMin(0), globalMax(0);
    AppGID *gidPtr = myGids_.get();  // for performance
    min = max = gidPtr[0];
    AppGID checkVal = min;
    bool consecutive = true;
  
    for (teuchos_size_t i=1; i < localNumberOfIds_; i++){
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
      Teuchos::reduceAll<int, AppGID>(*comm_, Teuchos::REDUCE_MIN, 2, 
        results.getRawPtr(), results.getRawPtr()+2);
    }
    catch (const std::exception &e) {
      Z2_THROW_OUTSIDE_ERROR(*env_, e);
    }

    if (results[2] != 1)       // min of consecutive flags
      consecutive=false;

    if (consecutive){
      globalMin = results[3];
      try{
        Teuchos::reduceAll<int, AppGID>(*comm_, Teuchos::REDUCE_MAX, 1, &max, 
          &globalMax);
      }
      catch (const std::exception &e) {
        Z2_THROW_OUTSIDE_ERROR(*env_, e);
      }

      if (globalMax - globalMin + 1 != static_cast<AppGID>(globalNumberOfIds_))
        consecutive = false;   // there are gaps in the gids
  
      if (consecutive){
        AppGID myStart = myGids_[0];

        gnoDist_ = ArrayRCP<GNO>(numProcs_ + 1);

        AppGID *startGID = static_cast<AppGID *>(gnoDist_.getRawPtr());
      
        try{
          Teuchos::gatherAll<int, AppGID>(*comm_, 1, &myStart, numProcs_, 
            startGID);
        }
        catch (const std::exception &e) {
          Z2_THROW_OUTSIDE_ERROR(*env_, e);
        }
      
        for (int p=1; p < numProcs_; p++){
          if (startGID[p] < startGID[p-1]){
            consecutive = false;  // gids do not increase with process rank
            break;
          }
        }
        if (consecutive){
          startGID[numProcs_] = globalMax + 1;
        }
      }
    }
  }
  else{

    // AppGIDs are not Ordinals.  We map them to consecutive 
    // global numbers starting with 0. 

    try{
      gnoDist_ = ArrayRCP<GNO>(numProcs_ + 1, 0);
    }
    catch (const std::exception &e) {
      Z2_THROW_OUTSIDE_ERROR(*env_, e);
    }

    GNO myNum = static_cast<GNO>(localNumberOfIds_);

    try{
      Teuchos::scan<int, GNO>(*comm_, Teuchos::REDUCE_SUM, 1, &myNum, 
        gnoDist_.getRawPtr() + 1);
    }
    catch (const std::exception &e) {
      Z2_THROW_OUTSIDE_ERROR(*env_, e);
    }
  }

  if (gnoDist_.size() == 0){

    // We need a hash table mapping the application global ID
    // to its index in myGids_.

    id2index_hash_t *p = NULL;
    if (localNumberOfIds_){
      try{
        p = new id2index_hash_t(localNumberOfIds_);
      }
      catch (const std::exception &e) 
        Z2_LOCAL_MEMORY_ASSERTION(*comm_, *env_, localNumberOfIds_, false); 
    }

    AppGID *gidPtr = myGids_.get();  // for performance

    for (teuchos_size_t i=0; i < localNumberOfIds_; i++){
      try{
        p->put(IdentifierTraits<AppGID>::key(gidPtr[i]), i);
      }
      catch (const std::exception &e) {
        Z2_THROW_OUTSIDE_ERROR(*env_, e);
      }
    }
    gidHash_ = RCP<id2index_hash_t>(p);
  }
}

}   // end namespace Zoltan2

#endif /* _ZOLTAN2_IDENTIFIERMAP_HPP_ */
