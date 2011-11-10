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
#include <Zoltan2_AlltoAll.hpp>

namespace Zoltan2
{
/*! \brief Identifier translations available from IdentifierMap.

  TRANSLATE_APP_TO_LIB  Translate the applications identifiers to
                          Zoltan's internal identifiers

  TRANSLATE_LIB_TO_APP  Translate Zoltan's internal identifiers back
                          to the application's original identifiers
*/

enum TranslationType {
  TRANSLATE_APP_TO_LIB,
  TRANSLATE_LIB_TO_APP
};

/*! Z2::IdentifierMap
    \brief An IdentifierMap manages a global space of unique object identifiers.

    TODO - trim down comments and code once we get it all straight
           replace new/delete with memory wrappers
           exception handling
           use Kokkos node 

    LID  is the data type used by application for local IDs, which are optional.
    GID  is the data type used by application for globls IDs
    LNO  is the data type used by Zoltan2 for local counts and indexes.
    GNO  is the integral data type used by Zoltan2 for global counts.
*/

////////////////////////////////////////////////////////////////////
// Declarations
////////////////////////////////////////////////////////////////////

template<typename LID, typename GID, typename LNO, typename GNO>
    class IdentifierMap{

public:

  /*! Constructor - Must be called by all processes 
   *
   * \param comm the problem communicator
   * \param env  the problem and library environment
   * \param gids  the application global IDs
   * \param lids  the application local IDs, if any
   *    TODO - provide some way to indicate that local IDs
   *              are consecutive ordinals, and omit lids array
   * \param gidsMustBeConsecutive  set to true if the algorithm
   *           or third party library requires consective ids
   *           If necessary the IdentifierMap will map the application's
   *           global IDs to integral IDs beginning at zero.
   */

  typedef LNO    lno_t;
  typedef GNO    gno_t;
  typedef LID    lid_t;
  typedef GID    gid_t;

  explicit IdentifierMap( const RCP<const Comm<int> > &comm, 
                          const RCP<Environment > &env, 
                          const ArrayRCP<GID> &gids, 
                          const ArrayRCP<LID> &lids,
                          bool gidsMustBeConsecutive=false);

  /*! Destructor */
  ~IdentifierMap();

  /*! Copy Constructor */
  IdentifierMap(const IdentifierMap &id);

  /*! Assignment operator */
  IdentifierMap &operator=(const IdentifierMap &id);

  /*! Return true if we are using the application global IDs 
   *  for our internal global numbers 
   */
  bool gnosAreGids() const;

  /*! Return true if our internal global numbers are consecutive.
   */
  bool gnosAreConsecutive() const;

  /*! Return true if consecutive Gids are required.
   */
  bool consecutiveGnosAreRequired() const;

  /*! Return the minimum Zoltan2 global Id across all processes.
   */
  GNO getMinimumGlobalId() const;

  /*! Return the maximum Zoltan2 global Id across all processes.
   */
  GNO getMaximumGlobalId() const;

  /*! Map the application global IDs to internal global numbers or vice versa.

      \param gid an array of caller's global IDs
      \param gno an array of Zoltan2 global numbers
      \param tt should be TRANSLATE_APP_TO_LIB or TRANSLATE_LIB_TO_APP

      This is a local call.  If gid is a vector of application global IDs, then
      gno will be set to the corresponding internal global numbers.  If the
      gno vector contains internal global numbers, they will be translated
      to application global IDs.  The application global IDs must be from
      those supplied by this process.

      TODO: Why does code fail to compile if I pass gid and gno by reference?
   */
  void gidTranslate(ArrayView<GID> gid, 
                    ArrayView<GNO> gno,
                    TranslationType tt) const;

  /*! Map application local IDs to internal global numbers or vice versa.

      \param lid an array of caller's local IDs
      \param gno an array of Zoltan2 global numbers
      \param tt should be TRANSLATE_APP_TO_LIB or TRANSLATE_LIB_TO_APP

      This is a local call.  If lid is a vector of application local IDs, then
      gno will be set to the corresponding internal global numbers.  If the
      gno vector contains internal global numbers, they will be translated
      to application local IDs.  The application local IDs must be from
      those supplied by this process.
   */
  void lidTranslate(ArrayView<LID> lid, 
                    ArrayView<GNO> gno,
                    TranslationType tt) const;

  /*! Map application global IDs to internal global numbers.

      \param in_gid input, an array of the global IDs
      \param out_gno output, an optional array of the corresponding 
          global numbers used by Zoltan2.  If out_gno.size() is zero,
          we assume global numbers are not needed.
      \param out_proc output, an array of the process ranks corresponding with
                         the in_gid and out_gno, out_proc[i] is the process
                         that supplied in_gid[i] to Zoltan2.

      All processes must call this.  The global IDs 
      supplied may belong to another process.  
   */
  void gidGlobalTranslate( ArrayView<const GID> in_gid,
                           ArrayView<GNO> out_gno,
                           ArrayView<int> out_proc) const;

private:

  // Input communicator

  const RCP<const Comm<int> > comm_;

  // Problem parameters, library configuration.

  const RCP<Environment> env_;

  // Application global and local IDs

  const ArrayRCP<GID> myGids_; 
  const ArrayRCP<LID> myLids_;

  // Zoltan2 GNOs will be consecutive if the application GIDs
  // were mapped to GNOs, or if the application GIDs happen
  // to be consecutive ordinals.  In either case, gnoDist_[p]
  // is the first GNO on process p.  gnoDist_[numProcs_] is the
  // global total of GNOs.

  ArrayRCP<GNO> gnoDist_;

  // If application GIDs are not consecutive ordinals, the gidHash_
  // table maps a key (a double generated from the GID) to the
  // location of the GID in the myGids_ list.

  RCP<Teuchos::Hashtable<double, LNO> >  gidHash_;

  // A hash table from application local ID key to our local index.
  // TODO: just save min LID if LIDs are consecutive ordinals

  RCP<Teuchos::Hashtable<double, LNO> >  lidHash_;

  typename Array<GNO>::size_type globalNumberOfIds_;
  typename Array<GNO>::size_type localNumberOfIds_;
  bool haveLocalIds_;
  int myRank_;
  int numProcs_;

  bool userGidsAreTeuchosOrdinal_;
  bool userGidsAreConsecutive_;
  bool userGidsAreZoltan2Gids_;
  bool zoltan2GidsAreConsecutive_;
  // TODO it should be possible to specify the base too
  bool consecutiveGidsAreRequired_;

  GNO minGlobalGno_;
  GNO maxGlobalGno_;

  void setupMap();
};

////////////////////////////////////////////////////////////////////
// Definitions
////////////////////////////////////////////////////////////////////

template<typename LID, typename GID, typename LNO, typename GNO> 
  IdentifierMap<LID,GID,LNO,GNO>::IdentifierMap(
    const RCP<const Comm<int> > &incomm, const RCP<Environment> &env,
    const ArrayRCP<GID> &gids, const ArrayRCP<LID> &lids,
    bool idsMustBeConsecutive) 
         : comm_(incomm),  env_(env), myGids_(gids), myLids_(lids),
           gnoDist_(), gidHash_(), lidHash_(), 
           globalNumberOfIds_(0), localNumberOfIds_(0), haveLocalIds_(false),
           myRank_(incomm->getRank()), numProcs_(incomm->getSize()), 
           userGidsAreTeuchosOrdinal_(false), userGidsAreConsecutive_(false), 
           userGidsAreZoltan2Gids_(false), zoltan2GidsAreConsecutive_(false), 
           consecutiveGidsAreRequired_(idsMustBeConsecutive),
           minGlobalGno_(0), maxGlobalGno_(0)
{
  setupMap();
}

  /*! Destructor */
template<typename LID, typename GID, typename LNO, typename GNO>
  IdentifierMap<LID,GID,LNO,GNO>::~IdentifierMap() 
  {
  }

#if 0
  /*! Copy Constructor */
template<typename LID, typename GID, typename LNO, typename GNO>
  IdentifierMap<LID,GID,LNO,GNO::IdentifierMap(const IdentifierMap &id)
{
    //TODO    default should work, shouldn't it?
}

  /*! Assignment operator */
template<typename LID, typename GID, typename LNO, typename GNO>
  IdentifierMap<LID,GID,LNO,GNO &
    IdentifierMapLID,GID,LNO,GNO>::operator=(const IdentifierMap<User> &id)
{
    //TODO
}
#endif

template<typename LID, typename GID, typename LNO, typename GNO>
  bool IdentifierMap<LID,GID,LNO,GNO>::gnosAreGids() const
{
  return userGidsAreZoltan2Gids_;
}

template<typename LID, typename GID, typename LNO, typename GNO>
  bool IdentifierMap<LID,GID,LNO,GNO>::gnosAreConsecutive() const
{
  return zoltan2GidsAreConsecutive_;
}

template<typename LID, typename GID, typename LNO, typename GNO>
  bool IdentifierMap<LID,GID,LNO,GNO>::consecutiveGnosAreRequired() const
{
  return consecutiveGidsAreRequired_;
}

template<typename LID, typename GID, typename LNO, typename GNO>
  GNO IdentifierMap<LID,GID,LNO,GNO>::getMinimumGlobalId() const
{
  return minGlobalGno_;
}

template<typename LID, typename GID, typename LNO, typename GNO>
  GNO IdentifierMap<LID,GID,LNO,GNO>::getMaximumGlobalId() const
{
  return maxGlobalGno_;
}

template<typename LID, typename GID, typename LNO, typename GNO>
  void IdentifierMap<LID,GID,LNO,GNO>::gidTranslate(
    ArrayView<GID> gid, 
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
    (tt==TRANSLATE_APP_TO_LIB) || (tt==TRANSLATE_LIB_TO_APP), 
    BASIC_ASSERTION);

  Z2_LOCAL_INPUT_ASSERTION(*comm_, *env_, 
    "Destination array is too small",
    ((tt==TRANSLATE_LIB_TO_APP) && (gid.size() >= gno.size())) || 
     ((tt==TRANSLATE_APP_TO_LIB) && (gno.size() >= gid.size())),
    BASIC_ASSERTION);

  if (userGidsAreZoltan2Gids_){   // our gnos are the app gids
    if (tt == TRANSLATE_LIB_TO_APP)
      for (teuchos_size_t i=0; i < len; i++)
        IdentifierTraits<GNO>::castTo(gno[i], gid[i]);
    else
      for (teuchos_size_t i=0; i < len; i++)
        IdentifierTraits<GID>::castTo(gid[i], gno[i]);
  }
  else{              // we mapped gids to consecutive gnos
    GNO firstGno = gnoDist_[myRank_];
    GNO endGno = gnoDist_[myRank_ + 1];

    if (tt == TRANSLATE_LIB_TO_APP){
      for (teuchos_size_t i=0; i < len; i++){

        Z2_LOCAL_INPUT_ASSERTION(*comm_, *env_, "invalid global number", 
        (gno[i] >= firstGno) && (gno[i] < endGno), BASIC_ASSERTION);

        gid[i] = myGids_[gno[i] - firstGno];
      }
    }
    else{
      LNO idx;
      for (teuchos_size_t i=0; i < len; i++){
        try{
          double key = Zoltan2::IdentifierTraits<GID>::key(gid[i]);
          idx = gidHash_->get(key);
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

template<typename LID, typename GID, typename LNO, typename GNO>
  void IdentifierMap<LID,GID,LNO,GNO>::lidTranslate(
    ArrayView<LID> lid, 
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
    (tt==TRANSLATE_LIB_TO_APP) || (tt==TRANSLATE_APP_TO_LIB), 
    BASIC_ASSERTION);

  Z2_LOCAL_INPUT_ASSERTION(*comm_, *env_, 
    "Destination array is too small",
    ((tt==TRANSLATE_LIB_TO_APP) && (lid.size() >= gno.size())) || ((tt==TRANSLATE_APP_TO_LIB) && (gno.size() >= lid.size())),
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
  
  if (tt == TRANSLATE_LIB_TO_APP){
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
          GID keyArg;
          IdentifierTraits<GNO>::castTo(gno[i], keyArg);
          idx = gidHash_->get(IdentifierTraits<GID>::key(keyArg));
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
        idx = lidHash_->get(IdentifierTraits<LID>::key(lid[i]));
      }
      catch (const std::exception &e) {
        Z2_THROW_OUTSIDE_ERROR(*env_, e);
      }

      if (gnoDist_.size() > 0)  // gnos are consecutive
        gno[i] = firstGno + idx;
      else                     // gnos must be the app gids
        IdentifierTraits<GID>::castTo(myGids_[idx], gno[i]);
    }
  }
}

template<typename LID, typename GID, typename LNO, typename GNO>
  void IdentifierMap<LID,GID,LNO,GNO>::gidGlobalTranslate(
    ArrayView<const GID> in_gid,
    ArrayView<GNO> out_gno,
    ArrayView<int> out_proc) const
{
  typedef typename Array<GNO>::size_type teuchos_size_t;
  typedef typename Teuchos::Hashtable<double, LNO> id2index_hash_t;

  teuchos_size_t len=in_gid.size();

  if (len == 0){
    return;
  }

  bool skipGno = (out_gno.size() == 0);

  Z2_LOCAL_INPUT_ASSERTION(*comm_, *env_, 
    "Destination array is too small", 
    (out_proc.size() >= len) && (skipGno || (out_gno.size() >= len)),
    BASIC_ASSERTION);

  if (userGidsAreZoltan2Gids_ && (gnoDist_.size() > 0)){

    // Easy case - communication is not needed.
    // Global numbers are the application global IDs and
    // they are increasing consecutively with rank.
 
    typename std::map<GNO, int> firstGnoToProc;
    typename std::map<GNO, int>::iterator pos;

    for (int p=0; p <= numProcs_; p++){
      firstGnoToProc[gnoDist_[p]] = p;
    }

    for (teuchos_size_t i=0; i < len; i++){
      GNO globalNumber;
      IdentifierTraits<GID>::castTo(in_gid[i], globalNumber);
      if (!skipGno)
        out_gno[i] = globalNumber;
      pos = firstGnoToProc.upper_bound(globalNumber);
      out_proc[i] = pos->second - 1;
    }

    return;
  }

  bool needGnoInfo = !userGidsAreZoltan2Gids_;

  ///////////////////////////////////////////////////////////////////////
  // First: Hash each of my GIDs to a process that will answer
  // for it.  Send my GIDs (and the Gnos if they are different)
  // to their assigned processes.  Build a search structure for
  // the GIDs that were assigned to me, so I can reply with
  // with the process owning them (and their Gnos if they are different).
  ///////////////////////////////////////////////////////////////////////

  Array<int> hashProc;
  Array<GID> gidOutBuf;
  Array<GNO> gnoOutBuf;
  Array<LNO> countOutBuf(numProcs_, 0);
  Array<LNO> offsetBuf(numProcs_ + 1, 0);

  ArrayRCP<GID> gidInBuf;
  ArrayRCP<GNO> gnoInBuf;
  ArrayRCP<LNO> countInBuf;

  if (localNumberOfIds_ > 0){

    try{ 
      hashProc.resize(localNumberOfIds_, 0);
    }
    catch(...)
      Z2_LOCAL_MEMORY_ASSERTION(*comm_, *env_, localNumberOfIds_, false); 

    try{ 
      gidOutBuf.resize(localNumberOfIds_); 
    }
    catch(...)
      Z2_LOCAL_MEMORY_ASSERTION(*comm_, *env_, localNumberOfIds_, false); 

    for (teuchos_size_t i=0; i < localNumberOfIds_; i++){
      hashProc[i] = IdentifierTraits<GID>::hashCode(myGids_[i]) % numProcs_;
      countOutBuf[hashProc[i]]++;
    }
  
    for (int p=1; p <= numProcs_; p++){
      offsetBuf[p] = offsetBuf[p-1] + countOutBuf[p-1];
    }
  
    if (needGnoInfo){   
      // The gnos are not the gids, which also implies that
      // gnos are consecutive numbers given by gnoDist_.
      gnoOutBuf.resize(localNumberOfIds_, 0);
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

  // Teuchos comment #1: The Array::() operator returns an ArrayView.
  // Teuchos comment #2: GID may not be a Teuchos Packet type,
  //                     so we wrote our own AlltoAllv.
  // Z2::AlltoAllv comment: Buffers are in process rank contiguous order.

  try{
    ArrayView<const GID> gidView = gidOutBuf();
    ArrayView<const LNO> countView = countOutBuf();
    AlltoAllv(*comm_, *env_, gidView, countView, gidInBuf, countInBuf);
  }
  catch (const std::exception &e)
    Z2_THROW_ZOLTAN2_ERROR(*env_, e);

  gidOutBuf.clear();
  
  if (needGnoInfo){
    countInBuf.release();
    ArrayView<const GNO> gnoView = gnoOutBuf();
    ArrayView<const LNO> countView = countOutBuf();
    try{
      AlltoAllv(*comm_, *env_, gnoView, countView, gnoInBuf, countInBuf);
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
  typename std::map<LNO, int>::iterator indexProcCurr;
  LNO indexTotal = 0;

  for (int p=0; p < numProcs_; p++){
    firstIndexToProc[indexTotal] = p;
    indexTotal += countInBuf[p];
  }

  firstIndexToProc[indexTotal] = numProcs_;

  id2index_hash_t gidToIndex(indexTotal);

  LNO total = 0;
  for (int p=0; p < numProcs_; p++){
    for (LNO i=0; i < countInBuf[p]; i++, total++){
      try{
        gidToIndex.put(IdentifierTraits<GID>::key(gidInBuf[total]), total);
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

  
  Array<double> uniqueGidQueries;
  Array<Array<LNO> > uniqueGidQueryIndices;
  teuchos_size_t numberOfUniqueGids = 0;
  Array<LNO> gidLocation;

  countOutBuf.resize(numProcs_, 0);

  std::map<double, std::vector<LNO> > gidIndices;
  typename std::map<double, std::vector<LNO> >::iterator next; 

  if (len > 0){
    for (LNO i=0; i < len; i++){
      double uniqueKey(IdentifierTraits<GID>::key(in_gid[i]));
      next = gidIndices.find(uniqueKey);
      if (next == gidIndices.end()){
        std::vector<LNO> v(1, i);
        gidIndices[uniqueKey] = v;
      }
      else{
        std::vector<LNO> &v = next->second;
        v.push_back(i);
      }
    }
  
    numberOfUniqueGids = gidIndices.size();
  
    try{ 
      gidOutBuf.resize(numberOfUniqueGids); 
    }
    catch(...)
      Z2_LOCAL_MEMORY_ASSERTION(*comm_, *env_, numberOfUniqueGids, false); 

    try{ 
      hashProc.resize(numberOfUniqueGids, 0);
    }
    catch(...)
      Z2_LOCAL_MEMORY_ASSERTION(*comm_, *env_, numberOfUniqueGids, false); 
  
    LNO idx = 0;
    for (next = gidIndices.begin(); next != gidIndices.end(); ++next, ++idx){
      double key = next->first;
      GID gid = IdentifierTraits<GID>::keyToGid(key);
      hashProc[idx] = IdentifierTraits<GID>::hashCode(gid) % numProcs_;
      countOutBuf[hashProc[idx]]++;
    }
  
    offsetBuf[0] = 0;
  
    for (int p=0; p < numProcs_; p++){
      offsetBuf[p+1] = offsetBuf[p] + countOutBuf[p];
    }
  
    try{ 
      gidLocation.resize(numberOfUniqueGids, 0);
    }
    catch(...)
      Z2_LOCAL_MEMORY_ASSERTION(*comm_, *env_, numberOfUniqueGids, false); 
  
    idx = 0;
    for (next = gidIndices.begin(); next != gidIndices.end(); ++next, ++idx){
      double key = next->first;
      GID gid = IdentifierTraits<GID>::keyToGid(key);
      gidLocation[idx] = offsetBuf[hashProc[idx]];
      gidOutBuf[gidLocation[idx]] = gid;
      offsetBuf[hashProc[idx]] = gidLocation[idx] + 1;
    }

    hashProc.clear();
  }

  try{
    ArrayView<const GID> gidView = gidOutBuf();
    ArrayView<const LNO> countView = countOutBuf();
    AlltoAllv(*comm_, *env_, gidView, countView, gidInBuf, countInBuf);
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
  ArrayRCP<int> procInBuf;

  if (needGnoInfo){
    try{ 
      gnoOutBuf.resize(total, 0);
    }
    catch(...)
      Z2_LOCAL_MEMORY_ASSERTION(*comm_, *env_, total, false); 
  }

  if (total > 0){
  
    total=0;
  
    for (int p=0; p < numProcs_; p++){
      for (LNO i=0; i < countInBuf[p]; i++, total++){
        double k(IdentifierTraits<GID>::key(gidInBuf[total]));
        LNO index(0);
        try{
          index = gidToIndex.get(k);
        }
        catch (const std::exception &e) {
          Z2_THROW_OUTSIDE_ERROR(*env_, e);
        }

        Z2_LOCAL_BUG_ASSERTION(*comm_, *env_, "gidToIndex table", 
          (index >= 0)&&(index<=indexTotal), BASIC_ASSERTION);
        
        indexProcCurr = firstIndexToProc.upper_bound(index);
        procOutBuf[total] = indexProcCurr->second-1;
  
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
    ArrayView<const int> procView = procOutBuf();
    ArrayView<const LNO> countView = countOutBuf();
    AlltoAllv(*comm_, *env_, procView, countView, procInBuf, countInBuf);
  }
  catch (const std::exception &e)
    Z2_THROW_ZOLTAN2_ERROR(*env_, e);

  procOutBuf.clear();

  if (needGnoInfo){
    try{
      ArrayView<const GNO> gnoView = gnoOutBuf();
      ArrayView<const LNO> countView = countOutBuf();
      AlltoAllv(*comm_, *env_, gnoView, countView, gnoInBuf, countInBuf);
    }
    catch (const std::exception &e)
      Z2_THROW_ZOLTAN2_ERROR(*env_, e);

    gnoOutBuf.clear();
  }

  countOutBuf.clear();

  ///////////////////////////////////////////////////////////////////////
  // Done.  Process the replies to my queries
  ///////////////////////////////////////////////////////////////////////

  LNO idx = 0;
  for (next = gidIndices.begin(); next != gidIndices.end(); ++next, ++idx){
    double key = next->first;
    std::vector<LNO> &v = next->second;
    int gidProc = procInBuf[gidLocation[idx]];

    GNO gno;
    if (needGnoInfo){
      gno = gnoInBuf[gidLocation[idx]];
    }
    else{
      GID gid = IdentifierTraits<GID>::keyToGid(key);
      IdentifierTraits<GID>::castTo(gid, gno);
    }

    for (teuchos_size_t j=0; j < v.size(); j++){
      out_proc[v[j]] = gidProc;
      if (!skipGno)
        out_gno[v[j]] = gno;
    }
  }
}


template<typename LID, typename GID, typename LNO, typename GNO> 
  void IdentifierMap<LID,GID,LNO,GNO>::setupMap(void)
{
  numProcs_ = comm_->getSize(); 
  myRank_ = comm_->getRank(); 

  Z2_GLOBAL_INPUT_ASSERTION( *comm_, *env_, 
           "application global ID type is not supported yet",
           IdentifierTraits<GID>::is_valid_id_type() == true,
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
    // TODO test to see if LIDs are consecutive, to avoid need for
    // hash for LID to its location in array.
    id2index_hash_t *p = NULL;
    if (localNumberOfIds_){
      try{
        p = new id2index_hash_t(localNumberOfIds_);
      }
      catch (const std::exception &e) 
        Z2_LOCAL_MEMORY_ASSERTION(*comm_, *env_, localNumberOfIds_, false); 
    }

    LID *lidPtr = myLids_.get();  // for performance

    for (teuchos_size_t i=0; i < localNumberOfIds_; i++){
      try{
        p->put(IdentifierTraits<LID>::key(lidPtr[i]), i);
      }
      catch (const std::exception &e) 
        Z2_THROW_OUTSIDE_ERROR(*env_, e);
    }

    lidHash_ = RCP<id2index_hash_t>(p);
  }

  // Determine whether the application's global ID data type (GID) is a 
  // Teuchos Ordinal, and if so, whether they are consecutive, and what the
  // GID base is.
  //
  // If they are non-consecutive Teuchos Ordinals, and consecutive IDs are
  // not required, or if they are consecutive Teuchos Ordinals, we can use 
  // the user's global IDs.  Otherwise we map their IDs to valid global numbers.

  userGidsAreTeuchosOrdinal_ = false;
  userGidsAreConsecutive_ = false;
  userGidsAreZoltan2Gids_ = false;

  if (IdentifierTraits<GID>::isGlobalOrdinal()){

    userGidsAreTeuchosOrdinal_ = true;
    GID *gidPtr = myGids_.get();

    GID min = gidPtr[0];
    GID max = gidPtr[localNumberOfIds_-1];

    std::pair<GID, GID> globalMinMax = 
      IdentifierTraits<GID>::globalMinMax(min, max, *comm_);
  
    GNO nMin, nMax;
    IdentifierTraits<GID>::castTo(globalMinMax.first, nMin);
    IdentifierTraits<GID>::castTo(globalMinMax.second, nMax);

    minGlobalGno_ = nMin;  // We'll overwrite these if
    maxGlobalGno_ = nMax;  // not userGidsAreZoltan2Gids_

    // Are the gids consecutive and increasing with process rank? 
    // If so GID/proc lookups can be optimized.

    bool locallyConsecutiveIncreasing = 
      IdentifierTraits<GID>::areConsecutive(gidPtr, localNumberOfIds_);
  
    int localFlag = (locallyConsecutiveIncreasing ? 1 : 0);
    int globalFlag = 0;

    try{
      Teuchos::reduceAll<int, int>(*comm_, Teuchos::REDUCE_MIN, 1, 
        &localFlag, &globalFlag);
    }
    catch (const std::exception &e) {
      Z2_THROW_OUTSIDE_ERROR(*env_, e);
    }

    if (globalFlag == 1){

      if (nMax - nMin + 1 == globalNumberOfIds_){

        Array<GID> sendBuf(numProcs_, min);
        ArrayRCP<GID> recvBuf;
  
        // We use Zoltan2::AlltoAll because GID may 
        // not be a Teuchos Packet type at compile time.
  
        try{
          AlltoAll<GID, LNO>(*comm_, *env_, sendBuf, LNO(1), recvBuf);
        }
        catch(std::exception &e){
          Z2_THROW_ZOLTAN2_ERROR(*env_, e);
        }
    
        userGidsAreConsecutive_ = true;

        for (int i=1; i < numProcs_; i++){
          if (IdentifierTraits<GID>::lessThan(recvBuf[i-1], recvBuf[i]))
            continue;
          userGidsAreConsecutive_ = false;
          break; 
        }

        if (userGidsAreConsecutive_){
          gnoDist_.resize(numProcs_ + 1, 0);
          for (int i=0; i < numProcs_; i++){
            IdentifierTraits<GID>::castTo(recvBuf[i], gnoDist_[i]);
          }
          gnoDist_[numProcs_] = globalNumberOfIds_ + gnoDist_[0];
        }
      }
    }
  }

  if (userGidsAreTeuchosOrdinal_ &&
      ((consecutiveGidsAreRequired_ && userGidsAreConsecutive_) ||
      !consecutiveGidsAreRequired_))
  {
    // We can use the application's global IDs
    userGidsAreZoltan2Gids_ = true;
    zoltan2GidsAreConsecutive_ = userGidsAreConsecutive_;
  }

  if (!userGidsAreZoltan2Gids_){
    // We map application gids to consecutive global numbers starting with 0. 

    zoltan2GidsAreConsecutive_ = true;

    try{
      gnoDist_.resize(numProcs_ + 1, 0);
    }
    catch (const std::exception &e) {
      Z2_THROW_OUTSIDE_ERROR(*env_, e);
    }

    GNO myNum = static_cast<GNO>(localNumberOfIds_);

    try{
      GNO *p = gnoDist_.getRawPtr();
      Teuchos::gatherAll<int, GNO>(*comm_, 1, &myNum, numProcs_, p+1);
    }
    catch (const std::exception &e) {
      Z2_THROW_OUTSIDE_ERROR(*env_, e);
    }

    for (int i=2; i <= numProcs_; i++){
      gnoDist_[i] += gnoDist_[i-1];
    }
  }

  if (!userGidsAreConsecutive_){

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

    GID *gidPtr = myGids_.get();  // for performance

    for (teuchos_size_t i=0; i < localNumberOfIds_; i++){
      try{
        // CRASH here when GIDs are std::pair<int,int> TODO
        p->put(IdentifierTraits<GID>::key(gidPtr[i]), i);
      }
      catch (const std::exception &e) {
        Z2_THROW_OUTSIDE_ERROR(*env_, e);
      }
    }
    gidHash_ = RCP<id2index_hash_t>(p);
  }

  if (gnoDist_.size() > 0){
    minGlobalGno_ = gnoDist_[0];
    maxGlobalGno_ = gnoDist_[numProcs_]-1;
  }
}

}   // end namespace Zoltan2

#endif /* _ZOLTAN2_IDENTIFIERMAP_HPP_ */
