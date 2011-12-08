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
#include <Zoltan2_IdentifierTraits.hpp>
#include <Zoltan2_AlltoAll.hpp>

#include <vector>
#include <map>

#include <Teuchos_as.hpp>
#include <Teuchos_Hashtable.hpp>

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

/*! Zoltan2::IdentifierMap
    \brief An IdentifierMap manages a global space of unique object identifiers.

    GID  is the data type used by application for globls IDs
    LNO  is the data type used by Zoltan2 for local counts and indexes.
    GNO  is the integral data type used by Zoltan2 for global counts.
*/

////////////////////////////////////////////////////////////////////
// Declarations
////////////////////////////////////////////////////////////////////

template<typename GID, typename LNO, typename GNO>
    class IdentifierMap{

public:

  /*! Constructor - Must be called by all processes 
   *
   * \param comm the problem communicator
   * \param env  the problem and library environment
   * \param gids  the application global IDs
   * \param gidsMustBeConsecutive  set to true if the algorithm
   *           or third party library requires consective ids
   *           If necessary the IdentifierMap will map the application's
   *           global IDs to consecutive integral IDs beginning at zero.
   */

  typedef LNO    lno_t;
  typedef GNO    gno_t;
  typedef GID    gid_t;

  explicit IdentifierMap( const RCP<const Environment > &env, 
                          const ArrayRCP<const GID> &gids, 
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
   */
  void gidTranslate(ArrayView<GID> gid, 
                    ArrayView<GNO> gno,
                    TranslationType tt) const;

  /*! Map application indices to internal global numbers or vice versa.

      \param lno an array of indices
      \param gno an array of Zoltan2 global numbers
      \param tt should be TRANSLATE_APP_TO_LIB or TRANSLATE_LIB_TO_APP

      This is a local call. 

      If gno contains a list of Zoltan2 internal
      global numbers, then lno[i] on return will be the index at which 
      the application's global ID associated with gno[i] appeared in the 
      input adapter's "get" method list. In this case, lno should be
      pre-allocated by the caller, and tt should be TRANSLATE_LIB_TO_APP.

      Similarly, if lno contains a list of indices ranging from 0 to N-1
      (where N would be the local number of objects), then gno[i] on
      return will be the internal global number associated with the
      application's global ID that appeared in location lno[i] in the 
      input adapter's "get" method list. In this case, gno should be
      pre-allocated by the caller and tt should be TRANSLATE_APP_TO_LIB.
   */
  void lnoTranslate(ArrayView<LNO> lno, 
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

  const RCP<const Environment> env_;

  // Application global IDs

  const ArrayRCP<const GID> myGids_; 

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

  global_size_t globalNumberOfIds_;
  size_t localNumberOfIds_;
  int myRank_;
  int numProcs_;

  // By "Consecutive" we mean globally consecutive increasing
  // with process rank.

  bool userGidsAreTeuchosOrdinal_;
  bool userGidsAreConsecutive_;
  bool userGidsAreZoltan2Gnos_;
  bool zoltan2GnosAreConsecutive_;
  
  bool consecutiveGidsAreRequired_;

  GNO minGlobalGno_;
  GNO maxGlobalGno_;

  void setupMap();
};

////////////////////////////////////////////////////////////////////
// Definitions
////////////////////////////////////////////////////////////////////

template<typename GID, typename LNO, typename GNO> 
  IdentifierMap<GID,LNO,GNO>::IdentifierMap( const RCP<const Environment> &env,
    const ArrayRCP<const GID> &gids, bool idsMustBeConsecutive) 
         : comm_(env->comm_),  env_(env), myGids_(gids), gnoDist_(), gidHash_(),
           globalNumberOfIds_(0), localNumberOfIds_(0),
           myRank_(0), numProcs_(1),
           userGidsAreTeuchosOrdinal_(false), userGidsAreConsecutive_(false), 
           userGidsAreZoltan2Gnos_(false), zoltan2GnosAreConsecutive_(false), 
           consecutiveGidsAreRequired_(idsMustBeConsecutive),
           minGlobalGno_(0), maxGlobalGno_(0)
{
  myRank_ = comm_->getRank();
  numProcs_ = comm_->getSize();
  setupMap();
}

  /*! Destructor */
template<typename GID, typename LNO, typename GNO>
  IdentifierMap<GID,LNO,GNO>::~IdentifierMap() 
  {
  }

template< typename GID, typename LNO, typename GNO>
  bool IdentifierMap<GID,LNO,GNO>::gnosAreGids() const
{
  return userGidsAreZoltan2Gnos_;
}

template< typename GID, typename LNO, typename GNO>
  bool IdentifierMap<GID,LNO,GNO>::gnosAreConsecutive() const
{
  return zoltan2GnosAreConsecutive_;
}

template< typename GID, typename LNO, typename GNO>
  bool IdentifierMap<GID,LNO,GNO>::consecutiveGnosAreRequired() const
{
  return consecutiveGidsAreRequired_;
}

template< typename GID, typename LNO, typename GNO>
  GNO IdentifierMap<GID,LNO,GNO>::getMinimumGlobalId() const
{
  return minGlobalGno_;
}

template< typename GID, typename LNO, typename GNO>
  GNO IdentifierMap<GID,LNO,GNO>::getMaximumGlobalId() const
{
  return maxGlobalGno_;
}

template< typename GID, typename LNO, typename GNO>
  void IdentifierMap<GID,LNO,GNO>::gidTranslate(
    ArrayView<GID> gid, 
    ArrayView<GNO> gno,
    TranslationType tt) const
{
  size_t len=gid.size();

  if (len == 0){
    return;
  }

  Z2_LOCAL_INPUT_ASSERTION(*env_, "invalid TranslationType", 
    (tt==TRANSLATE_APP_TO_LIB) || (tt==TRANSLATE_LIB_TO_APP), 
    BASIC_ASSERTION);

  Z2_LOCAL_INPUT_ASSERTION(*env_, "Destination array is too small",
    ((tt==TRANSLATE_LIB_TO_APP) && (gid.size() >= gno.size())) || 
     ((tt==TRANSLATE_APP_TO_LIB) && (gno.size() >= gid.size())),
    BASIC_ASSERTION);

  if (userGidsAreZoltan2Gnos_){   // our gnos are the app gids
    if (tt == TRANSLATE_LIB_TO_APP)
      for (size_t i=0; i < len; i++)
        gid[i] = Teuchos::as<GID>(gno[i]);
    else
      for (size_t i=0; i < len; i++)
        gno[i] = Teuchos::as<GNO>(gid[i]);
  }
  else{              // we mapped gids to consecutive gnos
    GNO firstGno = gnoDist_[myRank_];
    GNO endGno = gnoDist_[myRank_ + 1];

    if (tt == TRANSLATE_LIB_TO_APP){
      for (size_t i=0; i < len; i++){

        Z2_LOCAL_INPUT_ASSERTION(*env_, "invalid global number", 
        (gno[i] >= firstGno) && (gno[i] < endGno), BASIC_ASSERTION);

        gid[i] = myGids_[gno[i] - firstGno];
      }
    }
    else{
      LNO idx;
      if (userGidsAreConsecutive_){
        for (size_t i=0; i < len; i++){
          gno[i] = firstGno + IdentifierTraits<GID>::difference(
            myGids_[0], gid[i]);
          Z2_LOCAL_INPUT_ASSERTION(*env_, "invalid global id", 
            (gno[i] >= firstGno) && (gno[i] < endGno), BASIC_ASSERTION);
        }
      }
      else{
        for (size_t i=0; i < len; i++){
          try{
            double key = Zoltan2::IdentifierTraits<GID>::key(gid[i]);
            idx = gidHash_->get(key);
          }
          catch (const std::exception &e) {
            Z2_LOCAL_INPUT_ASSERTION(*env_, "invalid global id", 
              false, BASIC_ASSERTION);
          }
          
          gno[i] = firstGno + idx;
        }
      }
    }
  }
  return;
}

template< typename GID, typename LNO, typename GNO>
  void IdentifierMap<GID,LNO,GNO>::lnoTranslate(
    ArrayView<LNO> lno, 
    ArrayView<GNO> gno, 
    TranslationType tt) const
{
  size_t len=lno.size();

  if (len == 0){
    return;
  }
  Z2_LOCAL_INPUT_ASSERTION(*env_, "invalid TranslationType", 
    (tt==TRANSLATE_LIB_TO_APP) || (tt==TRANSLATE_APP_TO_LIB), 
    BASIC_ASSERTION);

  Z2_LOCAL_INPUT_ASSERTION(*env_, "Destination array is too small",
    ((tt==TRANSLATE_LIB_TO_APP) && (lno.size() >= gno.size())) || 
    ((tt==TRANSLATE_APP_TO_LIB) && (gno.size() >= lno.size())),
    BASIC_ASSERTION);

  GNO firstGno(0), endGno(0);
  if (gnoDist_.size() > 0){
    firstGno = gnoDist_[myRank_];
    endGno = gnoDist_[myRank_+1];
  }
  
  if (tt == TRANSLATE_LIB_TO_APP){
    if (gnoDist_.size() > 0) {   // gnos are consecutive
      for (size_t i=0; i < len; i++){
        Z2_LOCAL_INPUT_ASSERTION(*env_, "invalid global number", 
          (gno[i] >= firstGno) && (gno[i] < endGno), BASIC_ASSERTION);
        lno[i] = gno[i] - firstGno;
      }
    }
    else {                    // gnos must be the app gids
      if (userGidsAreConsecutive_){
        for (size_t i=0; i < len; i++){ 
          GID tmp = Teuchos::as<GID>(gno[i]);
          lno[i] = IdentifierTraits<GID>::difference(myGids_[0], tmp);
          Z2_LOCAL_INPUT_ASSERTION(*env_, "invalid global number",
            (lno[i] >= 0) && (lno[i] < localNumberOfIds_), BASIC_ASSERTION);
        }
      }
      else{
        for (size_t i=0; i < len; i++){ 
          try{
            GID keyArg = Teuchos::as<GID>(gno[i]);
            lno[i] = gidHash_->get(IdentifierTraits<GID>::key(keyArg));
          }
          catch (const std::exception &e) {
            Z2_LOCAL_INPUT_ASSERTION(*env_, "invalid global number",
              false, BASIC_ASSERTION);
          }
        }
      }
    }
  }
  else{                           // TRANSLATE_APP_TO_LIB
    for (size_t i=0; i < len; i++){
      LNO idx = lno[i];

      if (gnoDist_.size() > 0)  // gnos are consecutive
        gno[i] = firstGno + idx;
      else                     // gnos must be the app gids
        gno[i] = Teuchos::as<GNO>(myGids_[idx]);
    }
  }
}

template< typename GID, typename LNO, typename GNO>
  void IdentifierMap<GID,LNO,GNO>::gidGlobalTranslate(
    ArrayView<const GID> in_gid,
    ArrayView<GNO> out_gno,
    ArrayView<int> out_proc) const
{
  typedef typename Teuchos::Hashtable<double, LNO> id2index_hash_t;

  size_t len=in_gid.size();

  if (len == 0){
    return;
  }

  bool skipGno = (out_gno.size() == 0);

  Z2_LOCAL_INPUT_ASSERTION(*env_, "Destination array is too small", 
    (out_proc.size() >= len) && (skipGno || (out_gno.size() >= len)),
    BASIC_ASSERTION);

  if (userGidsAreZoltan2Gnos_ && (gnoDist_.size() > 0)){

    // Easy case - communication is not needed.
    // Global numbers are the application global IDs and
    // they are increasing consecutively with rank.
 
    typename std::map<GNO, int> firstGnoToProc;
    typename std::map<GNO, int>::iterator pos;

    for (int p=0; p <= numProcs_; p++){
      firstGnoToProc[gnoDist_[p]] = p;
    }

    for (size_t i=0; i < len; i++){
      GNO globalNumber = Teuchos::as<GNO>(in_gid[i]);;
      if (!skipGno)
        out_gno[i] = globalNumber;
      pos = firstGnoToProc.upper_bound(globalNumber);
      out_proc[i] = pos->second - 1;
    }

    return;
  }

  bool needGnoInfo = !userGidsAreZoltan2Gnos_;

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
      Z2_LOCAL_MEMORY_ASSERTION(*env_, localNumberOfIds_, false); 

    try{ 
      gidOutBuf.resize(localNumberOfIds_); 
    }
    catch(...)
      Z2_LOCAL_MEMORY_ASSERTION(*env_, localNumberOfIds_, false); 

    for (size_t i=0; i < localNumberOfIds_; i++){
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
  
    for (size_t i=0; i < localNumberOfIds_; i++){
      LNO offset = offsetBuf[hashProc[i]];
      gidOutBuf[offset] = myGids_[i];
      if (needGnoInfo)
        gnoOutBuf[offset] = gnoDist_[myRank_] + i;
      offsetBuf[hashProc[i]] = offset + 1;
    }
    hashProc.clear();
  }

  // Z2::AlltoAllv comment: Buffers are in process rank contiguous order.

  try{
    ArrayView<const GID> gidView = gidOutBuf();
    ArrayView<const LNO> countView = countOutBuf();
    AlltoAllv<GID, LNO>(*comm_, *env_, gidView, countView, gidInBuf, countInBuf);
  }
  Z2_FORWARD_EXCEPTIONS;

  gidOutBuf.clear();
  
  if (needGnoInfo){
    countInBuf.release();
    ArrayView<const GNO> gnoView = gnoOutBuf();
    ArrayView<const LNO> countView = countOutBuf();
    try{
      AlltoAllv<GNO, LNO>(*comm_, *env_, gnoView, countView, gnoInBuf, countInBuf);
    }
    Z2_FORWARD_EXCEPTIONS;
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
  size_t numberOfUniqueGids = 0;
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
      Z2_LOCAL_MEMORY_ASSERTION(*env_, numberOfUniqueGids, false); 

    try{ 
      hashProc.resize(numberOfUniqueGids, 0);
    }
    catch(...)
      Z2_LOCAL_MEMORY_ASSERTION(*env_, numberOfUniqueGids, false); 
  
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
      Z2_LOCAL_MEMORY_ASSERTION(*env_, numberOfUniqueGids, false); 
  
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
    AlltoAllv<GID,LNO>(*comm_, *env_, gidView, countView, gidInBuf, countInBuf);
  }
  Z2_FORWARD_EXCEPTIONS;

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
      Z2_LOCAL_MEMORY_ASSERTION(*env_, total, false); 
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

        Z2_LOCAL_BUG_ASSERTION(*env_, "gidToIndex table", 
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
    AlltoAllv<int,LNO>(*comm_, *env_, procView, countView, procInBuf, countInBuf);
  }
  Z2_FORWARD_EXCEPTIONS;

  procOutBuf.clear();

  if (needGnoInfo){
    try{
      ArrayView<const GNO> gnoView = gnoOutBuf();
      ArrayView<const LNO> countView = countOutBuf();
      AlltoAllv<GNO,LNO>(*comm_, *env_, gnoView, countView, gnoInBuf, countInBuf);
    }
    Z2_FORWARD_EXCEPTIONS;

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
      gno = Teuchos::as<GNO>(gid);
    }

    for (size_t j=0; j < v.size(); j++){
      out_proc[v[j]] = gidProc;
      if (!skipGno)
        out_gno[v[j]] = gno;
    }
  }
}


template< typename GID, typename LNO, typename GNO> 
  void IdentifierMap<GID,LNO,GNO>::setupMap(void)
{
  numProcs_ = comm_->getSize(); 
  myRank_ = comm_->getRank(); 

  Z2_GLOBAL_INPUT_ASSERTION(*env_, 
           "application global ID type is not supported yet",
           IdentifierTraits<GID>::is_valid_id_type() == true, BASIC_ASSERTION);

  localNumberOfIds_ = myGids_.size();

  typedef typename Teuchos::Hashtable<double, LNO> id2index_hash_t;

  global_size_t gtmp = localNumberOfIds_;

  try{
    reduceAll<int, global_size_t>(*comm_, Teuchos::REDUCE_SUM, 
      1, &gtmp, &globalNumberOfIds_);
  } 
  catch (const std::exception &e) {
    Z2_THROW_OUTSIDE_ERROR(*env_, e);
  }

  // Determine whether the application's global ID data type (GID) is a 
  // Teuchos Ordinal, and if so, whether they are consecutive, and what the
  // GID base is.
  //
  // This determines whether we can use the application's IDs internally.

  userGidsAreTeuchosOrdinal_ = false;
  userGidsAreConsecutive_ = false;
  bool baseZeroConsecutiveIds = false;

  const GID *gidPtr = myGids_.get();

  if (IdentifierTraits<GID>::isGlobalOrdinal()){

    userGidsAreTeuchosOrdinal_ = true;

    ArrayRCP<GID> tmpDist(numProcs_+1);

    userGidsAreConsecutive_= globallyConsecutiveOrdinals<GID>(gidPtr,
      localNumberOfIds_,  globalNumberOfIds_, *(env_->comm_), *env_,
      tmpDist(0, numProcs_+1));

    if (userGidsAreConsecutive_){    
      // A GNO is large enough to hold GIDs, but may be a different type.
      if (sizeof(GID) == sizeof(GNO)) {
        gnoDist_ = arcp_reinterpret_cast<GNO>(tmpDist);
      }
      else{
        gnoDist_.resize(numProcs_ + 1);
        for (LNO i=0; i <= numProcs_; i++){
          gnoDist_[i] = static_cast<GNO>(tmpDist[i]);
        }
      }

      if (gnoDist_[0] == 0){
        baseZeroConsecutiveIds = true;
      }
    }
  }

  // If user global IDs are not consecutive ordinals, create a hash table
  // mapping the global IDs to their location in myGids_.

  if (!userGidsAreConsecutive_){
    id2index_hash_t *p = NULL;
    if (localNumberOfIds_){
      try{
        p = new id2index_hash_t(localNumberOfIds_);
      }
      catch (const std::exception &e)
        Z2_LOCAL_MEMORY_ASSERTION(*env_, localNumberOfIds_, false);
    }

    for (size_t i=0; i < localNumberOfIds_; i++){
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

  // Can we use the user's global IDs as our internal global numbers?
  // If so, we can save memory.

  userGidsAreZoltan2Gnos_ = false;

  if (userGidsAreTeuchosOrdinal_ &&
      ((consecutiveGidsAreRequired_ && baseZeroConsecutiveIds) ||
      !consecutiveGidsAreRequired_))
  {
    // We can use the application's global IDs
    userGidsAreZoltan2Gnos_ = true;
  }

  if (userGidsAreZoltan2Gnos_){

    zoltan2GnosAreConsecutive_ = userGidsAreConsecutive_;

    if (userGidsAreConsecutive_){
      minGlobalGno_ = gnoDist_[0];
      maxGlobalGno_ = gnoDist_[numProcs_]-1;
    }
    else{
      intmax_t localMin, globalMin, localMax, globalMax;

      if (localNumberOfIds_ > 0){
        std::pair<GID, GID> minMax =
          IdentifierTraits<GID>::minMax(gidPtr, localNumberOfIds_);
        localMin = Teuchos::as<intmax_t>(minMax.first);
        localMax = Teuchos::as<intmax_t>(minMax.second);
      }
      else{
        localMin = INTMAX_MAX;
        localMax = INTMAX_MIN;
      }

      reduceAll<int, intmax_t>(*comm_, Teuchos::REDUCE_MAX, 
          1, &localMax, &globalMax);

      reduceAll<int, intmax_t>(*comm_, Teuchos::REDUCE_MIN, 
          1, &localMin, &globalMin);

      minGlobalGno_ = static_cast<GNO>(globalMin);
      maxGlobalGno_ = static_cast<GNO>(globalMax);
    }
  } else{
    // We map application gids to consecutive global numbers starting with 0.

    zoltan2GnosAreConsecutive_ = true;

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

    minGlobalGno_ = gnoDist_[0];
    maxGlobalGno_ = gnoDist_[numProcs_]-1;
  }
}

}   // end namespace Zoltan2

#endif /* _ZOLTAN2_IDENTIFIERMAP_HPP_ */
