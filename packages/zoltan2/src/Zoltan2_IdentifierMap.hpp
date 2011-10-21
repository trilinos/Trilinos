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

    lid_t  is the data type used by application for local IDs, which are optional.
    gid_t  is the data type used by application for globls IDs
    lno_t  is the data type used by Zoltan2 for local counts and indexes.
    gno_t  is the integral data type used by Zoltan2 for global counts.
*/

////////////////////////////////////////////////////////////////////
// Declarations
////////////////////////////////////////////////////////////////////

template<typename lid_t, typename gid_t, typename lno_t, typename gno_t>
    class IdentifierMap{

public:

  /*! Constructor - Must be called by all processes 
   *
   * \param comm the problem communicator
   * \param env  the problem and library environment
   * \param gids  the application global IDs
   * \param lids  the application local IDs, if any
   * \param gidsMustBeConsecutive  set to true if the algorithm
   *           or third party library requires consective ids
   *           If necessary the IdentifierMap will map the application's
   *           global IDs to integral IDs beginning at zero.
   */

  explicit IdentifierMap( const RCP<const Comm<int> > &comm, 
                          const RCP<Environment > &env, 
                          const ArrayRCP<gid_t> &gids, 
                          const ArrayRCP<lid_t> &lids,
                          bool gidsMustBeConsecutive=false);

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

  /*! Initialize object if not done in the constructor 
   *
   * \param comm the problem communicator
   * \param env  the problem and library environment
   * \param gids  the application global IDs
   * \param lids  the application local IDs, if any
   * \param gidsMustBeConsecutive  set to true if the algorithm
   *           or third party library requires consective ids
   *           If necessary the IdentifierMap will map the application's
   *           global IDs to integral IDs beginning at zero.
   *
   *  If lids is an empty array, and lid_t is integral, we assume local
   *  ids are consecutive and begin a zero.
   */
  void initialize(const RCP<const Comm<int> > &incomm_,
                  const RCP<Environment > &env,
                  const ArrayRCP<gid_t> &gids,
                  const ArrayRCP<lid_t> &lids,
                  bool gidsMustBeConsecutive=false);

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

  /*! Return the lowest global number if Gids are consecutive.
   */
  void getConsecutiveGidBase(gno_t &base) const;

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
  void gidTranslate(ArrayView<gid_t> gid, 
                    ArrayView<gno_t> gno,
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
  void lidTranslate(ArrayView<lid_t> lid, 
                    ArrayView<gno_t> gno,
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
  void gidGlobalTranslate( ArrayView<const gid_t> in_gid,
                           ArrayView<gno_t> out_gno,
                           ArrayView<int> out_proc) const;

private:

  // Input communicator

  RCP<const Comm<int> > comm_;

  // Problem parameters, library configuration.

  RCP<Environment> env_;

  // Application global and local IDs

  ArrayRCP<gid_t> myGids_; 
  ArrayRCP<lid_t> myLids_;

  // In the case of consecutive ordinal application global IDs,
  // gnoDist[p] is the first global number on process p, and
  // we don't need the gidHash_.

  ArrayRCP<gno_t> gnoDist_;

  // A hash table from application global ID key to our local index.

  RCP<Teuchos::Hashtable<double, lno_t> >  gidHash_;

  // A hash table from application local ID key to our local index.

  RCP<Teuchos::Hashtable<double, lno_t> >  lidHash_;

  typename Array<gno_t>::size_type globalNumberOfIds_;
  typename Array<gno_t>::size_type localNumberOfIds_;
  bool haveLocalIds_;
  int myRank_;
  int numProcs_;

  bool userGidsAreTeuchosOrdinal_;
  bool userGidsAreConsecutive_;
  bool userGidsAreZoltan2Gids_;
  bool zoltan2GidsAreConsecutive_;
  bool consecutiveGidsAreRequired_;

  void setupMap();

};

////////////////////////////////////////////////////////////////////
// Definitions
////////////////////////////////////////////////////////////////////

template<typename lid_t, typename gid_t, typename lno_t, typename gno_t> 
  IdentifierMap<lid_t,gid_t,lno_t,gno_t>::IdentifierMap(
    const RCP<const Comm<int> > &incomm_, const RCP<Zoltan2::Environment> &env,
    const ArrayRCP<gid_t> &gids, const ArrayRCP<lid_t> &lids,
    bool idsMustBeConsecutive) 
         : comm_(incomm_),  env_(env), myGids_(gids), myLids_(lids),
           globalNumberOfIds_(0), localNumberOfIds_(0), haveLocalIds_(false),
           myRank_(0), numProcs_(0), userGidsAreTeuchosOrdinal_(false),
           userGidsAreConsecutive_(false), userGidsAreZoltan2Gids_(false),
           zoltan2GidsAreConsecutive_(false), 
           consecutiveGidsAreRequired_(idsMustBeConsecutive)
{
  setupMap();
}

  /*! Constructor */
template<typename lid_t, typename gid_t, typename lno_t, typename gno_t>
  IdentifierMap<lid_t,gid_t,lno_t,gno_t>::IdentifierMap()  
         : comm_(), env_(), myGids_(), myLids_(),
           globalNumberOfIds_(0), localNumberOfIds_(0), haveLocalIds_(false),
           myRank_(0), numProcs_(0), userGidsAreTeuchosOrdinal_(false),
           userGidsAreConsecutive_(false), userGidsAreZoltan2Gids_(false),
           zoltan2GidsAreConsecutive_(false)
{
}

  /*! Destructor */
template<typename lid_t, typename gid_t, typename lno_t, typename gno_t>
  IdentifierMap<lid_t,gid_t,lno_t,gno_t>::~IdentifierMap() 
  {
  }

#if 0
  /*! Copy Constructor */
template<typename lid_t, typename gid_t, typename lno_t, typename gno_t>
  IdentifierMap<lid_t,gid_t,lno_t,gno_t::IdentifierMap(const IdentifierMap &id)
{
    //TODO    default should work, shouldn't it?
}

  /*! Assignment operator */
template<typename lid_t, typename gid_t, typename lno_t, typename gno_t>
  IdentifierMap<lid_t,gid_t,lno_t,gno_t &
    IdentifierMaplid_t,gid_t,lno_t,gno_t>::operator=(const IdentifierMap<User> &id)
{
    //TODO
}
#endif

template<typename lid_t, typename gid_t, typename lno_t, typename gno_t>
  void IdentifierMap<lid_t,gid_t,lno_t,gno_t>::initialize(
    const RCP<const Comm<int> > &incomm_, 
    const RCP<Zoltan2::Environment> &env,
    const ArrayRCP<gid_t> &gids, const ArrayRCP<lid_t> &lids,
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
  userGidsAreTeuchosOrdinal_ = false;
  userGidsAreConsecutive_ = false ;
  userGidsAreZoltan2Gids_ = false;
  zoltan2GidsAreConsecutive_ = false;
  consecutiveGidsAreRequired_ = idsMustBeConsecutive;

  setupMap();
}

template<typename lid_t, typename gid_t, typename lno_t, typename gno_t>
  bool IdentifierMap<lid_t,gid_t,lno_t,gno_t>::gnosAreGids() const
{
  return userGidsAreZoltan2Gids_;
}

template<typename lid_t, typename gid_t, typename lno_t, typename gno_t>
  bool IdentifierMap<lid_t,gid_t,lno_t,gno_t>::gnosAreConsecutive() const
{
  return zoltan2GidsAreConsecutive_;
}

template<typename lid_t, typename gid_t, typename lno_t, typename gno_t>
  bool IdentifierMap<lid_t,gid_t,lno_t,gno_t>::consecutiveGnosAreRequired() const
{
  return consecutiveGidsAreRequired_;
}

template<typename lid_t, typename gid_t, typename lno_t, typename gno_t>
  void IdentifierMap<lid_t,gid_t,lno_t,gno_t>::getConsecutiveGidBase(gno_t &base) const
{
  if (gnoDist_.size() > 0) 
    base = gnoDist_[0];
  else 
    base = 0;   // error or not?
}

template<typename lid_t, typename gid_t, typename lno_t, typename gno_t>
  void IdentifierMap<lid_t,gid_t,lno_t,gno_t>::gidTranslate(
    ArrayView<gid_t> gid, 
    ArrayView<gno_t> gno,
    TranslationType tt) const
{
  typedef typename Array<gno_t>::size_type teuchos_size_t;
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

  if (IdentifierTraits<gid_t>::isGlobalOrdinalType()){   
                              // our gnos are the app gids
    if (tt == TRANSLATE_LIB_TO_APP)
      for (teuchos_size_t i=0; i < len; i++)
        gid[i] = static_cast<gid_t>(gno[i]);
    else
      for (teuchos_size_t i=0; i < len; i++)
        gno[i] = static_cast<gno_t>(gid[i]);
  }
  else{              // we mapped gids to consecutive gnos
  
    gno_t firstGno = gnoDist_[myRank_];
    gno_t endGno = gnoDist_[myRank_ + 1];

    if (tt == TRANSLATE_LIB_TO_APP){
      for (teuchos_size_t i=0; i < len; i++){

        Z2_LOCAL_INPUT_ASSERTION(*comm_, *env_, "invalid global number", 
        (gno[i] < firstGno) || (gno[i] >= endGno), BASIC_ASSERTION);

        gid[i] = myGids_[gno[i] - firstGno];
      }
    }
    else{
      lno_t idx;
      for (teuchos_size_t i=0; i < len; i++){
        try{
          idx = gidHash_->get(Zoltan2::IdentifierTraits<gid_t>::key(gid[i]));
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

template<typename lid_t, typename gid_t, typename lno_t, typename gno_t>
  void IdentifierMap<lid_t,gid_t,lno_t,gno_t>::lidTranslate(
    ArrayView<lid_t> lid, 
    ArrayView<gno_t> gno, 
    TranslationType tt) const
{
  typedef typename Array<gno_t>::size_type teuchos_size_t;
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

  gno_t firstGno(0), endGno(0);
  if (gnoDist_.size() > 0){
    firstGno = gnoDist_[myRank_];
    endGno = gnoDist_[myRank_+1];
  }
  
  if (tt == TRANSLATE_LIB_TO_APP){
    for (teuchos_size_t i=0; i < len; i++){
      lno_t idx = 0;
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
            IdentifierTraits<gid_t>::key(static_cast<gid_t>(gno[i])));
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
      lno_t idx(0);
      try{
        idx = lidHash_->get(IdentifierTraits<lid_t>::key(lid[i]));
      }
      catch (const std::exception &e) {
        Z2_THROW_OUTSIDE_ERROR(*env_, e);
      }

      if (gnoDist_.size() > 0)  // gnos are consecutive
        gno[i] = firstGno + idx;
      else                     // gnos must be the app gids
        gno[i] = static_cast<gno_t>(myGids_[idx]);
    }
  }
}

template<typename lid_t, typename gid_t, typename lno_t, typename gno_t>
  void IdentifierMap<lid_t,gid_t,lno_t,gno_t>::gidGlobalTranslate(
    ArrayView<const gid_t> in_gid,
    ArrayView<gno_t> out_gno,
    ArrayView<int> out_proc) const
{
  typedef typename Array<gno_t>::size_type teuchos_size_t;
  typedef typename Teuchos::Hashtable<double, lno_t> id2index_hash_t;
  typedef typename Teuchos::Hashtable<double, Array<lno_t> > 
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

  if (IdentifierTraits<gid_t>::isGlobalOrdinalType() && (gnoDist_.size() > 0)){

    // Easy case - communication is not needed.
    // Global numbers are the application global IDs and
    // they are increasing consecutively with rank.
 
    typename std::map<gno_t, int> firstGnoToProc;
    typename std::map<gno_t, int>::iterator pos;

    for (int p=0; p <= numProcs_; p++){
      firstGnoToProc[gnoDist_[p]] = p;
    }

    for (teuchos_size_t i=0; i < len; i++){
      gno_t globalNumber = static_cast<gno_t>(in_gid[i]);
      if (!skipGno)
        out_gno[i] = globalNumber;
      pos = firstGnoToProc.upper_bound(globalNumber);
      out_proc[i] = pos->first - 1;
    }

    return;
  }

  bool needGnoInfo = !IdentifierTraits<gid_t>::isGlobalOrdinalType();

  ///////////////////////////////////////////////////////////////////////
  // First: Hash each of my gid_ts to a process that will answer
  // for it.  Send my gid_ts (and the Gnos if they are different)
  // to their assigned processes.  Build a search structure for
  // the gid_ts that were assigned to me, so I can reply with
  // with the process owning them (and their Gnos if they are different).
  ///////////////////////////////////////////////////////////////////////

  Array<int> hashProc(0);
  Array<gid_t> gidOutBuf(0);
  Array<gno_t> gnoOutBuf(0);
  Array<lno_t> countOutBuf(numProcs_, 0);
  Array<lno_t> offsetBuf(numProcs_ + 1, 0);

  ArrayRCP<gid_t> gidInBuf();
  ArrayRCP<gno_t> gnoInBuf();
  ArrayRCP<lno_t> countInBuf();

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
      hashProc[i] = IdentifierTraits<gid_t>::hashCode(myGids_[i]) % numProcs_;
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
      lno_t offset = offsetBuf[hashProc[i]];
      gidOutBuf[offset] = myGids_[i];
      if (needGnoInfo)
        gnoOutBuf[offset] = gnoDist_[myRank_] + i;
      offsetBuf[hashProc[i]] = offset + 1;
    }
    hashProc.clear();
  }

  // Teuchos comment #1: An Array can be passed for an ArrayView parameter.
  // Teuchos comment #2: gid_t need not be a Teuchos Packet type,
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

  std::map<lno_t, int> firstIndexToProc;
  lno_t total = 0;

  for (int p=0; p < numProcs_; p++){
    firstIndexToProc[total] = p;
    total += countInBuf[p];
  }

  firstIndexToProc[total] = numProcs_;

  id2index_hash_t gidToIndex(total);

  total = 0;
  for (int p=0; p < numProcs_; p++){
    for (lno_t i=countInBuf[p]; i < countInBuf[p+1]; i++, total++){
      try{
        gidToIndex.put(IdentifierTraits<gid_t>::key(gidInBuf[total]), total);
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
  Array<Array<lno_t> > uniqueGidQueryIndices(0);
  teuchos_size_t numberOfUniqueGids = 0;
  Array<lno_t> gidLocation(0);

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

    for (lno_t i=0; i < len; i++){

      double uniqueKey(IdentifierTraits<gid_t>::key(in_gid[i]));

      if (gidIndices->containsKey(uniqueKey)){
        Array<lno_t> v;
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
        Array<lno_t> v(sizeChunk);
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
      gid_t gid = IdentifierTraits<gid_t>::keyToGid(uniqueGidQueries[i]);
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
      for (lno_t i=0; i < countInBuf[p]; i++, total++){
        double k(IdentifierTraits<gid_t>::key(gidInBuf[total]));
        lno_t index(0);
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
    Array<lno_t> v(uniqueGidQueryIndices[i]);

    int gidProc = procInBuf[gidLocation[i]];

    lno_t gno;
    if (needGnoInfo){
      gno = gnoInBuf[gidLocation[i]];
    }
    else{
      gno = static_cast<gno_t>(IdentifierTraits<gid_t>::keyToGid(key));
    }

    for (teuchos_size_t j=0; j < v.size(); j++){
      out_proc[v[j]] = gidProc;
      if (!skipGno)
        out_gno[v[j]] = gno;
    }
  }
}


template<typename lid_t, typename gid_t, typename lno_t, typename gno_t> 
  void IdentifierMap<lid_t,gid_t,lno_t,gno_t>::setupMap(void)
{
  numProcs_ = comm_->getSize(); 
  myRank_ = comm_->getRank(); 

  Z2_GLOBAL_INPUT_ASSERTION( *comm_, *env_, 
           "application global ID type is not supported",
           IdentifierTraits<gid_t>::is_valid_id_type() == true,
           BASIC_ASSERTION);

  localNumberOfIds_ = myGids_.size();

  typedef typename Array<gno_t>::size_type teuchos_size_t;
  typedef typename Teuchos::Hashtable<double, lno_t> id2index_hash_t;

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

    lid_t *lidPtr = myLids_.get();  // for performance

    for (teuchos_size_t i=0; i < localNumberOfIds_; i++){
      try{
        p->put(IdentifierTraits<lid_t>::key(lidPtr[i]), i);
      }
      catch (const std::exception &e) 
        Z2_THROW_OUTSIDE_ERROR(*env_, e);
    }

    lidHash_ = RCP<id2index_hash_t>(p);
  }

  // Determine whether the application's global ID data type (gid_t) is a 
  // Teuchos Ordinal, and if so, whether they are consecutive, and what the
  // gid_t base is.
  //
  // If they are non-consecutive Teuchos Ordinals, and consecutive IDs are
  // not required, or if they are consecutive Teuchos Ordinals, we can use 
  // the user's global IDs.  Otherwise we map their IDs to valid global numbers.

  userGidsAreTeuchosOrdinal_ = false;
  userGidsAreConsecutive_ = false;
  userGidsAreZoltan2Gids_ = false;

  if (IdentifierTraits<gid_t>::isGlobalOrdinalType()){

    // Are the gids consecutive and increasing with process rank? 
    // If so GID/proc lookups can be optimized.

    userGidsAreTeuchosOrdinal_ = true;
    gid_t min(0), max(0), globalMin(0), globalMax(0);
    gid_t *gidPtr = myGids_.get();  // for performance
    min = max = gidPtr[0];
    gid_t checkVal = min;
  
    for (teuchos_size_t i=1; i < localNumberOfIds_; i++){
      if (userGidsAreConsecutive_ && (gidPtr[i] != ++checkVal)){
        userGidsAreConsecutive_=false;
        break;
      }
      if (gidPtr[i] < min)
        min = gidPtr[i];
      else if (gidPtr[i] > max)
        max = gidPtr[i];
    }

    Teuchos::Tuple<gid_t, 4> results;

    results[0] = static_cast<gid_t>(userGidsAreConsecutive_? 1 : 0);
    results[1] = min;
    results[2] = results[3] = 0;

    try{
      Teuchos::reduceAll<int, gid_t>(*comm_, Teuchos::REDUCE_MIN, 2, 
        results.getRawPtr(), results.getRawPtr()+2);
    }
    catch (const std::exception &e) {
      Z2_THROW_OUTSIDE_ERROR(*env_, e);
    }

    if (results[2] != 1)       // min of consecutive flags
      userGidsAreConsecutive_=false;

    if (userGidsAreConsecutive_){
      globalMin = results[3];
      try{
        Teuchos::reduceAll<int, gid_t>(*comm_, Teuchos::REDUCE_MAX, 1, &max, 
          &globalMax);
      }
      catch (const std::exception &e) {
        Z2_THROW_OUTSIDE_ERROR(*env_, e);
      }

      if (globalMax - globalMin + 1 != static_cast<gid_t>(globalNumberOfIds_))
        userGidsAreConsecutive_= false;   // there are gaps in the gids
  
      if (userGidsAreConsecutive_){
        gid_t myStart = myGids_[0];

        gnoDist_ = ArrayRCP<gno_t>(numProcs_ + 1);

        gid_t *startGID = static_cast<gid_t *>(gnoDist_.getRawPtr());
      
        try{
          Teuchos::gatherAll<int, gid_t>(*comm_, 1, &myStart, numProcs_, 
            startGID);
        }
        catch (const std::exception &e) {
          Z2_THROW_OUTSIDE_ERROR(*env_, e);
        }
      
        for (int p=1; p < numProcs_; p++){
          if (startGID[p] < startGID[p-1]){
            gnoDist_.clear();
            userGidsAreConsecutive_= false;  // gids do not increase with process rank
            break;
          }
        }
        if (userGidsAreConsecutive_){
          startGID[numProcs_] = globalMax + 1;
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
      gnoDist_ = ArrayRCP<gno_t>(numProcs_ + 1, 0);
    }
    catch (const std::exception &e) {
      Z2_THROW_OUTSIDE_ERROR(*env_, e);
    }

    gno_t myNum = static_cast<gno_t>(localNumberOfIds_);

    try{
      Teuchos::scan<int, gno_t>(*comm_, Teuchos::REDUCE_SUM, 1, &myNum, 
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

    gid_t *gidPtr = myGids_.get();  // for performance

    for (teuchos_size_t i=0; i < localNumberOfIds_; i++){
      try{
        p->put(IdentifierTraits<gid_t>::key(gidPtr[i]), i);
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
