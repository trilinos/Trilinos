// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/*! \file Zoltan2_IdentifierMap.hpp
    \brief Defines the IdentifierMap class.
*/

#ifndef _ZOLTAN2_IDENTIFIERMAP_HPP_
#define _ZOLTAN2_IDENTIFIERMAP_HPP_

#include <Zoltan2_IdentifierTraits.hpp>
#include <Zoltan2_InputTraits.hpp>
#include <Zoltan2_AlltoAll.hpp>
#include <Zoltan2_GidLookupHelper.hpp>

#include <vector>
#include <algorithm>
#include <map>
#include <iostream>

#include <Teuchos_as.hpp>
#include <Teuchos_CommHelpers.hpp>

namespace Zoltan2
{
/*! \brief Identifier translations available from IdentifierMap.
*/

enum TranslationType {
  TRANSLATE_APP_TO_LIB,  /*!< \brief convert user ids to internal ids */
  TRANSLATE_LIB_TO_APP   /*!< \brief convert internal ids to user ids */
};

/*!  \brief An IdentifierMap manages a global space of object identifiers.

    Data types:
    \li \c lno_t    local indices and local counts
    \li \c gno_t    global indices and global counts
    \li \c gid_t    application global Ids

    The template parameter \c User is a user-defined data type
    which, through a traits mechanism, provides the actual data types
    with which the Zoltan2 library will be compiled.
    \c User may be the actual class or structure used by application to
    represent a vector, or it may be the helper class BasicUserTypes.
    See InputTraits for more information.

      \todo test for global IDs that are std::pair<T1, T2>

  \todo we require that user's gid_ts have base zero if we are to
     use them as consecutive gno_ts.  This can be fixed if needed.
     We are not getting the efficiency
     advantage of using the user's gids simply because
     those IDs do not have base 0.  We can add a flag
     about base 0 being required, and then map only
     if base 0 is required.
*/

////////////////////////////////////////////////////////////////////
// Declarations
////////////////////////////////////////////////////////////////////

template<typename User>
    class IdentifierMap{

public:

  /*! \brief Constructor - Must be called by all processes
   *
   * \param env  the problem and library environment
   * \param comm the problem communicator
   * \param gids  the application global IDs
   * \param gidsMustBeConsecutive  set to true if the algorithm
   *           or third party library requires consective ids
   *           If necessary the IdentifierMap will map the application's
   *           global IDs to consecutive integral IDs beginning at zero.
   */

  typedef typename InputTraits<User>::lno_t lno_t;
  typedef typename InputTraits<User>::gno_t gno_t;
  typedef typename InputTraits<User>::gid_t gid_t;

  explicit IdentifierMap( const RCP<const Environment > &env,
                          const RCP<const Comm<int> > &comm,
                          const ArrayRCP<const gid_t> &gids,
                          bool gidsMustBeConsecutive=false);

  /*! \brief Destructor */
  ~IdentifierMap() {};

  /*! \brief Copy Constructor */
  IdentifierMap(const IdentifierMap &id);

  /*! \brief Assignment operator */
  IdentifierMap &operator=(const IdentifierMap &id);

  /*! \brief Return true if we are using the application global IDs
   *  for our internal global numbers
   */
  bool gnosAreGids() const;

  /*! \brief Return the global number of identifiers.
   */
  gno_t getGlobalNumberOfIds() const { return globalNumberOfIds_;}

  /*! \brief Return the local number of identifiers.
   */
  gno_t getLocalNumberOfIds() const { return localNumberOfIds_;}

  /*! \brief Return the minimum and maximum values of the internal
   *  global numbers
   */
  void getGnoRange(gno_t &min, gno_t &max) const;

  /*! \brief Return true if our internal global numbers are consecutive.
   */
  bool gnosAreConsecutive() const;

  /*! \brief Return true if consecutive Gids are required.
   */
  bool consecutiveGnosAreRequired() const;

  /*! \brief Return the minimum Zoltan2 global Id across all processes.
   */
  gno_t getMinimumGlobalId() const;

  /*! \brief Return the maximum Zoltan2 global Id across all processes.
   */
  gno_t getMaximumGlobalId() const;

  /*! \brief Map the application global IDs to internal global numbers or vice versa.

      \param gid an array of caller's global IDs
      \param gno an array of Zoltan2 global numbers
      \param tt should be TRANSLATE_APP_TO_LIB or TRANSLATE_LIB_TO_APP

      Both gid and gno are preallocated by the caller.

      This is a local call.  If gid is a vector of application global IDs, then
      gno will be set to the corresponding internal global numbers.  If the
      gno vector contains internal global numbers, they will be translated
      to application global IDs.  The application global IDs must be from
      those supplied by this process.

      \todo If the gnos need to be translated to gids, check first to
      see if the gnos are in their original order.  If so, we can
      save memory by returning an arrayview of the gids.
   */
  void gidTranslate(ArrayView<gid_t> gid,
                    ArrayView<gno_t> gno,
                    TranslationType tt) const;

  /*! \brief Map application indices to internal global numbers or vice versa.

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
  void lnoTranslate(ArrayView<lno_t> lno,
                    ArrayView<gno_t> gno,
                    TranslationType tt) const;

  /*! \brief Map application global IDs to internal global numbers
                 and/or process owners.

      \param in_gid input, an array of the global IDs
      \param out_gno output, an optional array of the corresponding
          global numbers used by Zoltan2.  If out_gno.size() is zero,
          we assume global numbers are not needed.
      \param out_proc output, an optional array of the corresponding
          process ranks owning the global Ids.  If out_proc.size() is zero,
          we assume process ranks are not needed.

      If in_gid[i] is not one of the global Ids supplied in this
      problem, then
      out_proc[i] will be -1, and out_gno[i] (if requested) is undefined.
      This behavior supports graph subsetting.

      All processes must call this.  The global IDs
      supplied may belong to another process.
   */
  void gidGlobalTranslate( ArrayView<const gid_t> in_gid,
                           ArrayView<gno_t> out_gno,
                           ArrayView<int> out_proc) const;

  /*! \brief Map library internal global numbers to application GIDs,
                        and their process owners

      \param in_gno input, an array of the global numbers.
      \param out_gid output, an optional array of the corresponding
          user global Ids.  If out_gid.size() is zero,
          we assume global Ids are not desired.
      \param out_proc output, an optional array of the corresponding
          process ranks owning the global numbers.  If out_proc.size() is zero,
          we assume process ranks are not needed.

      Global numbers are always consecutive numbers beginning with
      zero.  If the global number is out of range, a std::runtime_error()
      is thrown.

      All processes must call this.  The global numbers
      supplied may belong to another process.
   */
  void gnoGlobalTranslate( ArrayView<const gno_t> in_gno,
                           ArrayView<gid_t> out_gid,
                           ArrayView<int> out_proc) const;
private:

  // Problem parameters, library configuration.

  const RCP<const Environment> env_;

  // Problem communicator

  const RCP<const Comm<int> > comm_;

  // Application global IDs

  const ArrayRCP<const gid_t> myGids_;

  // Zoltan2 gno_ts will be consecutive if the application gid_ts
  // were mapped to gno_ts, or if the application gid_ts happen
  // to be consecutive ordinals.  In either case, gnoDist_[p]
  // is the first gno_t on process p.  gnoDist_[numProcs_] is the
  // global total of gno_ts.

  ArrayRCP<gno_t> gnoDist_;

  // If application gids are not consecutive ordinals, the gidLookup_
  // object allows lookup of the location of the gid_t in the myGids_ list.

  RCP<GidLookupHelper<gid_t, lno_t > > gidLookup_;

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

  gno_t minGlobalGno_;
  gno_t maxGlobalGno_;

  void setupMap();
};

////////////////////////////////////////////////////////////////////
// Definitions
////////////////////////////////////////////////////////////////////

template<typename User>
  IdentifierMap<User>::IdentifierMap( const RCP<const Environment> &env,
    const RCP<const Comm<int> > &comm,
    const ArrayRCP<const gid_t> &gids, bool idsMustBeConsecutive) :
      env_(env), comm_(comm), myGids_(gids), gnoDist_(), gidLookup_(),
      globalNumberOfIds_(0), localNumberOfIds_(0),
      myRank_(comm_->getRank()), numProcs_(comm_->getSize()),
      userGidsAreTeuchosOrdinal_(false), userGidsAreConsecutive_(false),
      userGidsAreZoltan2Gnos_(false), zoltan2GnosAreConsecutive_(false),
      consecutiveGidsAreRequired_(idsMustBeConsecutive),
      minGlobalGno_(0), maxGlobalGno_(0)
{
  env_->memory("Initial memory in use");

  setupMap();

  env_->memory("After user IDs mapped");
}

// TODO many of these should be inline

template< typename User>
  bool IdentifierMap<User>::gnosAreGids() const
{
  return userGidsAreZoltan2Gnos_;
}

template <typename User>
  void IdentifierMap<User>::getGnoRange(gno_t &min, gno_t &max) const
{
  min = minGlobalGno_;
  max = maxGlobalGno_;
}

template< typename User>
  bool IdentifierMap<User>::gnosAreConsecutive() const
{
  return zoltan2GnosAreConsecutive_;
}

template< typename User>
  bool IdentifierMap<User>::consecutiveGnosAreRequired() const
{
  return consecutiveGidsAreRequired_;
}

template< typename User>
  typename InputTraits<User>::gno_t IdentifierMap<User>::getMinimumGlobalId() const
{
  return minGlobalGno_;
}

template< typename User>
  typename InputTraits<User>::gno_t IdentifierMap<User>::getMaximumGlobalId() const
{
  return maxGlobalGno_;
}

template< typename User>
  void IdentifierMap<User>::gidTranslate(
    ArrayView<gid_t> gid,
    ArrayView<gno_t> gno,
    TranslationType tt) const
{
  typename ArrayView<gid_t>::size_type inLen=gid.size();

  if (inLen == 0){
    return;
  }

  env_->localInputAssertion(__FILE__, __LINE__, "invalid TranslationType",
    (tt==TRANSLATE_APP_TO_LIB) || (tt==TRANSLATE_LIB_TO_APP),
    BASIC_ASSERTION);

  env_->localInputAssertion(__FILE__, __LINE__,
    "Destination array is too small",
    ((tt==TRANSLATE_LIB_TO_APP) && (gid.size() >= gno.size())) ||
     ((tt==TRANSLATE_APP_TO_LIB) && (gno.size() >= gid.size())),
    BASIC_ASSERTION);

  if (userGidsAreZoltan2Gnos_){   // our gnos are the app gids
    if (tt == TRANSLATE_LIB_TO_APP)
      for (typename ArrayView<gid_t>::size_type i=0; i < inLen; i++)
        gid[i] = static_cast<gid_t>(gno[i]);
    else
      for (typename ArrayView<gid_t>::size_type i=0; i < inLen; i++)
        gno[i] = static_cast<gno_t>(gid[i]);
  }
  else{              // we mapped gids to consecutive gnos
    gno_t firstGno = gnoDist_[myRank_];
    gno_t endGno = gnoDist_[myRank_ + 1];

    if (tt == TRANSLATE_LIB_TO_APP){
      for (typename ArrayView<gid_t>::size_type i=0; i < inLen; i++){

        env_->localInputAssertion(__FILE__, __LINE__, "invalid global number",
        (gno[i] >= firstGno) && (gno[i] < endGno), BASIC_ASSERTION);

        gid[i] = myGids_[gno[i] - firstGno];
      }
    }
    else{
      lno_t idx=0;
      if (userGidsAreConsecutive_){
        for (typename ArrayView<gid_t>::size_type i=0; i < inLen; i++){
          gno[i] = firstGno + IdentifierTraits<gid_t>::difference(
            myGids_[0], gid[i]);
          env_->localInputAssertion(__FILE__, __LINE__, "invalid global id",
            (gno[i] >= firstGno) && (gno[i] < endGno), BASIC_ASSERTION);
        }
      }
      else{
        for (typename ArrayView<gid_t>::size_type i=0; i < inLen; i++){
          try{
            idx = gidLookup_->lookup(gid[i]);
          }
          catch (const std::exception &e) {
            env_->localInputAssertion(__FILE__, __LINE__, "invalid global id",
              false, BASIC_ASSERTION);
          }

          gno[i] = firstGno + idx;
        }
      }
    }
  }
  return;
}

template< typename User>
  void IdentifierMap<User>::lnoTranslate(
    ArrayView<lno_t> lno,
    ArrayView<gno_t> gno,
    TranslationType tt) const
{
  typename ArrayView<lno_t>::size_type inLen=lno.size();

  if (inLen == 0){
    return;
  }
  env_->localInputAssertion(__FILE__, __LINE__, "invalid TranslationType",
    (tt==TRANSLATE_LIB_TO_APP) || (tt==TRANSLATE_APP_TO_LIB),
    BASIC_ASSERTION);

  env_->localInputAssertion(__FILE__, __LINE__,
    "Destination array is too small",
    ((tt==TRANSLATE_LIB_TO_APP) && (lno.size() >= gno.size())) ||
    ((tt==TRANSLATE_APP_TO_LIB) && (gno.size() >= lno.size())),
    BASIC_ASSERTION);

  gno_t firstGno(0), endGno(0);
  if (gnoDist_.size() > 0){
    firstGno = gnoDist_[myRank_];
    endGno = gnoDist_[myRank_+1];
  }

  if (tt == TRANSLATE_LIB_TO_APP){
    if (gnoDist_.size() > 0) {   // gnos are consecutive
      for (typename ArrayView<lno_t>::size_type i=0; i < inLen; i++){
        env_->localInputAssertion(__FILE__, __LINE__, "invalid global number",
          (gno[i] >= firstGno) && (gno[i] < endGno), BASIC_ASSERTION);
        lno[i] = gno[i] - firstGno;
      }
    }
    else {                    // gnos must be the app gids
      if (userGidsAreConsecutive_){
        for (typename ArrayView<lno_t>::size_type i=0; i < inLen; i++){
          gid_t tmp = static_cast<gid_t>(gno[i]);
          lno[i] = IdentifierTraits<gid_t>::difference(myGids_[0], tmp);
          env_->localInputAssertion(__FILE__, __LINE__, "invalid global number",
            (lno[i] >= 0) && (lno[i] < lno_t(localNumberOfIds_)),
            BASIC_ASSERTION);
        }
      }
      else{
        for (typename ArrayView<lno_t>::size_type i=0; i < inLen; i++){
          try{
            gid_t keyArg = static_cast<gid_t>(gno[i]);
            lno[i] = gidLookup_->lookup(keyArg);
          }
          catch (const std::exception &e) {
            env_->localInputAssertion(__FILE__, __LINE__,
              "invalid global number", false, BASIC_ASSERTION);
          }
        }
      }
    }
  }
  else{                           // TRANSLATE_APP_TO_LIB
    for (typename ArrayView<lno_t>::size_type i=0; i < inLen; i++){
      lno_t idx = lno[i];

      if (gnoDist_.size() > 0)  // gnos are consecutive
        gno[i] = firstGno + idx;
      else                     // gnos must be the app gids
        gno[i] = static_cast<gno_t>(myGids_[idx]);
    }
  }
}

template< typename User>
  void IdentifierMap<User>::gidGlobalTranslate(
    ArrayView<const gid_t> in_gid,
    ArrayView<gno_t> out_gno,
    ArrayView<int> out_proc) const
{
  typename ArrayView<const gid_t>::size_type inLen = in_gid.size();

  // It's faster in loops to user raw pointers.
  const gid_t *gids = in_gid.getRawPtr();
  gno_t *gnos = out_gno.getRawPtr();
  int *pids = out_proc.getRawPtr();

  bool needGnos = false;

  if (out_gno.size() > 0){
    env_->localInputAssertion(__FILE__, __LINE__, "array too small",
      out_gno.size()>=static_cast<typename ArrayView<gno_t>::size_type>(inLen),
      BASIC_ASSERTION);

    if (userGidsAreZoltan2Gnos_){
      // Global numbers are the application global IDs
      for (typename ArrayView<const gid_t>::size_type i=0; i < inLen; i++)
        gnos[i] = static_cast<gno_t>(gids[i]);
    }
    else
      needGnos = true;
  }

  bool needProcs = false;

  if (out_proc.size() > 0){
    env_->localInputAssertion(__FILE__, __LINE__, "array too small",
      out_proc.size()>=static_cast<typename ArrayView<int>::size_type>(inLen),
      BASIC_ASSERTION);

    if (userGidsAreZoltan2Gnos_ && (gnoDist_.size() > 0)){

      // The gids are the gnos, which are consecutive
      // increasing with rank.

      gno_t *gnoDist = gnoDist_.getRawPtr();
      gno_t *gnoEnd = gnoDist + numProcs_ + 1;

      for (typename ArrayView<const gid_t>::size_type i=0; i < inLen; i++){
        gno_t gno = static_cast<gno_t>(gids[i]);
        // pointer to first entry greater than gno
        gno_t *ub = std::upper_bound(gnoDist, gnoEnd, gno);

        if (ub == gnoDist || ub == gnoEnd)
          pids[i] = -1;   // globally not one of our gids
        else
          pids[i] = (ub - gnoDist - 1);
      }
    }
    else
      needProcs = true;
  }

  int haveWork = ((needProcs || needGnos) ? 1 : 0);
  int someHaveWork = 0;

  Teuchos::reduceAll<int, int>(*comm_, Teuchos::REDUCE_MAX, 1,
    &haveWork, &someHaveWork);

  if (!someHaveWork)
    return;

  ///////////////////////////////////////////////////////////////////////
  // First: Hash each of my local gids to a process that will answer
  // for it.  Send my gids (and the Gnos if they are different)
  // to their assigned processes.  Build a lookup object for
  // the gids that were assigned to me, so I can reply with
  // with the process owning them (and their Gnos if they are different).
  ///////////////////////////////////////////////////////////////////////

  Array<gid_t> gidOutBuf;
  Array<gno_t> gnoOutBuf;
  Array<int> countOutBuf(numProcs_, 0);

  ArrayRCP<gid_t> gidInBuf;
  ArrayRCP<gno_t> gnoInBuf;
  ArrayRCP<int> countInBuf;

  Array<gno_t> offsetBuf(numProcs_ + 1, 0);

  if (localNumberOfIds_ > 0){

    try{
      gidOutBuf.resize(localNumberOfIds_);
    }
    catch(...){
      env_->localMemoryAssertion(__FILE__, __LINE__, localNumberOfIds_, false);
    }

    for (size_t i=0; i < localNumberOfIds_; i++){
      int hashProc = IdentifierTraits<gid_t>::hashCode(myGids_[i]) % numProcs_;
      countOutBuf[hashProc]++;
    }

    for (int p=1; p <= numProcs_; p++){
      offsetBuf[p] = offsetBuf[p-1] + countOutBuf[p-1];
    }

    if (needGnos){
      // The gnos are not the gids, which also implies that
      // gnos are consecutive numbers given by gnoDist_.
      gnoOutBuf.resize(localNumberOfIds_, 0);
    }

    for (size_t i=0; i < localNumberOfIds_; i++){
      int hashProc = IdentifierTraits<gid_t>::hashCode(myGids_[i]) % numProcs_;
      gno_t offset = offsetBuf[hashProc];
      gidOutBuf[offset] = myGids_[i];
      if (needGnos)
        gnoOutBuf[offset] = gnoDist_[myRank_] + i;
      offsetBuf[hashProc] = offset + 1;
    }
  }

  // Z2::AlltoAllv comment: Buffers are in process rank contiguous order.

  ArrayView<const gid_t> gidView = gidOutBuf();
  ArrayView<const int> countView = countOutBuf();
  try{
    AlltoAllv<gid_t>(*comm_, *env_, gidView, countView, gidInBuf, countInBuf);
  }
  Z2_FORWARD_EXCEPTIONS;

  gidOutBuf.clear();

  if (needGnos){
    ArrayView<const gno_t> gnoView = gnoOutBuf();
    try{
      AlltoAllv<gno_t>(*comm_, *env_, gnoView, countView, gnoInBuf, countInBuf);
    }
    Z2_FORWARD_EXCEPTIONS;
  }

  gnoOutBuf.clear();
  countOutBuf.clear();

  // Save the information that was hashed to me so I can do lookups.
  // This is a lookup from gid to its position in gidInBuf.
  // (The list of gids hashed to me is unique.)

  typedef GidLookupHelper<gid_t, lno_t> lookup_t;

  RCP<lookup_t> lookupMine = rcp(new lookup_t(env_, gidInBuf));

  // Use a vector to find process that sent gid.

  std::vector<gno_t> firstIndex;
  std::vector<int> sendProc;
  gno_t indexTotal = 0;

  if (needProcs){
    for (int p=0; p < numProcs_; p++){
      int bufLen = countInBuf[p];
      if (bufLen > 0){
        firstIndex.push_back(indexTotal);
        sendProc.push_back(p);
        indexTotal += bufLen;
      }
    }
    firstIndex.push_back(indexTotal);
  }

  ///////////////////////////////////////////////////////////////////////
  // Send a request for information to the "answer process" for each
  // of the unique gids in in_gid.
  //
  // We may be called by a model that wants to find out what
  // process owns neighbors or columns.  Such a list of gids is
  // likely to have many duplicates.
  //
  // It is possible that some of the gids do not belong to any process.
  // (This happens in graph subsetting.)
  ///////////////////////////////////////////////////////////////////////

  ArrayRCP<const gid_t> in_gidArray =
    arcp(in_gid.getRawPtr(), 0, inLen, false);
  RCP<lookup_t> lookupRequested = rcp(new lookup_t(env_, in_gidArray));

  size_t numberOfUniqueGids = lookupRequested->size();

  std::map<gid_t, lno_t> answerMap;// map from gid to position of answer

  ArrayRCP<lno_t> uniqueIndices;
  lookupRequested->getIndices(uniqueIndices);

  countOutBuf.resize(numProcs_, 0);

  if (numberOfUniqueGids > 0){
    try{
      gidOutBuf.resize(numberOfUniqueGids);
    }
    catch(...){
      env_->localMemoryAssertion(__FILE__, __LINE__, numberOfUniqueGids,
        false);
    }

    for (size_t i=0; i < numberOfUniqueGids; i++){
      gid_t gid = in_gid[uniqueIndices[i]];
      int hashProc = IdentifierTraits<gid_t>::hashCode(gid) % numProcs_;
      countOutBuf[hashProc]++;
    }

    offsetBuf[0] = 0;

    for (int p=0; p < numProcs_; p++){
      offsetBuf[p+1] = offsetBuf[p] + countOutBuf[p];
    }

    for (size_t i=0; i < numberOfUniqueGids; i++){
      gid_t gid = in_gid[uniqueIndices[i]];
      int hashProc = IdentifierTraits<gid_t>::hashCode(gid) % numProcs_;
      gno_t loc = offsetBuf[hashProc];

      answerMap[gid] = loc;

      gidOutBuf[loc] = gid;
      offsetBuf[hashProc] = loc + 1;
    }
  }

  try{
    ArrayView<const gid_t> gidView = gidOutBuf();
    ArrayView<const int> countView = countOutBuf();
    AlltoAllv<gid_t>(*comm_, *env_, gidView, countView,
      gidInBuf, countInBuf);
  }
  Z2_FORWARD_EXCEPTIONS;

  gidOutBuf.clear();

  ///////////////////////////////////////////////////////////////////////
  // Create and send answers to the processes that made requests of me.
  ///////////////////////////////////////////////////////////////////////

  gno_t total = 0;
  for (int p=0; p < numProcs_; p++){
    countOutBuf[p] = countInBuf[p];
    total += countOutBuf[p];
  }

  ArrayRCP<int> procOutBuf;
  ArrayRCP<int> procInBuf;

  if (needProcs){
    int *tmp = new int [total];
    env_->localMemoryAssertion(__FILE__, __LINE__, total, tmp);
    procOutBuf = arcp(tmp, 0, total, true);

  }

  if (needGnos){
    try{
      gnoOutBuf.resize(total, 0);
    }
    catch(...){
      env_->localMemoryAssertion(__FILE__, __LINE__, total, false);
    }
  }

  if (total > 0){

    total=0;
    typename std::vector<gno_t>::iterator indexFound;

    for (int p=0; p < numProcs_; p++){
      for (int i=0; i < countInBuf[p]; i++, total++){

        lno_t index(0);
        bool badGid = false;

        try{
          index = lookupMine->lookup(gidInBuf[total]);
        }
        catch (std::exception &e){
          badGid = true;
        }

        env_->localBugAssertion(__FILE__, __LINE__, "gidToIndex table",
          badGid || ((index >= 0)&&(index<=indexTotal)), BASIC_ASSERTION);

        if (!badGid){

          if (needProcs){
            indexFound =
              upper_bound(firstIndex.begin(), firstIndex.end(), index);
            int sendingProc = indexFound - firstIndex.begin() - 1;
            procOutBuf[total] = sendProc[sendingProc];
          }

          if (needGnos){
            gnoOutBuf[total] = gnoInBuf[index];
          }
        }
        else if (needProcs){
          // globally not one of our gids, can happen in subsetting
          procOutBuf[total] = -1;
        }
      }
    }
  }

  gidInBuf.clear();
  gnoInBuf.clear();
  lookupMine.~RCP<lookup_t>();

  if (needProcs){
    try{
      ArrayView<const int> procView = procOutBuf.view(0, total);
      ArrayView<const int> countView = countOutBuf();
      AlltoAllv<int>(*comm_, *env_, procView, countView, procInBuf, countInBuf);
    }
    Z2_FORWARD_EXCEPTIONS;

    procOutBuf.clear();
  }

  if (needGnos){
    try{
      ArrayView<const gno_t> gnoView = gnoOutBuf();
      ArrayView<const int> countView = countOutBuf();
      AlltoAllv<gno_t>(*comm_, *env_, gnoView, countView, gnoInBuf, countInBuf);
    }
    Z2_FORWARD_EXCEPTIONS;

    gnoOutBuf.clear();
  }

  countOutBuf.clear();

  ///////////////////////////////////////////////////////////////////////
  // Done.  Process the replies to my queries
  ///////////////////////////////////////////////////////////////////////

  for (typename ArrayView<const gid_t>::size_type i=0; i < inLen; i++){
    lno_t loc = answerMap[in_gid[i]];
    if (needProcs)
      out_proc[i] = procInBuf[loc];
    if (needGnos)
      out_gno[i] = gnoInBuf[loc];
  }
}

template< typename User>
  void IdentifierMap<User>::gnoGlobalTranslate(
    ArrayView<const gno_t> in_gno,
    ArrayView<gid_t> out_gid,
    ArrayView<int> out_proc) const
{
  typename ArrayView<const gno_t>::size_type inLen = in_gno.size();

  if (userGidsAreZoltan2Gnos_){

    // Easy case - use gidGlobalTranslate.

    const gno_t *gnos = in_gno.getRawPtr();
    ArrayView<const gid_t> gids(static_cast<const gid_t *>(gnos), inLen);
    ArrayView<gno_t> noGnos;

    try{
      gidGlobalTranslate(gids, noGnos, out_proc);
    }
    Z2_FORWARD_EXCEPTIONS;

    if (out_gid.size() ==
                static_cast<typename ArrayView<gid_t>::size_type>(inLen))
      for (typename ArrayView<const gno_t>::size_type i=0; i < inLen; i++)
        out_gid[i] = gids[i];

    return;
  }

  bool needProcs = out_proc.size() > 0;
  bool needGids = out_gid.size() > 0;

  int haveWork = ((needProcs || needGids) ? 1 : 0);
  int someHaveWork = 0;

  Teuchos::reduceAll<int, int>(*comm_, Teuchos::REDUCE_MAX, 1,
    &haveWork, &someHaveWork);

  if (!someHaveWork)
    return;

  env_->localInputAssertion(__FILE__, __LINE__, "array too small",
    (!needProcs ||
      out_proc.size()>=static_cast<typename ArrayView<int>::size_type>(inLen))
    &&
    (!needGids ||
      out_gid.size()>=static_cast<typename ArrayView<gid_t>::size_type>(inLen)),
     BASIC_ASSERTION);

  // Get process owning each gno.  We need it to get the gid.

  ArrayRCP<int> procArray;

  if (out_proc.size() == static_cast<typename ArrayView<int>::size_type>(inLen))
    procArray = arcp(out_proc.getRawPtr(), 0, inLen, false);
  else{
    int *tmp = new int [inLen];
    env_->localMemoryAssertion(__FILE__, __LINE__, inLen, tmp);
    procArray = arcp(tmp, 0, inLen, true);
  }

  // The fact that global Ids and global Numbers are different
  // means that the user's global Ids were mapped to consecutive increasing
  // global numbers, so we know the process owning each global number.

  gno_t *gnos = gnoDist_.getRawPtr();
  gno_t *final = gnos + numProcs_ + 1;
  bool remoteGno = false;
  int *pids = procArray.getRawPtr();  // in loops it's faster to use raw ptr

  for (typename ArrayView<const gno_t>::size_type i=0; i < inLen; i++){

    gno_t *ub = std::upper_bound(gnos, final, in_gno[i]);

    pids[i] = (ub - gnos - 1);

    env_->localInputAssertion(__FILE__, __LINE__, "invalid global number",
      (pids[i] >= 0) && (pids[i] < env_->numProcs_), BASIC_ASSERTION);

    if (!remoteGno && (pids[i] != myRank_))
      remoteGno = true;
  }

  if (!needGids)
    return;

  if (!remoteGno){

    // Make a local call to get the gids

    const gno_t *gnos = in_gno.getRawPtr();
    ArrayView<gno_t> gnoList(const_cast<gno_t *>(gnos), inLen);

    try{
      gidTranslate(out_gid, gnoList, TRANSLATE_LIB_TO_APP);
    }
    Z2_FORWARD_EXCEPTIONS;

    return;
  }

  // Get the global ID from the owner.
  // TODO: In gidGlobalTranslate, we go to much trouble to
  //   avoid duplicate requests.  (Because that may be used
  //   to identify owners of neighbors, which would likely
  //   include many duplicates.)  Any reason to go to this
  //   trouble here?  Could the list of gnos have many duplicates?

  ArrayRCP<gno_t> gnoOutBuf;
  ArrayRCP<gno_t> offsetBuf;
  gno_t *tmpGno = NULL;

  if (inLen){
    tmpGno = new gno_t [inLen];
    env_->localMemoryAssertion(__FILE__, __LINE__, inLen, tmpGno);
    gnoOutBuf = arcp(tmpGno, 0, inLen, true);
  }

  int *tmpCount = new int [numProcs_];
  env_->localMemoryAssertion(__FILE__, __LINE__, numProcs_, tmpCount);
  memset(tmpCount, 0, sizeof(int) * numProcs_);
  ArrayRCP<int> countOutBuf(tmpCount, 0, numProcs_, true);

  gno_t *tmpOff = new gno_t [numProcs_+1];
  env_->localMemoryAssertion(__FILE__, __LINE__, numProcs_+1, tmpOff);
  offsetBuf = arcp(tmpOff, 0, numProcs_+1, true);

  for (typename ArrayView<const gno_t>::size_type i=0; i < inLen; i++)
    tmpCount[pids[i]]++;

  offsetBuf[0] = 0;
  for (int i=0; i < numProcs_; i++)
    tmpOff[i+1] = tmpOff[i] + tmpCount[i];

  for (typename ArrayView<const gno_t>::size_type i=0; i < inLen; i++){
    int p = pids[i];
    gno_t off = tmpOff[p];
    tmpGno[off] = in_gno[i];
    tmpOff[p]++;
  }

  // Reset offsetBuf so we can find replies later.

  tmpOff[0] = 0;
  for (int i=0; i < numProcs_; i++)
    tmpOff[i+1] = tmpOff[i] + tmpCount[i];

  ArrayRCP<gno_t> gnoInBuf;
  ArrayRCP<int> countInBuf;

  try{
    AlltoAllv<gno_t>(*comm_, *env_,
      gnoOutBuf.view(0, inLen), countOutBuf.view(0, numProcs_),
      gnoInBuf, countInBuf);
  }
  Z2_FORWARD_EXCEPTIONS;

  gnoOutBuf.clear();
  countOutBuf.clear();

  gno_t numRequests = 0;
  for (int i=0; i < numProcs_; i++)
    numRequests += countInBuf[i];

  ArrayRCP<gid_t> gidQueryBuf;

  if (numRequests){
    gid_t *tmpGid = new gid_t [numRequests];
    env_->localMemoryAssertion(__FILE__, __LINE__, numRequests, tmpGid);
    gidQueryBuf = arcp(tmpGid, 0, numRequests);
  }

  try{
    gidTranslate(gidQueryBuf.view(0, numRequests),
                 gnoInBuf.view(0, numRequests),
                 TRANSLATE_LIB_TO_APP);
  }
  Z2_FORWARD_EXCEPTIONS;

  ArrayRCP<gid_t> gidInBuf;
  ArrayRCP<lno_t> newCountInBuf;

  try{
    AlltoAllv<gid_t>(*comm_, *env_,
      gidQueryBuf.view(0, numRequests), countInBuf.view(0, numProcs_),
      gidInBuf, newCountInBuf);
  }
  Z2_FORWARD_EXCEPTIONS;

  gidQueryBuf.clear();

  // copy the replies into the output gid list

  gid_t *gidInfo = gidInBuf.getRawPtr();  // faster in a loop

  for (typename ArrayView<const gno_t>::size_type i=0; i < inLen; i++){
    int p = pids[i];
    gno_t off = tmpOff[p];
    out_gid[i] = gidInfo[off];
    tmpOff[p]++;
  }
}

template< typename User>
  void IdentifierMap<User>::setupMap(void)
{
  env_->globalInputAssertion(__FILE__, __LINE__,
       "application global ID type is not supported yet",
       IdentifierTraits<gid_t>::is_valid_id_type() == true, BASIC_ASSERTION,
       comm_);

  if (IdentifierTraits<gid_t>::isGlobalOrdinal())
    userGidsAreTeuchosOrdinal_ = true;

  localNumberOfIds_ = myGids_.size();
  const gid_t *gidPtr = myGids_.get();
  ArrayRCP<gid_t> tmpDist;
  gid_t mingid=0, maxgid=0;

  userGidsAreConsecutive_ = globallyConsecutiveOrdinals<gid_t>(
    *comm_, *env_, gidPtr, localNumberOfIds_,      // input
    tmpDist, globalNumberOfIds_);                  // output

  bool baseZeroConsecutiveIds = false;

  if (userGidsAreTeuchosOrdinal_){
    if (!userGidsAreConsecutive_){
      mingid = tmpDist[0];
      maxgid = tmpDist[1];
    }
    else{
      // A gno_t is large enough to hold gid_ts, but may be a different type.
      if (sizeof(gid_t) == sizeof(gno_t)) {
        gnoDist_ = arcp_reinterpret_cast<gno_t>(tmpDist);
      }
      else{
        gnoDist_.resize(numProcs_ + 1);
        for (lno_t i=0; i <= numProcs_; i++){
          gnoDist_[i] = static_cast<gno_t>(tmpDist[i]);
        }
      }

      if (gnoDist_[0] == 0)
        baseZeroConsecutiveIds = true;
    }
  }

  // If user global IDs are not consecutive ordinals, create a lookup table
  // mapping the global IDs to their location in myGids_.

  typedef GidLookupHelper<gid_t, lno_t> lookup_t;

  if (!userGidsAreConsecutive_){
    try{
      gidLookup_ = rcp(new lookup_t(env_, myGids_));
    }
    Z2_FORWARD_EXCEPTIONS;
  }

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
      minGlobalGno_ = static_cast<gno_t>(mingid);
      maxGlobalGno_ = static_cast<gno_t>(maxgid);
    }
  } else{
    // We map application gids to consecutive global numbers starting with 0.

    zoltan2GnosAreConsecutive_ = true;

    try{
      gnoDist_.resize(numProcs_ + 1, 0);
    }
    Z2_THROW_OUTSIDE_ERROR(*env_);

    gno_t myNum = static_cast<gno_t>(localNumberOfIds_);

    try{
      gno_t *p = gnoDist_.getRawPtr();
      Teuchos::gatherAll<int, gno_t>(*comm_, 1, &myNum, numProcs_, p+1);
    }
    Z2_THROW_OUTSIDE_ERROR(*env_);

    for (int i=2; i <= numProcs_; i++){
      gnoDist_[i] += gnoDist_[i-1];
    }

    minGlobalGno_ = gnoDist_[0];
    maxGlobalGno_ = gnoDist_[numProcs_]-1;
  }
}

}   // end namespace Zoltan2

#endif /* _ZOLTAN2_IDENTIFIERMAP_HPP_ */
