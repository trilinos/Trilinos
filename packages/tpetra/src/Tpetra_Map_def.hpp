// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER

/// \file Tpetra_Map_def.hpp
///
/// Implementation of the methods of Tpetra::Map, and of related
/// nonmember constructors for Tpetra::Map.

#ifndef TPETRA_MAP_DEF_HPP
#define TPETRA_MAP_DEF_HPP

#include <Teuchos_as.hpp>
#include "Tpetra_Directory.hpp" // must include for implicit instantiation to work
#include "Tpetra_Util.hpp"
#include <stdexcept>

#ifdef DOXYGEN_USE_ONLY
  #include "Tpetra_Map_decl.hpp"
#endif

namespace {

  //! Whether localValue is true on all processes in the communicator.
  bool
  allProcsTrue (const Teuchos::Comm<int>& comm,
                const bool localValue)
  {
    using Teuchos::outArg;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;

    // reduceAll not yet implemented for bool, so we convert to char.
    char localChar = localValue ? 1 : 0;
    char globalChar = 1;
    reduceAll (comm, REDUCE_MIN, localChar, outArg (globalChar));
    return (globalChar == 1);
  }

  //! Whether the pair of integers (prevGid, curGid) is contiguous (repeats allowed).
  template<class GO>
  bool
  contiguous (const GO prevGid, const GO curGid) {
    return curGid == prevGid || curGid == prevGid + Teuchos::as<GO>(1);
  }

  // \fn areMyGidsLocallyContiguous
  // \brief Are my global IDs (GIDs) contiguous? and other useful info.
  //
  // "Contiguous" allows repeated entries.  Checking for contiguity
  // requires iterating through the GID list, so we save some passes
  // over the list by computing other useful things at the same time.
  //
  // This is a local function.  It does not involve distributed-memory
  // communication and need not be called collectively.
  //
  // \param myMinGid [out] My process' minimum GID.  Only meaningful
  //   if myGids.size() > 0.
  //
  // \param myMaxGid [out] My process' maximum GID.  Only meaningful
  //   if myGids.size() > 0.
  //
  // \param myLastInitContigArrayIndex [out] My process' last initial
  //   contiguous array index in the input view.  That is, if the
  //   elements myGids[0], myGids[1], ..., myGids[k] are contiguous
  //   and myGids[k+1] != myGids[k], then set this to k.  Only
  //   meaningful if myGids.size() > 0.
  //
  // \param myGids [in] My process' GIDs.  They need not be sorted,
  //   and there may be duplicates.
  //
  // \param node [in/out] My process' Kokkos Node instance.  At some
  //   point, we might use this to turn the check for contiguity into
  //   a parallel kernel.
  //
  // \return Whether my GIDs are contiguous.
  template<class Map>
  bool
  areMyGidsLocallyContiguous (typename Map::global_ordinal_type& myMinGid,
                              typename Map::global_ordinal_type& myMaxGid,
                              typename Teuchos::ArrayView<const typename Map::global_ordinal_type>::size_type& myLastInitContigArrayIndex,
                              const Teuchos::ArrayView<const typename Map::global_ordinal_type>& myGids,
                              const Teuchos::Ptr<typename Map::node_type>& node)
  {
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::as;
    typedef typename Map::global_ordinal_type GO;
    typedef typename ArrayView<const GO>::size_type size_type;

    // We're not using the Kokkos Node for now, but it wouldn't be
    // hard to write a parallel kernel to compute all the things we
    // want to compute.  (Checking contiguity is a scan, for example.)
    (void) node;

    bool locallyContiguous = true;
    GO minGid = as<GO> (0); // not meaningful if there are no GIDs
    GO maxGid = as<GO> (0); // ditto
    size_type lastContigInd = as<size_type> (0); // ditto

    if (myGids.size() > 0) { // trivially contiguous if size is zero
      GO prevGid = myGids[0];
      minGid = prevGid;
      maxGid = prevGid;

      for (size_type k = 1; k < myGids.size(); ++k) {
        const GO curGid = myGids[k];
        // We've phrased the comparison so that GO may be either
        // signed or unsigned.  Remember that GIDs need not be sorted,
        // and may have repeated entries (which are still counted as
        // contiguous).
        if (! contiguous (prevGid, curGid)) {
          locallyContiguous = false;
        }
        else {
          lastContigInd = k;
        }
        prevGid = curGid;

        // Update the min and max GID.
        if (curGid < minGid) {
          minGid = curGid;
        }
        if (curGid > maxGid) {
          maxGid = curGid;
        }
      }
    }
    myMinGid = minGid;
    myMaxGid = maxGid;
    myLastInitContigArrayIndex = lastContigInd;
    return locallyContiguous;
  }

  // \fn areGidsGloballyContiguous
  // \brief Are the GIDs globally contiguous? and if so, tell me each process' min GID (and other info).
  //
  // Contiguous Maps are much faster and take much less memory than
  // noncontiguous Maps, so it's worthwhile to check for contiguity,
  // even when using the most general Map constructor (that accepts a
  // list on each process of the GIDs that process owns).
  //
  // Global contiguity requires all of the following properties:
  // 1. The GIDs on each process are locally contiguous.
  // 2. For all process ranks p in the communicator with p > 0,
  //    process p's minimum GID is one plus process p-1's maximum GID.
  //
  // We can check #2 by computing an array containing the minimum GID
  // of all processes, replicated on all processes.  This array
  // minimum GIDs on each process is useful anyway for constructing
  // the Map's Directory, so we return the array so the Directory
  // doesn't need to compute it again.
  //
  // This function is a collective over the given communicator.
  //
  // \param allMinGids [out] If the GIDs are globally contiguous, then
  //   this is a globally replicated array of size comm->getSize()+1,
  //   for which the p-th entry is the min GID on process p, and the
  //   last entry is one plus the max GID over all processes (this is
  //   useful as an iteration guard, so that you can write loops over
  //   GIDs with the traditional strict less-than bound).  Otherwise,
  //   this array is undefined.
  // \param allMaxGids [out] If the GIDs are globally contiguous, then
  //   this is a globally replicated array of size comm->getSize(),
  //   for which the p-th entry is the max GID on process p.
  //   Otherwise, this array is undefined.
  // \param myMinGid [out] My process' minimum GID.  Only meaningful
  //   if myGids.size() > 0.
  // \param myMaxGid [out] My process' maximum GID.  Only meaningful
  //   if myGids.size() > 0.
  // \param myLastInitContigArrayIndex [out] My process' the last
  //   initial contiguous array index in the input view.  That is, if
  //   the elements myGids[0], myGids[1], ..., myGids[k] are
  //   contiguous and myGids[k+1] != myGids[k], then set this to k.
  //   Only meaningful if myGids.size() > 0.
  // \param myGids [in] My process' GIDs.  They need not be sorted,
  //   and there may be duplicates.
  // \param comm [in] The communicator over which to check global
  //   contiguity.
  // \param node [in/out] My process' Kokkos Node instance.  At some
  //   point, we might use this to turn the check for contiguity into
  //   a parallel kernel.
  //
  // \return Whether all processes' GIDs are contiguous.
  //
  // \note In the globally contiguous case, the globally replicated
  //   allMinGids array functions as the "directory" that maps any GID
  //   to its owning process without communication.
  //
  // \note On implementation: We could just check for global
  //   contiguity on Proc 0, with a gather instead of an all-gather.
  //   However, the Directory uses the globally replicated array of
  //   min GIDs and would need to compute it anyway in the contiguous
  //   Map case, so we just compute it here.
  template<class Map>
  bool
  areGidsGloballyContiguous (Teuchos::ArrayRCP<typename Map::global_ordinal_type>& allMinGids,
                             Teuchos::ArrayRCP<typename Map::global_ordinal_type>& allMaxGids,
                             typename Map::global_ordinal_type& myMinGid,
                             typename Map::global_ordinal_type& myMaxGid,
                             typename Teuchos::ArrayView<const typename Map::global_ordinal_type>::size_type& myLastInitContigArrayIndex,
                             const Teuchos::ArrayView<const typename Map::global_ordinal_type>& myGids,
                             const Teuchos::Ptr<const Teuchos::Comm<int> >& comm,
                             const Teuchos::Ptr<typename Map::node_type>& node)
  {
    using Teuchos::arcp;
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::as;
    using Teuchos::gatherAll;
    using Teuchos::null;
    using Teuchos::outArg;
    using Teuchos::ptr;
    using Teuchos::REDUCE_MAX;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    typedef typename Map::global_ordinal_type GO;
    typedef typename ArrayView<const GO>::size_type size_type;

    allMinGids = null;
    allMaxGids = null;

    // Are my GIDs (locally) contiguous?  Also, get other info.
    const bool locallyContiguous =
      areMyGidsLocallyContiguous (myMinGid, myMaxGid,
                                  myLastInitContigArrayIndex, myGids, node);

    // Are all processes' GIDs locally contiguous?  That is a
    // necessary but not sufficient condition for them to be globally
    // contiguous.
    const bool allLocallyContiguous = allProcsTrue (*comm, locallyContiguous);
    if (! allLocallyContiguous) {
      return false;
    }
    //
    // Gather each process' min and max GID onto all processes.  Do so
    // via a temporary array that compresses two all-gathers into one.
    //
    const int numProcs = comm->getSize();
    // Leave an extra space at the end to put one plus the global max GID.
    allMinGids = arcp<GO> (as<size_type> (numProcs+1));
    allMaxGids = arcp<GO> (as<size_type> (numProcs));
    {
      Array<GO> allMinMax (as<size_type> (numProcs * 2));
      GO myMinMax[2];
      myMinMax[0] = myMinGid;
      myMinMax[1] = myMaxGid;
      gatherAll<int,GO> (*comm, 2, &myMinMax, 2, &allMinMax[0]);
      // Unpack into separate mins and maxes arrays.
      for (size_type k = 0; k < numProcs; ++k) {
        allMinGids[k] = allMinMax[2*k];
        allMaxGids[k] = allMinMax[2*k+1];
      }
      // Set the iteration guard (see public documentation).
      allMinGids[numProcs] = allMaxGids[numProcs-1] + as<GO> (1);
    }
    //
    // The GIDs are globally contiguous when:
    // 1. Each process' GIDs are contiguous
    // 2. The previous max GID and the current min GID are contiguous
    //
    bool globallyContiguous = true;
    if (! allLocallyContiguous) {
      globallyContiguous = false;
    }
    else {
      // We know there is at least one process in the communicator.
      GO prevMaxGid = allMaxGids[0];
      for (size_type k = 1; k < allMinGids.size(); ++k) {
        const GO curMinGid = allMinGids[k];
        const GO curMaxGid = allMaxGids[k];
        if (prevMaxGid != curMinGid && prevMaxGid != curMinGid + as<GO>(1)) {
          globallyContiguous = false;
          break;
        }
        prevMaxGid = curMaxGid;
      }
    }
    return globallyContiguous;
  }
} // namespace (anonymous)

//
// compute isDistributed. it will be global.
// min/max GID are always computed (from indexBase and numGlobal), and are therefore globally coherent as long as deps are.
// indexBase and numGlobal must always be verified.
// isContiguous is true for the "easy" constructors, assume false for the expert constructor
//
// so we explicitly check    : isCont, numGlobal, indexBase
// then we explicitly compute: MinAllGID, MaxAllGID
// Data explicitly checks    : isDistributed

namespace Tpetra {

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Map<LocalOrdinal,GlobalOrdinal,Node>::
  Map (global_size_t numGlobalElements_in,
       GlobalOrdinal indexBase_in,
       const Teuchos::RCP<const Teuchos::Comm<int> > &comm_in,
       LocalGlobal lOrG,
       const Teuchos::RCP<Node> &node_in)
  : comm_(comm_in)
  , node_(node_in)
  {
    using Teuchos::as;
    using Teuchos::broadcast;
    using Teuchos::OrdinalTraits;
    using Teuchos::outArg;
    using Teuchos::typeName;
    typedef GlobalOrdinal GO;

    // Distribute the elements across the processes in the given
    // communicator so that global IDs (GIDs) are
    //
    // - Nonoverlapping (only one process owns each GID)
    // - Contiguous (the sequence of GIDs is nondecreasing, and no two
    //   adjacent GIDs differ by more than one)
    // - As evenly distributed as possible (the numbers of GIDs on two
    //   different processes do not differ by more than one)

    const global_size_t GSTI = OrdinalTraits<global_size_t>::invalid();

    std::string errPrefix = typeName (*this) +
      " constructor (numGlobalElements, indexBase, comm, lOrG, node): ";

    if (lOrG == GloballyDistributed) {
      const int numImages = comm_->getSize();
      const int myImageID = comm_->getRank();

      // This particular Map constructor requires that the given
      // number of global elements (numGlobalElements) be valid; it
      // does not compute that number, as other Map constructors do.
      // All processes in the given communicator must provide the same
      // numGlobalElements and indexBase values.  (We check this by
      // broadcasting Rank 0's values and comparing with the local
      // values.)
      global_size_t rootNGE = numGlobalElements_in;
      GlobalOrdinal rootIB  = indexBase_in;
      Teuchos::broadcast<int,global_size_t> (*comm_, 0, &rootNGE);
      Teuchos::broadcast<int,GlobalOrdinal> (*comm_, 0, &rootIB);
      int localChecks[2], globalChecks[2];
      localChecks[0] = -1;   // fail or pass
      localChecks[1] = 0;    // fail reason
      if (numGlobalElements_in != rootNGE) {
        localChecks[0] = myImageID;
        localChecks[1] = 1;
      }
      else if (indexBase_in != rootIB) {
        localChecks[0] = myImageID;
        localChecks[1] = 2;
      }
      // REDUCE_MAX will give us the rank ("image ID") of the
      // highest-rank process that DID NOT pass, as well as the
      // reason.  These will be -1 resp. 0 if all processes passed.
      Teuchos::reduceAll<int,int>(*comm_,Teuchos::REDUCE_MAX,2,localChecks,globalChecks);
      if (globalChecks[0] != -1) {
        if (globalChecks[1] == 1) {
          TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
              errPrefix << "numGlobal must be the same on all nodes (examine node " << globalChecks[0] << ").");
        }
        else if (globalChecks[1] == 2) {
          TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
              errPrefix << "indexBase must be the same on all nodes (examine node " << globalChecks[0] << ").");
        }
        else {
          // logic error on our part
          TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
              errPrefix << "logic error. Please contact the Tpetra team.");
        }
      }
      // All processes have the same numGlobalElements, but we still
      // need to check that it is valid.  numGlobalElements must be
      // positive and not the "invalid" value (GSTI).
      //
      // This comparison looks funny, but it avoids compiler warnings
      // for comparing unsigned integers (numGlobalElements_in is a
      // global_size_t, which is unsigned).
      TEUCHOS_TEST_FOR_EXCEPTION((numGlobalElements_in < 1 && numGlobalElements_in != 0) || numGlobalElements_in == GSTI, std::invalid_argument,
          errPrefix << "numGlobalElements (== " << rootNGE << ") must be >= 0.");

      indexBase_ = rootIB;
      numGlobalElements_ = rootNGE;

      /* compute numLocalElements
         We can write numGlobalElements as follows:
         numGlobalElements == numImages * B + remainder
         Each image is allocated elements as follows:
         [ B+1    iff myImageID <  remainder
         numLocalElements = [
         [ B      iff myImageID >= remainder
         In the case that remainder == 0, then all images fall into the
         latter case: numLocalElements == B == numGlobalElements / numImages
         It can then be shown that
         numImages
         \Sum      numLocalElements_i  == numGlobalElements
         i=0
         This strategy is simple, requires no communication, and is optimal vis-a-vis
         uniform distribution of elements.
         This strategy is valid for any value of numGlobalElements and numImages,
         including the following border cases:
         - numImages == 1         -> remainder == 0 && numGlobalElements == numLocalElements
         - numelements < numImages -> remainder == numGlobalElements && numLocalElements \in [0,1]
       */
      numLocalElements_ = as<size_t>(numGlobalElements_ / as<global_size_t>(numImages));
      int remainder = as<int>(numGlobalElements_ % as<global_size_t>(numImages));
#ifdef HAVE_TEUCHOS_DEBUG
      // The above code assumes truncation. Is that safe?
      SHARED_TEST_FOR_EXCEPTION(numLocalElements_ * numImages + remainder != numGlobalElements_,
          std::logic_error, "Tpetra::Map::constructor(numGlobalElements,indexBase,comm,localOrGlobal,node): GlobalOrdinal does not implement division with truncation."
          << " Please contact Tpetra team.",*comm_);
#endif
      GlobalOrdinal start_index;
      if (myImageID < remainder) {
        ++numLocalElements_;
        /* the myImageID images before were each allocated
           numGlobalElements/numImages+1
           ergo, my offset is:
           myImageID * (numGlobalElements/numImages+1)
           == myImageID * numLocalElements
         */
        start_index = as<GlobalOrdinal>(myImageID) * as<GlobalOrdinal>(numLocalElements_);
      }
      else {
        /* a quantity (equal to remainder) of the images before me
           were each allocated
           numGlobalElements/numImages+1
           elements. a quantity (equal to myImageID-remainder) of the remaining
           images before me were each allocated
           numGlobalElements/numImages
           elements. ergo, my offset is:
           remainder*(numGlobalElements/numImages+1) + (myImageID-remainder)*numGlobalElements/numImages
           == remainder*numLocalElements + remainder + myImageID*numLocalElements - remainder*numLocalElements
           == myImageID*numLocalElements + remainder
         */
        start_index = as<GlobalOrdinal>(myImageID)*as<GlobalOrdinal>(numLocalElements_) + as<GlobalOrdinal>(remainder);
      }

      // compute the min/max global IDs
      minMyGID_  = start_index + indexBase_;
      maxMyGID_  = minMyGID_ + numLocalElements_ - 1;
      minAllGID_ = indexBase_;
      maxAllGID_ = indexBase_ + numGlobalElements_ - 1;
      contiguous_ = true;
      distributed_ = (numImages > 1 ? true : false);
    }
    else {  // lOrG == LocallyReplicated
      // compute the min/max global IDs
      indexBase_ = indexBase_in;
      numGlobalElements_ = numGlobalElements_in;
      numLocalElements_  = as<size_t>(numGlobalElements_in);
      minAllGID_ = indexBase_;
      maxAllGID_ = indexBase_ + numGlobalElements_ - 1;
      minMyGID_  = minAllGID_;
      maxMyGID_  = maxAllGID_;
      contiguous_ = true;
      distributed_ = false;
    }
    setupDirectory();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Map<LocalOrdinal,GlobalOrdinal,Node>::
  Map (global_size_t numGlobalElements_in,
       size_t numLocalElements_in,
       GlobalOrdinal indexBase_in,
       const Teuchos::RCP<const Teuchos::Comm<int> > &comm_in,
       const Teuchos::RCP<Node> &node_in)
  : comm_(comm_in)
  , node_(node_in)
  {
    using Teuchos::as;
    using Teuchos::outArg;
    using Teuchos::OrdinalTraits;
    using Teuchos::REDUCE_MAX;
    using Teuchos::REDUCE_SUM;
    using Teuchos::reduceAll;

    // Distribute the elements across the nodes so that they are
    // - non-overlapping
    // - contiguous
    // This differs from Map(Ord,Ord,Plat) in that the user has
    // specified the number of elements per node, so that they are not
    // (necessarily) evenly distributed.

    const global_size_t GSTI = OrdinalTraits<global_size_t>::invalid();

    std::string errPrefix;
    errPrefix = Teuchos::typeName(*this) + "::constructor(numGlobal,numLocal,indexBase,platform): ";

    const int myImageID = comm_->getRank();

    { // begin scoping block
      // for communicating failures
      int localChecks[2], globalChecks[2];
      /* compute the global size
         we are computing the number of global elements because exactly ONE of the following is true:
         - the user didn't specify it, and we need it
         - the user did specify it, but we need to
           + validate it against the sum of the local sizes, and
           + ensure that it is the same on all nodes
       */
      global_size_t global_sum;
      reduceAll<int,global_size_t> (*comm_,
                                    REDUCE_SUM,
                                    as<global_size_t> (numLocalElements_in),
                                    outArg (global_sum));
      /* there are three errors we should be detecting:
         - numGlobalElements != invalid() and it is incorrect/invalid
         - numLocalElements invalid (<0)
      */
      localChecks[0] = -1;
      localChecks[1] = 0;
      if (numLocalElements_in < 1 && numLocalElements_in != 0) {
        // invalid
        localChecks[0] = myImageID;
        localChecks[1] = 1;
      }
      else if (numGlobalElements_in < 1 && numGlobalElements_in != 0 && numGlobalElements_in != GSTI) {
        // invalid
        localChecks[0] = myImageID;
        localChecks[1] = 2;
      }
      else if (numGlobalElements_in != GSTI && numGlobalElements_in != global_sum) {
        // incorrect
        localChecks[0] = myImageID;
        localChecks[1] = 3;
      }
      // now check that indexBase is equivalent across images
      GlobalOrdinal rootIB = indexBase_in;
      Teuchos::broadcast<int,GlobalOrdinal>(*comm_,0,&rootIB);   // broadcast one ordinal from node 0
      if (indexBase_in != rootIB) {
        localChecks[0] = myImageID;
        localChecks[1] = 4;
      }
      // REDUCE_MAX will give us the image ID of the highest rank proc that DID NOT pass
      // this will be -1 if all procs passed
      reduceAll<int,int> (*comm_, REDUCE_MAX, 2, localChecks, globalChecks);
      if (globalChecks[0] != -1) {
        if (globalChecks[1] == 1) {
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
            errPrefix << "numLocal is not valid on at least one node (possibly node "
            << globalChecks[0] << ").");
        }
        else if (globalChecks[1] == 2) {
          TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
            errPrefix << "numGlobal is not valid on at least one node (possibly node "
            << globalChecks[0] << ").");
        }
        else if (globalChecks[1] == 3) {
          TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
            errPrefix << "Input global number of elements = " << numGlobalElements_in
            << " doesn't match sum of numLocal (== " << global_sum << ") on at least "
            "one node (possibly node " << globalChecks[0] << ").");
        }
        else if (globalChecks[1] == 4) {
          TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
            errPrefix << "indexBase is not the same on all nodes (examine node "
            << globalChecks[0] << ").");
        }
        else {
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, errPrefix
            << "Should never get here!  localChecks = " << localChecks[0] << ","
            << localChecks[1] << " and globalChecks = " << globalChecks[0]
            << "," << globalChecks[1] << ".  Please report this bug to the "
            "Tpetra developers.");
        }
      }
      // set numGlobalElements
      if (numGlobalElements_in == GSTI) {
        numGlobalElements_ = global_sum;
      }
      else {
        numGlobalElements_ = numGlobalElements_in;
      }
      numLocalElements_ = numLocalElements_in;
      indexBase_ = indexBase_in;
    } // end of scoping block

    // compute my local offset
    GlobalOrdinal start_index;
    Teuchos::scan<int,GlobalOrdinal>(*comm_,Teuchos::REDUCE_SUM,numLocalElements_,Teuchos::outArg(start_index));
    start_index -= numLocalElements_;

    minAllGID_ = indexBase_;
    maxAllGID_ = indexBase_ + numGlobalElements_ - 1;
    minMyGID_ = start_index + indexBase_;
    maxMyGID_ = minMyGID_ + numLocalElements_ - 1;
    contiguous_ = true;
    distributed_ = checkIsDist();
    setupDirectory();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Map<LocalOrdinal,GlobalOrdinal,Node>::
  Map (global_size_t numGlobalElements_in,
       const Teuchos::ArrayView<const GlobalOrdinal> &entryList,
       GlobalOrdinal indexBase_in,
       const Teuchos::RCP<const Teuchos::Comm<int> > &comm_in,
       const Teuchos::RCP<Node> &node_in)
  : comm_(comm_in)
  , node_(node_in)
  {
    using Teuchos::as;
    using Teuchos::broadcast;
    using Teuchos::outArg;
    using Teuchos::ptr;
    using Teuchos::REDUCE_MAX;
    using Teuchos::REDUCE_MIN;
    using Teuchos::REDUCE_SUM;
    using Teuchos::reduceAll;
    using Teuchos::typeName;

    // The user has specified the distribution of elements over the
    // nodes, via entryList.  The distribution is not necessarily
    // contiguous or equally shared over the nodes.

    const global_size_t GSTI = Teuchos::OrdinalTraits<global_size_t>::invalid();

    // The length of entryList on this node is the number of local
    // elements (on this node), even though entryList contains global
    // indices.  We assume that the number of local elements can be
    // stored in a size_t; numLocalElements_ is a size_t, so this
    // variable and that should have the same type.
    const size_t numLocalElements_in = as<size_t> (entryList.size ());

    const std::string errPrefix = typeName (*this) +
      "::Map(numGlobal,entryList,indexBase,comm,node): ";

    const int myImageID = comm_->getRank();
    { // Begin scoping block for communicating failures.
      int localChecks[2], globalChecks[2];

      // Compute the global number of elements.
      //
      // We are doing this because exactly ONE of the following is true:
      // * the user didn't specify it, and we need it
      // * the user _did_ specify it, but we need to
      // ** validate it against the sum of the local sizes, and
      // ** ensure that it is the same on all nodes
      //
      global_size_t global_sum; // Global number of elements
      reduceAll (*comm_, REDUCE_SUM, as<global_size_t> (numLocalElements_in),
                 outArg (global_sum));
      // localChecks[0] == -1 means that the calling process did not
      // detect an error.  If the calling process did detect an error,
      // it sets localChecks[0] to its rank, and localChecks[1] to the
      // type of error.  Later, we will do a max reduce-all over both
      // elements of localChecks, to find out if at least one process
      // reported an error.
      localChecks[0] = -1;
      localChecks[1] = 0;
      // If the user supplied the number of global elements (i.e., if
      // it's not invalid (== GSTI)), then make sure that there is at
      // least one global element.  The first two clauses of the test
      // are apparently redundant, but help avoid compiler warnings
      // about comparing signed and unsigned integers.
      if (numGlobalElements_in < 1 &&
          numGlobalElements_in != 0 &&
          numGlobalElements_in != GSTI) {
        // Number of global elements is not the "invalid" value
        // (meaning "let the constructor compute it"), and is
        // negative.  That's not a valid input.
        localChecks[0] = myImageID;
        localChecks[1] = 1; // "Negative number of global elements given"
      }
      else if (numGlobalElements_in != GSTI && numGlobalElements_in != global_sum) {
        // If the user specifies a number of global elements not equal
        // to the "invalid" value, then we assume that this equals the
        // sum of all the local counts of elements (including overlap)
        // on all processes in the communicator.  If not, then we
        // report an error.
        localChecks[0] = myImageID;
        localChecks[1] = 2; // "Incorrect number of global elements given"
      }
      // Check that all nodes have the same indexBase value.  We do so
      // by broadcasting the indexBase value from Proc 0 to all the
      // processes, and then checking locally on each process whether
      // it's the same as indexBase_in.  (This does about half as much
      // communication as an all-reduce.)
      GlobalOrdinal rootIndexBase = indexBase_in;
      const int rootRank = 0;
      broadcast (*comm_, rootRank, outArg(rootIndexBase));

      if (indexBase_in != rootIndexBase) {
        localChecks[0] = myImageID;
        localChecks[1] = 3; // "indexBase values not the same on all processes"
      }
      // After the reduceAll below, globalChecks[0] will be -1 if all
      // processes passed their tests, else it will be the rank
      // ("image ID") of the highest-rank process that did NOT pass.
      // In the latter case, globalChecks[1] will be an error code, though
      // not necessarily the error code experienced by that node.
      reduceAll (*comm_, REDUCE_MAX, 2, localChecks, globalChecks);
      if (globalChecks[0] != -1) {
        if (globalChecks[1] == 1) {
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
            errPrefix << "You specified a negative number of global elements "
            "(numGlobalElements_in argument to the Map constructor) on at least "
            "one of the processes in the input communicator (including process "
            << globalChecks[0] << ").  If you want Map's constructor to "
            "compute the global number of elements, you must set the "
            "numGlobaElements_in argument to Teuchos::Ordinal_Traits"
            "<global_size_t>::invalid() on all processes in the communicator.");
        }
        else if (globalChecks[1] == 2) {
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
            errPrefix << "On at least one process in the input communicator "
            "(including process " << globalChecks[0] << ", the given number of "
            "global elements (numGlobalElements_in argument to the Map "
            "constructor) does not match the sum " << global_sum << " of the "
            "number of elements on each process.  The latter is the sum of "
            "entryList.size() over all processes in the communicator; elements "
            "that overlap over multiple processes or that are duplicated on "
            "the same process are counted multiple times.");
        }
        else if (globalChecks[1] == 3) {
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
            errPrefix << "The given values for the index base (indexBase_in "
            "argument to the Map constructor) do not match on all the processes."
            "  This includes process " << globalChecks[0] << ".");
        }
        else {
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
            errPrefix << "Should never get here!  globalChecks[0] == "
            << globalChecks[0] << " and globalChecks[1] == " << globalChecks[1]
            << ".  Please report this bug to the Tpetra developers.");
        }
      }
      //
      // We've successfully validated or computed the number of global
      // elements, and validated the index base.
      //
      if (numGlobalElements_in == GSTI) {
        numGlobalElements_ = global_sum;
      }
      else {
        numGlobalElements_ = numGlobalElements_in;
      }
      numLocalElements_ = numLocalElements_in;
      indexBase_ = indexBase_in;
    } // end scoping block

    // Assume for now that there are numLocalElements (there may be
    // less, if some GIDs are duplicated in entryList).
    // NOTE: cgb: This could result in an incorrecdt minMyGID_ if indexBase_ is less than all of the GIDs
    minMyGID_ = indexBase_;
    maxMyGID_ = indexBase_;
    //
    // Create the GID to LID map.  Do not assume that the GIDs in
    // entryList are distinct.  In the case that a GID is duplicated,
    // each duplication gets a new LID.
    //
    glMap_ = rcp(new global_to_local_table_type());
    if (numLocalElements_ > 0) {
      lgMap_ = Teuchos::arcp<GlobalOrdinal>(numLocalElements_);
      // While iterating through entryList, we compute its (local)
      // min and max elements.
      minMyGID_ = entryList[0];
      maxMyGID_ = entryList[0];
      for (size_t i=0; i < numLocalElements_; i++) {
        lgMap_[i] = entryList[i];      // lgMap_:  LID to GID
        (*glMap_)[entryList[i]] = i;   // glMap_: GID to LID
        if (entryList[i] < minMyGID_) {
          minMyGID_ = entryList[i];
        }
        if (entryList[i] > maxMyGID_) {
          maxMyGID_ = entryList[i];
        }
      }
    }

    // Compute the min and max of all processes' global IDs.
    reduceAll (*comm_, REDUCE_MIN, minMyGID_, outArg (minAllGID_));
    reduceAll (*comm_, REDUCE_MAX, maxMyGID_, outArg (maxAllGID_));
    contiguous_  = false;
    distributed_ = checkIsDist();

    using std::endl;
    TEUCHOS_TEST_FOR_EXCEPTION(
      minAllGID_ < indexBase_,
      std::invalid_argument,
      errPrefix << std::endl << "Minimum global ID (== " << minAllGID_
      << ") over all process(es) is less than the given indexBase (== "
      << indexBase_ << ").");
#if 0
    using std::cerr;
    cerr << "Tpetra::Map: Just for your information:"
         << endl << "- minMyGID_ = " << minMyGID_
         << endl << "- minAllGID_ = " << minAllGID_
         << endl << "- maxMyGID_ = " << maxMyGID_
         << endl << "- maxAllGID_ = " << maxAllGID_
         << endl << "- numLocalElements_ = " << numLocalElements_
         << endl << "- numGlobalElements_ = " << numGlobalElements_
         << endl << "- entryList = " << toString (entryList)
         << endl;
#endif // 0

    setupDirectory();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Map<LocalOrdinal,GlobalOrdinal,Node>::~Map ()
  {}

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  LocalOrdinal Map<LocalOrdinal,GlobalOrdinal,Node>::getLocalElement(GlobalOrdinal globalIndex) const
  {
    if (isContiguous()) {
      if (globalIndex < getMinGlobalIndex() || globalIndex > getMaxGlobalIndex()) {
        return Teuchos::OrdinalTraits<LocalOrdinal>::invalid();
      }
      return Teuchos::as<LocalOrdinal>(globalIndex - getMinGlobalIndex());
    }
    else {
      typedef typename global_to_local_table_type::const_iterator iter_type;
      iter_type i = glMap_->find(globalIndex);
      if (i == glMap_->end()) {
        return Teuchos::OrdinalTraits<LocalOrdinal>::invalid();
      }
      return i->second;
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  GlobalOrdinal Map<LocalOrdinal,GlobalOrdinal,Node>::getGlobalElement(LocalOrdinal localIndex) const {
    if (localIndex < getMinLocalIndex() || localIndex > getMaxLocalIndex()) {
      return Teuchos::OrdinalTraits<GlobalOrdinal>::invalid();
    }
    if (isContiguous()) {
      return getMinGlobalIndex() + localIndex;
    }
    else {
      return lgMap_[localIndex];
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool Map<LocalOrdinal,GlobalOrdinal,Node>::isNodeLocalElement(LocalOrdinal localIndex) const {
    if (localIndex < getMinLocalIndex() || localIndex > getMaxLocalIndex()) {
      return false;
    }
    return true;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool Map<LocalOrdinal,GlobalOrdinal,Node>::isNodeGlobalElement(GlobalOrdinal globalIndex) const {
    if (isContiguous()) {
      return (getMinGlobalIndex() <= globalIndex) && (globalIndex <= getMaxGlobalIndex());
    }
    else {
      typedef typename global_to_local_table_type::const_iterator iter_type;
      iter_type i = glMap_->find(globalIndex);
      return (i != glMap_->end());
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool Map<LocalOrdinal,GlobalOrdinal,Node>::isContiguous() const {
    return contiguous_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool Map<LocalOrdinal,GlobalOrdinal,Node>::isCompatible (const Map<LocalOrdinal,GlobalOrdinal,Node> &map) const {
    using Teuchos::outArg;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;

    // Do both Maps have the same number of elements, both globally
    // and on the calling process?
    char locallyCompat = 0;
    if (getGlobalNumElements() != map.getGlobalNumElements() ||
          getNodeNumElements() != map.getNodeNumElements()) {
      locallyCompat = 0; // NOT compatible on this process
    }
    else {
      locallyCompat = 1; // compatible on this process
    }

    char globallyCompat = 0;
    reduceAll (*comm_, REDUCE_MIN, locallyCompat, outArg (globallyCompat));
    return (globallyCompat == 1);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  Map<LocalOrdinal,GlobalOrdinal,Node>::
  isSameAs (const Map<LocalOrdinal,GlobalOrdinal,Node> &map) const
  {
    using Teuchos::outArg;
    if (this == &map) {
      // we assume that this is globally coherent
      // if they are the same object, then they are equivalent maps
      return true;
    }

    // check all other globally coherent properties
    // if they do not share each of these properties, then they cannot be
    // equivalent maps
    if ( (getMinAllGlobalIndex()   != map.getMinAllGlobalIndex())   ||
         (getMaxAllGlobalIndex()   != map.getMaxAllGlobalIndex())   ||
         (getGlobalNumElements() != map.getGlobalNumElements()) ||
         (isDistributed()       != map.isDistributed())       ||
         (getIndexBase()        != map.getIndexBase())         )  {
      return false;
    }

    // If we get this far, we need to check local properties and the
    // communicate same-ness across all nodes
    // we prefer local work over communication, ergo, we will perform all
    // comparisons and conduct a single communication
    char isSame_lcl = 1;

    // check number of entries owned by this node
    if (getNodeNumElements() != map.getNodeNumElements()) {
      isSame_lcl = 0;
    }

    // check the identity of the entries owned by this node
    // only do this if we haven't already determined not-same-ness
    if (isSame_lcl == 1) {
      // if they are contiguous, we can check the ranges easily
      // if they are not contiguous, we must check the individual LID -> GID mappings
      // the latter approach is valid in either case, but the former is faster
      if (isContiguous() && map.isContiguous()) {
        if ( (getMinGlobalIndex() != map.getMinGlobalIndex()) ||
             (getMaxGlobalIndex() != map.getMaxGlobalIndex()) ){
          isSame_lcl = 0;
        }
      }
      else {
        /* Note: std::equal requires that the latter range is as large as the former.
         * We know the ranges have equal length, because they have the same number of
         * local entries.
         * However, we only know that lgMap_ has been filled for the Map that is not
         * contiguous (one is potentially contiguous.) Calling getNodeElementList()
         * will create it. */
        Teuchos::ArrayView<const GlobalOrdinal> ge1, ge2;
        ge1 =     getNodeElementList();
        ge2 = map.getNodeElementList();
        if (!std::equal(ge1.begin(),ge1.end(),ge2.begin())) {
          isSame_lcl = 0;
        }
      }
    }

    // now, determine if we detected not-same-ness on any node
    char isSame_gbl;
    Teuchos::reduceAll<int,char>(*comm_,Teuchos::REDUCE_MIN,isSame_lcl,outArg(isSame_gbl));
    return (isSame_gbl == 1);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayView<const GlobalOrdinal>
  Map<LocalOrdinal,GlobalOrdinal,Node>::getNodeElementList() const {
    // If the local-to-global mapping doesn't exist yet, and if we
    // have local entries, then create and fill the local-to-global
    // mapping.
    if (lgMap_.is_null() && numLocalElements_ > 0) {
#ifdef HAVE_TEUCHOS_DEBUG
      // The local-to-global mapping should have been set up already
      // for a noncontiguous map.
      TEUCHOS_TEST_FOR_EXCEPTION( ! isContiguous(), std::logic_error,
        "Tpetra::Map::getNodeElementList: The local-to-global mapping (lgMap_) "
        "should have been set up already for a noncontiguous Map.  Please report"
        " this bug to the Tpetra team.");
#endif // HAVE_TEUCHOS_DEBUG
      lgMap_ = Teuchos::arcp<GlobalOrdinal>(numLocalElements_);
      Teuchos::ArrayRCP<GlobalOrdinal> lgptr = lgMap_;
      for (GlobalOrdinal gid=minMyGID_; gid <= maxMyGID_; ++gid) {
        *(lgptr++) = gid;
      }
    }
    return lgMap_();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool Map<LocalOrdinal,GlobalOrdinal,Node>::isDistributed() const {
    return distributed_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  std::string Map<LocalOrdinal,GlobalOrdinal,Node>::description() const {
    std::ostringstream oss;
    oss << Teuchos::Describable::description();
    oss << "{getGlobalNumElements() = " << getGlobalNumElements()
        << ", getNodeNumElements() = " << getNodeNumElements()
        << ", isContiguous() = " << isContiguous()
        << ", isDistributed() = " << isDistributed()
        << "}";
    return oss.str();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void Map<LocalOrdinal,GlobalOrdinal,Node>::describe( Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const {
    using std::endl;
    using std::setw;
    using Teuchos::VERB_DEFAULT;
    using Teuchos::VERB_NONE;
    using Teuchos::VERB_LOW;
    using Teuchos::VERB_MEDIUM;
    using Teuchos::VERB_HIGH;
    using Teuchos::VERB_EXTREME;

    const size_t nME = getNodeNumElements();
    Teuchos::ArrayView<const GlobalOrdinal> myEntries = getNodeElementList();
    int myImageID = comm_->getRank();
    int numImages = comm_->getSize();

    Teuchos::EVerbosityLevel vl = verbLevel;
    if (vl == VERB_DEFAULT) vl = VERB_LOW;

    size_t width = 1;
    for (size_t dec=10; dec<getGlobalNumElements(); dec *= 10) {
      ++width;
    }
    width = std::max<size_t>(width,Teuchos::as<size_t>(12)) + 2;

    Teuchos::OSTab tab(out);

    if (vl == VERB_NONE) {
      // do nothing
    }
    else if (vl == VERB_LOW) {
      out << this->description() << endl;
    }
    else {  // MEDIUM, HIGH or EXTREME
      for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
        if (myImageID == imageCtr) {
          if (myImageID == 0) { // this is the root node (only output this info once)
            out << endl
                << "Number of Global Entries = " << getGlobalNumElements()  << endl
                << "Maximum of all GIDs      = " << getMaxAllGlobalIndex() << endl
                << "Minimum of all GIDs      = " << getMinAllGlobalIndex() << endl
                << "Index Base               = " << getIndexBase()         << endl;
          }
          out << endl;
          if (vl == VERB_HIGH || vl == VERB_EXTREME) {
            out << "Number of Local Elements   = " << nME           << endl
                << "Maximum of my GIDs         = " << getMaxGlobalIndex() << endl
                << "Minimum of my GIDs         = " << getMinGlobalIndex() << endl;
            out << endl;
          }
          if (vl == VERB_EXTREME) {
            out << std::setw(width) << "Node ID"
                << std::setw(width) << "Local Index"
                << std::setw(width) << "Global Index"
                << endl;
            for (size_t i=0; i < nME; i++) {
              out << std::setw(width) << myImageID
                  << std::setw(width) << i
                  << std::setw(width) << myEntries[i]
                  << endl;
            }
            out << std::flush;
          }
        }
        // Do a few global ops to give I/O a chance to complete
        comm_->barrier();
        comm_->barrier();
        comm_->barrier();
      }
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void Map<LocalOrdinal,GlobalOrdinal,Node>::setupDirectory() {
    if (directory_ == Teuchos::null) {
      directory_ = Teuchos::rcp( new Directory<LocalOrdinal,GlobalOrdinal,Node>(Teuchos::rcp(this,false)) );
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  LookupStatus Map<LocalOrdinal,GlobalOrdinal,Node>::getRemoteIndexList(
                    const Teuchos::ArrayView<const GlobalOrdinal> & GIDList,
                    const Teuchos::ArrayView<int> & imageIDList,
                    const Teuchos::ArrayView<LocalOrdinal> & LIDList) const {
    TEUCHOS_TEST_FOR_EXCEPTION(GIDList.size() > 0 && getGlobalNumElements() == 0, std::runtime_error,
        Teuchos::typeName(*this) << "::getRemoteIndexList(): getRemoteIndexList() cannot be called, zero entries in Map.");
    return directory_->getDirectoryEntries(GIDList, imageIDList, LIDList);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  LookupStatus Map<LocalOrdinal,GlobalOrdinal,Node>::getRemoteIndexList(
                    const Teuchos::ArrayView<const GlobalOrdinal> & GIDList,
                    const Teuchos::ArrayView<int> & imageIDList) const {
    TEUCHOS_TEST_FOR_EXCEPTION(GIDList.size() > 0 && getGlobalNumElements() == 0, std::runtime_error,
        Teuchos::typeName(*this) << "::getRemoteIndexList(): getRemoteIndexList() cannot be called, zero entries in Map.");
    return directory_->getDirectoryEntries(GIDList, imageIDList);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  const Teuchos::RCP<const Teuchos::Comm<int> > &
  Map<LocalOrdinal,GlobalOrdinal,Node>::getComm() const {
    return comm_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  const Teuchos::RCP<Node> &
  Map<LocalOrdinal,GlobalOrdinal,Node>::getNode() const {
    return node_;
  }

  template <class LocalOrdinal,class GlobalOrdinal, class Node>
  bool Map<LocalOrdinal,GlobalOrdinal,Node>::checkIsDist() const {
    using Teuchos::outArg;
    bool global = false;
    if(comm_->getSize() > 1) {
      // The communicator has more than one process, but that doesn't
      // necessarily mean the Map is distributed.
      char localRep = 0;
      if (numGlobalElements_ == Teuchos::as<global_size_t>(numLocalElements_)) {
        // The number of local elements on this process equals the
        // number of global elements.
        //
        // NOTE (mfh 22 Nov 2011) Does this still work if there were
        // duplicates in the global ID list on input (the third Map
        // constructor), so that the number of local elements (which
        // are not duplicated) on this process could be less than the
        // number of global elements, even if this process owns all
        // the elements?
        localRep = 1;
      }
      char allLocalRep;
      Teuchos::reduceAll<int>(*comm_,Teuchos::REDUCE_MIN,localRep,outArg(allLocalRep));
      if (allLocalRep != 1) {
        // At least one process does not own all the elements.
        // This makes the Map a distributed Map.
        global = true;
      }
    }
    // If the communicator has only one process, then the Map is not
    // distributed.
    return global;
  }

} // Tpetra namespace

template <class LocalOrdinal, class GlobalOrdinal>
Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Kokkos::DefaultNode::DefaultNodeType> >
Tpetra::createLocalMap(size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm) {
  return createLocalMapWithNode<LocalOrdinal,GlobalOrdinal,Kokkos::DefaultNode::DefaultNodeType>(numElements, comm, Kokkos::DefaultNode::getDefaultNode());
}

template <class LocalOrdinal, class GlobalOrdinal>
Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Kokkos::DefaultNode::DefaultNodeType> >
Tpetra::createUniformContigMap(global_size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm) {
  return createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal,Kokkos::DefaultNode::DefaultNodeType>(numElements, comm, Kokkos::DefaultNode::getDefaultNode());
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
Tpetra::createUniformContigMapWithNode(global_size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP<Node> &node) {
  Teuchos::RCP< Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > map;
  map = Teuchos::rcp( new Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>(numElements,                           // num elements, global and local
                                                              Teuchos::OrdinalTraits<GlobalOrdinal>::zero(),  // index base is zero
                                                              comm, GloballyDistributed, node) );
  return map.getConst();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
Tpetra::createLocalMapWithNode(size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node) {
  Teuchos::RCP< Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > map;
  map = Teuchos::rcp( new Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>((Tpetra::global_size_t)numElements,    // num elements, global and local
                                                              Teuchos::OrdinalTraits<GlobalOrdinal>::zero(),  // index base is zero
                                                              comm, LocallyReplicated, node) );
  return map.getConst();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
Tpetra::createContigMapWithNode(Tpetra::global_size_t numElements, size_t localNumElements,
                                const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node) {
  Teuchos::RCP< Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > map;
  map = Teuchos::rcp( new Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>(numElements,localNumElements,
                                                              Teuchos::OrdinalTraits<GlobalOrdinal>::zero(),  // index base is zero
                                                              comm, node) );
  return map.getConst();
}

template <class LocalOrdinal, class GlobalOrdinal>
Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Kokkos::DefaultNode::DefaultNodeType> >
Tpetra::createContigMap(Tpetra::global_size_t numElements, size_t localNumElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm) {
  return Tpetra::createContigMapWithNode<LocalOrdinal,GlobalOrdinal,Kokkos::DefaultNode::DefaultNodeType>(numElements, localNumElements, comm, Kokkos::DefaultNode::getDefaultNode() );
}


template <class LocalOrdinal, class GlobalOrdinal>
Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Kokkos::DefaultNode::DefaultNodeType> >
Tpetra::createNonContigMap(const Teuchos::ArrayView<const GlobalOrdinal> &elementList,
                           const Teuchos::RCP<const Teuchos::Comm<int> > &comm)
{
  return Tpetra::createNonContigMapWithNode<LocalOrdinal,GlobalOrdinal,Kokkos::DefaultNode::DefaultNodeType>(elementList, comm, Kokkos::DefaultNode::getDefaultNode() );
}


template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
Tpetra::createNonContigMapWithNode(const Teuchos::ArrayView<const GlobalOrdinal> &elementList,
                                   const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
                                   const Teuchos::RCP<Node> &node)
{
  Teuchos::RCP< Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > map;
  map = rcp(new Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>(
                      Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
                      elementList,
                      Teuchos::OrdinalTraits<Tpetra::global_size_t>::zero(),
                      comm, node ) );
  return map;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
Tpetra::createWeightedContigMapWithNode(int myWeight, Tpetra::global_size_t numElements,
                                        const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node) {
  Teuchos::RCP< Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > map;
  int sumOfWeights, elemsLeft, localNumElements;
  const int numImages = comm->getSize(),
            myImageID = comm->getRank();
  Teuchos::reduceAll<int>(*comm,Teuchos::REDUCE_SUM,myWeight,Teuchos::outArg(sumOfWeights));
  const double myShare = ((double)myWeight) / ((double)sumOfWeights);
  localNumElements = (int)std::floor( myShare * ((double)numElements) );
  // std::cout << "numElements: " << numElements << "  myWeight: " << myWeight << "  sumOfWeights: " << sumOfWeights << "  myShare: " << myShare << std::endl;
  Teuchos::reduceAll<int>(*comm,Teuchos::REDUCE_SUM,localNumElements,Teuchos::outArg(elemsLeft));
  elemsLeft = numElements - elemsLeft;
  // std::cout << "(before) localNumElements: " << localNumElements << "  elemsLeft: " << elemsLeft << std::endl;
  // i think this is true. just test it for now.
  TEUCHOS_TEST_FOR_EXCEPT(elemsLeft < -numImages || numImages < elemsLeft);
  if (elemsLeft < 0) {
    // last elemsLeft nodes lose an element
    if (myImageID >= numImages-elemsLeft) --localNumElements;
  }
  else if (elemsLeft > 0) {
    // first elemsLeft nodes gain an element
    if (myImageID < elemsLeft) ++localNumElements;
  }
  // std::cout << "(after) localNumElements: " << localNumElements << std::endl;
  return createContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(numElements,localNumElements,comm,node);

}

template<class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
Tpetra::createOneToOne (Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > &M)
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef Tpetra::Map<LO,GO,Node> map_type;
  int myID = M->getComm()->getRank();

  //Based off Epetra's one to one.

  Tpetra::Directory<LO, GO, Node> directory(M);

  size_t numMyElems = M->getNodeNumElements();

  ArrayView<const GlobalOrdinal> myElems = M->getNodeElementList();

  Array<int> owner_procs_vec (numMyElems);

  directory.getDirectoryEntries(myElems, owner_procs_vec ());

  Array<GO> myOwned_vec (numMyElems);
  size_t numMyOwnedElems = 0;

  for(size_t i=0; i<numMyElems; ++i)
  {
    GO GID = myElems[i];
    int owner = owner_procs_vec[i];

    if (myID == owner)
    {
      myOwned_vec[numMyOwnedElems++]=GID;
    }
  }
  myOwned_vec.resize (numMyOwnedElems);

  RCP< const Tpetra::Map<LO,GO,Node> > one_to_one_map =
    rcp (new map_type (Teuchos::OrdinalTraits<GO>::invalid (), myOwned_vec (),
                       M->getIndexBase (), M->getComm (), M->getNode()));

  return(one_to_one_map);

}
//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

//! Explicit instantiation macro supporting the Map class. Instantiates the class and the non-member constructors.
#define TPETRA_MAP_INSTANT(LO,GO,NODE) \
  \
  template class Map< LO , GO , NODE >; \
  \
  template Teuchos::RCP< const Map<LO,GO,NODE> > \
  createLocalMapWithNode<LO,GO,NODE>(size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< NODE > &node); \
  \
  template Teuchos::RCP< const Map<LO,GO,NODE> > \
  createContigMapWithNode<LO,GO,NODE>(global_size_t numElements, size_t localNumElements, \
                                      const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< NODE > &node); \
  \
  template Teuchos::RCP< const Map<LO,GO,NODE> > \
  createNonContigMapWithNode(const Teuchos::ArrayView<const GO> &elementList, \
                             const RCP<const Teuchos::Comm<int> > &comm,      \
                             const RCP<NODE> &node);                          \
  template Teuchos::RCP< const Map<LO,GO,NODE> > \
  createUniformContigMapWithNode<LO,GO,NODE>(global_size_t numElements, \
                                             const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< NODE > &node); \
  \
  template Teuchos::RCP< const Map<LO,GO,NODE> > \
  createWeightedContigMapWithNode<LO,GO,NODE>(int thisNodeWeight, global_size_t numElements, \
                                              const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< NODE > &node); \
  \
  template Teuchos::RCP<const Map<LO,GO,NODE> > \
  createOneToOne (Teuchos::RCP<const Map<LO,GO,NODE> > &M);


//! Explicit instantiation macro supporting the Map class, on the default node for specified ordinals.
#define TPETRA_MAP_INSTANT_DEFAULTNODE(LO,GO) \
  template Teuchos::RCP< const Map<LO,GO> > \
  createLocalMap<LO,GO>( size_t, const Teuchos::RCP< const Teuchos::Comm< int > > &); \
  \
  template Teuchos::RCP< const Map<LO,GO> > \
  createContigMap<LO,GO>( global_size_t, size_t, \
                          const Teuchos::RCP< const Teuchos::Comm< int > > &); \
  \
  template Teuchos::RCP< const Map<LO,GO> >  \
  createNonContigMap(const Teuchos::ArrayView<const GO> &,          \
                     const RCP<const Teuchos::Comm<int> > &); \
  \
  template Teuchos::RCP< const Map<LO,GO> >  \
  createUniformContigMap<LO,GO>( global_size_t, \
                                 const Teuchos::RCP< const Teuchos::Comm< int > > &); \

#endif // TPETRA_MAP_DEF_HPP
