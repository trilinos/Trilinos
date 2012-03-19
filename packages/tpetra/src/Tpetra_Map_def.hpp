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

// TODO: make sure that Ordinal values in constructors aren't invalid()

#ifndef TPETRA_MAP_DEF_HPP
#define TPETRA_MAP_DEF_HPP

#include <Teuchos_as.hpp>
#include "Tpetra_Directory.hpp" // must include for implicit instantiation to work
#include "Tpetra_Util.hpp"
#include <stdexcept>

#ifdef DOXYGEN_USE_ONLY
  #includee "Tpetra_Map_decl.hpp"
#endif

namespace {
  // \fn areMyGidsLocallyContiguous
  // \brief Are my global IDs (GIDs) contiguous, and if so, what is my process' minimum GID?
  //
  // If the first return value is true, then the second return value
  // is my process' minimum GID.  Otherwise, the second return value
  // is invalid.
  //
  // We merge the min GID computation into the contiguous computation
  // to avoid extra passes over all the GIDs.  This makes computing
  // the min "free" if myGids.size() is sufficiently large.
  template<class Map>
  std::pair<bool, typename Map::global_ordinal_type>
  areMyGidsLocallyContiguous (const Teuchos::ArrayView<const typename Map::global_ordinal_type>& myGids,
			      const Teuchos::Ptr<typename Map::node_type>& /* node */)
  {
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::as;
    typedef typename Map::global_ordinal_type GO;
    typedef typename ArrayView<const GO>::size_type size_type;

    bool locallyContiguous = true; 
    GO minGid = as<GO> (0); // not meaningful if there are no GIDs
    if (myGids.size() > 0) { // trivially true if size is zero
      GO prevGid = myGids[0];
      minGid = prevGid;
      for (size_type k = 1; k < myGids.size(); ++k) {
	const GO curGid = myGids[k];
	// Allow GO to be unsigned.  Remember that GIDs need not be
	// sorted, and may have repeated entries (which are still
	// counted as contiguous).
	if (curGid != prevGid && curGid != prevGid + as<GO>(1)) {
	  locallyContiguous = false;
	  break;
	}
	prevGid = curGid;
	if (curGid > prevGid) {
	  minGid = curGid;
	}
      }
    }
    return std::make_pair (locallyContiguous, minGid);
  }

  // \fn areGidsGloballyContiguous
  // \brief Are all processes' global IDs (GIDs) contiguous? and if so, tell me the min GID on each process.
  //
  // Contiguous Maps are much faster and take much less memory than
  // noncontiguous Maps, so it's worthwhile to check for contiguity,
  // even when using the most general Map constructor (that accepts a
  // list on each process of the GIDs that process owns).  One way to
  // check global contiguity is to compute each process' minimum GID.
  // If the GIDs do turn out to be contiguous, the list of min GIDs on
  // each process is useful for constructing the Map's Directory, so
  // we return that so the Directory doesn't need to compute it again.
  //
  // First return value: whether all processes' GIDs are contiguous.
  // If false, then the second return value is invalid (null).  If
  // true, then the second return value is a globally replicated array
  // of size comm->getSize(), for which the p-th entry is the min GID
  // on process p.  In the globally contiguous case, this globally
  // replicated array functions as the "directory" that maps any GID
  // to its owning process without communication.
  //
  // Implementation note: We could just check for global contiguity on
  // Proc 0, with a gather instead of an all-gather.  However, the
  // Directory uses the globally replicated array of min GIDs and
  // would need to compute it anyway in the contiguous Map case, so we
  // just compute it here.
  template<class Map>
  std::pair<bool, Teuchos::ArrayRCP<typename Map::global_ordinal_type> >
  areGidsGloballyContiguous (const Teuchos::ArrayView<const typename Map::global_ordinal_type>& myGids,
			     const Teuchos::Ptr<const Teuchos::Comm<int> >& comm,
			     const Teuchos::Ptr<typename Map::node_type>& node)
  {
    using Teuchos::arcp;
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::as;
    using Teuchos::ptr;
    using Teuchos::REDUCE_MAX;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    typedef typename Map::global_ordinal_type GO;
    typedef typename ArrayView<const GO>::size_type size_type;

    ArrayRCP<GO> allMinGids;

    // Are my GIDs (locally) contiguous?  Also get my min GID, which
    // is useful for checking for global contiguity.
    std::pair<bool, GO> localResult = areMyGidsLocallyContiguous (myGids, node);
    const bool locallyContiguous = localResult.first;
    const GO myMinGid = localResult.second; 

    // Are all processes' GIDs locally contiguous?  That is a
    // necessary but not sufficient condition for them to be globally
    // contiguous.
    bool allLocallyContiguous = false;
    { // reduceAll not yet implemented for bool.
      int locallyContig = locallyContiguous ? 1 : 0;
      int allLocallyContig = 1;
      reduceAll (*comm, REDUCE_MIN, locallyContig, ptr (&allLocallyContig));
      allLocallyContiguous = (allLocallyContig == 1);
    }
    if (! allLocallyContiguous) {
      return std::make_pair (false, allMinGids);
    }

    // Compute each process' min GID.  This is useful if the GIDs are
    // globally contiguous, so we return it as well.  (In the globally
    // contiguous case, this globally replicated array functions as
    // the "directory" that maps any GID to its owning process without
    // communication.)
    const int numProcs = comm->getSize();
    const int myRank   = comm->getRank();
    allMinGids = arcp<GO> (as<typename ArrayRCP<GO>::size_type> (numProcs));
    Teuchos::gatherAll (*comm, 1, &myMinGid, numProcs, &allMinGids.front());

    // If the min GIDs are nondecreasing, and each process' GIDs are
    // contiguous, then the GIDs are globally contiguous.  We know
    // there is at least one process in the communicator.
    bool globallyContiguous = true; // trivially, if numProcs == 1.
    GO prevMinGid = allMinGids[0];
    for (typename ArrayRCP<GO>::size_type k = 1; k < allMinGids.size(); ++k) {
      const GO curMinGid = allMinGids[k];
      if (curMinGid < prevMinGid) {
	globallyContiguous = false;
	break;
      }
    }
    return std::make_pair (globallyContiguous, allMinGids);
  }
} // namespace (anonymous)

/** \file Tpetra_Map_def.hpp 

    The implementations for the members of Tpetra::Map and related non-member constructors.
 */

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
    // The user wants us to distribute the elements across the nodes
    // so that they are
    // - non-overlapping
    // - contiguous
    // - as evenly distributed as possible
    using Teuchos::as;
    using Teuchos::outArg;
    const global_size_t GST0 = Teuchos::OrdinalTraits<global_size_t>::zero();
    const global_size_t GST1 = Teuchos::OrdinalTraits<global_size_t>::one();
    const global_size_t GSTI = Teuchos::OrdinalTraits<global_size_t>::invalid();
    const GlobalOrdinal G1 = Teuchos::OrdinalTraits<GlobalOrdinal>::one();

    std::string errPrefix = Teuchos::typeName(*this) + 
      "::constructor(numGlobal,indexBase,comm,lOrG): ";

    if (lOrG == GloballyDistributed) {
      const int numImages = comm_->getSize();
      const int myImageID = comm_->getRank();

      // This constructor requires that the given number of global
      // elements (numGlobalElements) be valid, since this constructor
      // doesn't compute it.  All processes in the given communicator
      // must provide the same numGlobalElements and indexBase values.
      // (We check this by broadcasting Rank 0's values and comparing
      // with the local values.)
      global_size_t rootNGE = numGlobalElements_in;
      GlobalOrdinal rootIB  = indexBase_in;
      Teuchos::broadcast<int,global_size_t>(*comm_,0,&rootNGE);
      Teuchos::broadcast<int,GlobalOrdinal>(*comm_,0,&rootIB);
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
      TEUCHOS_TEST_FOR_EXCEPTION((numGlobalElements_in < GST1 && numGlobalElements_in != GST0) || numGlobalElements_in == GSTI, std::invalid_argument,
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
      maxMyGID_  = minMyGID_ + numLocalElements_ - G1;
      minAllGID_ = indexBase_;
      maxAllGID_ = indexBase_ + numGlobalElements_ - G1;
      contiguous_ = true;
      distributed_ = (numImages > 1 ? true : false);
    }
    else {  // lOrG == LocallyReplicated
      // compute the min/max global IDs
      indexBase_ = indexBase_in;
      numGlobalElements_ = numGlobalElements_in;
      numLocalElements_  = as<size_t>(numGlobalElements_in);
      minAllGID_ = indexBase_;
      maxAllGID_ = indexBase_ + numGlobalElements_ - G1;
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
    // Distribute the elements across the nodes so that they are 
    // - non-overlapping
    // - contiguous
    // This differs from Map(Ord,Ord,Plat) in that the user has
    // specified the number of elements per node, so that they are not
    // (necessarily) evenly distributed.

    using Teuchos::outArg;

    const size_t  L0 = Teuchos::OrdinalTraits<size_t>::zero();
    const size_t  L1 = Teuchos::OrdinalTraits<size_t>::one();
    const global_size_t GST0 = Teuchos::OrdinalTraits<global_size_t>::zero();
    const global_size_t GST1 = Teuchos::OrdinalTraits<global_size_t>::one();
    const global_size_t GSTI = Teuchos::OrdinalTraits<global_size_t>::invalid();
    const GlobalOrdinal G1 = Teuchos::OrdinalTraits<GlobalOrdinal>::one();

    std::string errPrefix;
    errPrefix = Teuchos::typeName(*this) + "::constructor(numGlobal,numLocal,indexBase,platform): ";

    // get a internodal communicator from the Platform
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
      Teuchos::reduceAll<int,global_size_t>(*comm_,Teuchos::REDUCE_SUM,
        Teuchos::as<global_size_t>(numLocalElements_in),outArg(global_sum));
      /* there are three errors we should be detecting:
         - numGlobalElements != invalid() and it is incorrect/invalid
         - numLocalElements invalid (<0)
      */
      localChecks[0] = -1;
      localChecks[1] = 0;
      if (numLocalElements_in < L1 && numLocalElements_in != L0) {
        // invalid
        localChecks[0] = myImageID;
        localChecks[1] = 1;
      }
      else if (numGlobalElements_in < GST1 && numGlobalElements_in != GST0 && numGlobalElements_in != GSTI) {
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
      Teuchos::reduceAll<int,int>(*comm_,Teuchos::REDUCE_MAX,2,localChecks,globalChecks);
      if (globalChecks[0] != -1) {
        if (globalChecks[1] == 1) {
          TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
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
              errPrefix << "numGlobal doesn't match sum of numLocal (== " 
              << global_sum << ") on at least one node (possibly node " 
              << globalChecks[0] << ").");
        }
        else if (globalChecks[1] == 4) {
          TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
              errPrefix << "indexBase is not the same on all nodes (examine node " 
              << globalChecks[0] << ").");
        }
        else {
          // logic error on my part
          TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
              errPrefix << "logic error. Please contact the Tpetra team.");
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
    maxAllGID_ = indexBase_ + numGlobalElements_ - G1;
    minMyGID_ = start_index + indexBase_;
    maxMyGID_ = minMyGID_ + numLocalElements_ - G1;
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
    // The user has specified the distribution of elements over the
    // nodes, via entryList.  The distribution is not necessarily
    // contiguous or equally shared over the nodes.
    using Teuchos::as;
    using Teuchos::outArg;
    const size_t  L0 = Teuchos::OrdinalTraits<size_t>::zero();
    const global_size_t GST0 = Teuchos::OrdinalTraits<global_size_t>::zero();
    const global_size_t GST1 = Teuchos::OrdinalTraits<global_size_t>::one();
    const global_size_t GSTI = Teuchos::OrdinalTraits<global_size_t>::invalid();

    // The length of entryList on this node is the number of local
    // elements (on this node), even though entryList contains global
    // indices.  We assume that the number of local elements can be
    // stored in a LocalOrdinal.
    LocalOrdinal numLocalElements_in = Teuchos::as<LocalOrdinal>(entryList.size());

    std::string errPrefix = Teuchos::typeName(*this) + 
      "::constructor(numGlobal,entryList,indexBase,platform): ";

    const int myImageID = comm_->getRank();
    { // begin scoping block
      // for communicating failures 
      int localChecks[2], globalChecks[2];

      /* Compute the global number of elements.
	 We are doing this because exactly ONE of the following is true:
         - the user didn't specify it, and we need it
         - the user did specify it, but we need to 
           + validate it against the sum of the local sizes, and
           + ensure that it is the same on all nodes
       */
      global_size_t global_sum; // Global number of elements
      Teuchos::reduceAll<int,global_size_t>(*comm_, Teuchos::REDUCE_SUM,
        as<global_size_t>(numLocalElements_in), outArg(global_sum));
      localChecks[0] = -1;
      localChecks[1] = 0;
      // If the user supplied the number of global elements (i.e., if
      // it's not invalid (== GSTI)), then make sure that there is at
      // least one global element.  The first two clauses of the test
      // are apparently redundant, but help avoid compiler warnings
      // about comparing signed and unsigned integers.
      if (numGlobalElements_in < GST1 && numGlobalElements_in != GST0 && numGlobalElements_in != GSTI) {
        // Number of global elements is not the "invalid" value, and
        // is nonpositive.
        localChecks[0] = myImageID;
        localChecks[1] = 1;
      }
      else if (numGlobalElements_in != GSTI && numGlobalElements_in != global_sum) {
        // Number of global elements is not the "invalid" value, and
        // is not equal to the global sum.
        localChecks[0] = myImageID;
        localChecks[1] = 2;
      }
      //
      // Now check that all nodes have the same indexBase value.
      //
      GlobalOrdinal rootIB = indexBase_in;
      Teuchos::broadcast<int,GlobalOrdinal>(*comm_,0,&rootIB);   // broadcast one ordinal from node 0
      if (indexBase_in != rootIB) {
        localChecks[0] = myImageID;
        localChecks[1] = 3;
      }
      // REDUCE_MAX will give us the rank ("image ID") of the
      // highest-rank process that DID NOT pass.  This will be -1 if
      // all processes passed.
      Teuchos::reduceAll<int,int>(*comm_,Teuchos::REDUCE_MAX,2,localChecks,globalChecks);
      if (globalChecks[0] != -1) {
        if (globalChecks[1] == 1) {
          TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
              errPrefix << "numGlobal is not valid on at least one node (possibly node "
              << globalChecks[0] << ").");
        }
        else if (globalChecks[1] == 2) {
          TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
              errPrefix << "numGlobal doesn't match sum of numLocal (" 
              << global_sum << ") on at least one node (possibly node "
              << globalChecks[0] << ").");
        }
        else if (globalChecks[1] == 3) {
          TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
              errPrefix << "indexBase is not the same on all nodes (possibly node "
              << globalChecks[0] << ").");
        }
        else {
          // logic error on my part
          TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
              errPrefix << "logic error. Please contact the Tpetra team.");
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
    minMyGID_ = indexBase_;
    maxMyGID_ = indexBase_;
    // Create the GID to LID map.  Do not assume that the GIDs in
    // entryList are distinct.  In the case that a GID is duplicated,
    // use the same LID for all duplicates.  This is necessary so that
    // the LIDs are in [0,numLocal).
    size_t numUniqueGIDs = 0;
    if (numLocalElements_ > L0) {
      lgMap_ = Teuchos::arcp<GlobalOrdinal>(numLocalElements_);
      for (size_t i=0; i < numLocalElements_; i++) {
        lgMap_[numUniqueGIDs] = entryList[i];   // lgMap_:  LID to GID
        glMap_[entryList[i]] = numUniqueGIDs;   // glMap_: GID to LID
        numUniqueGIDs++;
      }

      // shrink lgMap appropriately
      if (numLocalElements_ != numUniqueGIDs) {
        numLocalElements_ = numUniqueGIDs;
        lgMap_ = lgMap_.persistingView(0,numLocalElements_);
      }
      minMyGID_ = *std::min_element(lgMap_.begin(), lgMap_.end());
      maxMyGID_ = *std::max_element(lgMap_.begin(), lgMap_.end());
    }

    // Compute the min and max of all processes' global IDs.
    Teuchos::reduceAll<int,GlobalOrdinal>(*comm_,Teuchos::REDUCE_MIN,minMyGID_,Teuchos::outArg(minAllGID_));
    Teuchos::reduceAll<int,GlobalOrdinal>(*comm_,Teuchos::REDUCE_MAX,maxMyGID_,Teuchos::outArg(maxAllGID_));
    contiguous_  = false;
    distributed_ = checkIsDist();
    TEUCHOS_TEST_FOR_EXCEPTION(minAllGID_ < indexBase_, std::invalid_argument,
      errPrefix << "Minimum global ID (== " << minAllGID_ << ") over all process"
      "(es) is less than the given indexBase (== " << indexBase_ << ")");
    setupDirectory();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Map<LocalOrdinal,GlobalOrdinal,Node>::~Map () 
  {}

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  LocalOrdinal Map<LocalOrdinal,GlobalOrdinal,Node>::getLocalElement(GlobalOrdinal globalIndex) const 
  {
    if (contiguous_) {
      if (globalIndex < getMinGlobalIndex() || globalIndex > getMaxGlobalIndex()) {
        return Teuchos::OrdinalTraits<LocalOrdinal>::invalid();
      }
      return Teuchos::as<LocalOrdinal>(globalIndex - getMinGlobalIndex());
    }
    else {
      typename std::map<GlobalOrdinal,LocalOrdinal>::const_iterator i;
      i = glMap_.find(globalIndex);
      if (i == glMap_.end()) {
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
    if (contiguous_) {
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
    if (contiguous_) {
      return (getMinGlobalIndex() <= globalIndex) && (globalIndex <= getMaxGlobalIndex());
    }
    else {
      typename std::map<GlobalOrdinal,LocalOrdinal>::const_iterator i;
      i = glMap_.find(globalIndex);
      return (i != glMap_.end());
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool Map<LocalOrdinal,GlobalOrdinal,Node>::isContiguous() const {
    return contiguous_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool Map<LocalOrdinal,GlobalOrdinal,Node>::isCompatible (const Map<LocalOrdinal,GlobalOrdinal,Node> &map) const {
    using Teuchos::outArg;
    // check to make sure distribution is the same
    char iscompat_lcl;
    if (getGlobalNumElements() != map.getGlobalNumElements() ||
          getNodeNumElements() != map.getNodeNumElements()) {
      // NOT compat on this node
      iscompat_lcl = 0;
    }
    else {
      // compat on this node
      iscompat_lcl = 1;
    }
    char iscompat_gbl;
    Teuchos::reduceAll<int,char>(*comm_,Teuchos::REDUCE_MIN,iscompat_lcl,
      outArg(iscompat_gbl));
    return (iscompat_gbl == 1);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool Map<LocalOrdinal,GlobalOrdinal,Node>::isSameAs(const Map<LocalOrdinal,GlobalOrdinal,Node> &map) const {
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
      if (contiguous_ == true && map.contiguous_ == true) {
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
    if (lgMap_ == Teuchos::null && numLocalElements_ > 0) {
#ifdef HAVE_TEUCHOS_DEBUG
      // The local-to-global mapping should have been set up already
      // for a noncontiguous map.
      TEUCHOS_TEST_FOR_EXCEPTION(contiguous_ != true, std::logic_error,
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
  map = Teuchos::rcp( new Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>(numElements,                                    // num elements, global and local
                                                              Teuchos::OrdinalTraits<GlobalOrdinal>::zero(),  // index base is zero
                                                              comm, GloballyDistributed, node) );
  return map.getConst();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
Tpetra::createLocalMapWithNode(size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node) {
  Teuchos::RCP< Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > map;
  map = Teuchos::rcp( new Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>((Tpetra::global_size_t)numElements,                     // num elements, global and local
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


#endif // TPETRA_MAP_DEF_HPP
