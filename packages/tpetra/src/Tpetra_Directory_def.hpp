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

#ifndef TPETRA_DIRECTORY_HPP
#define TPETRA_DIRECTORY_HPP

#include <Teuchos_as.hpp>
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_Distributor.hpp"
#include "Tpetra_Map.hpp"

#ifdef DOXYGEN_USE_ONLY
  #include "Tpetra_Directory_decl.hpp"
#endif

namespace Tpetra {

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Directory<LocalOrdinal,GlobalOrdinal,Node>::Directory(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map_in)
  : map_(map_in) 
  {
    // initialize Comm instance
    comm_ = map_->getComm();

    // A directory is not necessary for a non-global ES.
    if (map_->isDistributed()) {
      // If map_ is contiguously allocated, we can construct the 
      // directory from the minMyGID value from each image.
      if (map_->isContiguous()) {
        // Make room for the min GID on each proc, plus one entry at
        // the end for the max cap.
        allMinGIDs_.resize(comm_->getSize() + 1);
        // Get my process' min GID.
        GlobalOrdinal minMyGID = map_->getMinGlobalIndex();
        // Gather all of the min GIDs into the first getSize() entries
        // of allMinGIDs_.
        Teuchos::gatherAll<int,GlobalOrdinal>(*comm_,1,&minMyGID,comm_->getSize(),&allMinGIDs_.front());
        // Put the max cap at the end.  Adding one lets us write loops
        // over GIDs with the traditional strict less-than bound.
        allMinGIDs_.back() = map_->getMaxAllGlobalIndex() + Teuchos::OrdinalTraits<GlobalOrdinal>::one();
      }
      // Otherwise we have to generate the directory using MPI calls.
      else {
        generateDirectory();
      }
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Directory<LocalOrdinal,GlobalOrdinal,Node>::~Directory() {}

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  LookupStatus Directory<LocalOrdinal,GlobalOrdinal,Node>::getDirectoryEntries(
              const Teuchos::ArrayView<const GlobalOrdinal> &globalIDs, 
              const Teuchos::ArrayView<int> &nodeIDs) const {
    const bool computeLIDs = false;
    return getEntries(globalIDs, nodeIDs, Teuchos::null, computeLIDs);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  LookupStatus 
  Directory<LocalOrdinal,GlobalOrdinal,Node>::
  getDirectoryEntries (const Teuchos::ArrayView<const GlobalOrdinal> &globalIDs, 
		       const Teuchos::ArrayView<int> &nodeIDs, 
		       const Teuchos::ArrayView<LocalOrdinal> &localIDs) const 
  {
    const bool computeLIDs = true;
    return getEntries (globalIDs, nodeIDs, localIDs, computeLIDs);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  LookupStatus 
  Directory<LocalOrdinal,GlobalOrdinal,Node>::
  getEntries (const Teuchos::ArrayView<const GlobalOrdinal> &globalIDs, 
              const Teuchos::ArrayView<int> &nodeIDs, 
              const Teuchos::ArrayView<LocalOrdinal> &localIDs, 
              bool computeLIDs) const 
  {
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::as;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;

    const LO LINVALID = Teuchos::OrdinalTraits<LO>::invalid();
    LookupStatus res = AllIDsPresent;

    // Ensure that globalIDs, nodeIDs, and localIDs (if applicable)
    // all have the same size, before modifying any output arguments.
    TEUCHOS_TEST_FOR_EXCEPTION(nodeIDs.size() != globalIDs.size(), 
      std::invalid_argument, Teuchos::typeName(*this) << "::getEntries(): Output" 
      " buffers are not allocated properly.  nodeIDs.size() = " << nodeIDs.size()
      << " != globalIDs.size() = " << globalIDs.size() << ".");
    if (computeLIDs) {
      TEUCHOS_TEST_FOR_EXCEPTION(localIDs.size() != globalIDs.size(), 
        std::invalid_argument, Teuchos::typeName(*this) << "::getEntries(): "
        "Output buffers are not allocated properly.  localIDs.size() = " 
        << localIDs.size() << " != globalIDs.size() = " << globalIDs.size() 
        << ".");
    }

    // Initially, fill nodeIDs and localIDs (if applicable) with
    // invalid values.  The "invalid" process ID is -1; the invalid
    // LO comes from OrdinalTraits.
    std::fill (nodeIDs.begin(), nodeIDs.end(), -1);
    if (computeLIDs) {
      std::fill (localIDs.begin(), localIDs.end(), LINVALID);
    }

    const int numImages  = comm_->getSize();
    const int myImageID  = comm_->getRank();
    const size_t numEntries = globalIDs.size();
    const global_size_t nOverP = map_->getGlobalNumElements() / numImages;

    if (! map_->isDistributed()) {
      // Easiest case: Map is serial or locally replicated. 
      typename ArrayView<int>::iterator procIter = nodeIDs.begin();
      typename ArrayView<LO>::iterator lidIter = localIDs.begin();
      typename ArrayView<const GO>::iterator gidIter;
      for (gidIter = globalIDs.begin(); gidIter != globalIDs.end(); ++gidIter) {
        if (map_->isNodeGlobalElement (*gidIter)) {
          *procIter++ = myImageID;
          if (computeLIDs) {
            *lidIter++ = map_->getLocalElement (*gidIter);
          }
        }
        else {
          // Advance the pointers, leaving these values set to invalid
          procIter++;
          if (computeLIDs) {
            lidIter++;
          }
          res = IDNotPresent;
        }
      }
    }
    else if (map_->isContiguous()) {
      // Next easiest case: Map is distributed but allocated contiguously.
      typename ArrayView<int>::iterator imgptr = nodeIDs.begin();
      typename ArrayView<LO>::iterator lidptr = localIDs.begin();
      typename ArrayView<const GO>::iterator gid;
      for (gid = globalIDs.begin(); gid != globalIDs.end(); ++gid) {
        LO LID = LINVALID; // Assume not found until proven otherwise
        int image = -1;
        GO GID = *gid;
        // Guess uniform distribution and start a little above it
        // TODO: replace by a binary search
	int curimg;
	{ // We go through all this trouble to avoid overflow and
	  // signed / unsigned casting mistakes (that were made in
	  // previous versions of this code).
	  const GO one = Teuchos::OrdinalTraits<GO>::one();
	  const GO two = one + one;
	  const GO nOverP_GID = static_cast<GO> (nOverP);
	  const GO lowerBound = GID / std::max(nOverP_GID, one) + two;
	  // It's probably not OK to cast this to int in general.  It
	  // works as long as |GID| <= the global number of entries
	  // and nOverP is appropriately sized for int.  Trouble may
	  // ensue if the index base has an exotic value.
	  const int lowerBound_int = static_cast<int> (lowerBound);
	  curimg = std::min(lowerBound_int, numImages - 1);
	}
        bool found = false;
        while (curimg >= 0 && curimg < numImages) {
          if (allMinGIDs_[curimg] <= GID) {
            if (GID < allMinGIDs_[curimg + 1]) {
              found = true;
              break;
            }
            else {
              curimg++;
            }
          }
          else {
            curimg--;
          }
        }
        if (found) {
          image = curimg;
          LID = Teuchos::as<LO>(GID - allMinGIDs_[image]);
        }
        else {
          res = IDNotPresent;
        }
        *imgptr++ = image;
        if (computeLIDs) {
          *lidptr++ = LID;
        }
      }
    }
    else {
      // General Case: Map is distributed and allocated arbitrarily
      // Here we need to set up an actual directory structure
      int packetSize = 2;
      if (computeLIDs) {
        packetSize = 3;
      }

      Distributor distor(comm_);

      // Get directory locations for the requested list of entries
      Array<int> dirImages(numEntries);
      res = directoryMap_->getRemoteIndexList(globalIDs, dirImages());
      // Check for unfound globalIDs and set corresponding nodeIDs to -1
      size_t numMissing = 0;
      if (res == IDNotPresent) {
        for (size_t i=0; i < numEntries; ++i) {
          if (dirImages[i] == -1) {
            nodeIDs[i] = -1;
            if (computeLIDs) {
              localIDs[i] = LINVALID;
            }
            numMissing++;
          }
        }
      }

      ArrayRCP<GO> sendGIDs; 
      ArrayRCP<int> sendImages;
      distor.createFromRecvs(globalIDs, dirImages(), sendGIDs, sendImages);
      size_t numSends = sendGIDs.size();

      //    global_size_t >= GO
      //    global_size_t >= size_t >= int
      //    global_size_t >= size_t >= LO
      // Therefore, we can safely stored all of these in a global_size_t
      Array<global_size_t> exports(packetSize*numSends);
      {
        LO curLID;
        typename Array<global_size_t>::iterator ptr = exports.begin();
        typename ArrayRCP<GO>::const_iterator gidptr;
        for (gidptr = sendGIDs.begin(); gidptr != sendGIDs.end(); ++gidptr) {
          *ptr++ = as<global_size_t>(*gidptr);
          curLID = directoryMap_->getLocalElement(*gidptr);
          TEUCHOS_TEST_FOR_EXCEPTION(curLID == LINVALID, std::logic_error,
              Teuchos::typeName(*this) << "::getEntries(): Internal logic error. Please contact Tpetra team.");
          *ptr++ = as<global_size_t>(nodeIDs_[curLID]);
          if (computeLIDs) {
            *ptr++ = as<global_size_t>(LIDs_[curLID]);
          }
        }
      }

      Array<global_size_t> imports(packetSize*distor.getTotalReceiveLength());
      distor.doPostsAndWaits(exports().getConst(), packetSize, imports());

      typename Array<global_size_t>::iterator ptr = imports.begin();
      const size_t numRecv = numEntries - numMissing;

      Array<GO> sortedIDs(globalIDs);
      ArrayRCP<GO> offset = arcp<GO>(numEntries);
      GO ii=0;
      for (typename ArrayRCP<GO>::iterator oo = offset.begin(); oo != offset.end(); ++oo,++ii)
        *oo = ii;
      sort2(sortedIDs.begin(),sortedIDs.begin()+numEntries,offset.begin());

      typedef typename Array<GO>::iterator IT;
      // we know these conversions are in range, because we loaded this data
      for (size_t i = 0; i < numRecv; ++i) {
        GO curGID = as<GO> (*ptr++);
        std::pair<IT, IT> p1 = std::equal_range (sortedIDs.begin(), sortedIDs.end(), curGID);
        if (p1.first != p1.second) {
          //found it
          size_t j = p1.first - sortedIDs.begin();
          nodeIDs[offset[j]] = as<int>(*ptr++);
          if (computeLIDs) {
            localIDs[offset[j]] = as<LO>(*ptr++);
          }
          if (nodeIDs[offset[j]] == -1) res = IDNotPresent;
        }
      }
    }
    return res;
  }


  // directory setup for non-contiguous Map
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void 
  Directory<LocalOrdinal,GlobalOrdinal,Node>::generateDirectory() 
  {
    using Teuchos::as;
    using Teuchos::rcp;
    using Teuchos::Array;
    using Teuchos::ArrayView;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef Map<LO, GO, Node> map_type;

    const LO  LINVALID = Teuchos::OrdinalTraits<LO>::invalid();
    const GO minAllGID = map_->getMinAllGlobalIndex();
    const GO maxAllGID = map_->getMaxAllGlobalIndex();

    // DirectoryMap will have a range of elements from the minimum to the maximum
    // GID of the user Map, and an indexBase of minAllGID from the user Map
    const global_size_t numGlobalEntries = maxAllGID - minAllGID + 1;

    // We can't afford to store the whole directory on each node, so
    // create a uniform contiguous Map that describes how we will
    // distribute the directory over processes.
    directoryMap_ = rcp (new map_type (numGlobalEntries, minAllGID, comm_, GloballyDistributed, map_->getNode()));

    // The number of Directory elements that my process owns.
    const size_t dir_numMyEntries = directoryMap_->getNodeNumElements();

    // NOTE (mfh 21 Mar 2012): I rephrased the comment below so that
    // it made a little more sense, but I don't fully understand what
    // it means yet.
    //
    // Allocate process ID List and LID list.  Initialize to invalid
    // values, in case the user global element list does fill all IDs
    // from minAllGID to maxAllGID (e.g., allows global indices to be
    // all even integers).
    nodeIDs_.resize(dir_numMyEntries, -1);
    LIDs_.resize(dir_numMyEntries, LINVALID);

    // Get list of process IDs that own the directory entries for the
    // Map GIDs.  These will be the targets of the sends that the
    // Distributor will do.
    const int myImageID = comm_->getRank();
    const size_t numMyEntries = map_->getNodeNumElements();
    Array<int> sendImageIDs(numMyEntries);
    ArrayView<const GO> myGlobalEntries = map_->getNodeElementList();
    // An ID not present in this lookup indicates that it lies outside
    // of the range [minAllGID,maxAllGID] (from map_).  this means
    // something is wrong with map_, our fault.
    TEUCHOS_TEST_FOR_EXCEPTION( directoryMap_->getRemoteIndexList(myGlobalEntries, sendImageIDs) == IDNotPresent, 
      std::logic_error, Teuchos::typeName(*this) << "::generateDirectory(): the "
      "Directory Map could not find out where one or more of my Map's elements "
      "should go.  This probably means there is a bug in Map.  Please report "
      "this bug to the Tpetra developers.");

    // Initialize the distributor using the list of process IDs to
    // which to send.  We'll use the distributor to send out triples
    // of (GID, process ID, LID).  We're sending the entries to the
    // processes that the Directory Map says should own them, which is
    // why we called directoryMap_->getRemoteIndexList() above.
    Distributor distor (comm_);      
    const size_t numReceives = distor.createFromSends (sendImageIDs);

    // NOTE (mfh 21 Mar 2012) The following code assumes that 
    // sizeof(GO) >= sizeof(int) and sizeof(GO) >= sizeof(LO).
    //
    // Create and fill buffer of (GID, process ID, LID) triples to
    // send out.  We pack the (GID, process ID, LID) triples into a
    // single Array of GO, casting the process ID from int to GO and
    // the LID from LO to GO as we do so.
    int packetSize = 3; // We're sending triples, so packet size is 3.
    Array<GO> exportEntries (packetSize * numMyEntries); // data to send out
    {
      typename Array<GO>::iterator ptr = exportEntries.begin();
      for (size_t i=0; i < numMyEntries; ++i) {
        *ptr++ = myGlobalEntries[i];
        *ptr++ = as<GO>(myImageID);
        *ptr++ = as<GO>(i);
      }
    }
    // Buffer of data to receive.  The Distributor figured out for us
    // how many packets we're receiving, when we called its
    // createFromSends() method to set up the distribution plan.
    Array<GO> importElements (packetSize * distor.getTotalReceiveLength()); 

    // Distribute the triples of (GID, process ID, LID).
    distor.doPostsAndWaits (exportEntries().getConst(), packetSize, importElements());

    // Unpack the redistributed data.
    {
      typename Array<GO>::iterator ptr = importElements.begin();
      for (size_t i = 0; i < numReceives; ++i) {
	// Each "packet" (contiguous chunk of importElements) contains
	// a triple: (GID, process ID, LID).
	//
	// Convert incoming GID to Directory LID.
        const LO currLID = directoryMap_->getLocalElement (*ptr++); 
        TEUCHOS_TEST_FOR_EXCEPTION(currLID == LINVALID, std::logic_error,
            Teuchos::typeName(*this) << "::generateDirectory(): logic error. Please notify the Tpetra team.");
        nodeIDs_[currLID] = *ptr++;
        LIDs_[currLID] = *ptr++;
      }
    }
  } // end generateDirectory()
    
} // namespace Tpetra

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_DIRECTORY_INSTANT(LO,GO,NODE) \
  \
  template class Directory< LO , GO , NODE >; \

#endif // TPETRA_DIRECTORY_HPP
