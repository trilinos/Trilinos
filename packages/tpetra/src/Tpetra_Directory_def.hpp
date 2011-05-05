// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
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
  : map_(map_in) {
    // initialize Comm instance
    comm_ = map_->getComm();

    // A directory is not necessary for a non-global ES.
    if (map_->isDistributed()) {
      // If map_ is contiguously allocated, we can construct the 
      // directory from the minMyGID value from each image.
      if (map_->isContiguous()) {
        // make room for the min on each proc, plus one entry at the end for the max cap
        allMinGIDs_.resize(comm_->getSize() + 1);
        // get my min
        GlobalOrdinal minMyGID = map_->getMinGlobalIndex();
        // gather all of the mins into the first getSize() entries of allMinDIGs_
        Teuchos::gatherAll<int,GlobalOrdinal>(*comm_,1,&minMyGID,comm_->getSize(),&allMinGIDs_.front());
        // put the max cap at the end
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
  LookupStatus Directory<LocalOrdinal,GlobalOrdinal,Node>::getDirectoryEntries(
              const Teuchos::ArrayView<const GlobalOrdinal> &globalIDs, 
              const Teuchos::ArrayView<int> &nodeIDs, 
              const Teuchos::ArrayView<LocalOrdinal> &localIDs) const {
    const bool computeLIDs = true;
    return getEntries(globalIDs, nodeIDs, localIDs, computeLIDs);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  LookupStatus Directory<LocalOrdinal,GlobalOrdinal,Node>::getEntries(
              const Teuchos::ArrayView<const GlobalOrdinal> &globalIDs, 
              const Teuchos::ArrayView<int> &nodeIDs, 
              const Teuchos::ArrayView<LocalOrdinal> &localIDs, 
              bool computeLIDs) const {
    const LocalOrdinal LINVALID = Teuchos::OrdinalTraits<LocalOrdinal>::invalid();

    LookupStatus res = AllIDsPresent;

    // fill nodeIDs and localIDs with -1s
    TEST_FOR_EXCEPTION(nodeIDs.size() != globalIDs.size(), std::runtime_error,
        Teuchos::typeName(*this) << "::getEntries(): Output buffers are not allocated properly.");
    std::fill(nodeIDs.begin(),nodeIDs.end(), -1);
    if (computeLIDs) {
      TEST_FOR_EXCEPTION(localIDs.size() != globalIDs.size(), std::runtime_error,
          Teuchos::typeName(*this) << "::getEntries(): Output buffers are not allocated properly.");
      std::fill(localIDs.begin(),localIDs.end(), LINVALID);
    }

    const int numImages  = comm_->getSize();
    const int myImageID  = comm_->getRank();
    const size_t numEntries = globalIDs.size();
    const global_size_t nOverP = map_->getGlobalNumElements() / numImages;

    if (map_->isDistributed() == false) {
      // Easiest case: Map is serial or locally-replicated
      typename Teuchos::ArrayView<int>::iterator imgptr = nodeIDs.begin();
      typename Teuchos::ArrayView<LocalOrdinal>::iterator lidptr = localIDs.begin();
      typename Teuchos::ArrayView<const GlobalOrdinal>::iterator gid;
      for (gid = globalIDs.begin(); gid != globalIDs.end(); ++gid) 
      {
        if (map_->isNodeGlobalElement(*gid)) {
          *imgptr++ = myImageID;
          if (computeLIDs) {
            *lidptr++ = map_->getLocalElement(*gid);
          }
        }
        else {
          // advance the pointers, leaving these values set to invalid
          imgptr++;
          if (computeLIDs) {
            lidptr++;
          }
          res = IDNotPresent;
        }
      }
    }
    else if (map_->isContiguous()) {
      // Next Easiest Case: Map is distributed but allocated contiguously
      typename Teuchos::ArrayView<int>::iterator imgptr = nodeIDs.begin();
      typename Teuchos::ArrayView<LocalOrdinal>::iterator lidptr = localIDs.begin();
      typename Teuchos::ArrayView<const GlobalOrdinal>::iterator gid;
      for (gid = globalIDs.begin(); gid != globalIDs.end(); ++gid) 
      {
        LocalOrdinal LID = LINVALID; // Assume not found
        int image = -1;
        GlobalOrdinal GID = *gid;
        // Guess uniform distribution and start a little above it
        // TODO: replace by a binary search
        int curimg = TEUCHOS_MIN((int)(GID / TEUCHOS_MAX(nOverP, 1)) + 2, numImages - 1);
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
          LID = Teuchos::as<LocalOrdinal>(GID - allMinGIDs_[image]);
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
      using Teuchos::as;
      int packetSize = 2;
      if (computeLIDs) {
        packetSize = 3;
      }

      Distributor distor(comm_);

      // Get directory locations for the requested list of entries
      Teuchos::Array<int> dirImages(numEntries);
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

      Teuchos::ArrayRCP<GlobalOrdinal> sendGIDs; 
      Teuchos::ArrayRCP<int> sendImages;
      distor.createFromRecvs(globalIDs, dirImages(), sendGIDs, sendImages);
      size_t numSends = sendGIDs.size();

      //    global_size_t >= GlobalOrdinal
      //    global_size_t >= size_t >= int
      //    global_size_t >= size_t >= LocalOrdinal
      // Therefore, we can safely stored all of these in a global_size_t
      Teuchos::Array<global_size_t> exports(packetSize*numSends);
      {
        LocalOrdinal curLID;
        typename Teuchos::Array<global_size_t>::iterator ptr = exports.begin();
        typename Teuchos::ArrayRCP<GlobalOrdinal>::const_iterator gidptr;
        for (gidptr = sendGIDs.begin(); gidptr != sendGIDs.end(); ++gidptr) {
          *ptr++ = as<global_size_t>(*gidptr);
          curLID = directoryMap_->getLocalElement(*gidptr);
          TEST_FOR_EXCEPTION(curLID == LINVALID, std::logic_error,
              Teuchos::typeName(*this) << "::getEntries(): Internal logic error. Please contact Tpetra team.");
          *ptr++ = as<global_size_t>(nodeIDs_[curLID]);
          if (computeLIDs) {
            *ptr++ = as<global_size_t>(LIDs_[curLID]);
          }
        }
      }

      Teuchos::Array<global_size_t> imports(packetSize*distor.getTotalReceiveLength());
      distor.doPostsAndWaits(exports().getConst(), packetSize, imports());

      typename Teuchos::Array<global_size_t>::iterator ptr = imports.begin();
      const size_t numRecv = numEntries - numMissing;
      // we know these conversions are in range, because we loaded this data
      for (size_t i = 0; i < numRecv; ++i) {
        GlobalOrdinal curGID = as<GlobalOrdinal>(*ptr++);
        for (size_t j = 0; j < numEntries; ++j) {
          if (curGID == globalIDs[j]) {
            nodeIDs[j] = as<int>(*ptr++);
            if (computeLIDs) {
              localIDs[j] = as<LocalOrdinal>(*ptr++);
            }
            if (nodeIDs[j] == -1) res = IDNotPresent;
            break;
          }
        }
      }
    }
    return res;
  }


  // directory setup for non-contiguous Map
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void Directory<LocalOrdinal,GlobalOrdinal,Node>::generateDirectory() {
    using Teuchos::as;
    const LocalOrdinal  LINVALID = Teuchos::OrdinalTraits<LocalOrdinal>::invalid();
          
    const GlobalOrdinal minAllGID = map_->getMinAllGlobalIndex();
    const GlobalOrdinal maxAllGID = map_->getMaxAllGlobalIndex();

    // DirectoryMap will have a range of elements from the minimum to the maximum
    // GID of the user Map, and an indexBase of minAllGID from the user Map
    const global_size_t numGlobalEntries = maxAllGID - minAllGID + 1;

    // Obviously, we can't afford to store the whole directory on each node
    // Create a uniform linear map to contain the directory to split up the storage among all nodes
    directoryMap_ = Teuchos::rcp(new Map<LocalOrdinal,GlobalOrdinal,Node>(numGlobalEntries, minAllGID, comm_, GloballyDistributed, map_->getNode() ));

    const size_t dir_numMyEntries = directoryMap_->getNodeNumElements();

    // Allocate imageID List and LID List.  Initialize to -1s.
    // Initialize values to -1 in case the user global element list does
    // fill all IDs from minAllGID to maxAllGID (e.g., allows global indices to be 
    // all even integers).
    nodeIDs_.resize(dir_numMyEntries, -1);
    LIDs_.resize(dir_numMyEntries, LINVALID);

    // Get list of nodeIDs owning the directory entries for the Map GIDs
    const int myImageID = comm_->getRank();
    const size_t numMyEntries = map_->getNodeNumElements();
    Teuchos::Array<int> sendImageIDs(numMyEntries);
    Teuchos::ArrayView<const GlobalOrdinal> myGlobalEntries = map_->getNodeElementList();
    // an ID not present in this lookup indicates that it lies outside of the range [minAllGID,maxAllGID] (from map_). 
    // this means something is wrong with map_, our fault.
    TEST_FOR_EXCEPTION( directoryMap_->getRemoteIndexList(myGlobalEntries, sendImageIDs) == IDNotPresent, std::logic_error,
        Teuchos::typeName(*this) << "::generateDirectory(): logic error. Please contact Tpetra team.");

    // Create distributor & call createFromSends
    size_t numReceives = 0;
    Distributor distor(comm_);      
    numReceives = distor.createFromSends(sendImageIDs);

    // Execute distributor plan
    // Transfer GIDs, ImageIDs, and LIDs that we own to all nodeIDs
    // End result is all nodeIDs have list of all GIDs and corresponding ImageIDs and LIDs
    int packetSize = 3; // We will send GIDs, ImageIDs, and LIDs.
    // GlobalOrdinal >= LocalOrdinal and int, so this is safe
    Teuchos::Array<GlobalOrdinal> exportEntries(packetSize*numMyEntries);
    {
      typename Teuchos::Array<GlobalOrdinal>::iterator ptr = exportEntries.begin();
      for (size_t i=0; i < numMyEntries; ++i) {
        *ptr++ = myGlobalEntries[i];
        *ptr++ = as<GlobalOrdinal>(myImageID);
        *ptr++ = as<GlobalOrdinal>(i);
      }
    }

    Teuchos::Array<GlobalOrdinal> importElements(packetSize*distor.getTotalReceiveLength());
    distor.doPostsAndWaits(exportEntries().getConst(), packetSize, importElements());

    {// begin scoping block
      typename Teuchos::Array<GlobalOrdinal>::iterator ptr = importElements.begin();
      for (size_t i = 0; i < numReceives; ++i) {
        LocalOrdinal currLID = directoryMap_->getLocalElement(*ptr++); // Convert incoming GID to Directory LID
        TEST_FOR_EXCEPTION(currLID == LINVALID, std::logic_error,
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
