// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2004) Sandia Corporation
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

#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_as.hpp>
#include "Tpetra_DirectoryDecl.hpp"
#include "Tpetra_Distributor.hpp"
#include "Tpetra_Map.hpp"

namespace Tpetra {

  template<typename Ordinal>
  Directory<Ordinal>::Directory(const Map<Ordinal> & map)  
    : Teuchos::Object("Tpetra::Directory") 
    , map_(map) 
  {
    // initialize Comm instance
    comm_ = map_.getComm();

    // A directory is not necessary for a non-global ES.
    if(map.isDistributed()) {
      // If map is contiguously allocated, we can construct the 
      // directory from the minMyGID value from each image.
      if(map.isContiguous()) {
        const Ordinal one = Teuchos::OrdinalTraits<Ordinal>::one();
        // make room for the min on each proc, plus one entry at the end for the max cap
        allMinGIDs_.resize(comm_->getSize() + one);
        // get my min
        Ordinal minMyGID = map.getMinGlobalIndex();
        // gather all of the mins into the first getSize() entries of allMinDIGs_
        Teuchos::gatherAll(*comm_,one,&minMyGID,(Ordinal)(map.getComm()->getSize()),&allMinGIDs_.front());
        // put the max cap at the end
        allMinGIDs_.back() = map.getMaxAllGlobalIndex() + one; // FINISH: is this right?
      }
      // Otherwise we have to generate the directory using MPI calls.
      else {
        generateDirectory();
      }
    }
  }

  template<typename Ordinal>
  Directory<Ordinal>::Directory(const Directory<Ordinal> & directory)
    : Teuchos::Object(directory.label()) 
    , map_(directory.map_) 
    , comm_(directory.comm_)
    , allMinGIDs_(directory.allMinGIDs_)
    , imageIDs_(directory.imageIDs_)
    , LIDs_(directory.LIDs_)
  {}

  template<typename Ordinal>
  Directory<Ordinal>::~Directory() {}

  template<typename Ordinal>
  bool Directory<Ordinal>::getDirectoryEntries(
      const Teuchos::ArrayView<const Ordinal> &globalEntries, 
      const Teuchos::ArrayView<Ordinal> &images) const 
  {
    return getEntries(globalEntries, images, Teuchos::ArrayView<Ordinal>(Teuchos::null), false);
  }

  template<typename Ordinal>
  bool Directory<Ordinal>::getDirectoryEntries(
      const Teuchos::ArrayView<const Ordinal> &globalEntries, 
      const Teuchos::ArrayView<Ordinal> &images, 
      const Teuchos::ArrayView<Ordinal> &localEntries) const 
  {
    return getEntries(globalEntries, images, localEntries, true);
  }

  template<typename Ordinal>
  bool Directory<Ordinal>::getEntries(
      const Teuchos::ArrayView<const Ordinal> &globalEntries, 
      const Teuchos::ArrayView<Ordinal> &images, 
      const Teuchos::ArrayView<Ordinal> &localEntries, 
            bool computeLIDs) const 
  {
    const Ordinal ZERO = Teuchos::OrdinalTraits<Ordinal>::zero();
    const Ordinal ONE  = Teuchos::OrdinalTraits<Ordinal>::one();
    const Ordinal NEGONE = ZERO - ONE;

    bool invalidGIDs = false;

    // fill images and localEntries with -1s
    TEST_FOR_EXCEPTION(images.size() != globalEntries.size(), std::runtime_error,
        "Tpetra::Directory<" << Teuchos::TypeNameTraits<Ordinal>::name() << 
        ">::getEntries(): Output buffers are not allocated properly.");
    std::fill(images.begin(),images.end(),NEGONE);
    if (computeLIDs) {
      TEST_FOR_EXCEPTION(localEntries.size() != globalEntries.size(), std::runtime_error,
          "Tpetra::Directory<" << Teuchos::TypeNameTraits<Ordinal>::name() << 
          ">::getEntries(): Output buffers are not allocated properly.");
      std::fill(localEntries.begin(),localEntries.end(),NEGONE);
    }

    const Ordinal numImages  = comm_->getSize();
    const Ordinal myImageID  = comm_->getRank();
    const Ordinal numEntries = globalEntries.size();
    const Ordinal nOverP     = map_.getNumGlobalEntries() / numImages;

    if (map_.isDistributed() == false) {
      // Easiest case: Map is serial or locally-replicated
      typename Teuchos::ArrayView<Ordinal>::iterator imgptr = images.begin(),    
                                                     lidptr = localEntries.begin();
      for (typename Teuchos::ArrayView<const Ordinal>::iterator gid = globalEntries.begin(); gid != globalEntries.end(); ++gid) {
        if (map_.isMyGlobalIndex(*gid)) {
          *imgptr++ = myImageID;
          if (computeLIDs) {
            *lidptr++ = map_.getLocalIndex(*gid);
          }
        }
        else {
          imgptr++;
          if (computeLIDs) {
            lidptr++;
          }
          invalidGIDs = true;
        }
      }
    }
    else if (map_.isContiguous()) {
      // Next Easiest Case: Map is distributed but allocated contiguously
      typename Teuchos::ArrayView<Ordinal>::iterator imgptr = images.begin(),    
                                                     lidptr = localEntries.begin();
      for (typename Teuchos::ArrayView<const Ordinal>::iterator gid = globalEntries.begin(); gid != globalEntries.end(); ++gid) {
        Ordinal LID = NEGONE; // Assume not found
        Ordinal image = NEGONE;
        Ordinal GID = *gid;
        // Guess uniform distribution and start a little above it
        Ordinal image1 = TEUCHOS_MIN((GID / TEUCHOS_MAX(nOverP, ONE)) + Teuchos::as<Ordinal>(2), numImages - ONE);
        bool found = false;
        while (image1 >= ZERO && image1 < numImages) {
          if (allMinGIDs_[image1] <= GID) {
            if (GID < allMinGIDs_[image1 + ONE]) {
              found = true;
              break;
            }
            else {
              image1++;
            }
          }
          else {
            image1--;
          }
        }
        if (found) {
          image = image1;
          LID = GID - allMinGIDs_[image];
        }
        else {
          invalidGIDs = true;
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
      Ordinal packetSize = Teuchos::as<Ordinal>(2);
      if (computeLIDs) {
        ++packetSize;
      }

      Distributor<Ordinal> distor(comm_);

      // Get directory locations for the requested list of entries
      Teuchos::Array<Ordinal> dirImages(numEntries);
      invalidGIDs = directoryMap_->getRemoteIndexList(globalEntries, dirImages());
      // Check for unfound globalEntries and set corresponding images to -1
      Ordinal numMissing = ZERO;
      if (invalidGIDs) {
        for (Ordinal i = ZERO; i < numEntries; ++i) {
          if (dirImages[i] == NEGONE) {
            images[i] = NEGONE;
            if (computeLIDs) {
              localEntries[i] = NEGONE;
            }
            numMissing++;
          }
        }
      }

      Teuchos::ArrayRCP<Ordinal> sendGIDs, sendImages;
      distor.createFromRecvs(globalEntries, dirImages(), sendGIDs, sendImages);
      Ordinal numSends = Teuchos::as<Ordinal>(sendGIDs.size());

      Ordinal currLID;
      Teuchos::Array<Ordinal> exports(packetSize*numSends);
      {
        typename Teuchos::Array<Ordinal>::iterator ptr = exports.begin();
        for(Ordinal i = ZERO; i < numSends; i++) {
          Ordinal currGID = sendGIDs[i];
          *ptr++ = currGID;
          currLID = directoryMap_->getLocalIndex(currGID);
          assert(currLID != NEGONE); // Internal error
          *ptr++ = imageIDs_[currLID];
          if(computeLIDs) {
            *ptr++ = LIDs_[currLID];
          }
        }
      }

      Teuchos::Array<Ordinal> imports(packetSize*distor.getTotalReceiveLength());
      distor.doPostsAndWaits(exports().getConst(), packetSize, imports());

      typename Teuchos::Array<Ordinal>::iterator ptr = imports.begin();
      const Ordinal numRecv = numEntries - numMissing;
      for(Ordinal i = ZERO; i < numRecv; i++) {
        currLID = *ptr++;
        for(Ordinal j = ZERO; j < numEntries; j++) {
          if(currLID == globalEntries[j]) {
            images[j] = *ptr++;
            if(computeLIDs) {
              localEntries[j] = *ptr++;
            }
            break;
          }
        }
      }
    }
    return invalidGIDs;
  }


  // directory setup for non-contiguous Map
  template<typename Ordinal>
  void Directory<Ordinal>::generateDirectory() 
  {
    const Ordinal ONE = Teuchos::OrdinalTraits<Ordinal>::one();
    const Ordinal ZERO = Teuchos::OrdinalTraits<Ordinal>::zero();
    const Ordinal NEGONE = ZERO - ONE;
          
    const Ordinal minAllGID = map_.getMinAllGlobalIndex();
    const Ordinal maxAllGID = map_.getMaxAllGlobalIndex();

    comm_ = map_.getComm();

    // DirectoryMap will have a range of elements from the minimum to the maximum
    // GID of the user Map, and an indexBase of minAllGID from the user Map
    Ordinal numGlobalEntries = maxAllGID - minAllGID + ONE;

    // Create a uniform linear map to contain the directory
    directoryMap_ = Teuchos::rcp(new Map<Ordinal>(numGlobalEntries, minAllGID, *map_.getPlatform()));

    Ordinal dir_numMyEntries = directoryMap_->getNumMyEntries();

    // Allocate imageID List and LID List.  Initialize to -1s.
    // Initialize values to -1 in case the user global element list does
    // fill all IDs from minAllGID to maxAllGID (e.g., allows global indices to be 
    // all even integers).
    imageIDs_.resize(dir_numMyEntries, NEGONE);
    LIDs_.resize(dir_numMyEntries, NEGONE);

    // Get list of images owning the directory entries for the Map GIDs
    Ordinal myImageID = comm_->getRank();
    Ordinal numMyEntries = map_.getNumMyEntries();
    std::vector<Ordinal> sendImageIDs(numMyEntries);
    Teuchos::ArrayView<const Ordinal> myGlobalEntries = map_.getMyGlobalEntries();
    // a "true" return here indicates that one of myGlobalEntries (from map_) is not on the map directoryMap_, indicating that 
    // it lies outside of the range [minAllGID,maxAllGID] (from map_). this means something is wrong with map_.
    TEST_FOR_EXCEPTION( directoryMap_->getRemoteIndexList(myGlobalEntries, sendImageIDs) == true, std::logic_error,
        "Tpetra::Directory::generateDirectory(): logic error. Please contact Tpetra team.");

    // Create distributor & call createFromSends
    Ordinal numReceives = ZERO;
    Distributor<Ordinal> distor(comm_);      
    distor.createFromSends(sendImageIDs, numReceives);

    // Execute distributor plan
    // Transfer GIDs, ImageIDs, and LIDs that we own to all images
    // End result is all images have list of all GIDs and corresponding ImageIDs and LIDs
    Ordinal packetSize = Teuchos::as<Ordinal>(3); // We will send GIDs, ImageIDs, and LIDs.
    Teuchos::Array<Ordinal> exportEntries(packetSize*numMyEntries);
    {
      typename Teuchos::Array<Ordinal>::iterator ptr = exportEntries.begin();
      for(Ordinal i = ZERO; i < numMyEntries; ++i) {
        *ptr++ = myGlobalEntries[i];
        *ptr++ = myImageID;
        *ptr++ = i;
      }
    }

    Teuchos::Array<Ordinal> importElements(packetSize*distor.getTotalReceiveLength());
    distor.doPostsAndWaits(exportEntries().getConst(), packetSize, importElements());

    {
      typename Teuchos::Array<Ordinal>::iterator ptr = importElements.begin();
      for(Ordinal i = ZERO; i < numReceives; i++) {
        Ordinal currLID = directoryMap_->getLocalIndex(*ptr++); // Convert incoming GID to Directory LID
        assert(currLID != NEGONE); // Internal error
        TEST_FOR_EXCEPTION(currLID == NEGONE, std::logic_error,
            "Tpetra::Directory<" << Teuchos::OrdinalTraits<Ordinal>::name() << ">::generateDirectory(): logic error. Please notify the Tpetra team.");
        imageIDs_[currLID] = *ptr++;
        LIDs_[currLID] = *ptr++;
      }
    }
  }
    
} // namespace Tpetra

#endif // TPETRA_DIRECTORY_HPP
