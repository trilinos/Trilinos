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
  void Directory<Ordinal>::getDirectoryEntries(
      std::vector<Ordinal> const& globalEntries, 
      std::vector<Ordinal>& images) const 
  {
    getEntries(globalEntries, images, images, false);
  }

  template<typename Ordinal>
  void Directory<Ordinal>::getDirectoryEntries(
      const std::vector<Ordinal> & globalEntries, 
      std::vector<Ordinal> & images, 
      std::vector<Ordinal> & localEntries) const 
  {
    getEntries(globalEntries, images, localEntries, true);
  }

  template<typename Ordinal>
  void Directory<Ordinal>::getEntries(
      const std::vector<Ordinal> & globalEntries, 
            std::vector<Ordinal> & images, 
            std::vector<Ordinal> & localEntries, 
            bool computeLIDs) const 
  {
    const Ordinal ZERO = Teuchos::OrdinalTraits<Ordinal>::zero();
    const Ordinal ONE  = Teuchos::OrdinalTraits<Ordinal>::one();
    const Ordinal NEGONE = ZERO - ONE;

    // allocate space in images and localEntries
    // resize to same length as globalEntries and fill with -1s.
    images.assign(globalEntries.size(), NEGONE);
    if(computeLIDs) {
      localEntries.assign(globalEntries.size(), NEGONE);
    }

    bool ierr = false;
    const Ordinal numImages  = comm_->getSize();
    const Ordinal myImageID  = comm_->getRank();
    const Ordinal numEntries = globalEntries.size();
    const Ordinal nOverP     = map_.getNumGlobalEntries() / numImages;

    if(!map_.isDistributed()) {
      // Easiest case: Map is serial or locally-replicated
      for(Ordinal i = ZERO; i < numEntries; ++i) {
        if(!map_.isMyGlobalIndex(globalEntries[i])) { 
          // This means something bad happened as there should be no non-local entries in a non-global ES
          ierr = true;                        
          break;
        }
        else {
          images[i] = myImageID;
          if(computeLIDs) {
            localEntries[i] = map_.getLocalIndex(globalEntries[i]);
          }
        }
      }
    }
    else if(map_.isContiguous()) {
      // Next Easiest Case: Map is distributed but allocated contiguously
      Ordinal minAllGID = map_.getMinAllGlobalIndex();
      Ordinal maxAllGID = map_.getMaxAllGlobalIndex();
      for(Ordinal i = ZERO; i < numEntries; i++) {
        Ordinal LID = NEGONE; // Assume not found
        Ordinal image = NEGONE;
        Ordinal GID = globalEntries[i];
        if (GID < minAllGID || GID > maxAllGID) {
          ierr = true;
          break;
        }
        else {
          // Guess uniform distribution and start a little above it
          Ordinal image1 = TPETRA_MIN((GID / TPETRA_MAX(nOverP, ONE)) + Teuchos::as<Ordinal>(2), numImages - ONE);
          bool found = false;
          while (image1 >= ZERO && image1 < numImages) {
            if(allMinGIDs_[image1] <= GID) {
              if(GID < allMinGIDs_[image1 + ONE]) {
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
          if(found) {
            image = image1;
            LID = GID - allMinGIDs_[image];
          }
        }
        images[i] = image;
        if(computeLIDs) {
          localEntries[i] = LID;
        }
      }
    }
    else {
      // General Case: Map is distributed and allocated arbitrarily
      // Here we need to set up an actual directory structure
      Ordinal packetSize = Teuchos::as<Ordinal>(2);
      if(computeLIDs) {
        packetSize++;
      }

      Distributor<Ordinal> distor(comm_);

      // Get directory locations for the requested list of entries
      std::vector<Ordinal> dirImages(numEntries);
      directoryMap_->getRemoteIndexList(globalEntries, dirImages);

      // Check for unfound globalEntries and set corresponding images to -1
      Ordinal numMissing = ZERO;
      for(Ordinal i = ZERO; i < numEntries; ++i) {
        if(dirImages[i] == NEGONE) {
          images[i] = NEGONE;
          if(computeLIDs) {
            localEntries[i] = NEGONE;
          }
          numMissing++;
        }
      }

      std::vector<Ordinal> sendGIDs;
      std::vector<Ordinal> sendImages;
      distor.createFromRecvs(globalEntries, dirImages, sendGIDs, sendImages);
      Ordinal numSends = Teuchos::as<Ordinal>(sendGIDs.size());

      Ordinal currLID;
      std::vector<Ordinal> exports;
      exports.reserve(packetSize * numSends);
      for(Ordinal i = ZERO; i < numSends; i++) {
        Ordinal currGID = sendGIDs[i];
        exports.push_back(currGID);
        currLID = directoryMap_->getLocalIndex(currGID);
        assert(currLID != NEGONE); // Internal error
        exports.push_back(imageIDs_[currLID]);
        if(computeLIDs) {
          exports.push_back(LIDs_[currLID]);
        }
      }

      std::vector<Ordinal> imports;
      distor.doPostsAndWaits(exports, packetSize, imports);

      typename std::vector<Ordinal>::iterator ptr = imports.begin();
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
    SHARED_TEST_FOR_EXCEPTION(ierr, std::invalid_argument, 
        "Tpetra::Directory<" << Teuchos::TypeNameTraits<Ordinal>::name() << 
        ">::getEntries(): Invalid GIDs given.",*comm_);
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

    // initialize Comm instance
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
    const std::vector<Ordinal> & myGlobalEntries = map_.getMyGlobalEntries();
    directoryMap_->getRemoteIndexList(myGlobalEntries, sendImageIDs);

    // Create distributor & call createFromSends
    Ordinal numReceives = ZERO;
    Distributor<Ordinal> distor(comm_);      
    distor.createFromSends(sendImageIDs, numReceives);

    // Execute distributor plan
    // Transfer GIDs, ImageIDs, and LIDs that we own to all images
    // End result is all images have list of all GIDs and corresponding ImageIDs and LIDs
    std::vector<Ordinal> exportEntries;
    Ordinal packetSize = Teuchos::as<Ordinal>(3); // We will send GIDs, ImageIDs, and LIDs.

    exportEntries.reserve(packetSize * numMyEntries);
    for(Ordinal i = ZERO; i < numMyEntries; i++) {
      exportEntries.push_back(myGlobalEntries[i]);
      exportEntries.push_back(myImageID);
      exportEntries.push_back(i);
    }

    std::vector<Ordinal> importElements;
    distor.doPostsAndWaits(exportEntries, packetSize, importElements);

    typename std::vector<Ordinal>::iterator ptr = importElements.begin();
    for(Ordinal i = ZERO; i < numReceives; i++) {
      Ordinal currLID = directoryMap_->getLocalIndex(*ptr++); // Convert incoming GID to Directory LID
      assert(currLID != NEGONE); // Internal error
      TEST_FOR_EXCEPTION(currLID == NEGONE, std::logic_error,
        "Tpetra::Directory<" << Teuchos::OrdinalTraits<Ordinal>::name() << ">::generateDirectory(): logic error. Please notify the Tpetra team.");
      imageIDs_[currLID] = *ptr++;
      LIDs_[currLID] = *ptr++;
    }
  }
    
} // namespace Tpetra

#endif // TPETRA_DIRECTORY_HPP
