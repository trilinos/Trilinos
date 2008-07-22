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
#include "Tpetra_DirectoryDecl.hpp"
#include "Tpetra_Distributor.hpp"

namespace Tpetra {

  template<typename OrdinalType>
  Directory<OrdinalType>::Directory(const Map<OrdinalType> & map)  
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
        const OrdinalType one = Teuchos::OrdinalTraits<OrdinalType>::one();
        // make room for the min on each proc, plus one element at the end for the max cap
        allMinGIDs_.resize(comm_->getSize() + one);
        // get my min
        OrdinalType minMyGID = map.getMinGlobalIndex();
        // gather all of the mins into the first getSize() elements of allMinDIGs_
        Teuchos::gatherAll(*comm_,one,&minMyGID,(OrdinalType)(map.getComm()->getSize()),&allMinGIDs_.front());
        // put the max cap at the end
        allMinGIDs_.back() = map.getMaxAllGlobalIndex() + one; // FINISH: is this right?
      }
      // Otherwise we have to generate the directory using MPI calls.
      else {
        generateDirectory();
      }
    }
  }

  template<typename OrdinalType>
  Directory<OrdinalType>::Directory(const Directory<OrdinalType> & directory)
    : Teuchos::Object(directory.label()) 
    , map_(directory.ElementSpace_) 
    , comm_(directory.comm_)
    , allMinGIDs_(directory.allMinGIDs_)
    , imageIDs_(directory.imageIDs_)
    , LIDs_(directory.LIDs_)
  {}

  template<typename OrdinalType>
  Directory<OrdinalType>::~Directory() {}
    
  template<typename OrdinalType>
  void Directory<OrdinalType>::getDirectoryEntries(
      std::vector<OrdinalType> const& globalEntries, 
      std::vector<OrdinalType>& images) const 
  {
    getEntries(globalEntries, images, images, false);
  }

  template<typename OrdinalType>
  void Directory<OrdinalType>::getDirectoryEntries(
      std::vector<OrdinalType> const& globalEntries, 
      std::vector<OrdinalType>& images, 
      std::vector<OrdinalType>& localEntries) const 
  {
    getEntries(globalEntries, images, localEntries, true);
  }

  template<typename OrdinalType>
  void Directory<OrdinalType>::getEntries(
      const std::vector<OrdinalType> & globalEntries, 
            std::vector<OrdinalType> & images, 
            std::vector<OrdinalType> & localEntries, 
             bool computeLIDs) const 
  {
    (void)globalEntries;
    (void)images;
    (void)localEntries;
    (void)computeLIDs;
  /*
    OrdinalType const zero = Teuchos::OrdinalTraits<OrdinalType>::zero();
    OrdinalType const one = Teuchos::OrdinalTraits<OrdinalType>::one();
    OrdinalType const negOne = zero - one;

    // allocate space in images and localEntries
    // resize to same length as globalEntries and fill with -1s.
    images.assign(globalEntries.size(), negOne);
    if(computeLIDs)
      localEntries.assign(globalEntries.size(), negOne);

    bool ierr = false;
    OrdinalType const myImageID = map_.getComm()->getRank();
    OrdinalType const numImages = map_.getComm()->getSize();
    OrdinalType const numEntries = globalEntries.size();
    OrdinalType const nOverP = map_.getNumGlobalElements() / numImages;

    // Easiest case: Map is serial or locally-replicated
    if(!map_.isGlobal()) {
      for(OrdinalType i = zero; i < numEntries; i++) {
        if(!map_.isMyGID(globalEntries[i])) { // This means something bad happened
          ierr = true;                        // As there should be no non-local elements in a non-global ES
        }
        else {
          images[i] = myImageID;
          if(computeLIDs)
            localEntries[i] = ElementSpace_.getLID(globalEntries[i]);
        }
      }
      if(ierr)
        throw reportError("Non-local GIDs given but this Map is not distributed globally", 1);
    }

    // Next Easiest Case: Map is distributed but allocated contiguously
    else if(map_.isContiguous()) {
      OrdinalType minAllGID = map_.getMinAllGlobalIndex();
      OrdinalType maxAllGID = map_.getMaxAllGlobalIndex();
      for(OrdinalType i = zero; i < numEntries; i++) {
        OrdinalType LID = negOne; // Assume not found
        OrdinalType image = negOne;
        OrdinalType GID = globalEntries[i];
        if(GID < minAllGID || GID > maxAllGID) {
          cerr << "ERROR (Image " << myImageID << ") GID " << GID 
            << " is outside the range of this ES (" << minAllGID 
            << " - " << maxAllGID << ")" << endl;
          ierr = true;
        }
        else {
          // Guess uniform distribution and start a little above it
          OrdinalType image1 = TPETRA_MIN((GID / TPETRA_MAX(nOverP, one)) + one + one, numImages - one);
          bool found = false;
          while(image1 >= zero && image1 < numImages) {
            if(allMinGIDs_[image1] <= GID) {
              if(GID < allMinGIDs_[image1 + one]) {
                found = true;
                break;
              }
              else
                image1++;
            }
            else
              image1--;
          }
          if(found) {
            image = image1;
            LID = GID - allMinGIDs_[image];
          }
        }
        images[i] = image;
        if(computeLIDs)
          localEntries[i] = LID;
      }
      //if(ierr)
      //throw reportError("Some GIDs specified were not found in this Map", 2);
      //cerr << "ERROR: BasicDirectory::getEntries - Some GIDs specified were not found in this Map" << endl;
    }

    // General Case: Map is distributed and allocated arbitrarily
    // Here we need to set up an actual directory structure
    else {
      OrdinalType packetSize = one + one; // We will send at least the GID and imageID. Might also send LID.
      if(computeLIDs)
        packetSize++;

      Distributor<OrdinalType> distor(comm_);

      // Get directory locations for the requested list of entries
      std::vector<OrdinalType> dirImages(numEntries);
      directoryES_->getRemoteIDList(globalEntries, dirImages);

      // Check for unfound globalEntries and set corresponding images to -1
      OrdinalType numMissing = zero;
      for(OrdinalType i = zero; i < numEntries; i++) {
        if(dirImages[i] == negOne) {
          images[i] = negOne;
          if(computeLIDs)
            localEntries[i] = negOne;
          numMissing++;
        }
      }

      OrdinalType numSends;
      std::vector<OrdinalType> sendGIDs;
      std::vector<OrdinalType> sendImages;
      distor.createFromRecvs(numEntries, globalEntries, dirImages, true, numSends, sendGIDs, sendImages);

      OrdinalType currLID;
      std::vector<OrdinalType> exports;
      exports.reserve(packetSize * numSends);
      for(OrdinalType i = zero; i < numSends; i++) {
        OrdinalType currGID = sendGIDs[i];
        exports.push_back(currGID);
        currLID = directoryES_->getLID(currGID);
        assert(currLID != negOne); // Internal error
        exports.push_back(imageIDs_[currLID]);
        if(computeLIDs)
          exports.push_back(LIDs_[currLID]);
      }

      std::vector<OrdinalType> imports;
      comm_->doPostsAndWaits(distor, exports, packetSize, imports);

      typename std::vector<OrdinalType>::iterator ptr = imports.begin();
      OrdinalType const numRecv = numEntries - numMissing;
      for(OrdinalType i = zero; i < numRecv; i++) {
        currLID = *ptr++;
        for(OrdinalType j = zero; j < numEntries; j++) {
          if(currLID == globalEntries[j]) {
            images[j] = *ptr++;
            if(computeLIDs)
              localEntries[j] = *ptr++;
            break;
          }
        }
      }
    }
  */
  }
    
    // directory setup for non-contiguous Map
  template<typename OrdinalType>
  void Directory<OrdinalType>::generateDirectory() 
  {
    const OrdinalType one = Teuchos::OrdinalTraits<OrdinalType>::one();
    const OrdinalType zero = Teuchos::OrdinalTraits<OrdinalType>::zero();
    const OrdinalType negOne = zero - one;
          
    const OrdinalType minAllGID = map_.getMinAllGlobalIndex();
    const OrdinalType maxAllGID = map_.getMaxAllGlobalIndex();

    // initialize Comm instance
    comm_ = map_.getComm();

    // DirectoryMap will have a range of elements from the minimum to the maximum
    // GID of the user Map, and an indexBase of minAllGID from the user Map
    OrdinalType numGlobalEntries = maxAllGID - minAllGID + one;

    // Create a uniform linear map to contain the directory
    directoryMap_ = Teuchos::rcp(new Map<OrdinalType>(numGlobalEntries, minAllGID, *map_.getPlatform()));

    OrdinalType dir_numMyEntries = directoryMap_->getNumMyEntries();

    // Allocate imageID List and LID List.  Initialize to -1s.
    // Initialize values to -1 in case the user global element list does
    // fill all IDs from minAllGID to maxAllGID (e.g., allows global indices to be 
    // all even integers).
    imageIDs_.resize(dir_numMyEntries, negOne);
    LIDs_.resize(dir_numMyEntries, negOne);


    // Get list of images owning the directory entries for the Map GIDs
    OrdinalType myImageID = comm_->getRank();
    OrdinalType numMyEntries = map_.getNumMyEntries();
    std::vector<OrdinalType> sendImageIDs(numMyEntries);
    const std::vector<OrdinalType> & myGlobalEntries = map_.getMyGlobalEntries();
    directoryMap_->getRemoteIndexList(myGlobalEntries, sendImageIDs);

    // Create distributor & call createFromSends
    OrdinalType numReceives = Teuchos::OrdinalTraits<OrdinalType>::zero();
    Distributor<OrdinalType> distor(comm_);      
    distor.createFromSends(sendImageIDs, numReceives);

    // Execute distributor plan
    // Transfer GIDs, ImageIDs, and LIDs that we own to all images
    // End result is all images have list of all GIDs and corresponding ImageIDs and LIDs
    std::vector<OrdinalType> exportEntries;
    OrdinalType packetSize = one + one + one; // We will send GIDs, ImageIDs, and LIDs.

    exportEntries.reserve(packetSize * numMyEntries);
    for(OrdinalType i = zero; i < numMyEntries; i++) {
      exportEntries.push_back(myGlobalEntries[i]);
      exportEntries.push_back(myImageID);
      exportEntries.push_back(i);
    }

    std::vector<OrdinalType> importElements;
    distor.doPostsAndWaits(exportEntries, packetSize, importElements);

    typename std::vector<OrdinalType>::iterator ptr = importElements.begin();
    for(OrdinalType i = zero; i < numReceives; i++) {
      OrdinalType currLID = directoryMap_->getLocalIndex(*ptr++); // Convert incoming GID to Directory LID
      assert(currLID != negOne); // Internal error
      TEST_FOR_EXCEPTION(currLID == negOne, std::logic_error,
        "Tpetra::Directory<" << Teuchos::OrdinalTraits<OrdinalType>::name() << ">::generateDirectory(): logic error. Contact Tpetra team.");
      imageIDs_[currLID] = *ptr++;
      LIDs_[currLID] = *ptr++;
    }
  }
    
} // namespace Tpetra

#endif // TPETRA_DIRECTORY_HPP
