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

#ifndef TPETRA_BASICDIRECTORY_HPP
#define TPETRA_BASICDIRECTORY_HPP

#include "../test/tpetra_test_util.hpp"
#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include "Tpetra_Object.hpp"
#include "Tpetra_Directory.hpp"
#include "Tpetra_ElementSpace.hpp"

namespace Tpetra {

  //! Tpetra::BasicDirectory - Default implementation of Tpetra::Directory
  
  /*! For ElementSpace objects, a Directory object must be created to allow referencing
      of non-local elements. BasicDirectory produces and contains a uniform linear
      ElementSpace and a list of imageIDs allowing non-local elements to be accessed
      by dereferencing throught the BasicDirectory.
      
      This class currently has one constructor, taking an ElementSpace object.
  */
  
  template<typename OrdinalType>
  class BasicDirectory : public Object, public virtual Directory<OrdinalType> {
  public:
    
    //@{ \name Constructors/Destructor.
    
    //! constructor
    BasicDirectory(ElementSpace<OrdinalType> const& ElementSpace)	
      : Object("Tpetra::BasicDirectory") 
      , ElementSpace_(ElementSpace) 
    {
      // A directory is not necessary for a non-global ES.
      if(ElementSpace.isGlobal()) {
        // If ElementSpace is contiguously allocated, we can construct the 
        // directory from the minMyGID value from each image.
        if(ElementSpace.isContiguous()) {
          OrdinalType const one = Teuchos::OrdinalTraits<OrdinalType>::one();
          allMinGIDs_.reserve(ElementSpace.platform().getNumImages() + one);
          OrdinalType minMyGID = ElementSpace.getMinMyGID();
          ElementSpace.comm().gatherAll(&minMyGID, &allMinGIDs_.front(), one);
          allMinGIDs_.back() = ElementSpace.getMaxAllGID() + one; // Set max cap
        }
        // Otherwise we have to generate the directory using MPI calls.
        else {
          generateDirectory();
        }
      }
    };
    
    //! copy constructor
    BasicDirectory(BasicDirectory<OrdinalType> const& Directory)
      : Object(Directory.label()) 
      , ElementSpace_(Directory.ElementSpace_) 
      , allMinGIDs_(Directory.allMinGIDs_)
      , imageIDs_(Directory.imageIDs_)
      , LIDs_(Directory.LIDs_)
    {};
    
    //! destructor.
    ~BasicDirectory() {};
    
    //@}
    
    //@{ \name Query methods.
    
    //! getDirectoryEntries : Returns image info for non-local ElementSpace entries
    /*! Given a list of Global Entry IDs, this function returns the list of
        IDs of the owning memory image that correspond to the list of entries.
      \param In
             globalEntries - List of Global IDs being passed in.
      \param Out
             images - On return contains list of Image IDs owning the Global IDs in question.
    */
    void getDirectoryEntries(std::vector<OrdinalType> const& globalEntries, std::vector<OrdinalType>& images) const {
      getEntries(globalEntries, images, images, false);
    };
    
    //! getDirectoryEntries : Returns image and local id info for non-local ElementSpace entries
    /*! Given a list of Global Entry IDs, this function returns the list of
        image IDs and local IDs on the owning memory image that correspond
        to the list of entries.  If LocalEntries is 0, then local IDs are 
        not returned.  If EntrySizes is nonzero, it will contain a list of corresponding 
        element sizes for the requested global entries.
      \param In
             globalEntries - List of Global IDs being passed in.
      \param Out
             images - On return contains list of Image IDs owning the Global IDs in question.
      \param Out
             localEntries - On return contains the local ID of the global on the owning image. 
    */
    void getDirectoryEntries(std::vector<OrdinalType> const& globalEntries, std::vector<OrdinalType>& images, std::vector<OrdinalType>& localEntries) const {
      getEntries(globalEntries, images, localEntries, true);
    };
    
    //@}
    
  private:
    ElementSpace<OrdinalType> const ElementSpace_;
    std::vector<OrdinalType> allMinGIDs_;
    std::vector<OrdinalType> imageIDs_;
    std::vector<OrdinalType> LIDs_;
    Teuchos::RefCountPtr< ElementSpace<OrdinalType> > directoryES_;
    
    //! Convenience function for accessing ElementSpace
    ElementSpace<OrdinalType> const& es() const {return(ElementSpace_);};

    //! Assignment operator (declared but not defined, do not use)
    BasicDirectory<OrdinalType>& operator = (BasicDirectory<OrdinalType> const& Source);

    // common code for both versions of getDirectoryEntries
    void getEntries(std::vector<OrdinalType> const& globalEntries, std::vector<OrdinalType>& images, std::vector<OrdinalType>& localEntries, bool computeLIDs) const {
      OrdinalType const zero = Teuchos::OrdinalTraits<OrdinalType>::zero();
      OrdinalType const one = Teuchos::OrdinalTraits<OrdinalType>::one();
      OrdinalType const negOne = zero - one;

      // allocate space in images and localEntries
      // resize to same length as globalEntries and fill with -1s.
      images.assign(globalEntries.size(), negOne);
      if(computeLIDs)
        localEntries.assign(globalEntries.size(), negOne);

      bool ierr = false;
      OrdinalType const myImageID = es().platform().getMyImageID();
      OrdinalType const numImages = es().platform().getNumImages();
      OrdinalType const numEntries = globalEntries.size();
      OrdinalType const nOverP = es().getNumGlobalElements() / numImages;

      // Easiest case: ElementSpace is serial or locally-replicated
      if(!es().isGlobal()) {
        for(OrdinalType i = zero; i < numEntries; i++) {
          if(!es().isMyGID(globalEntries[i])) { // This means something bad happened
            ierr = true;                        // As there should be no non-local elements in a non-global ES
          }
          else {
            images[i] = myImageID;
            if(computeLIDs)
              localEntries[i] = ElementSpace_.getLID(globalEntries[i]);
          }
        }
        if(ierr)
          throw reportError("Non-local GIDs given but this ElementSpace is not distributed globally", 1);
      }
      
      // Next Easiest Case: ElementSpace is distributed but allocated contiguously
      else if(es().isContiguous()) {
        OrdinalType minAllGID = es().getMinAllGID();
        OrdinalType maxAllGID = es().getMaxAllGID();
        for(OrdinalType i = zero; i < numEntries; i++) {
          OrdinalType LID = negOne; // Assume not found
          OrdinalType image = negOne;
          OrdinalType GID = globalEntries[i];
          if(GID < minAllGID || GID > maxAllGID) {
            cerr << "ERROR (Image " << myImageID << ") GID " << GID << " is outside the range of this ES (" << minAllGID << " - " << maxAllGID << ")" << endl;
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
          //throw reportError("Some GIDs specified were not found in this ElementSpace", 2);
          //cerr << "ERROR: BasicDirectory::getEntries - Some GIDs specified were not found in this ElementSpace" << endl;
      }
        
      // General Case: ElementSpace is distributed and allocated arbitrarily
      // Here we need to set up an actual directory structure
      else {
        OrdinalType packetSize = one + one; // We will send at least the GID and imageID. Might also send LID.
        if(computeLIDs)
          packetSize++;
        
        Teuchos::RefCountPtr< Distributor<OrdinalType> > distor = es().platform().createDistributor();

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
        distor->createFromRecvs(numEntries, globalEntries, dirImages, true, numSends, sendGIDs, sendImages);
        
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
        
        char* cImports = 0;
        OrdinalType* imports = 0;
        OrdinalType lenImports = zero;
        OrdinalType numRecv = numEntries - numMissing;
        distor->doPostsAndWaits(reinterpret_cast<char*>(&exports.front()), packetSize * sizeof(OrdinalType), lenImports, cImports);
        imports = reinterpret_cast<OrdinalType*>(cImports);

        OrdinalType* ptr = imports;
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
        delete[] cImports;
      }
    };
    
    // directory setup for non-contiguous ES
    void generateDirectory() {
      OrdinalType const one = Teuchos::OrdinalTraits<OrdinalType>::one();
      OrdinalType const zero = Teuchos::OrdinalTraits<OrdinalType>::zero();
      OrdinalType const negOne = zero - one;
      
      OrdinalType const minAllGID = es().getMinAllGID();
      OrdinalType const maxAllGID = es().getMaxAllGID();

      // DirectoryES will have a range of elements from the minimum to the maximum
      // GID of the user ES, and an indexBase of minAllGID from the user ES
      OrdinalType numGlobalElements = maxAllGID - minAllGID + one;
      
      // Create a uniform linear map to contain the directory
      directoryES_ = Teuchos::rcp(new ElementSpace<OrdinalType>(numGlobalElements, minAllGID, es().platform()));
      
      OrdinalType directory_numMyElements = directoryES_->getNumMyElements();

      // Allocate imageID List and LID List.  Initialize to -1s.
      // Initialize values to -1 in case the user global element list does
      // fill all IDs from minAllGID to maxAllGID (e.g., allows global indices to be 
      // all even integers).
      imageIDs_.resize(directory_numMyElements, negOne);
      LIDs_.resize(directory_numMyElements, negOne);
      

      // Get list of images owning the directory entries for the ElementSpace GIDs
      OrdinalType myImageID = es().platform().getMyImageID();
      OrdinalType numMyElements = es().getNumMyElements();
      std::vector<OrdinalType> sendImageIDs(numMyElements);
      std::vector<OrdinalType> myGlobalElements = es().getMyGlobalElements();
      directoryES_->getRemoteIDList(myGlobalElements, sendImageIDs);
      
      // Create distributor & call createFromSends
      OrdinalType numReceives = Teuchos::OrdinalTraits<OrdinalType>::zero();
      Teuchos::RefCountPtr< Distributor<OrdinalType> > distor = es().platform().createDistributor();      
      distor->createFromSends(es().getNumMyElements(), sendImageIDs, true, numReceives);

      // Execute distributor plan
      // Transfer GIDs, ImageIDs, and LIDs that we own to all images
      // End result is all images have list of all GIDs and corresponding ImageIDs and LIDs
      std::vector<OrdinalType> exportElements;
      char* cImportElements = 0;
      OrdinalType* importElements = 0;
      OrdinalType lenImportElements;

      OrdinalType packetSize = one + one + one; // We will send GIDs, ImageIDs, and LIDs.

      exportElements.reserve((one + one) * numMyElements);
      for(OrdinalType i = zero; i < numMyElements; i++) {
        exportElements.push_back(myGlobalElements[i]);
        exportElements.push_back(myImageID);
        exportElements.push_back(i);
      }

      distor->doPostsAndWaits(reinterpret_cast<char*>(&exportElements.front()), packetSize * sizeof(OrdinalType), lenImportElements, cImportElements);
      importElements = reinterpret_cast<OrdinalType*>(cImportElements);

      OrdinalType currLID;
      OrdinalType* ptr = importElements;
      for(OrdinalType i = zero; i < numReceives; i++) {
        currLID = directoryES_->getLID(*ptr++); // Convert incoming GID to Directory LID
        assert(currLID != negOne); // Internal error
        imageIDs_[currLID] = *ptr++;
        LIDs_[currLID] = *ptr++;
      }
        
      delete[] cImportElements;
    };
    
  }; // class MpiDirectory
  
} // namespace Tpetra

#endif // TPETRA_MPIDIRECTORY_HPP
