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

#ifndef TPETRA_ELEMENTSPACE_HPP
#define TPETRA_ELEMENTSPACE_HPP

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include "Tpetra_ConfigDefs.hpp" // for STL map and vector
#include "Tpetra_Object.hpp"
#include "Tpetra_Platform.hpp"
#include "Tpetra_Directory.hpp"
#include "Tpetra_Util.hpp" // for toString

namespace Tpetra {

  // forward declarations
  template<typename PacketType, typename OrdinalType> class Comm;
  template<typename OrdinalType> class ElementSpaceData;
  
  //! Tpetra::ElementSpace: A class for constructing and using template<ordinalType> ElementSpaces.
  /*! ElementSpace objects are defined to have an element size of 1. Variable element sizes are implemented 
      in Tpetra::BlockElementSpace. Some ElementSpace methods throw exceptions, and should be enclosed 
      in a try/catch block. All Tpetra_ElementSpace objects require a Tpetra_Platform object. 
      Local IDs (LIDs) are always in the range [0, numMyElements).
      
      ElementSpace error codes (positive for non-fatal, negative for fatal):
      <ol>
      <li> +1  Specified Global ID not found on this image.
      <li> +2  Specified Local ID not found on this image.
      <li> -1  numGlobalElements < -1.  Should be >= -1 (Should be >= 0 for a Tpetra-defined ElementSpace).
      <li> -2  numMyElements < 0.  Should be >= 0.
      <li> -3  Invalid numGlobalElements.  Should equal sum of myGlobalElements, or set to -1 to compute automatically.
      <li> -4  Minimum global element index is less than index base.
      <li> -99 Internal ElementSpace error.  Contact developer.
      </ol>*/
  
  
  template<typename OrdinalType>
  class ElementSpace : public Object {
    
  public:
    
    //@{ \name Constructor/Destructor Methods
    
    //! Tpetra::ElementSpace constructor with Tpetra-defined contiguous uniform distribution.
    ElementSpace(OrdinalType numGlobalElements, OrdinalType indexBase, 
                 Platform<OrdinalType, OrdinalType> const& Platform)
      : Object("Tpetra::ElementSpace")
      , ElementSpaceData_()
    {
      const OrdinalType one = Teuchos::OrdinalTraits<OrdinalType>::one();
      const OrdinalType zero = Teuchos::OrdinalTraits<OrdinalType>::zero();
      
      // initial throws
      if (numGlobalElements < zero)
        throw reportError("numGlobalElements = " + toString(numGlobalElements) + ".  Should be >= " + toString(zero) + ".", -1);
      
      // platform & comm setup
      Teuchos::RefCountPtr< Comm<OrdinalType, OrdinalType> > comm = Platform.createOrdinalComm();
      OrdinalType numImages = Platform.getNumImages();
      OrdinalType myImageID = Platform.getMyImageID();
      
      // compute numMyElements
      OrdinalType numMyElements = numGlobalElements / numImages;
      OrdinalType remainder = numGlobalElements % numImages;
      OrdinalType start_index = myImageID * (numMyElements + one);
      if (myImageID < remainder)
        numMyElements++;
      else
        start_index -= (myImageID - remainder);
      
      // setup lgmap & glmap
      map<OrdinalType, OrdinalType> lgMap;
      map<OrdinalType, OrdinalType> glMap;
      
      // setup min/maxs
      OrdinalType minAllGID = indexBase;
      OrdinalType maxAllGID = minAllGID + numGlobalElements - one;
      OrdinalType minMyGID = start_index + indexBase;
      OrdinalType maxMyGID = minMyGID + numMyElements - one;
      
      // call ESData constructor
      ElementSpaceData_ = Teuchos::rcp(new ElementSpaceData<OrdinalType>(indexBase, numGlobalElements, numMyElements, minAllGID, maxAllGID, 
                                                                         minMyGID, maxMyGID, lgMap, glMap, true, Platform.clone(), comm));
      
      // initialize directory
      directorySetup();
    };
    
    //! Tpetra::ElementSpace constructor with user-defined contiguous distribution.
    ElementSpace(OrdinalType numGlobalElements, OrdinalType numMyElements, 
                 OrdinalType indexBase, Platform<OrdinalType, OrdinalType> const& Platform)
      : Object("Tpetra::ElementSpace")
      , ElementSpaceData_()
    {
      const OrdinalType one = Teuchos::OrdinalTraits<OrdinalType>::one();
      const OrdinalType zero = Teuchos::OrdinalTraits<OrdinalType>::zero();
      const OrdinalType negOne = zero - one;
      
      // initial throws
      if(numGlobalElements < negOne) 
        throw reportError("numGlobalElements = " + toString(numGlobalElements) + ".  Should be >= " + toString(negOne) + ".", -1);
      if(numMyElements < zero) 
        throw reportError("numMyElements = " + toString(numMyElements) + ".  Should be >= " + toString(zero) + ".", -2);
      
      // platform & comm setup
      Teuchos::RefCountPtr< Comm<OrdinalType, OrdinalType> > comm = Platform.createOrdinalComm();
      
      // check for invalid numGlobalElements
      //   Sum up all local element counts to get global count, and then
      //   check to see if user's value for numGlobalElements is either -1 
      //   (in which case we use our computed value) or matches ours.
      OrdinalType global_sum;
      comm->sumAll(&numMyElements, &global_sum, one);
      if(numGlobalElements == negOne)
        numGlobalElements = global_sum;
      else if(numGlobalElements != global_sum) 
        throw reportError("Invalid numGlobalElements.  numGlobalElements = " + toString(numGlobalElements) + 
                          ".  Should equal " + toString(global_sum) + ", or be set to -1 to compute automatically", -3);
      
      // setup lgmap & glmap
      map<OrdinalType, OrdinalType> lgMap;
      map<OrdinalType, OrdinalType> glMap;
      
      // setup min/maxs
      OrdinalType minAllGID = indexBase;
      OrdinalType maxAllGID = minAllGID + numGlobalElements - one;
      OrdinalType start_index;
      comm->scanSum(&numMyElements, &start_index, one);
      start_index -= numMyElements;
      OrdinalType minMyGID = start_index + indexBase;
      OrdinalType maxMyGID = minMyGID + numMyElements - one;
      
      // call ESData constructor
      ElementSpaceData_ = Teuchos::rcp(new ElementSpaceData<OrdinalType>(indexBase, numGlobalElements, numMyElements, minAllGID, maxAllGID, 
                                                                         minMyGID, maxMyGID, lgMap, glMap, true, Platform.clone(), comm));
      
      // initialize directory
      directorySetup();
    }
    
    //! Tpetra::ElementSpace constructor with user-defined non-contiguous (arbitrary) distribution.
    ElementSpace(OrdinalType numGlobalElements, OrdinalType numMyElements, std::vector<OrdinalType> const& elementList, 
                 OrdinalType indexBase, Platform<OrdinalType, OrdinalType> const& Platform)
      : Object("Tpetra::ElementSpace")
      , ElementSpaceData_()
    {
      const OrdinalType one = Teuchos::OrdinalTraits<OrdinalType>::one();
      const OrdinalType zero = Teuchos::OrdinalTraits<OrdinalType>::zero();
      const OrdinalType negOne = zero - one;
      
      // initial throws
      if(numGlobalElements < negOne) 
        throw reportError("numGlobalElements = " + toString(numGlobalElements) + ".  Should be >= " + toString(negOne) + ".", -1);
      if(numMyElements < zero) 
        throw reportError("numMyElements = " + toString(numMyElements) + ".  Should be >= " + toString(zero) + ".", -2);
      
      // platform & comm setup
      Teuchos::RefCountPtr< Comm<OrdinalType, OrdinalType> > comm = Platform.createOrdinalComm();
      
      // check for invalid numGlobalElements
      //   Sum up all local element counts to get global count, and then
      //   check to see if user's value for numGlobalElements is either -1 
      //   (in which case we use our computed value) or matches ours.
      OrdinalType global_sum;
      comm->sumAll(&numMyElements, &global_sum, one);
      if(numGlobalElements == negOne)
        numGlobalElements = global_sum;
      else if(numGlobalElements != global_sum)
        throw reportError("Invalid numGlobalElements.  numGlobalElements = " + toString(numGlobalElements) + 
                          ".  Should equal " + toString(global_sum) + ", or be set to -1 to compute automatically", -3);
      
      // setup lgmap and glmap, and min/maxMyGIDs
      map<OrdinalType, OrdinalType> lgMap;
      map<OrdinalType, OrdinalType> glMap;
      OrdinalType minMyGID = indexBase;
      OrdinalType maxMyGID = indexBase;
      if(numMyElements > zero) {
        for(OrdinalType i = zero; i < numMyElements; i++) {
          lgMap[i + zero] = elementList[i]; // lgmap: LID=key, GID=mapped
          glMap[elementList[i]] = (i + zero); // glmap: GID=key, LID=mapped
        }
        minMyGID = *min_element(elementList.begin(), elementList.end());
        maxMyGID = *max_element(elementList.begin(), elementList.end());
      }
      
      // set min/maxAllGIDs
      OrdinalType minAllGID;
      OrdinalType maxAllGID;
      comm->minAll(&minMyGID, &minAllGID, one);
      comm->maxAll(&maxMyGID, &maxAllGID, one);
      if (minAllGID < indexBase)
        throw reportError("Minimum global element index = " + toString(minAllGID) + 
                          " is less than index base = " + toString(indexBase) +".", -4);
      
      // call ESData constructor
      ElementSpaceData_ = Teuchos::rcp(new ElementSpaceData<OrdinalType>(indexBase, numGlobalElements, numMyElements, minAllGID, maxAllGID, 
                                                                         minMyGID, maxMyGID, lgMap, glMap, false, Platform.clone(), comm));
      
      // initialize directory
      directorySetup();
    };
    
    //! Tpetra::ElementSpace copy constructor.
    ElementSpace(ElementSpace<OrdinalType> const& ElementSpace)
      : Object(ElementSpace.label())
      , ElementSpaceData_(ElementSpace.ElementSpaceData_)
    {};
    
    //! Tpetra::ElementSpace destructor.
    ~ElementSpace() {};
    
    //@}
    
    
    //@{ \name Local/Global ID Accessor Methods
    
    //! Returns the image IDs and corresponding local IDs for a given list of global IDs.
    void getRemoteIDList(std::vector<OrdinalType> const& GIDList, std::vector<OrdinalType>& imageIDList, std::vector<OrdinalType>& LIDList) const {
      data().Directory_->getDirectoryEntries(GIDList, imageIDList, LIDList);
    };
    //! Returns only the image IDs for a given list of global IDs.
    void getRemoteIDList(std::vector<OrdinalType> const& GIDList, std::vector<OrdinalType>& imageIDList) const {
      data().Directory_->getDirectoryEntries(GIDList, imageIDList);
    };
    
    //! Returns local ID of global ID passed in, throws exception 1 if not found on this image.
    OrdinalType getLID(OrdinalType GID) const {
      if(!isMyGID(GID)) 
        throw reportError("Global ID " + toString(GID) + " was not found on this image.", 1);
      else if(isContiguous()) 
        return(GID - getMinMyGID()); //compute with offset
      else {
        return((data().glMap_.find(GID))->second);
      }
    };
    
    //! Returns global ID of local ID passed in, throws exception 2 if not found on this image.
    OrdinalType getGID(OrdinalType LID) const {
      if(!isMyLID(LID))
        throw reportError("Local ID " + toString(LID) + " was not found on this image.", 2);
      else if(isContiguous()) 
        return(LID + getMinMyGID()); //compute with offset
      else {
        return((data().lgMap_.find(LID))->second);
      }
    };
    
    //! Returns true if global ID passed in belongs to the calling image, returns false if it doesn't.
    bool isMyGID(OrdinalType GID) const {
      if(GID < getMinMyGID() || GID > getMaxMyGID())
        return(false);
      else if(isContiguous())
        return(true);
      else {
        return(data().glMap_.find(GID) != data().glMap_.end());
      }
    };
    
    //! Returns true if the local ID passed in belongs to the calling image, returns false if it doesn't.
    bool isMyLID(OrdinalType LID) const {
      if(LID < getMinLID() || LID > getMaxLID())
        return(false);
      else if(isContiguous())
        return(true);
      else {
        return (data().lgMap_.find(LID) != data().lgMap_.end());
      }
    };
    
    //! Returns the minimum global ID in this ElementSpace.
    OrdinalType getMinAllGID() const {return(data().minAllGID_);};
    
    //! Returns the maximum global ID in this ElementSpace.
    OrdinalType getMaxAllGID() const {return(data().maxAllGID_);};
    
    //! Returns the minimum global ID owned by this image.
    OrdinalType getMinMyGID() const {return(data().minMyGID_);};
    
    //! Returns the maximum global ID owned by this image.
    OrdinalType getMaxMyGID() const {return(data().maxMyGID_);};
    
    //! Returns the minimum local ID on the calling image.
    OrdinalType getMinLID() const {return(data().minLID_);};
    
    //! Returns the maximum local ID on the calling image.
    OrdinalType getMaxLID() const {return(data().maxLID_);};
    
    //@}
    
    
    //@{ \name Size & Dimension Accessor Methods
    
    //! Returns the number of elements in this ElementSpace.
    OrdinalType getNumGlobalElements() const {return(data().numGlobalElements_);};
    
    //! Returns the number of elements belonging to the calling image.
    OrdinalType getNumMyElements() const {return(data().numMyElements_);};
    
    //! Returns the Index base for this ElementSpace. Normally 0 for C/C++ or 1 for Fortran, but can be anything.
    OrdinalType getIndexBase() const {return(data().indexBase_);};
    
    //@}
    
    
    //@{ \name Misc. Boolean Tests
    
    //! Returns true if this ElementSpace is distributed contiguously, returns false otherwise.
    bool isContiguous() const {return(data().contiguous_);};
    
    //! Returns true if this ElementSpace is distributed across more than one image, returns false otherwise.
    bool isGlobal() const {return(data().global_);};
    
    //! Returns true if the ElementSpace passed in is identical to this ElementSpace. Also implemented through the == and != operators.
    bool isSameAs(ElementSpace<OrdinalType> const& ElementSpace) const {
      // Quickest test: See if we share an inner data class
      if(ElementSpaceData_.shares_resource(ElementSpace.ElementSpaceData_))
        return(true);
      
      // Next check other global properties that are easy global attributes
      if(getMinAllGID() != ElementSpace.getMinAllGID() || 
         getMaxAllGID() != ElementSpace.getMaxAllGID() ||
         getNumGlobalElements() != ElementSpace.getNumGlobalElements() || 
         getIndexBase() != ElementSpace.getIndexBase() ||
         isGlobal() != ElementSpace.isGlobal() || 
         isContiguous() != ElementSpace.isContiguous()) 
        return(false);
      
      // If we get this far, we need to check local properties and then check across
      // all images to see if local properties are all true
      
      // check that maps have the same number of local elements
      int mySameSpace = 1;
      if(getNumMyElements() != ElementSpace.getNumMyElements()) 
        mySameSpace = 0;
      
      // then check that GIDs are the same, by checking that the maps match
      // *possible optimization* is this only necessary in a non-contiguous elementspace?
      if(mySameSpace == 1)
        if(data().lgMap_ != ElementSpace.data().lgMap_)
          mySameSpace = 0;
      
      // Now get min of mySameSpace across all images
      int globalSameSpace = 0;
      comm().minAll(&mySameSpace, &globalSameSpace, 1);
      
      // if globalSameSpace is 1, that means none of the images set it to 0,
      // so the ElementSpaces are identical on all images. If this is the case, then we should return true.
      return(globalSameSpace == 1);
    };

    bool operator==(ElementSpace<OrdinalType> const& ElementSpace) const {return(isSameAs(ElementSpace));};
    bool operator!=(ElementSpace<OrdinalType> const& ElementSpace) const {return(!isSameAs(ElementSpace));};
    
    //@}
    
    
    //@{ \name Array Accessor Functions
    
    //! Returns a reference to an internal array containing a list of the Global IDs assigned to the calling image.
    // If myGlobalElements_ is empty, and we have any elements, then we create it.
    std::vector<OrdinalType> const& getMyGlobalElements() const {
      OrdinalType const nME = getNumMyElements();
      OrdinalType const zero = Teuchos::OrdinalTraits<OrdinalType>::zero();
      if((data().myGlobalElements_.empty()) && (nME > zero)) {
        data().myGlobalElements_.resize(nME);
        if(isContiguous()) {
          OrdinalType const minMyGID = getMinMyGID();
          for(OrdinalType i = zero; i < nME; i++)
            data().myGlobalElements_[i] = minMyGID + i;
        }
        else { // not contiguous
          typename std::map<OrdinalType, OrdinalType>::const_iterator lgi = data().lgMap_.begin();
          typename std::map<OrdinalType, OrdinalType>::const_iterator lgmax = data().lgMap_.end();
          for(OrdinalType i = zero; lgi != lgmax; i++) {
            data().myGlobalElements_[i] = lgi->second;
            lgi++;
          }
        }
      }
      return(data().myGlobalElements_);
    };
    
    //! Puts list of global elements on this image into the user-provided array.
    void getMyGlobalElements(std::vector<OrdinalType>& elementList) const {
      elementList = getMyGlobalElements();
    };
    //@}
    

    //@{ \name Misc.
    
    //! Prints the ElementSpace object to the output stream.
    /*! An << operator is inherited from Tpetra::Object, which uses the print method.*/
    void print(ostream& os) const {
      OrdinalType const zero = Teuchos::OrdinalTraits<OrdinalType>::zero();
      OrdinalType const nME = getNumMyElements();
      
      getMyGlobalElements(); // throw away output, we call this to make sure list is generated
      int myImageID = platform().getMyImageID();
      int numImages = platform().getNumImages();
      
      for(int imageCtr = 0; imageCtr < numImages; imageCtr++) {
        if(myImageID == imageCtr) {
          if(myImageID == 0) { // this is the root node (only output this info once)
            os <<  "\nNumber of Global Elements  = "; os << getNumGlobalElements(); os << endl;
            os <<    "Maximum of all GIDs        = "; os << getMaxAllGID(); os << endl;
            os <<    "Minimum of all GIDs        = "; os << getMinAllGID(); os << endl;
            os <<    "Index Base                 = "; os << getIndexBase(); os << endl;
          }
          os << endl;
          
          os <<    "Number of Local Elements   = "; os << nME; os << endl;
          os <<    "Maximum of my GIDs         = "; os << getMaxMyGID(); os << endl;
          os <<    "Minimum of my GIDs         = "; os << getMinMyGID(); os << endl;
          os << endl;
          
          os.width(14);
          os <<  "   ImageID"; os << "    ";
          os.width(14);
          os <<  "       Local Index "; os << " ";
          os.width(14);
          os <<  "      Global Index "; os << " ";
          os << endl;
          
          for(OrdinalType i = zero; i < nME; i++) {
            os.width(14);
            os <<  myImageID; os << "    ";
            os.width(14);
            os << i; os << "    ";
            os.width(14);
            os << data().myGlobalElements_[i]; os << "    ";
            os << endl;
          }
          
          os << flush;
          
        }
        // Do a few global ops to give I/O a chance to complete
        comm().barrier();
        comm().barrier();
        comm().barrier();
      }
    };
      
    //! Access functions for the Tpetra::Comm and Tpetra::Platform communicators.
    Comm<OrdinalType, OrdinalType> const& comm() const {return(*data().Comm_);};
    Platform<OrdinalType, OrdinalType> const& platform() const {return(*(data().Platform_));};
    
    //! Assignment operator
    ElementSpace<OrdinalType>& operator = (ElementSpace<OrdinalType> const& Source) {
      ElementSpaceData_ = Source.ElementSpaceData_;
      return(*this);
    }
    
    //@}
    
  private:
    // private data members
    Teuchos::RefCountPtr< ElementSpaceData<OrdinalType> > ElementSpaceData_; // Teuchos smart pointer
    
    // private functions
    void directorySetup() {
      if(getNumGlobalElements() != Teuchos::OrdinalTraits<OrdinalType>::zero())
        if(data().haveDirectory_ == false) {
          data().Directory_ = platform().createDirectory(*this); // Make directory
          data().haveDirectory_ = true;
        }
    };

  	// convenience functions for returning inner data class, both const and nonconst versions.
  	ElementSpaceData<OrdinalType>& data() {return(*ElementSpaceData_);};
  	ElementSpaceData<OrdinalType> const& data() const {return(*ElementSpaceData_);};
    
    
  }; // ElementSpace class
  
} // Tpetra namespace

#include "Tpetra_ElementSpaceData.hpp"

#endif // TPETRA_ELEMENTSPACE_HPP
