#ifndef _TPETRA_ELEMENTSPACE_HPP_
#define _TPETRA_ELEMENTSPACE_HPP_

#include "Tpetra_Object.hpp"
#include "Tpetra_Directory.hpp"
#include "Tpetra_ConfigDefs.hpp" // for STL map
#include <Teuchos_RefCountPtr.hpp>
#include "Tpetra_ElementSpaceData.hpp"
#include "Tpetra_Platform.hpp"

namespace Tpetra {

// forward declarations
template<typename PacketType, typename OrdinalType> class Comm;

//! Tpetra::ElementSpace: A class for constructing and using template<ordinalType> ElementSpaces.
/*! ElementSpace objects are defined to have an element size of 1. Variable element sizes are implemented 
	in Tpetra::BlockElementSpace. Some ElementSpace methods throw exceptions, and should be enclosed 
	in a try/catch block. All Tpetra_ElementSpace objects require a Tpetra_Platform object. 
	Local IDs (LIDs) are always in the range 0 to numMyElements - 1.

	ElementSpace error codes (positive for non-fatal, negative for fatal):
  <ol>
  <li> +1  Specified Global ID not found on this image.
  <li> +2  Specified Local ID not found on this image.
  <li> +3  Pointer passed to getMyGlobalElements(ordinalType*) does not have child allocated. (Null pointer).
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
						 Platform<OrdinalType, OrdinalType> const& Platform);

//! Tpetra::ElementSpace constructor with user-defined contiguous distribution.
ElementSpace(OrdinalType numGlobalElements, OrdinalType numMyElements, 
						 OrdinalType indexBase, Platform<OrdinalType, OrdinalType> const& Platform);

//! Tpetra::ElementSpace constructor with user-defined non-contiguous (arbitrary) distribution.
ElementSpace(OrdinalType numGlobalElements, OrdinalType numMyElements, OrdinalType* elementList, 
						 OrdinalType indexBase, Platform<OrdinalType, OrdinalType> const& Platform);

//! Tpetra::ElementSpace copy constructor.
ElementSpace(ElementSpace<OrdinalType> const& ElementSpace);

//! Tpetra::ElementSpace destructor.
~ElementSpace();
  
//@}


//@{ \name Local/Global ID Accessor Methods

//! Returns the image IDs and corresponding local IDs for a given list of global IDs.
void getRemoteIDList(OrdinalType numIDs, OrdinalType* GIDList, OrdinalType* imageIDList, OrdinalType* LIDList) const 
	{ElementSpaceData_->Directory_->getDirectoryEntries(numIDs, GIDList, imageIDList, LIDList);};

//! Returns local ID of global ID passed in, throws exception -1 if not found on this image.
OrdinalType getLID(OrdinalType GID) const;

//! Returns global ID of local ID passed in, throws exception -1 if not found on this image.
OrdinalType getGID(OrdinalType LID) const;

//! Returns true if global ID passed in belongs to the calling image, returns false if it doesn't.
bool isMyGID(OrdinalType GID) const;

//! Returns true if the local ID passed in belongs to the calling image, returns false if it doesn't.
bool isMyLID(OrdinalType LID) const;

//! Returns the minimum global ID in this ElementSpace.
OrdinalType getMinAllGID() const {return(ElementSpaceData_->minAllGID_);};

//! Returns the maximum global ID in this ElementSpace.
OrdinalType getMaxAllGID() const {return(ElementSpaceData_->maxAllGID_);};

//! Returns the minimum global ID owned by this image.
OrdinalType getMinMyGID() const {return(ElementSpaceData_->minMyGID_);};

//! Returns the maximum global ID owned by this image.
OrdinalType getMaxMyGID() const {return(ElementSpaceData_->maxMyGID_);};

//! Returns the minimum local ID on the calling image.
OrdinalType getMinLID() const {return(ElementSpaceData_->minLID_);};

//! Returns the maximum local ID on the calling image.
OrdinalType getMaxLID() const {return(ElementSpaceData_->maxLID_);};

//@}


//@{ \name Size & Dimension Accessor Methods

//! Returns the number of elements in this ElementSpace.
OrdinalType getNumGlobalElements() const {return(ElementSpaceData_->numGlobalElements_);};

//! Returns the number of elements belonging to the calling image.
OrdinalType getNumMyElements() const {return(ElementSpaceData_->numMyElements_);};

//! Puts list of global elements on this processor into the user-provided array.
void getMyGlobalElements(OrdinalType* elementList) const;

//! Returns the Index base for this ElementSpace. Normally 0 for C/C++ or 1 for Fortran, but can be anything.
OrdinalType getIndexBase() const {return(ElementSpaceData_->indexBase_);};

//@}


//@{ \name Misc. Boolean Tests

//! Returns true if this ElementSpace is distributed contiguously, returns false otherwise.
bool isContiguous() const {return(ElementSpaceData_->contiguous_);};

//! Returns true if this ElementSpace is distributed across more than one image, returns false otherwise.
bool isGlobal() const {return(ElementSpaceData_->global_);};

//! Returns true if the ElementSpace passed in is identical to this ElementSpace. Also implemented through the == and != operators.
bool isSameAs(ElementSpace<OrdinalType> const& ElementSpace) const;
bool operator==(ElementSpace<OrdinalType> const& ElementSpace) const {return(isSameAs(ElementSpace));};
bool operator!=(ElementSpace<OrdinalType> const& ElementSpace) const {return(!isSameAs(ElementSpace));};

//@}


//@{ \name Array Accessor Functions

//! Returns a pointer to an internal array containing a list of the Global IDs assigned to the calling image.
OrdinalType* getMyGlobalElements() const;

//@}


//@{ \name Misc.

//! Prints the ElementSpace object to the output stream.
/*! An << operator is inherited from Tpetra::Object, which uses the print method.*/
void print(ostream& os) const;

//! Access function for the Tpetra::Comm and Tpetra::Platform communicators.
Comm<OrdinalType, OrdinalType> const& comm() const {return(*ElementSpaceData_->Comm_);};
Platform<OrdinalType, OrdinalType> const& platform() const {return(*ElementSpaceData_->Platform_);};

//@}

private:
// private data members
Teuchos::RefCountPtr< ElementSpaceData<OrdinalType> > ElementSpaceData_; // Teuchos smart pointer

// private functions
void directorySetup();
OrdinalType getZero() const {return(ElementSpaceData_->zero_);};

}; // ElementSpace class

// begin Tpetra_ElementSpace.cpp
//=============================================================================

// constructor #1, tpetra contig
//=======================================================================
template<typename OrdinalType>
ElementSpace<OrdinalType>::ElementSpace(OrdinalType numGlobalElements, OrdinalType indexBase, 
																				Platform<OrdinalType, OrdinalType> const& Platform)
  : Object("Tpetra::ElementSpace")
	, ElementSpaceData_()
{
	// initial throws
	if (numGlobalElements < 0)
    throw reportError("numGlobalElements = " + toString(numGlobalElements) + ".  Should be >= 0.", -1);

	// platform & comm setup
	Comm<OrdinalType, OrdinalType>* comm = Platform.createOrdinalComm();
  OrdinalType numImages = comm->getNumImages();
  OrdinalType myImageID = comm->getMyImageID();
  
	// compute numMyElements
  OrdinalType numMyElements = numGlobalElements / numImages;
  OrdinalType remainder = numGlobalElements % numImages;
  OrdinalType start_index = myImageID * (numMyElements + 1);
  if (myImageID < remainder)
    numMyElements++;
  else
    start_index -= (myImageID - remainder);

	// setup lgmap & glmap
	map<OrdinalType, OrdinalType> lgMap;
	map<OrdinalType, OrdinalType> glMap;
  
	// setup min/maxs
  OrdinalType minAllGID = indexBase;
  OrdinalType maxAllGID = minAllGID + numGlobalElements - 1;
  OrdinalType minMyGID = start_index + indexBase;
  OrdinalType maxMyGID = minMyGID + numMyElements - 1;
	
	// call ESData constructor
	ElementSpaceData_ = Teuchos::rcp(new ElementSpaceData<OrdinalType>(indexBase, numGlobalElements, numMyElements, minAllGID, maxAllGID, 
																																		 minMyGID, maxMyGID, lgMap, glMap, true, Platform, comm));
  
	// initialize directory
  directorySetup();
}

// constructor #2, user contig
//=======================================================================
template<typename OrdinalType>
ElementSpace<OrdinalType>::ElementSpace(OrdinalType numGlobalElements, OrdinalType numMyElements, OrdinalType indexBase, 
																				Platform<OrdinalType, OrdinalType> const& Platform)
	: Object("Tpetra::ElementSpace")
	, ElementSpaceData_()
{
	// initial throws
  if(numGlobalElements < -1) 
    throw reportError("numGlobalElements = " + toString(numGlobalElements) + ".  Should be >= -1.", -1);
  if(numMyElements < 0) 
    throw reportError("numMyElements = " + toString(numMyElements) + ".  Should be >= 0.", -2);

	// platform & comm setup
	Comm<OrdinalType, OrdinalType>* comm = Platform.createOrdinalComm();
  OrdinalType numImages = comm->getNumImages();
  OrdinalType myImageID = comm->getMyImageID();

	// check for invalid numGlobalElements
	//   Sum up all local element counts to get global count, and then
	//   check to see if user's value for numGlobalElements is either -1 
	//   (in which case we use our computed value) or matches ours.
  OrdinalType global_sum;
  comm->sumAll(&numMyElements, &global_sum, 1);
	if(numGlobalElements == -1)
		numGlobalElements = global_sum;
	else if(numGlobalElements != global_sum) 
		throw reportError("Invalid numGlobalElements.  numGlobalElements = " + toString(numGlobalElements) + 
											".  Should equal " + toString(global_sum) + ", or be set to -1 to compute automatically", -3);

	// setup lgmap & glmap
	map<OrdinalType, OrdinalType> lgMap;
	map<OrdinalType, OrdinalType> glMap;
	
	// setup min/maxs
  OrdinalType minAllGID = indexBase;
  OrdinalType maxAllGID = minAllGID + numGlobalElements - 1;
	OrdinalType start_index;
	comm->scanSum(&numMyElements, &start_index, 1);
	start_index -= numMyElements;
	OrdinalType minMyGID = start_index + indexBase;
	OrdinalType maxMyGID = minMyGID + numMyElements - 1;

	// call ESData constructor
	ElementSpaceData_ = Teuchos::rcp(new ElementSpaceData<OrdinalType>(indexBase, numGlobalElements, numMyElements, minAllGID, maxAllGID, 
																																		minMyGID, maxMyGID, lgMap, glMap, true, Platform, comm));
  
	// initialize directory
  directorySetup();
}

// constructor #3, user non-contig
//=======================================================================
template<typename OrdinalType>
ElementSpace<OrdinalType>::ElementSpace(OrdinalType numGlobalElements, OrdinalType numMyElements, OrdinalType* elementList, 
																				OrdinalType indexBase, Platform<OrdinalType, OrdinalType> const& Platform)
  : Object("Tpetra::ElementSpace")
	, ElementSpaceData_()
{
	// initial throws
  if(numGlobalElements < -1) 
    throw reportError("numGlobalElements = " + toString(numGlobalElements) + ".  Should be >= -1.", -1);
  if(numMyElements < 0) 
    throw reportError("numMyElements = " + toString(numMyElements) + ".  Should be >= 0.", -2);

	// platform & comm setup
	Comm<OrdinalType, OrdinalType>* comm = Platform.createOrdinalComm();
  OrdinalType numImages = comm->getNumImages();
  OrdinalType myImageID = comm->getMyImageID();

	// check for invalid numGlobalElements
  //   Sum up all local element counts to get global count, and then
	//   check to see if user's value for numGlobalElements is either -1 
	//   (in which case we use our computed value) or matches ours.
  OrdinalType global_sum;
  comm->sumAll(&numMyElements, &global_sum, 1);
	if(numGlobalElements == -1)
		numGlobalElements = global_sum;
	else if(numGlobalElements != global_sum)
		throw reportError("Invalid numGlobalElements.  numGlobalElements = " + toString(numGlobalElements) + 
											".  Should equal " + toString(global_sum) + ", or be set to -1 to compute automatically", -3);
  
	// setup lgmap and glmap, and min/maxMyGIDs
	map<OrdinalType, OrdinalType> lgMap;
	map<OrdinalType, OrdinalType> glMap;
	OrdinalType minMyGID = indexBase;
	OrdinalType maxMyGID = indexBase;
	const OrdinalType zero = Teuchos::OrdinalTraits<OrdinalType>::zero();
	if(numMyElements > 0) {
		for(OrdinalType i = 0; i < numMyElements; i++) {
			lgMap[i + zero] = elementList[i]; // lgmap: LID=key, GID=mapped
			glMap[elementList[i]] = (i + zero); // glmap: GID=key, LID=mapped
		}
    minMyGID = elementList[0];
    maxMyGID = elementList[numMyElements - 1];
	}

	// set min/maxAllGIDs
	OrdinalType minAllGID;
	OrdinalType maxAllGID;
	comm->minAll(&minMyGID, &minAllGID, 1);
	comm->maxAll(&maxMyGID, &maxAllGID, 1);
  if (minAllGID < indexBase)
    throw reportError("Minimum global element index = " + toString(minAllGID) + 
											" is less than index base = " + toString(indexBase) +".", -4);
	
	// call ESData constructor
	ElementSpaceData_ = Teuchos::rcp(new ElementSpaceData<OrdinalType>(indexBase, numGlobalElements, numMyElements, minAllGID, maxAllGID, 
																																		 minMyGID, maxMyGID, lgMap, glMap, false, Platform, comm));

	// initialize directory
  directorySetup();
}

// copy constructor
//=======================================================================
template<typename OrdinalType>
ElementSpace<OrdinalType>::ElementSpace (ElementSpace<OrdinalType> const& ElementSpace) 
  : Object(ElementSpace.label())
	, ElementSpaceData_(ElementSpace.ElementSpaceData_)
{}

//=======================================================================
template<typename OrdinalType>
ElementSpace<OrdinalType>::~ElementSpace() {}

//=======================================================================
template<typename OrdinalType>
OrdinalType ElementSpace<OrdinalType>::getLID (OrdinalType GID) const {
  if(!isMyGID(GID)) 
    throw reportError("Global ID " + toString(GID) + " was not found on this processor.", 1);
  else if(isContiguous()) 
    return(GID - getMinMyGID() + getZero()); //compute with offset
  else {
    return((ElementSpaceData_->glMap_.find(GID))->second);
	}
}

//=======================================================================
template<typename OrdinalType>
OrdinalType ElementSpace<OrdinalType>::getGID (OrdinalType LID) const {
  if(!isMyLID(LID))
    throw reportError("Local ID " + toString(LID) + " was not found on this processor.", 2);
  else if(isContiguous()) 
    return(LID + getMinMyGID()); //compute with offset
  else {
    return((ElementSpaceData_->lgMap_.find(LID))->second);
	}
}

//=======================================================================
template<typename OrdinalType>
bool ElementSpace<OrdinalType>::isMyGID (OrdinalType GID) const {
  if(GID < getMinMyGID() || GID > getMaxMyGID())
    return(false);
  else if(isContiguous())
    return(true);
  else {
    return (ElementSpaceData_->glMap_.find(GID) != ElementSpaceData_->glMap_.end());
	}
}

//=======================================================================
template<typename OrdinalType>
bool ElementSpace<OrdinalType>::isMyLID (OrdinalType LID) const {
  if(LID < getMinLID() || LID > getMaxLID())
    return(false);
  else if(isContiguous())
    return(true);
  else {
    return (ElementSpaceData_->lgMap_.find(LID) != ElementSpaceData_->lgMap_.end());
	}
}

//=======================================================================
template<typename OrdinalType>
void ElementSpace<OrdinalType>::getMyGlobalElements(OrdinalType* elementList) const {
  if(elementList == 0)
    throw reportError("Pointer does not have child allocated.", 3);
  else if(isContiguous()) {
		OrdinalType nME = getNumMyElements();
		OrdinalType minMyGID = getMinMyGID();
    for(OrdinalType i = 0; i < nME; i++)
      elementList[i] = minMyGID + i;
	}
  else { // not contiguous
    typename std::map<OrdinalType, OrdinalType>::iterator lgi = ElementSpaceData_->lgMap_.begin();
    typename std::map<OrdinalType, OrdinalType>::iterator lgmax = ElementSpaceData_->lgMap_.end();
    for(OrdinalType i = 0; lgi != lgmax; i++) {
      elementList[i] = lgi->second;
      lgi++;
    }
  }
}

//=======================================================================
template<typename OrdinalType>
OrdinalType* ElementSpace<OrdinalType>::getMyGlobalElements() const {
	OrdinalType nME = getNumMyElements();
  if((ElementSpaceData_->myGlobalElements_ == 0) && (nME > 0)) {
		ElementSpaceData_->myGlobalElements_ = new OrdinalType[nME];
    getMyGlobalElements(ElementSpaceData_->myGlobalElements_);
  }
  return(ElementSpaceData_->myGlobalElements_);
}

//=======================================================================
template<typename OrdinalType>
bool ElementSpace<OrdinalType>::isSameAs (ElementSpace<OrdinalType> const& ElementSpace) const {
  if (this == &ElementSpace) 
    return(true);
  if (getMinAllGID() != ElementSpace.getMinAllGID() || 
      getMaxAllGID() != ElementSpace.getMaxAllGID() ||
      getNumGlobalElements() != ElementSpace.getNumGlobalElements() || 
      getIndexBase() != ElementSpace.getIndexBase() ||
      isGlobal() != ElementSpace.isGlobal() || 
      isContiguous() != ElementSpace.isContiguous()) 
    return(false);

  // If we get this far, we need to check local properties and then check across
  // all processors to see if local properties are all true
	
  int mySameSpace = 1;
  if(getNumMyElements() != ElementSpace.getNumMyElements()) 
    mySameSpace=0;
	
  if(!isContiguous() && mySameSpace == 1)
    if(ElementSpaceData_->lgMap_ != ElementSpace.ElementSpaceData_->lgMap_)
      mySameSpace=0;

  // Now get min of mySameSpace across all processors
  int globalSameSpace = 0;
  comm().minAll(&mySameSpace, &globalSameSpace, 1);
  return(globalSameSpace==1);
}

//=======================================================================
template<typename OrdinalType>
void ElementSpace<OrdinalType>::print(ostream& os) const {
  OrdinalType* myGlobalElements1 = getMyGlobalElements();
  OrdinalType myImageID = comm().getMyImageID();
  OrdinalType numImages = comm().getNumImages();
	OrdinalType minLID = getMinLID();
	OrdinalType nME = getNumMyElements();
  
  for (int imageCtr = 0; imageCtr < numImages; imageCtr++) {
    if (myImageID == imageCtr) {
      if (myImageID == 0) {
				os <<  "\nNumber of Global Elements  = "; os << getNumGlobalElements(); os << endl;
				os <<    "Maximum of all GIDs        = "; os << getMaxAllGID(); os << endl;
				os <<    "Minimum of all GIDs        = "; os << getMinAllGID(); os << endl;
				os <<    "Index Base                 = "; os << getIndexBase(); os << endl;
      }
      os << endl;

      os <<    "Number of Local Elements   = "; os << getNumMyElements(); os << endl;
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
    
      for (OrdinalType i = 0, lid = minLID; i < nME; i++, lid++) {
				os.width(14);
				os <<  myImageID; os << "    ";
				os.width(14);
				os << lid; os << "    ";
				os.width(14);
				os <<  myGlobalElements1[i]; os << "    ";
				os << endl;
      }
      
      os << flush;
      
    }
    // Do a few global ops to give I/O a chance to complete
    comm().barrier();
    comm().barrier();
    comm().barrier();
  }
}

//=======================================================================
template<typename OrdinalType>
void ElementSpace<OrdinalType>::directorySetup() {
  if(getNumGlobalElements() != 0)
    if(ElementSpaceData_->Directory_ == 0)
      ElementSpaceData_->Directory_ = platform().createDirectory(*this); // Make directory
}

//=============================================================================
// end Tpetra_ElementSpace.cpp

} // Tpetra namespace

#endif // _TPETRA_ELEMENTSPACE_HPP_
