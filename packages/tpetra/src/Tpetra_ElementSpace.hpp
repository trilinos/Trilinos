/*Paul
11-June-2002 Initial writeup.
12-June-2002 Initial writeup finished.
15-June-2002 changed linear to contiguous, changed some ints to ordinalTypes.
18-June-2002 Added private variables
20-June-2002 Cosmetic name changes, other trivial changes.
21-June-2002 STL calls done.
01-July-2002 Everything finished, myGlobalElements() left commented out.
06-July-2002 myGlobalElements() implemented.
09-July-2002 Directory calls taken out. (cout instead). Untested.
24-July-2002 everything including directory. untemplated. gMGE & print still not const.
27-July-2002 gMGE & print const, templated for OrdinalType.
06-Aug-2002 Switched to images.
03-Sept-2002 Added == and != operators.
21-Sept-2002 Comm/Platform split.
07-Oct-2002 ElementSpaceData move started
12-Nov-2002 Updated to use createOrdinalComm() instead of createComm() (nothing changed)
*/

#ifndef _TPETRA_ELEMENTSPACE_HPP_
#define _TPETRA_ELEMENTSPACE_HPP_

#include "Tpetra_Object.hpp"
#include "Tpetra_Directory.hpp"
#include <map>
#include <boost/shared_ptr.hpp>
#include "Tpetra_ElementSpaceData.hpp"

namespace Tpetra {

// forward declarations
template<typename PacketType, typename OrdinalType> class Comm;
template<typename ScalarType, typename OrdinalType> class Platform;

//! Tpetra::ElementSpace: A class for constructing and using template<ordinalType> ElementSpaces.
/*! ElementSpace objects are defined to have an element size of 1. Variable element sizes are implemented in Tpetra::BlockElementSpace. Some ElementSpace methods throw exceptions, and should be enclosed in a try/catch block. All Tpetra_ElementSpace objects require a Tpetra_Comm object. Local IDs (LIDs) are always in the range indexBase to numMyElements - indexBase.

ElementSpace error codes (positive for non-fatal, negative for fatal):
  <ol>
  <li> +1  Specified Global ID not found on this processor.
  <li> +2  Specified Local ID not found on this processor.
  <li> +3  Pointer passed to getMyGlobalElements(ordinalType*) does not have child allocated. (Null pointer).
  <li> -1  numGlobalElements < -1.  Should be >= -1 (Should be >= 0 for a Tpetra-defined ElementSpace).
  <li> -2  numMyElements < indexBase.  Should be >= indexBase.
  <li> -3  Invalid numGlobalElements.  Should equal sum of myGlobalElements, or set to indexBase-1 to compute automatically.
  <li> -4  Minimum global element index is less than index base.
  <li> -99 Internal ElementSpace error.  Contact developer.
  </ol>*/


template<typename OrdinalType>
class ElementSpace : public Object {

public:
  
//@{ \name Constructor/Destructor Methods
  
//! Tpetra::ElementSpace constructor with Tpetra-defined contiguous uniform distribution.
ElementSpace(OrdinalType numGlobalElements, OrdinalType indexBase, const Platform<OrdinalType, OrdinalType>& Platform);

//! Tpetra::ElementSpace constructor with user-defined contiguous distribution.
ElementSpace(OrdinalType numGlobalElements, OrdinalType numMyElements, OrdinalType indexBase, 
						 const Platform<OrdinalType, OrdinalType>& Platform);

//! Tpetra::ElementSpace constructor with user-defined non-contiguous (arbitrary) distribution.
ElementSpace(OrdinalType numGlobalElements, OrdinalType numMyElements, OrdinalType* elementList, 
						 OrdinalType indexBase, const Platform<OrdinalType, OrdinalType>& Platform);

//! Tpetra::ElementSpace copy constructor.
ElementSpace(const ElementSpace<OrdinalType>& ElementSpace);

//! Tpetra::ElementSpace destructor.
~ElementSpace();
  
//@}


//@{ \name Local/Global ID Accessor Methods

//! Returns the image IDs and corresponding local index values for a given list of global indices.
void getRemoteIDList(OrdinalType numIDs, const OrdinalType* GIDList, OrdinalType* imageIDList, OrdinalType* LIDList) const 
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
bool isSameAs(const ElementSpace<OrdinalType>& ElementSpace) const;
bool operator==(const ElementSpace<OrdinalType>& ElementSpace) const {return(isSameAs(ElementSpace));};
bool operator!=(const ElementSpace<OrdinalType>& ElementSpace) const {return(!isSameAs(ElementSpace));};

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
const Comm<OrdinalType, OrdinalType>& comm() const {return(*ElementSpaceData_->Comm_);};
const Platform<OrdinalType, OrdinalType>& platform() const {return(*ElementSpaceData_->Platform_);};

//@}

private:
// private data members
boost::shared_ptr< ElementSpaceData<OrdinalType> > ElementSpaceData_; // Boost shared pointer

// private functions
void directorySetup();

}; // ElementSpace class
} // Tpetra namespace

#include "Tpetra_ElementSpace.cpp"

#endif // _TPETRA_ELEMENTSPACE_HPP_
