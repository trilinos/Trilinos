/*Paul
27-Jan-2003 Updated for .hpp and for new const syntax.
06-Feb-2003 Ordering of ES<->ESData mutuality updated.
*/

#ifndef _TPETRA_ELEMENTSPACE_HPP_
#define _TPETRA_ELEMENTSPACE_HPP_

#include "Tpetra_Object.hpp"
#include "Tpetra_Directory.hpp"
#include <map>
#include <boost/shared_ptr.hpp>

namespace Tpetra {

// forward declarations
template<typename PacketType, typename OrdinalType> class Comm;
template<typename ScalarType, typename OrdinalType> class Platform;
template<typename OrdinalType> class ElementSpaceData;

//! Tpetra::ElementSpace: A class for constructing and using template<ordinalType> ElementSpaces.
/*! ElementSpace objects are defined to have an element size of 1. Variable element sizes are implemented 
	in Tpetra::BlockElementSpace. Some ElementSpace methods throw exceptions, and should be enclosed 
	in a try/catch block. All Tpetra_ElementSpace objects require a Tpetra_Comm object. 
	Local IDs (LIDs) are always in the range indexBase to numMyElements - indexBase.

	ElementSpace error codes (positive for non-fatal, negative for fatal):
  <ol>
  <li> +1  Specified Global ID not found on this image.
  <li> +2  Specified Local ID not found on this image.
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
boost::shared_ptr< ElementSpaceData<OrdinalType> > ElementSpaceData_; // Boost shared pointer

// private functions
void directorySetup();

}; // ElementSpace class
} // Tpetra namespace

#include "Tpetra_ElementSpaceData.hpp"
#include "Tpetra_ElementSpace.cpp"

#endif // _TPETRA_ELEMENTSPACE_HPP_
