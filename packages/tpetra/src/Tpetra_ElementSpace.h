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
*/

#ifndef _TPETRA_ELEMENTSPACE_H_
#define _TPETRA_ELEMENTSPACE_H_

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

#include "Tpetra_Object.h"
#include "Tpetra_Directory.h"
#include <map>

namespace Tpetra {

// forward declarations
template<typename OrdinalType> class Comm;

template<typename OrdinalType>
class ElementSpace : public Tpetra::Object {

public:
  
//@{ \name Constructor/Destructor Methods
  
//! Tpetra::ElementSpace constructor with Tpetra-defined contiguous uniform distribution.
ElementSpace(OrdinalType numGlobalElements, OrdinalType indexBase, const Comm<OrdinalType>& Comm);

//! Tpetra::ElementSpace constructor with user-defined contiguous distribution.
ElementSpace(OrdinalType numGlobalElements, OrdinalType numMyElements, OrdinalType indexBase, const Comm<OrdinalType>& Comm);

//! Tpetra::ElementSpace constructor with user-defined non-contiguous (arbitrary) distribution.
ElementSpace(OrdinalType numGlobalElements, OrdinalType numMyElements, OrdinalType* elementList, OrdinalType indexBase, const Comm<OrdinalType>& Comm);

//! Tpetra::ElementSpace copy constructor.
ElementSpace(const ElementSpace<OrdinalType>& ElementSpace);

//! Tpetra::ElementSpace destructor.
~ElementSpace();
  
//@}


//@{ \name Local/Global ID Accessor Methods

//! Returns the processor IDs and corresponding local index values for a given list of global indices.
void getRemoteIDList(OrdinalType numIDs, const OrdinalType* GIDList, OrdinalType* PIDList, OrdinalType* LIDList) const {Directory_->getDirectoryEntries(numIDs, GIDList, PIDList, LIDList);};

//! Returns local ID of global ID passed in, throws exception -1 if not found on this processor.
OrdinalType getLID(OrdinalType GID) const;

//! Returns global ID of local ID passed in, throws exception -1 if not found on this processor.
OrdinalType getGID(OrdinalType LID) const;

//! Returns true if global ID passed in belongs to the calling processor, returns false if it doesn't.
bool isMyGID(OrdinalType GID) const;

//! Returns true if the local ID passed in belongs to the calling processor, returns false if it doesn't.
bool isMyLID(OrdinalType LID) const;

//! Returns the minimum global ID in this ElementSpace.
OrdinalType getMinAllGID() const {return(minAllGID_);};

//! Returns the maximum global ID in this ElementSpace.
OrdinalType getMaxAllGID() const {return(maxAllGID_);};

//! Returns the minimum global ID owned by this processor.
OrdinalType getMinMyGID() const {return(minMyGID_);};

//! Returns the maximum global ID owned by this processor.
OrdinalType getMaxMyGID() const {return(maxMyGID_);};

//! Returns the minimum local ID on the calling processor.
OrdinalType getMinLID() const {return(minLID_);};

//! Returns the maximum local ID on the calling processor.
OrdinalType getMaxLID() const {return(maxLID_);};

//@}


//@{ \name Size & Dimension Accessor Methods

//! Returns the number of elements in this ElementSpace.
OrdinalType getNumGlobalElements() const {return(numGlobalElements_);};

//! Returns the number of elements belonging to the calling processor.
OrdinalType getNumMyElements() const {return(numMyElements_);};

//! Puts list of global elements on this processor into the user-provided array.
void getMyGlobalElements(OrdinalType* elementList) const;

//! Returns the Index base for this ElementSpace. Normally 0 for C/C++ or 1 for Fortran, but can be anything.
OrdinalType getIndexBase() const {return(indexBase_);};

//@}


//@{ \name Misc. Boolean Tests

//! Returns true if this ElementSpace is distributed contiguously, returns false otherwise.
bool isContiguous() const {return(contiguous_);};

//! Returns true if this ElementSpace is distributed across more than one processor, returns false otherwise.
bool isGlobal() const {return(global_);};

//! Returns true if the ElementSpace passed in is identical to this ElementSpace. Also implemented as the == operator.
bool isSameAs(const ElementSpace<OrdinalType>& ElementSpace) const;

//@}


//@{ \name Array Accessor Functions

//! Returns a pointer to an internal array containing a list of the Global IDs assigned to the calling processor.
OrdinalType* getMyGlobalElements() const;

//@}


//@{ \name Misc.

//! Prints the ElementSpace object to the output stream.
/*! An << operator is inherited from Tpetra::Object, which uses the print method.*/
void print(ostream& os) const;

//! Access function for the Tpetra_Comm communicator.
const Comm<OrdinalType>& comm() const {return(*Comm_);};

//@}

private:
// private data members
mutable map<OrdinalType, OrdinalType> lgMap_; // mutable because STL is finicky about iterators
map<OrdinalType, OrdinalType> glMap_;         // and gMGE(OT*) is only logically const, not truely const
OrdinalType numGlobalElements_;
OrdinalType numMyElements_;
OrdinalType indexBase_;
OrdinalType minLID_;
OrdinalType maxLID_;
OrdinalType minMyGID_;
OrdinalType maxMyGID_;
OrdinalType minAllGID_;
OrdinalType maxAllGID_;
bool contiguous_;
bool global_;
mutable OrdinalType* myGlobalElements_; // mutable so that gMGE() can modify it

// private functions
bool checkGlobalness (OrdinalType numGlobalElements, OrdinalType numMyElements);
void directorySetup();

// private instances of other classes
const Comm<OrdinalType>* Comm_;
Directory<OrdinalType>* Directory_;

}; // ElementSpace class
} // Tpetra namespace

#include "Tpetra_ElementSpace.cpp"

#endif // _TPETRA_ELEMENTSPACE_H_
