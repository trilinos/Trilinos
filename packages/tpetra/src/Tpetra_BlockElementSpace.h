/*Paul
11-June-2002 Initial writeup.
12-June-2002 Initial writeup finished.
15-June-2002 Changed constructors to all take ElementSpace objects. Changed some ints to ordinalTypes.
01-August-2002 Real writeup starts. Minor changes.
06-August-2002 Switched to images.
*/

#ifndef _TPETRA_BLOCKELEMENTSPACE_H_
#define _TPETRA_BLOCKELEMENTSPACE_H_

#include "Tpetra_Object.h"
#include "Tpetra_ElementSpace.h"

namespace Tpetra {

//! Tpetra::BlockElementSpace: A class for constructing and using template<OrdinalType> BlockElementSpaces.
/*! BlockElementSpace objects can have variable element sizes. (If variable element sizes are not needed, an ElementSpace object should probably be used instead.) Some BlockElementSpace methods throw exceptions, and should be enclosed in a try/catch block. All BlockElementSpace objects require an ElementSpace object, which requires a Comm object.  

BlockElementSpace error codes (positive for non-fatal, negative for fatal):
  <ol>
  <li> +1  Specified Point ID not found on this image.
  <li> +2  Specified Local ID not found on this image.
  <li> +3  elementSize requested in a variable-sized BlockElementSpace.
  <li> -1  elementSize (or element in elementSizeList) <= 0.  Should be > 0.
  <li> -99 Internal BlockElementSpace error.  Contact developer.
  </ol>*/

template<class OrdinalType> 
class BlockElementSpace : public Object {

public:

//@{ \name Constructor/Destructor Methods.

//! Tpetra::BlockElementSpace constructor with constant element size.
BlockElementSpace(ElementSpace<OrdinalType>& ElementSpace, OrdinalType elementSize);

//! Tpetra::BlockElementSpace constructor with arbitrary element sizes.
BlockElementSpace(ElementSpace<OrdinalType>& ElementSpace, OrdinalType* elementSizeList);

//! Tpetra::BlockElementSpace copy constructor.
BlockElementSpace(BlockElementSpace<OrdinalType>& BlockElementSpace);

//! Tpetra::BlockElementSpace destructor.
~BlockElementSpace();

//@}


//@{ \name Local/Global ID Accessor Methods

//! Returns the image IDs, corresponding local index values, and element sizes for a given list of global indices.
/*! Theimage IDs, local index values, and element sizes are placed into arrays passed in by the user. The list of global indices used to create these is also passed in by the user. Exceptions might be thrown. */ 
void getRemoteIDList(OrdinalType numIDs, const OrdinalType* GIDList, OrdinalType* imageIDList, OrdinalType* LIDList, OrdinalType* elementSizeList) const;

//! Returns the local ID of the element that contains the given local Point ID, and the offset of the point in that element.
/*! The local ID and offset are placed in OrdinalType variables passed in by reference by the user. */
void getLocalElementID(OrdinalType pointID, OrdinalType& elementID, OrdinalType& elementOffset) const;

//@}


//@{ \name Size & Dimension Accessor Methods

//! Returns the size of elements in the BlockElementSpace. Throws an exception of +2 if not all elements are the same size.
OrdinalType getElementSize() const;

//! Returns the size of the element whose local ID is passed in. Throws an exception of +1 if the local ID is not found on the calling image.
OrdinalType getElementSize(OrdinalType LID) const;

//! Returns the number of global points in the BlockElementSpace; equals the sum of all element sizes across all images.
OrdinalType getNumGlobalPoints() const {return(numGlobalPoints_);};

//! Returns the number of global points on this image; equals the sum of all element sizes on the calling image.
OrdinalType getNumMyPoints() const {return(numMyPoints_);};

//! Returns the minimum element size on the calling image.
OrdinalType getMinMyElementSize() const {return(minMySize_);};

//! Returns the maximum element size on the calling image.
OrdinalType getMaxMyElementSize() const {return(maxMySize_);};

//! Returns the minimum element size in the BlockElementSpace.
OrdinalType getMinElementSize() const {return(minGlobalSize_);};

//! Returns the maximum element size in the BlockElementSpace.
OrdinalType getMaxElementSize() const {return(maxGlobalSize_);};

//@}


//@{ \name Misc. Boolean Tests

//! Returns true if all elements have a constant size, returns false otherwise.
bool isConstantElementSize() const {return(constantSize_);};

//! Returns true if this BlockElementSpace is identical to the one passed in, returns false otherwise. Also implemented as the == operator.
bool isSameAs(const BlockElementSpace<OrdinalType>& BlockElementSpace) const;

//@}


//@{ \name Array Accessor Functions. Each of these methods is implemented twice, one that returns a pointer, and one that copies the array into one passed in by the user.

//! Returns a pointer to the internal array of the mapping between the local elements, and the first local point number in each element.
OrdinalType* getFirstPointInElementList() const;
void getFirstPointInElementList(OrdinalType* firstPointInElementList) const;

//! Returns a pointer to array of the sizes of all the elements that belong to the calling image.
OrdinalType* getElementSizeList() const {return(elementSizeList_);};
void getElementSizeList(OrdinalType* elementSizeList) const;

//! Returns a pointer to an array that lists the LID of the element that each point belongs to.
OrdinalType* getPointToElementList() const;
void getPointToElementList(OrdinalType* pointToElementList) const;

//@}


//@{ Misc.

//! Prints the BlockElementSpace object to the output stream. (Used by the overloaded << operator inherited from Object).
void print(ostream& os) const;

//! Access function for ElementSpace object.
const ElementSpace<OrdinalType>& elementSpace() const {return(*ElementSpace_);};

//@}

private:

bool constantSize_;
OrdinalType elementSize_;
OrdinalType numMyPoints_;
OrdinalType numGlobalPoints_;
OrdinalType minMySize_;
OrdinalType maxMySize_;
OrdinalType minGlobalSize_;
OrdinalType maxGlobalSize_;
OrdinalType* elementSizeList_;
mutable OrdinalType* pointToElementList_;
mutable OrdinalType* firstPointList_;
ElementSpace<OrdinalType>* ElementSpace_;

}; // BlockElementSpace class

} // Tpetra namespace

#include "Tpetra_BlockElementSpace.cpp"

#endif // _TPETRA_BLOCKELEMENTSPACE_H_
