/*Paul
11-June-2002 Initial writeup.
12-June-2002 Initial writeup finished.
15-June-2002 Changed constructors to all take ElementSpace objects. Changed some ints to ordinalTypes.
01-August-2002 Real writeup starts. Minor changes.
06-August-2002 Switched to images.
03-Sept-2002 Added == and != operators
16-Oct-2002 Updated for BESData
12-Nov-2002 Updated to use createOrdinalComm() instead of createComm() (nothing changed)
14-Nov-2002 Changed slightly for switch to massive constructor call
*/

#ifndef _TPETRA_BLOCKELEMENTSPACE_HPP_
#define _TPETRA_BLOCKELEMENTSPACE_HPP_

#include "Tpetra_Object.hpp"
#include "Tpetra_ElementSpace.hpp"
#include <boost/shared_ptr.hpp>

namespace Tpetra {

template<typename OrdinalType> class BlockElementSpaceData;

//! Tpetra::BlockElementSpace: A class for constructing and using template<OrdinalType> BlockElementSpaces.
/*! BlockElementSpace objects can have variable element sizes. (If variable element sizes are not needed, an ElementSpace object should probably be used instead.) Some BlockElementSpace methods throw exceptions, and should be enclosed in a try/catch block. All BlockElementSpace objects require an ElementSpace object, which requires a Comm object.  

BlockElementSpace error codes (positive for non-fatal, negative for fatal):
  <ol>
  <li> +1  Specified Point ID not found on this image.
  <li> +2  Specified Local ID not found on this image.
  <li> +3  elementSize requested in a variable-sized BlockElementSpace.
  <li> +4  Pointer passed to getFirstPointInElementList, getElementSizeList, or getPointToElementList does not have child allocated. (Null pointer).
  <li> -1  elementSize (or element in elementSizeList) <= 0.  Should be > 0.
  <li> -99 Internal BlockElementSpace error.  Contact developer.
  </ol>*/

template<typename OrdinalType> 
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
 ~BlockElementSpace() { /*cout << "BES destructor called." << endl;*/};

//@}


//@{ \name Local/Global ID Accessor Methods

//! Returns the image IDs, corresponding local index values, and element sizes for a given list of global indices.
/*! Theimage IDs, local index values, and element sizes are placed into arrays passed in by the user. The list of global indices used to create these is also passed in by the user. Exceptions might be thrown. */ 
void getRemoteIDList(const OrdinalType numIDs, const OrdinalType* GIDList, OrdinalType* imageIDList, OrdinalType* LIDList, OrdinalType* elementSizeList) const;

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
OrdinalType getNumGlobalPoints() const {return(BlockElementSpaceData_->numGlobalPoints_);};

//! Returns the number of global points on this image; equals the sum of all element sizes on the calling image.
OrdinalType getNumMyPoints() const {return(BlockElementSpaceData_->numMyPoints_);};

//! Returns the minimum element size on the calling image.
OrdinalType getMinMyElementSize() const {return(BlockElementSpaceData_->minMySize_);};

//! Returns the maximum element size on the calling image.
OrdinalType getMaxMyElementSize() const {return(BlockElementSpaceData_->maxMySize_);};

//! Returns the minimum element size in the BlockElementSpace.
OrdinalType getMinElementSize() const {return(BlockElementSpaceData_->minGlobalSize_);};

//! Returns the maximum element size in the BlockElementSpace.
OrdinalType getMaxElementSize() const {return(BlockElementSpaceData_->maxGlobalSize_);};

//@}


//@{ \name Misc. Boolean Tests

//! Returns true if all elements have a constant size, returns false otherwise.
bool isConstantElementSize() const {return(BlockElementSpaceData_->constantSize_);};

//! Returns true if this BlockElementSpace is identical to the one passed in, returns false otherwise. Also implemented through the == and != operators.
bool isSameAs(const BlockElementSpace<OrdinalType>& BlockElementSpace) const;
bool operator==(const BlockElementSpace<OrdinalType>& BlockElementSpace) const {return(isSameAs(BlockElementSpace));};
bool operator!=(const BlockElementSpace<OrdinalType>& BlockElementSpace) const {return(!isSameAs(BlockElementSpace));};

//@}


//@{ \name Array Accessor Functions. Each of these methods is implemented twice, one that returns a pointer, and one that copies the array into one passed in by the user.

//! Returns a pointer to array of the sizes of all the elements that belong to the calling image.
const OrdinalType* getElementSizeList() const {return(BlockElementSpaceData_->elementSizeList_);};
void getElementSizeList(OrdinalType* elementSizeList) const;

//! Returns a pointer to the internal array of the mapping between the local elements, and the first local point number in each element.
const OrdinalType* getFirstPointInElementList() const;
void getFirstPointInElementList(OrdinalType* firstPointInElementList) const;

//! Returns a pointer to an array that lists the LID of the element that each point belongs to.
const OrdinalType* getPointToElementList() const;
void getPointToElementList(OrdinalType* pointToElementList) const;

//@}


//@{ Misc.

//! Prints the BlockElementSpace object to the output stream. (Used by the overloaded << operator inherited from Object).
void print(ostream& os) const;

//! Access function for ElementSpace object.
const ElementSpace<OrdinalType>& elementSpace() const {return(*BlockElementSpaceData_->ElementSpace_);};

//@}

private:
boost::shared_ptr<BlockElementSpaceData<OrdinalType> > BlockElementSpaceData_; // Boost shared pointer

}; // BlockElementSpace class

} // Tpetra namespace

#include "Tpetra_BlockElementSpaceData.hpp"
#include "Tpetra_BlockElementSpace.cpp"
#endif // _TPETRA_BLOCKELEMENTSPACE_HPP_
