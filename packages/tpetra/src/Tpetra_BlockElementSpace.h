// 11-June-2002 Initial writeup.
// 12-June-2002 Initial writeup finished.
// 15-June-2002 Changed constructors to all take ElementSpace objects. Changed some ints to ordinalTypes.

#ifndef _TPETRA_BLOCKELEMENTSPACE_H_
#define _TPETRA_BLOCKELEMENTSPACE_H_


#include "Tpetra_Object.h"
#include "Tpetra_ElementSpace.h"

namespace Tpetra {

//! Tpetra::BlockElementSpace: A class for constructing and using template<ordinalType> BlockElementSpaces.
/*! BlockElementSpace objects can have variable element sizes. (If variable element sizes are not needed, a Tpetra::ElementSpace object should probably be used instead. Some BlockElementSpace methods throw exceptions, and should be enclosed in a try/catch block. All BlockElementSpace objects require a Tpetra_ElementSpace object, as well as a Tpetra_Comm object. */
template<class ordinalType> 
class BlockElementSpace : public Object {

public:

//@{ \name Constructor/Destructor Methods.

//! Tpetra::BlockElementSpace constructor with constant element size.
Tpetra_BlockElementSpace (Tpetra_ElementSpace<ordinalType> & ElementSpace, ordinalType elementSize);

//! Tpetra::BlockElementSpace constructor with arbitrary element sizes.
Tpetra_BlockElementSpace (Tpetra_ElementSpace<ordinalType> & ElementSpace, ordinalType * elementSizeList);

//! Tpetra::BlockElementSpace copy constructor.
Tpetra_BlockElementSpace (Tpetra_BlockElementSpace<ordinalType> & BlockElementSpace);

//! Tpetra::BlockElementSpace destructor.
~Tpetra_BlockElementSpace ();

//@}


//@{ \name Local/Global ID Accessor Methods

//! Returns the processor IDs, corresponding local index values, and element sizes for a given list of global indices.
/*! The processor IDs, local index values, and element sizes are placed into arrays passed in by the user. The list of global indices used to create these is also passed in by the user. Exceptions might be thrown. */ 
void getRemoteIDList (ordinalType numIDs, const ordinalType * GIDList, ordinalType * PIDList, ordinalType * LIDList, ordinalType * elementSizeList) const;

//! Returns the local ID of the element that contains the given local Point ID, and the offset of the point in that element.
/*! The local ID and offset are placed in ordinalType variables passed in by reference by the user. */
void findLocalElementID (ordinalType pointID, ordinalType & elementID, ordinalType & elementOffset) const;

//@}


//@{ \name Size & Dimension Accessor Methods

//! Returns the size of elements in the BlockElementSpace. Throws an exception of -2 if not all elements are the same size.
ordinalType getElementSize () const;

//! Returns the size of the element whose local ID is passed in. Throws an exception of -1 if the local ID is not found on the calling processor.
ordinalType getElementSize (ordinalType LID) const;

//! Returns the number of global points in the BlockElementSpace; equals the sum of all element sizes across all processors.
ordinalType getNumGlobalPoints () const;

//! Returns the number of global points on this processor; equals the sum of all element sizes on the calling processor.
ordinalType getNumMyPoints () const;

//! Returns the minimum element size on the calling processor.
ordinalType getMinMyElementSize () const;

//! Returns the maximum element size on the calling processor.
ordinalType getMaxMyElementSize () const;

//! Returns the minimum element size in the BlockElementSpace.
ordinalType getMinElementSize () const;

//! Returns the maximum element size in the BlockElementSpace.
ordinalType getMaxElementSize () const;

//@}


//@{ \name Misc. Boolean Tests

//! Returns true if all elements have a constant size, returns false otherwise.
bool constantElementSize () const;

//! Returns true if this BlockElementSpace is identical to the one passed in, returns false otherwise. Also implemented as the == operator.
bool sameAs (const Tpetra_ElementSpace<ordinalType> & ElementSpace) const;

//@}


//@{ \name Array Accessor Functions. Each of these methods is implemented twice, one that returns a pointer, and one that copies the array into one passed in by the user.

//! Returns a pointer to the internal array of the mapping between the local elements, and the first local point number in each element.
ordinalType * getFirstPointInElementList () const;
void getFirstPointInElementList (ordinalType * firstPointInElementList) const;

//! Returns a pointer to array of the sizes of all the elements that belong to the calling processor.
ordinalType * getElementSizeList () const;
void getElementSizeList (ordinalType * elementSizeList) const;

//! Returns a pointer to an array that lists the LID of the element that each point belongs to.
ordinalType * getPointToElementList () const;
void getPointToElementList (ordinalType * pointToElementList) const;

//@}


//@{ Misc.

//! Prints the BlockElementSpace object to the output stream.
void print (ostream & os) const;

//! Access function for Tpetra_ElementSpace object.
const Tpetra_ElementSpace<ordinalType> & ElementSpace () const {return(*ElementSpace_);};

//@}

 private:

  Tpetra_ElementSpace<ordinalType> * ElementSpace_;

}; // Tpetra_BlockElementSpace class

} // Tpetra namespace

//@{ \name defined outside of Tpetra namespace

//! Overloaded == operator.
const bool operator == (Tpetra_BlockElementSpace<ordinalType>, Tpetra_BlockElementSpace<ordinalType>) const;

//! Overloaded << operator. Used by the BlockElementSpace::print method.
inline ostream & operator << (ostream & os, const & Tpetra_BlockElementSpace::obj);

//@}

#endif _TPETRA_BLOCKELEMENTSPACE_H_
