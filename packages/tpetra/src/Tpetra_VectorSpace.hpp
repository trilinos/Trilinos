/*Paul
26-Jan-2003 Initial writeup. Unimplemented.
05-Feb-2003
*/

#ifndef _TPETRA_VECTORSPACE_HPP_
#define _TPETRA_VECTORSPACE_HPP_

#include <boost/shared_ptr.hpp>
#include "Tpetra_Object.hpp"
//#include "Tpetra_Vector.hpp"

namespace Tpetra {

// forward declarations
template<typename OrdinalType> class ElementSpace;
template<typename OrdinalType> class BlockElementSpace;
template<typename OrdinalType, typename ScalarType> class Platform;
template<typename ScalarType, typename OrdinalType> class Comm;

//! Tpetra::VectorSpace
/*! VectorSpace serves two purposes. In addition to creating Tpetra::Vectors,
	  it acts as an "insulating" class between Vectors and ElementSpace/BlockElementSpace.
		Through this mechanism, Vectors can be created and manipulated using one nonambiguous
		set of vector indices, regardless of if it uses an ElementSpace or a BlockElementSpace
		for distribution.
*/

template<typename OrdinalType, typename ScalarType>
class VectorSpace : public Object {

public:
  
//@{ \name Constructor/Destructor Methods
  
//! Tpetra::VectorSpace constructor taking an ElementSpace object.
VectorSpace(ElementSpace<OrdinalType> const& elementSpace, Platform<OrdinalType, ScalarType> const& platform);

//! Tpetra::VectorSpace constructor taking a BlockElementSpace object.
//VectorSpace(BlockElementSpace<OrdinalType> const& blockElementSpace, Platform<OrdinalType, ScalarType> const& platform);

//! Tpetra::VectorSpace copy constructor.
//VectorSpace(VectorSpace<OrdinalType, ScalarType> const& vectorSpace);

//! Tpetra::VectorSpace destructor.
~VectorSpace();
  
//@}


//@{ \name VectorSpace Attribute Methods

//! Returns the number of entries in this VectorSpace.
OrdinalType getNumGlobalEntries() const {return(numGlobalEntries_);};

//! Returns the number of entries belonging to the calling image.
OrdinalType getNumMyEntries() const {return(numMyEntries_);};

//! Returns the index base for this VectorSpace.
OrdinalType getIndexBase() const {return(indexBase_);};

//! Min/Max Indices
OrdinalType getMinLocalIndex() const {return(getIndexBase());};
OrdinalType getMaxLocalIndex() const {return(getIndexBase() + getNumMyEntries());};
OrdinalType getMinGlobalIndex() const {return(getIndexBase());};
OrdinalType getMaxGlobalIndex() const {return(getIndexBase() + getNumGlobalEntries());};

//! Return the local index for a given global index
OrdinalType getLocalIndex(OrdinalType globalIndex) const;

//! Return the global index for a given local index
OrdinalType getGlobalIndex(OrdinalType localIndex) const;

//@}

//@{ \name Vector Creation

//! Creates a Tpetra::Vector, with all entries set to 0.
//Vector<OrdinalType, ScalarType>* createVector() const {
//	Vector<OrdinalType, ScalarType>* vector = new Vector<OrdinalType, ScalarType>(this);
//	return(vector);
//};

//@{ \name Boolean Tests

//! Returns true if the Vector passed in is compatible with this VectorSpace.
//bool compatibleVector(VectorSpace<ScalarType, OrdinalType> const& vector) const;

//! Returns true if the VectorSpace passed in is identical to this VectorSpace. Also implemented through the == and != operators.
//bool isSameAs(VectorSpace<ScalarType, OrdinalType> const& vectorSpace) const;
//bool operator==(VectorSpace<ScalarType, OrdinalType> const& vectorSpace) const {return(isSameAs(vectorSpace));};
//bool operator!=(VectorSpace<ScalarType, OrdinalType> const& vectorSpace) const {return(!isSameAs(vectorSpace));};

//@}

//@{ \name Misc.

//! Prints the VectorSpace object to the output stream.
/*! An << operator is inherited from Tpetra::Object, which uses the print method.*/
//void print(ostream& os) const;

//! Access function for the Tpetra::Platform and Tpetra::Comm communicators.
Platform<OrdinalType, ScalarType> const& platform() const {return(*Platform_);};
Comm<OrdinalType, ScalarType> const& comm() const {return(*Comm_);};

//@}

private:

ElementSpace<OrdinalType> const& elementSpace() const {return(*ElementSpace_);};
BlockElementSpace<OrdinalType> const& blockElementSpace() const {return(*BlockElementSpace_);};

bool const blockspace_;
OrdinalType const indexBase_;
OrdinalType const numMyEntries_;
OrdinalType const numGlobalEntries_;
boost::shared_ptr< ElementSpace<OrdinalType> const > const ElementSpace_;
boost::shared_ptr< BlockElementSpace<OrdinalType> const > const BlockElementSpace_;
boost::shared_ptr< Platform<OrdinalType, ScalarType> const > const Platform_;
boost::shared_ptr< Comm<ScalarType, OrdinalType> const > const Comm_;

}; // VectorSpace class
} // Tpetra namespace

#include "Tpetra_VectorSpace.cpp"

#endif // _TPETRA_VECTORSPACE_HPP_
