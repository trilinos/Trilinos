#ifndef _TPETRA_VECTORSPACE_HPP_
#define _TPETRA_VECTORSPACE_HPP_

#include <Teuchos_RefCountPtr.hpp>
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
Teuchos::RefCountPtr< ElementSpace<OrdinalType> const > const ElementSpace_;
Teuchos::RefCountPtr< BlockElementSpace<OrdinalType> const > const BlockElementSpace_;
Teuchos::RefCountPtr< Platform<OrdinalType, ScalarType> const > const Platform_;
Teuchos::RefCountPtr< Comm<ScalarType, OrdinalType> const > const Comm_;

}; // VectorSpace class
} // Tpetra namespace

// begin Tpetra_Vectorspace.cpp
//=======================================================================

#include "Tpetra_ElementSpace.hpp"
#include "Tpetra_BlockElementSpace.hpp"
#include "Tpetra_Platform.hpp"
#include "Tpetra_Comm.hpp"

namespace Tpetra {

// ElementSpace constructor
template<typename OrdinalType, typename ScalarType>
VectorSpace<OrdinalType, ScalarType>::VectorSpace(ElementSpace<OrdinalType> const& elementSpace, 
																									Platform<OrdinalType, ScalarType> const& platform)
	: Object("Tpetra::VectorSpace")
	, blockspace_(false)
	, indexBase_(elementSpace.getIndexBase())
	, numMyEntries_(elementSpace.getNumMyElements())
	, numGlobalEntries_(elementSpace.getNumGlobalElements())
	, ElementSpace_(0)
	, BlockElementSpace_(0)
	, Platform_(0)
	, Comm_(0)
{
	ElementSpace_ = Teuchos::rcp(&elementSpace);
	Platform_ = Teuchos::rcp(&platform);
	Comm_ = Teuchos::rcp(platform.createScalarComm());
}

// BlockElementSpace constructor
//template<typename OrdinalType, typename ScalarType>
//VectorSpace<OrdinalType, ScalarType>::VectorSpace(BlockElementSpace<OrdinalType> const& blockElementSpace, 
//																									Platform<OrdinalType, ScalarType> const& platform) {
//}

// copy constructor.
//template<typename OrdinalType, typename ScalarType>
//VectorSpace<OrdinalType, ScalarType>::VectorSpace(VectorSpace<OrdinalType, ScalarType> const& vectorSpace) {
//}

// destructor.
template<typename OrdinalType, typename ScalarType>
VectorSpace<OrdinalType, ScalarType>::~VectorSpace() {
}

// getLocalIndex
template<typename OrdinalType, typename ScalarType>
OrdinalType VectorSpace<OrdinalType, ScalarType>::getLocalIndex(OrdinalType globalIndex) const {
	if(!blockspace_)
		return(elementSpace().getLID(globalIndex));
}

// getGlobalIndex
template<typename OrdinalType, typename ScalarType>
OrdinalType VectorSpace<OrdinalType, ScalarType>::getGlobalIndex(OrdinalType localIndex) const {
	if(!blockspace_)
		return(elementSpace().getGID(localIndex));
}

// compatibleVector
//template<typename OrdinalType, typename ScalarType>
//bool VectorSpace<OrdinalType, ScalarType>::compatibleVector(VectorSpace<ScalarType, OrdinalType> const& vector) const {
//}

// isSameAs
//template<typename OrdinalType, typename ScalarType>
//bool VectorSpace<OrdinalType, ScalarType>::isSameAs(VectorSpace<ScalarType, OrdinalType> const& vectorSpace) const {
//	if(blockspace_)
//		return(blockElementSpace().isSameAs(vectorSpace.blockElementSpace())); // compare BlockElementSpaces
//	else
//		return(elementSpace().isSameAs(vectorSpace.elementSpace())); // compare ElementSpaces
//}

// print
/*template<typename OrdinalType, typename ScalarType>
void VectorSpace<OrdinalType, ScalarType>::print(ostream& os) const { 
  OrdinalType myImageID = comm().getMyImageID(); 
  OrdinalType numImages = comm().getNumImages(); 
  
  for (int imageCtr = 0; imageCtr < numImages; imageCtr++) { 
    if (myImageID == imageCtr) { 
      if (myImageID == 0) { 
				os << endl << "Number of Global Entries  = " << getNumGlobalEntries() << endl; 
				os <<         "Index Base                = " << getIndexBase() << endl; 
      } 

      os << endl <<   "Number of Local Entries   = " << getNumMyEntries() << endl; 
      os << endl; 
		} 
	} 
	if(blockspace_) 
		blockElementSpace.print(); 
	else 
		elementSpace().print(); 
		}*/

} // namespace Tpetra

//=======================================================================
// end Tpetra_Vectorspace.cpp

#endif // _TPETRA_VECTORSPACE_HPP_
