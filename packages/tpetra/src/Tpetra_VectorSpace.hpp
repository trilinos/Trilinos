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

#ifndef _TPETRA_VECTORSPACE_HPP_
#define _TPETRA_VECTORSPACE_HPP_

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include "Tpetra_Object.hpp"
#include "Tpetra_ElementSpace.hpp"
#include "Tpetra_BlockElementSpace.hpp"
#include "Tpetra_Platform.hpp"
#include "Tpetra_Comm.hpp"

namespace Tpetra {

// forward declaration of Vector and VectorSpaceData, needed to prevent circular inclusions
// actual #include statements are at the end of this file
template<typename OrdinalType, typename ScalarType> class Vector;
template<typename OrdinalType, typename ScalarType> class VectorSpaceData;

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
	VectorSpace(ElementSpace<OrdinalType> const& elementSpace, Platform<OrdinalType, ScalarType> const& platform)
		: Object("Tpetra::VectorSpace")
		, VectorSpaceData_()
	{
		VectorSpaceData_ = Teuchos::rcp(new VectorSpaceData<OrdinalType, ScalarType>(false, 
																																								 elementSpace.getIndexBase(), 
																																								 elementSpace.getNumMyElements(),
																																								 elementSpace.getNumGlobalElements(),
																																								 platform));
		
		VectorSpaceData_->ElementSpace_ = Teuchos::rcp(new ElementSpace<OrdinalType>(elementSpace));
	};
	
	//! Tpetra::VectorSpace constructor taking a BlockElementSpace object.
	VectorSpace(BlockElementSpace<OrdinalType> const& blockElementSpace, Platform<OrdinalType, ScalarType> const& platform)
		: Object("Tpetra::VectorSpace")
		, VectorSpaceData_()
	{
		VectorSpaceData_ = Teuchos::rcp(new VectorSpaceData<OrdinalType, ScalarType>(true,
																																								 blockElementSpace.elementSpace().getIndexBase(),
																																								 blockElementSpace.getNumMyPoints(),
																																								 blockElementSpace.getNumGlobalPoints(),
																																								 platform));
		
		VectorSpaceData_->BlockElementSpace_ = Teuchos::rcp(new BlockElementSpace<OrdinalType>(blockElementSpace));
		VectorSpaceData_->ElementSpace_ = Teuchos::rcp(blockElementSpace.generateCompatibleElementSpace());
	};
	
	//! Tpetra::VectorSpace copy constructor.
	VectorSpace(VectorSpace<OrdinalType, ScalarType> const& vectorSpace)
		: Object(vectorSpace.label())
		, VectorSpaceData_(vectorSpace.VectorSpaceData_)
	{}

	//! Tpetra::VectorSpace destructor.
	~VectorSpace() {};
  
	//@}
	
	
	//@{ \name VectorSpace Attribute Methods
	
	//! Returns the number of entries in this VectorSpace.
	OrdinalType getNumGlobalEntries() const {return(VectorSpaceData_->numGlobalEntries_);};
	
	//! Returns the number of entries belonging to the calling image.
	OrdinalType getNumMyEntries() const {return(VectorSpaceData_->numMyEntries_);};
	
	//! Returns the index base for this VectorSpace.
	OrdinalType getIndexBase() const {return(VectorSpaceData_->indexBase_);};
	
	//! Min/Max Indices
	OrdinalType getMinLocalIndex() const {return(Teuchos::OrdinalTraits<OrdinalType>::zero());};
	OrdinalType getMaxLocalIndex() const {return(Teuchos::OrdinalTraits<OrdinalType>::zero() + getNumMyEntries());};
	OrdinalType getMinGlobalIndex() const {return(getIndexBase());};
	OrdinalType getMaxGlobalIndex() const {return(getIndexBase() + getNumGlobalEntries());};
	
	//! Return the local index for a given global index
	/*! If this VectorSpace was created using a BlockElementSpace,
		  LIDs and GIDs from the compatible ElementSpace will be used.
	*/
	OrdinalType getLocalIndex(OrdinalType globalIndex) const {
		return(elementSpace().getLID(globalIndex));
	};
	
	//! Return the global index for a given local index
	/*! If this VectorSpace was created using a BlockElementSpace,
		  LIDs and GIDs from the compatible ElementSpace will be used.
	*/
	OrdinalType getGlobalIndex(OrdinalType localIndex) const {
		return(elementSpace().getGID(localIndex));
	};

  //! Returns true if the local index value passed in is found on the calling image, returns false if it doesn't.
  bool isMyLocalIndex(OrdinalType localIndex) const {
    return(elementSpace().isMyLID(localIndex));
  }

  //! Returns true if the global index value passed in is found the calling image, returns false if it doesn't.
  bool isMyGlobalIndex(OrdinalType globalIndex) const {
    return(elementSpace().isMyGID(globalIndex));
  }
	
	//@}
	
	//@{ \name Vector Creation
	
	//! Creates a Tpetra::Vector, with all entries set to 0.
	Vector<OrdinalType, ScalarType>* createVector() const {
		Vector<OrdinalType, ScalarType>* vector = new Vector<OrdinalType, ScalarType>(*this);
		return(vector);
	};
	
	//@{ \name Boolean Tests
	
	//! Returns true if the VectorSpace passed in is compatible with this VectorSpace.
	bool isCompatible(VectorSpace<OrdinalType, ScalarType> const& vectorSpace) const {
		// first check global length and local lenght on this image
		if(vectorSpace.getNumGlobalEntries() != getNumGlobalEntries() ||
			 vectorSpace.getNumMyEntries() != getNumMyEntries())
			return(false);
		
		// then check to make sure distribution is the same
		ScalarType sameNumLocal = Teuchos::ScalarTraits<ScalarType>::one(); // we already know the length matches on this image
		ScalarType sameNumGlobal = Teuchos::ScalarTraits<ScalarType>::zero();
		comm().minAll(&sameNumLocal, &sameNumGlobal, Teuchos::OrdinalTraits<OrdinalType>::one());
		return(sameNumGlobal == Teuchos::ScalarTraits<ScalarType>::one());
	};
	
	//! Returns true if the VectorSpace passed in is identical to this VectorSpace. Also implemented through the == and != operators.
	bool isSameAs(VectorSpace<OrdinalType, ScalarType> const& vectorSpace) const {
		if(isBlockSpace())
			return(blockElementSpace().isSameAs(vectorSpace.blockElementSpace())); // compare BlockElementSpaces
		else
			return(elementSpace().isSameAs(vectorSpace.elementSpace())); // compare ElementSpaces
	};
	bool operator==(VectorSpace<OrdinalType, ScalarType> const& vectorSpace) const {return(isSameAs(vectorSpace));};
	bool operator!=(VectorSpace<OrdinalType, ScalarType> const& vectorSpace) const {return(!isSameAs(vectorSpace));};
	
	//@}
	
	//@{ \name Misc.
	
	//! Prints the VectorSpace object to the output stream.
	/*! An << operator is inherited from Tpetra::Object, which uses the print method.*/
	void print(ostream& os) const {
		OrdinalType myImageID = comm().getMyImageID();
		OrdinalType numImages = comm().getNumImages();
		OrdinalType const ordinalZero = Teuchos::OrdinalTraits<OrdinalType>::zero();
		
		for (int imageCtr = ordinalZero; imageCtr < numImages; imageCtr++) {
			if (myImageID == imageCtr) {
				if (myImageID == ordinalZero) {
					os << endl << "Number of Global Entries  = " << getNumGlobalEntries() << endl;
					os <<         "Index Base                = " << getIndexBase() << endl;
				}
				os << endl <<   "ImageID = " << myImageID << endl;
				os <<           "Number of Local Entries   = " << getNumMyEntries() << endl;
				os << endl;
			}
		}
		if(isBlockSpace()) {
			os << "Built on a BlockElementSpace" << endl;
			blockElementSpace().print(os);
			os << "Compatible ElementSpace:" << endl;
			elementSpace().print(os);
		}
		else {
			os << "Built on an ElementSpace" << endl;
			elementSpace().print(os);
		}
	};
	
	
	//! Access functions for the Tpetra::Platform and Tpetra::Comm communicators.
	Platform<OrdinalType, ScalarType> const& platform() const {return(*VectorSpaceData_->Platform_);};
	Comm<ScalarType, OrdinalType> const& comm() const {return(*VectorSpaceData_->Comm_);}; // Comm is <ST, OT> because ST represents PT

	//! Access function for the ElementSpace used by this VectorSpace
	ElementSpace<OrdinalType> const& elementSpace() const {return(*VectorSpaceData_->ElementSpace_);};
	
	//! Assignment operator
	VectorSpace<OrdinalType, ScalarType>& operator = (VectorSpace<OrdinalType, ScalarType> const& Source) {
    VectorSpaceData_ = Source.VectorSpaceData_;
    return(*this);
  }

	//@}
	
private:
	
	BlockElementSpace<OrdinalType> const& blockElementSpace() const {return(*VectorSpaceData_->BlockElementSpace_);};
	bool isBlockSpace() const {return(VectorSpaceData_->blockspace_);};
    
	Teuchos::RefCountPtr< VectorSpaceData<OrdinalType, ScalarType> > VectorSpaceData_;

}; // VectorSpace class

} // Tpetra namespace

#include "Tpetra_Vector.hpp"
#include "Tpetra_VectorSpaceData.hpp"

#endif // _TPETRA_VECTORSPACE_HPP_
