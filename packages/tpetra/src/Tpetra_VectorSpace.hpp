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
#include "Tpetra_Object.hpp"
#include "Tpetra_ElementSpace.hpp"
#include "Tpetra_BlockElementSpace.hpp"
#include "Tpetra_Platform.hpp"
#include "Tpetra_Comm.hpp"

namespace Tpetra {

// forward declaration of Vector and VectorSpaceData, needed to prevent circular inclusions
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
        //, blockspace_(false)
        , VectorSpaceData_()
        , zero_(Teuchos::OrdinalTraits<OrdinalType>::zero())
        , one_(Teuchos::OrdinalTraits<OrdinalType>::one())
        //, indexBase_(elementSpace.getIndexBase())
        //, numMyEntries_(elementSpace.getNumMyElements())
        //, numGlobalEntries_(elementSpace.getNumGlobalElements())
        , ElementSpace_()
        , BlockElementSpace_()
        //, Platform_()
        //, Comm_()
	{
        VectorSpaceData_ = Teuchos::rcp(new VectorSpaceData<OrdinalType, ScalarType>(false, 
                                                                                     elementSpace.getIndexBase(), 
                                                                                     elementSpace.getNumMyElements(),
                                                                                     elementSpace.getNumGlobalElements(),
                                                                                     platform));
            
		ElementSpace_ = Teuchos::rcp(new ElementSpace<OrdinalType>(elementSpace));
		//Platform_ = Teuchos::rcp(platform.clone());
		//Comm_ = Teuchos::rcp(platform.createScalarComm());
	};
	
	//! Tpetra::VectorSpace constructor taking a BlockElementSpace object.
	VectorSpace(BlockElementSpace<OrdinalType> const& blockElementSpace, Platform<OrdinalType, ScalarType> const& platform)
		: Object("Tpetra::VectorSpace")
        //, blockspace_(true)
        , VectorSpaceData_()
        , zero_(Teuchos::OrdinalTraits<OrdinalType>::zero())
        , one_(Teuchos::OrdinalTraits<OrdinalType>::one())
        //, indexBase_(blockElementSpace.elementSpace().getIndexBase())
        //, numMyEntries_(blockElementSpace.getNumMyPoints())
        //, numGlobalEntries_(blockElementSpace.getNumGlobalPoints())
        , ElementSpace_()
        , BlockElementSpace_()
        //, Platform_()
        //, Comm_()
	{
        VectorSpaceData_ = Teuchos::rcp(new VectorSpaceData<OrdinalType, ScalarType>(true,
                                                                                     blockElementSpace.elementSpace().getIndexBase(),
                                                                                     blockElementSpace.getNumMyPoints(),
                                                                                     blockElementSpace.getNumGlobalPoints(),
                                                                                     platform));
        
		BlockElementSpace_ = Teuchos::rcp(new BlockElementSpace<OrdinalType>(blockElementSpace));
		ElementSpace_ = Teuchos::rcp(blockElementSpace.generateCompatibleElementSpace());
		//Platform_ = Teuchos::rcp(platform.clone());
		//Comm_ = Teuchos::rcp(platform.createScalarComm());
	};
	
	//! Tpetra::VectorSpace copy constructor.
	VectorSpace(VectorSpace<OrdinalType, ScalarType> const& vectorSpace)
		: Object(vectorSpace.label())
        , VectorSpaceData_(vectorSpace.VectorSpaceData_)
		//, blockspace_(vectorSpace.blockspace_)
		, zero_(vectorSpace.zero_)
		, one_(vectorSpace.one_)
		//, indexBase_(vectorSpace.indexBase_)
		//, numMyEntries_(vectorSpace.numMyEntries_)
		//, numGlobalEntries_(vectorSpace.numGlobalEntries_)
		, ElementSpace_(vectorSpace.ElementSpace_)
		, BlockElementSpace_(vectorSpace.BlockElementSpace_)
		//, Platform_(vectorSpace.Platform_)
		//, Comm_(vectorSpace.Comm_)
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
	OrdinalType getMinLocalIndex() const {return(zero_);};
	OrdinalType getMaxLocalIndex() const {return(zero_ + getNumMyEntries());};
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
		ScalarType sameNumLocal = 1.0; // we already know the length matches on this image
		ScalarType sameNumGlobal = 0.0;
		if(vectorSpace.getNumMyEntries() == getNumMyEntries())
			sameNumLocal = 1.0;
		comm().minAll(&sameNumLocal, &sameNumGlobal, 1);
		return(sameNumGlobal == 1.0);
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
		
		for (int imageCtr = 0; imageCtr < numImages; imageCtr++) {
			if (myImageID == imageCtr) {
				if (myImageID == 0) {
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
	
	
	//! Access function for the Tpetra::Platform and Tpetra::Comm communicators.
	Platform<OrdinalType, ScalarType> const& platform() const {return(*VectorSpaceData_->Platform_);};
	Comm<ScalarType, OrdinalType> const& comm() const {return(*VectorSpaceData_->Comm_);}; // Comm is <ST, OT> because ST represents PT
	
	//! Assignment operator (declared but not defined, do not use)
	VectorSpace<OrdinalType, ScalarType>& operator = (VectorSpace<OrdinalType, ScalarType> const& Source);

	//@}
	
private:
	
	ElementSpace<OrdinalType> const& elementSpace() const {return(*ElementSpace_);};
	BlockElementSpace<OrdinalType> const& blockElementSpace() const {return(*BlockElementSpace_);};
    bool isBlockSpace() const {return(VectorSpaceData_->blockspace_);};
    
    Teuchos::RefCountPtr< VectorSpaceData<OrdinalType, ScalarType> > VectorSpaceData_;
	
	//bool const blockspace_;
	OrdinalType const zero_;
	OrdinalType const one_;
	//OrdinalType const indexBase_;
	//OrdinalType const numMyEntries_;
	//OrdinalType const numGlobalEntries_;
	Teuchos::RefCountPtr< ElementSpace<OrdinalType> const > ElementSpace_;
	Teuchos::RefCountPtr< BlockElementSpace<OrdinalType> const > BlockElementSpace_;
	//Teuchos::RefCountPtr< Platform<OrdinalType, ScalarType> const > Platform_;
	//Teuchos::RefCountPtr< Comm<ScalarType, OrdinalType> const > Comm_; // Comm is <ST, OT> because ST represents PT
	
}; // VectorSpace class

} // Tpetra namespace

#include "Tpetra_Vector.hpp"
#include "Tpetra_VectorSpaceData.hpp"

#endif // _TPETRA_VECTORSPACE_HPP_
