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

#ifndef TPETRA_ELEMENTSPACEDATA_HPP
#define TPETRA_ELEMENTSPACEDATA_HPP

#include <Teuchos_OrdinalTraits.hpp>

namespace Tpetra {

template<typename OrdinalType>
class ElementSpaceData : public Object {
	friend class ElementSpace<OrdinalType>;
 public:
	ElementSpaceData(OrdinalType const indexBase, 
									 OrdinalType const numGlobalElements,
									 OrdinalType const numMyElements,
									 OrdinalType const minAllGID,
									 OrdinalType const maxAllGID,
									 OrdinalType const minMyGID,
									 OrdinalType const maxMyGID,
									 map<OrdinalType, OrdinalType> const lgMap,
									 map<OrdinalType, OrdinalType> const glMap,
									 bool const contiguous,
									 Platform<OrdinalType, OrdinalType> const& platform,
									 Comm<OrdinalType, OrdinalType> const* comm)
		: Object("Tpetra::ElementSpaceData")
		, Platform_(&platform) 
		, Comm_(comm) 
		, numGlobalElements_(numGlobalElements)
		, numMyElements_(numMyElements)
		, indexBase_(indexBase)
		, minLID_(Teuchos::OrdinalTraits<OrdinalType>::zero())
		, maxLID_(minLID_ + numMyElements_ - Teuchos::OrdinalTraits<OrdinalType>::one())
		, minMyGID_(minMyGID)
		, maxMyGID_(maxMyGID)
		, minAllGID_(minAllGID)
		, maxAllGID_(maxAllGID)
		, contiguous_(contiguous)
		, global_(checkGlobalness())
		, lgMap_(lgMap)
		, glMap_(glMap)
		, myGlobalElements_(0)
		, Directory_(0)
		{};

	~ElementSpaceData() {
		if(Directory_ != 0) {
			delete Directory_;
			Directory_ = 0;
		}
		if(myGlobalElements_ != 0) {
			delete [] myGlobalElements_;
			myGlobalElements_ = 0;
		}
		if(Comm_ != 0) {
			delete Comm_;
			Comm_ = 0;
		}
	};

 protected:
	Platform<OrdinalType, OrdinalType> const* Platform_;
	Comm<OrdinalType, OrdinalType> const* Comm_;
	OrdinalType const numGlobalElements_;
	OrdinalType const numMyElements_;
	OrdinalType const indexBase_;
	OrdinalType const minLID_;
	OrdinalType const maxLID_;
	OrdinalType const minMyGID_;
	OrdinalType const maxMyGID_;
	OrdinalType const minAllGID_;
	OrdinalType const maxAllGID_;
	bool const contiguous_;
	bool const global_;
  map<OrdinalType, OrdinalType> lgMap_;
  map<OrdinalType, OrdinalType> const glMap_;
	OrdinalType* myGlobalElements_;
	Directory<OrdinalType>* Directory_;

 private:
	bool checkGlobalness() {
		bool global = false;
		if(Comm_->getNumImages() > 1) {
			int localRep = 0;
			int allLocalRep;
			if(numGlobalElements_ == numMyElements_)
				localRep = 1;
			Comm_->minAll(&localRep, &allLocalRep, 1);
			if(allLocalRep != 1)
				global = true;
		}
		return(global);
	}

	//! Copy constructor (declared but not defined, do not use)
	ElementSpaceData(ElementSpaceData<OrdinalType> const& Source);
	//! Assignment operator (declared but not defined, do not use)
	ElementSpaceData<OrdinalType>& operator = (ElementSpaceData<OrdinalType> const& Source);

}; // class ElementSpaceData

} // namespace Tpetra

#endif // TPETRA_ELEMENTSPACEDATA_HPP
