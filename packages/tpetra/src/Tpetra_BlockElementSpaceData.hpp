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

#ifndef TPETRA_BLOCKELEMENTSPACEDATA_HPP
#define TPETRA_BLOCKELEMENTSPACEDATA_HPP

namespace Tpetra {

template<typename OrdinalType>
class BlockElementSpaceData : public Object {
	friend class BlockElementSpace<OrdinalType>;
 public:
	BlockElementSpaceData(ElementSpace<OrdinalType> const& ElementSpace, 
												bool const constantSize, 
												OrdinalType const elementSize,
												OrdinalType const numMyPoints,
												OrdinalType const numGlobalPoints,
												OrdinalType const minMySize,
												OrdinalType const maxMySize,
												OrdinalType const minGlobalSize,
												OrdinalType const maxGlobalSize,
												OrdinalType const* elementSizeList) 
		: Object("Tpetra::BlockElementSpaceData")
		, ElementSpace_(&ElementSpace) 
		, constantSize_(constantSize)
		, elementSize_(elementSize)
		, numMyPoints_(numMyPoints)
		, numGlobalPoints_(numGlobalPoints)
		, minMySize_(minMySize)
		, maxMySize_(maxMySize)
		, minGlobalSize_(minGlobalSize)
		, maxGlobalSize_(maxGlobalSize)
		, elementSizeList_(elementSizeList)
		, pointToElementList_(0)
		, firstPointList_(0) 
		{};

	~BlockElementSpaceData() {
		if(firstPointList_ != 0) {
			delete[] firstPointList_;
			firstPointList_ = 0;
		}
		if(pointToElementList_ != 0) {
			delete[] pointToElementList_;
			pointToElementList_ = 0;
		}
		if(elementSizeList_ != 0) {
			delete[] elementSizeList_;
			elementSizeList_ = 0;
		}
	};

 protected:
	ElementSpace<OrdinalType> const* ElementSpace_;
	bool const constantSize_;
	OrdinalType const elementSize_;
	OrdinalType const numMyPoints_;
	OrdinalType const numGlobalPoints_;
	OrdinalType const minMySize_;
	OrdinalType const maxMySize_;
	OrdinalType const minGlobalSize_;
	OrdinalType const maxGlobalSize_;
	OrdinalType const* elementSizeList_;
	OrdinalType const* pointToElementList_;
	OrdinalType const* firstPointList_;
	
 private:
	//! Copy constructor (declared but not defined, do not use)
	BlockElementSpaceData(BlockElementSpaceData<OrdinalType> const& Source);
	//! Assignment operator (declared but not defined, do not use)
	BlockElementSpaceData<OrdinalType>& operator = (BlockElementSpaceData<OrdinalType> const& Source);

}; // class BlockElementSpaceData

} // namespace Tpetra

#endif // TPETRA_BLOCKELEMENTSPACEDATA_HPP
