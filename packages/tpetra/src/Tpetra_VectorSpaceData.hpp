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

#ifndef _TPETRA_VECTORSPACEDATA_HPP_
#define _TPETRA_VECTORSPACEDATA_HPP_

namespace Tpetra {

template<typename OrdinalType, typename ScalarType>
class VectorSpaceData : public Object {
	friend class VectorSpace<OrdinalType, ScalarType>;
 public:
    // default constructor
	VectorSpaceData(bool blockspace, OrdinalType indexBase, OrdinalType numMyEntries, OrdinalType numGlobalEntries, Platform<OrdinalType, ScalarType> const& platform) 
    : Object("Tpetra::VectorSpaceData")
    , blockspace_(blockspace)
    , indexBase_(indexBase)
    , numMyEntries_(numMyEntries)
    , numGlobalEntries_(numGlobalEntries)
    , Platform_()
    , Comm_()
    {
        Platform_ = Teuchos::rcp(platform.clone());
        Comm_ = Teuchos::rcp(platform.createScalarComm());
    };

    // destructor. no heap-data, so no need to override
	~VectorSpaceData() {};

 protected:
    bool blockspace_;
    OrdinalType const indexBase_;
	OrdinalType const numMyEntries_;
	OrdinalType const numGlobalEntries_;
	Teuchos::RefCountPtr< Platform<OrdinalType, ScalarType> const > Platform_;
	Teuchos::RefCountPtr< Comm<ScalarType, OrdinalType> const > Comm_; // Comm is <ST, OT> because ST represents PT

 private:
	//! Copy constructor (declared but not defined, do not use)
	VectorSpaceData(VectorSpaceData<OrdinalType, ScalarType> const& Source);
	//! Assignment operator (declared but not defined, do not use)
	VectorSpaceData<OrdinalType, ScalarType>& operator = (VectorSpaceData<OrdinalType, ScalarType> const& Source);

}; // class VectorSpaceData

} // namespace Tpetra

#endif // _TPETRA_VECTORSPACEDATA_HPP_
