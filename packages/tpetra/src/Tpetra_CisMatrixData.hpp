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

#ifndef TPETRA_CISMATRIXDATA_HPP
#define TPETRA_CISMATRIXDATA_HPP

namespace Tpetra {

template<typename OrdinalType, typename ScalarType>
class CisMatrixData : public Object {
	friend class CisMatrix<OrdinalType, ScalarType>;
 public:
    // default constructor
	CisMatrixData(VectorSpace<OrdinalType, ScalarType> primaryDist, bool rowOriented)
	  : Object("Tpetra::CisMatrixData")
    , primary_(primaryDist)
    , secondary_(primaryDist)
    , domain_(primaryDist)
    , range_(primaryDist)
    , haveSecondary_(false) // initially assume we don't have any of these
    , haveDomain_(false)    // then set to true as we assign them
    , haveRange_(false)
    , haveRow_(rowOriented)  // if rowOriented is true, haveRow is true and haveCol is false
    , haveCol_(!rowOriented) // if rowOriented is false, haveRow is false and haveCol is true
	  , rowOriented_(rowOriented)
    , fillCompleted_(false)
    , numMyNonzeros_(Teuchos::OrdinalTraits<OrdinalType>::zero())
  {
	  platform_ = Teuchos::rcp(primaryDist.platform().clone());
	  comm_ = Teuchos::rcp(primaryDist.platform().createScalarComm());
    ordinalComm_ = Teuchos::rcp(primaryDist.platform().createOrdinalComm());
	}

	CisMatrixData(VectorSpace<OrdinalType, ScalarType> primaryDist, VectorSpace<OrdinalType, ScalarType> secondaryDist, bool rowOriented)
	  : Object("Tpetra::CisMatrixData")
    , primary_(primaryDist)
    , secondary_(secondaryDist)
    , domain_(primaryDist)
    , range_(primaryDist)
    , haveSecondary_(true)
    , haveDomain_(false) // initially assume we don't have any of these
    , haveRange_(false)  // then set to true as we assign them
    , haveRow_(true)
    , haveCol_(true)
	  , rowOriented_(rowOriented)
    , fillCompleted_(false)
    , numMyNonzeros_(Teuchos::OrdinalTraits<OrdinalType>::zero())
  {
	  platform_ = Teuchos::rcp(primaryDist.platform().clone());
	  comm_ = Teuchos::rcp(primaryDist.platform().createScalarComm());
    ordinalComm_ = Teuchos::rcp(primaryDist.platform().createOrdinalComm());
	}

    // destructor. no heap-data, so no need to override
	~CisMatrixData() {}

 protected:
  // map of maps that stores indices and values
	map< OrdinalType, map<OrdinalType, ScalarType> > indicesAndValues_;
  // vectors used by fillComplete to create HbMatrix
  std::vector<OrdinalType> pntr_;
  std::vector<OrdinalType> indx_;
  std::vector<ScalarType> values_;
  Kokkos::HbMatrix<OrdinalType, ScalarType> HbMatrix_;
  Kokkos::BaseSparseMultiply<OrdinalType, ScalarType> axy_;
  Kokkos::DenseVector<OrdinalType, ScalarType> kx_;
  Kokkos::DenseVector<OrdinalType, ScalarType> ky_;
  

	// VectorSpaces
	VectorSpace<OrdinalType, ScalarType> const primary_;
	VectorSpace<OrdinalType, ScalarType> secondary_;
	VectorSpace<OrdinalType, ScalarType> domain_;
	VectorSpace<OrdinalType, ScalarType> range_;

  // booleans
  bool haveSecondary_;
  bool haveDomain_;
  bool haveRange_;
  bool haveRow_;
  bool haveCol_;
	bool rowOriented_;
  bool fillCompleted_;

  // other variables used prior to fillComplete
  OrdinalType numMyNonzeros_;

  // Platform & Comm
	Teuchos::RefCountPtr<Platform<OrdinalType, ScalarType> const> platform_;
	Teuchos::RefCountPtr<Comm<ScalarType, OrdinalType> const> comm_;
  Teuchos::RefCountPtr<Comm<OrdinalType, OrdinalType> const> ordinalComm_;

 private:
	//! Copy constructor (declared but not defined, do not use)
	CisMatrixData(CisMatrixData<OrdinalType, ScalarType> const& rhs);
	//! Assignment operator (declared but not defined, do not use)
	CisMatrixData<OrdinalType, ScalarType>& operator = (CisMatrixData<OrdinalType, ScalarType> const& rhs);

}; // class CisMatrixData

} // namespace Tpetra

#endif // TPETRA_CISMATRIXDATA_HPP
