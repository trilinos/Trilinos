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

#ifndef _TPETRA_SPARSEMATRIX_HPP_
#define _TPETRA_SPARSEMATRIX_HPP_

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_CompObject.hpp>
#include "Tpetra_Object.hpp"

namespace Tpetra {

// forward declaration of SparseMatrixData, needed to prevent circular inclusions
// actual #include statement at the end of this file
template<typename OrdinalType, typename ScalarType> class SparseMatrixData;

//! Tpetra::SparseMatrix
/*! 
*/

template<typename OrdinalType, typename ScalarType>
class SparseMatrix : public Object, public Teuchos::CompObject {

public:
  
	//@{ \name Constructor/Destructor Methods
  
	//! default constructor.
	SparseMatrix()
		: Object("Tpetra::SparseMatrix")
		, SparseMatrixData_()
	{
		SparseMatrixData_ = Teuchos::rcp(new SparseMatrixData<OrdinalType, ScalarType>());
	}
	
	//! copy constructor.
	SparseMatrix(SparseMatrix<OrdinalType, ScalarType> const& rhs)
		: Object(rhs.label())
		, SparseMatrixData_(rhs.SparseMatrixData_)
	{}

	//! destructor.
	~SparseMatrix() {}
  
	//@}
	
	//! Assignment operator (declared but not defined, do not use)
	SparseMatrix<OrdinalType, ScalarType>& operator = (SparseMatrix<OrdinalType, ScalarType> const& rhs);

	//@}
	
private:

	Teuchos::RefCountPtr< SparseMatrixData<OrdinalType, ScalarType> > SparseMatrixData_;

}; // SparseMatrix class

} // Tpetra namespace

#include "Tpetra_SparseMatrixData.hpp"

#endif // _TPETRA_SPARSEMATRIX_HPP_
