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
#include "Tpetra_CombineMode.hpp"
#include "Tpetra_VectorSpace.hpp"

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
  
	//! Constructor with variable number of indices per row/column.
	SparseMatrix(VectorSpace<OrdinalType, ScalarType> const& rowMap, OrdinalType* numEntriesPerIndex)
		: Object("Tpetra::SparseMatrix")
		, SparseMatrixData_()
	{
		SparseMatrixData_ = Teuchos::rcp(new SparseMatrixData<OrdinalType, ScalarType>());
	}
  
  //! Constructor with fixed number of indices per row/column.
	SparseMatrix(VectorSpace<OrdinalType, ScalarType> const& rowMap, OrdinalType numEntriesPerIndex)
		: Object("Tpetra::SparseMatrix")
		, SparseMatrixData_()
  {
    SparseMatrixData_ = Teuchos::rcp(new SparseMatrixData<OrdinalType, ScalarType>());
  }
  
  //! Constructor with variable number of indices per row/column.
	SparseMatrix(VectorSpace<OrdinalType, ScalarType> const& rowMap, 
               VectorSpace<OrdinalType, ScalarType> const& colMap, 
               OrdinalType* numEntriesPerIndex)
		: Object("Tpetra::SparseMatrix")
		, SparseMatrixData_()
  {
    SparseMatrixData_ = Teuchos::rcp(new SparseMatrixData<OrdinalType, ScalarType>());
  }
  
  //! Constructor with fixed number of indices per row/column.
	SparseMatrix(VectorSpace<OrdinalType, ScalarType> const& rowMap, 
               VectorSpace<OrdinalType, ScalarType> const& colMap, 
               OrdinalType numEntriesPerIndex)
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
  
  //@{ \name Post-Construction Modification Routines
  
  //! Set all matrix entries equal to scalarThis
  void setAllToScalar(ScalarType scalarThis);
  
  //! Scale the current values of a matrix, \e this = scalarThis*\e this.
  void scale(ScalarType scalarThis);
  
  //! Submit entries using global IDs. Behavoir is defined by the CombineMode passed in.
  void submitGlobalValues(CombineMode CM, OrdinalType numEntries, ScalarType* values, OrdinalType* indices);
  
  //! Submit entries using local IDs. Behavoir is defined by the CombineMode passed in.
  void submitLocalValues(CombineMode CM, OrdinalType numEntries, ScalarType* values, OrdinalType* indices);
  
  //! Signal that data entry is complete. Matrix data is converted into a more optimized form.
  void fillComplete();
  
  //! Signal that data entry is complete. Matrix data is converted into a more optimized form.
  void fillComplete(VectorSpace<OrdinalType, ScalarType> const& domainMap, VectorSpace<OrdinalType, ScalarType> const& rangeMap);
  
  //@}
  
  //@{ \name Computational Methods
  
  //! Computes the matrix-vector multiplication y = Ax
  void apply(Vector<OrdinalType, ScalarType> const& x, Vector<OrdinalType, ScalarType>& y) const;
  
  //! Returns the global one-norm of the matrix
  ScalarType normOne() const;
  
  //! Returns the global infinity norm of the matrix
  ScalarType normInf() const;
  
  //@}
  
  //@{ \name Attribute Access Methods
  
	//! Returns the number of nonzero entries in the global matrix.
	OrdinalType getNumGlobalNonzeros() const;
  
  //! Returns the number of nonzero entries in the calling image's portion of the matrix.
	OrdinalType getNumMyNonzeros() const;
	
	//! Returns the number of global matrix rows.
	OrdinalType getNumGlobalRows() const;
	
	//! Returns the number of global matrix columns.
	OrdinalType getNumGlobalCols();
  
  //! Returns the number of matrix rows owned by the calling image.
	OrdinalType getNumMyRows() const;
	
	//! Returns the number of matrix columns owned by the calling image.
	OrdinalType getNumMyCols() const;
	
	//! Returns the number of global nonzero diagonal entries, based on global row/column index comparisons.
	OrdinalType getNumGlobalDiagonals() const;
	
	//! Returns the number of local nonzero diagonal entries, based on global row/column index comparisons.
	OrdinalType getNumMyDiagonals() const;
	
	//! Returns the current number of nonzero entries in specified global index on this image.
	OrdinalType getNumGlobalEntries(OrdinalType index) const;
  
  //! Returns the current number of nonzero entries in specified local index on this image.
	OrdinalType getNumMyEntries(OrdinalType index) const;
	
	//! Returns the allocated number of nonzero entries in specified global index on this image.
	OrdinalType getNumAllocatedGlobalEntries(OrdinalType index);
  
  //! Returns the allocated number of nonzero entries in specified local index on this image.
	OrdinalType getNumAllocatedMyEntries(OrdinalType index) const;
	
	//! Returns the maximum number of nonzero entries across all rows on all images.
	OrdinalType getGlobalMaxNumEntries() const;
  
	//! Returns the maximum number of nonzero entries across all rows on this image.
	OrdinalType getMyMaxNumEntries() const;
	
	//! Returns the index base for global indices for this matrix.
	OrdinalType getIndexBase() const;
  
  //! Returns false if this matrix shares data with another SparseMatrix instance (or instances).
  bool isSoleOwner() const;
  
  //! Returns the VectorSpace that describes the row distribution in this matrix.
  VectorSpace<OrdinalType, ScalarType> const& getRowMap() const;
  
  //! Returns the VectorSpace that describes the column distribution in this matrix.
  VectorSpace<OrdinalType, ScalarType> const& getColMap() const;
  
  //! Returns the VectorSpace associated with the domain of this matrix.
  VectorSpace<OrdinalType, ScalarType> const& getDomainMap() const;
  
  //! Returns the VectorSpace associated with the range of this matrix.
  VectorSpace<OrdinalType, ScalarType> const& getRangeMap() const;
  
  //! Returns the Platform object used by this matrix
  Platform<OrdinalType, ScalarType> const& platform() const;
  
  //@}
  
  //@{ \name I/O Methods
  
  // Print method, used by the overloaded << operator
  void print(ostream& os) const {
    os << "SparseMatrix print function not implemented yet.";
  }
  
  //@}
  
  //@{ \name Miscellaneous Methods
  
  //! Returns true if this is a row-oriented matrix, false if it this is a column-oriented matrix.
  bool isRowOriented() const;
  
  //! Set to true to specify that this is a row-oriented matrix, set to false if this is a column-oriented matrix.
  void setRowOriented(bool rowOrColumn);
  
  //! Inlined bracket operator, non-const version.
  ScalarType* operator[] (OrdinalType index);
  
  //! Inlined bracket operator, const version.
  ScalarType const* operator[] (OrdinalType index) const;
  
	//! Assignment operator (declared but not defined, do not use)
	SparseMatrix<OrdinalType, ScalarType>& operator = (SparseMatrix<OrdinalType, ScalarType> const& rhs);
  
  //@}

	
private:

	Teuchos::RefCountPtr< SparseMatrixData<OrdinalType, ScalarType> > SparseMatrixData_;

}; // SparseMatrix class

} // Tpetra namespace

#include "Tpetra_SparseMatrixData.hpp"

#endif // _TPETRA_SPARSEMATRIX_HPP_
