
//@HEADER
/*
************************************************************************

              Tpetra: Linear Algebra Services Package 
                Copyright (2001) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.
 
This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact Michael A. Heroux (maherou@sandia.gov) 

************************************************************************
*/
//@HEADER

#ifndef TPETRA_ROWMATRIXTRANSPOSER_HPP
#define TPETRA_ROWMATRIXTRANSPOSER_HPP
#include <Teuchos_RCP.hpp>
#include <Kokkos_DefaultNode.hpp>
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_ConfigDefs.hpp"

//! Tpetra_RowMatrixTransposer: A class for transposing an Tpetra_RowMatrix object.

namespace Tpetra {
/*! This class provides capabilities to construct a transpose matrix of an existing Tpetra_RowMatrix
	  object and (optionally) redistribute it across a parallel distributed memory machine.
*/

template <class Scalar, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
class RowMatrixTransposer {
    
  public:

    //! @name Constructors/destructors
  //@{ 
  //! Primary Tpetra_RowMatrixTransposer constructor.
  /*!
    \param Matrix (In) An existing Tpetra_RowMatrix object.  The Tpetra_RowMatrix, the LHS and RHS pointers
		       do not need to be defined before this constructor is called.

    \return Pointer to a Tpetra_RowMatrixTransposer object.

  */ 
  RowMatrixTransposer(const Teuchos::RCP<const RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > origMatrix);

  //! Tpetra_RowMatrixTransposer copy constructor.
  
//  RowMatrixTransposer(const RowMatrixTransposer& Source);
  
  //! Tpetra_RowMatrixTransposer destructor.
  
  virtual ~RowMatrixTransposer();
  //@}
  
  //! @name Forward transformation methods
  //@{ 
  
  //! Generate a new Tpetra_CrsMatrix as the transpose of an Tpetra_RowMatrix passed into the constructor.
  /*! Constructs a new Tpetra_CrsMatrix that is a copy of the Tpetra_RowMatrix passed in to the constructor.
		
		\param MakeDataContiguous (In) Causes the output matrix, LHS and RHS to be stored in a form compatible with
		       Fortran-style solvers.  The output matrix will be compatible with the Harwell-Boeing compressed
					 column format.  The RHS and LHS will be stored such that the last value in column j of the 
					 multivector is stored next to the first value in column j+1.
		\param TransposeRowMap (Optional/In) If this argument is defined, the transpose matrix will be distributed
		       using this map as the row map for the transpose.  If it is set to zero, the transpose matrix will use
					 the OrigMatrix->RowMatrixDomainMap as the row map.

		\return Integer error code, 0 if no errors.  Negative if some fatal error occured.
					 
  */
  int createTranspose(const OptimizeOption optimizeTranspose, Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > &transposeMatrix);

	
 private: 
	const Teuchos::RCP<const RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > origMatrix_;
	Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > transposeMatrix_;
	bool transposeCreated_;
	OptimizeOption optimizeTranspose_;	
};


}

#endif /* TPETRA_ROWMATRIXTRANSPOSER_DECL_HPP */
