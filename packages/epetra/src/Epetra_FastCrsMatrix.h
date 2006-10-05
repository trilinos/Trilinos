
//@HEADER
/*
************************************************************************

              Epetra: Linear Algebra Services Package 
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

#ifndef EPETRA_FASTCRSMATRIX_H
#define EPETRA_FASTCRSMATRIX_H

#include "Epetra_Operator.h"
#include "Epetra_CrsMatrix.h"
class Epetra_Vector;
class Epetra_MultiVector;

//! Epetra_FastCrsOperator: A class for constructing matrix objects optimized for common kernels.

/*! The Epetra_FastCrsOperator class takes an existing Epetra_CrsMatrix ojbect, analyzes it and build
    upon it for the purposes of obtaining the best possible performance on basic operations.
    
*/    

class Epetra_FastCrsOperator: public Epetra_CompObject, public virtual Epetra_Operator  {
      
 public:

   //! @name Constructors/Destructor
  //@{ 
  //! Epetra_FastCrsOperator constuctor.
  Epetra_FastCrsOperator(const Epetra_CrsMatrix & Matrix, bool UseFloats = false);

  //! Epetra_FastCrsOperator Destructor
  virtual ~Epetra_FastCrsOperator();
  //@}
  
  //! @name Post-construction modifications
  //@{ 
  //! Update values using a matrix with identical structure.
  /* Updates the values only using a matrix that has exactly the same structure as
     the matrix used to construct this Epetra_FastCrsOperator object.
  */ 
  int UpdateValues(const Epetra_CrsMatrix & Matrix);
  //@}
  
  //! @name Additional methods required to support the Epetra_Operator interface
  //@{ 

    //! Returns a character string describing the operator
    char * Label() const {return(CrsMatrix_.Label());};
    
    //! If set true, transpose of this operator will be applied.
    /*! This flag allows the transpose of the given operator to be used implicitly.  Setting this flag
        affects only the Apply() and ApplyInverse() methods.  If the implementation of this interface 
	does not support transpose use, this method should return a value of -1.
      
    \param In
	   UseTranspose -If true, multiply by the transpose of operator, otherwise just use operator.

    \return Always returns 0.
  */
  int SetUseTranspose(bool UseTranspose) {return(CrsMatrix_.SetUseTranspose());};

  //! Returns the infinity norm of the global matrix.
  /* Returns the quantity \f$ \| A \|_\infty\f$ such that
     \f[\| A \|_\infty = \max_{1\lei\lem} \sum_{j=1}^n |a_{ij}| \f].
  */ 
  double NormInf() const {return(CrsMatrix_.NormInf());};
  
  //! Returns a pointer to the Epetra_Comm communicator associated with this matrix.
  const Epetra_Comm & Comm() const {return(CrsMatrix_.Comm());};
  
  //! Returns the result of a Epetra_Operator applied to a Epetra_MultiVector X in Y.
  /*! 
    \param X (In) - A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Y (Out) - A Epetra_MultiVector of dimension NumVectors containing result.
    
    \return Integer error code, set to 0 if successful.
  */
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    //! Returns the result of a Epetra_Operator inverse applied to an Epetra_MultiVector X in Y.
    /*! 

    \param X (In) - A Epetra_MultiVector of dimension NumVectors to solve for.
    \param Y (Out) - A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
  */
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! Returns true because this class can compute an Inf-norm.
  bool HasNormInf() const {return(CrsMatrix_.HasNormInf());};
  
  //! Returns the current UseTranspose setting.
  bool UseTranspose() const {return(CrsMatrix_.UseTranspose());};
  
  //! Returns the Epetra_Map object associated with the domain of this matrix operator.
  const Epetra_Map & OperatorDomainMap() const {return(CrsMatrix_.OperatorDomainMap());};
  
  //! Returns the Epetra_Map object associated with the range of this matrix operator.
  const Epetra_Map & OperatorRangeMap() const {return(CrsMatrix_.OperatorRangeMap());};

  //@}


 protected:
  int Allocate(bool UseFloats);
  const Epetra_CrsMatrix & CrsMatrix_;
  int NumMyRows_;
  int NumMyNonzeros_;
  double *Values_;
  float * FloatValues_;
  int * Indices_;
  unsigned short * ShortIndices_;

  bool UsingFloats_;
  bool UsingShorts_;
  bool ValuesAllocated_;
    

  mutable Epetra_MultiVector * ImportVector_;
  mutable Epetra_MultiVector * ExportVector_;

};
#endif /* EPETRA_FASTCRSMATRIX_H */
