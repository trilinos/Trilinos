
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

#ifndef _EPETRA_FASTCRSMATRIX_H_
#define _EPETRA_FASTCRSMATRIX_H_

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

  //@{ \name Constructors/Destructor.
  //! Epetra_FastCrsOperator constuctor.
  Epetra_FastCrsOperator(const Epetra_CrsMatrix & Matrix, bool UseFloats = false);

  //! Epetra_FastCrsOperator Destructor
  virtual ~Epetra_FastCrsOperator();
  //@}
  
  //@{ \name Post-construction modifications.
  //! Update values using a matrix with identical structure.
  /* Updates the values only using a matrix that has exactly the same structure as
     the matrix used to construct this Epetra_FastCrsOperator object.
  */ 
  int UpdateValues(const Epetra_CrsMatrix & Matrix);
  //@}
  
  //@{ \name Additional methods required to support the Epetra_Operator interface.

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
#endif /* _EPETRA_FASTCRSMATRIX_H_ */
