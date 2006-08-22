
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

#ifndef EPETRA_JADOPERATOR_H
#define EPETRA_JADOPERATOR_H

#include "Epetra_Operator.h"
#include "Epetra_Object.h"
#include "Epetra_CompObject.h"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_IntSerialDenseVector.h"


class Epetra_Vector;
class Epetra_MultiVector;
class Epetra_RowMatrix;
class Epetra_Import;
class Epetra_Export;

//! Epetra_JadOperator: A class for constructing matrix objects optimized for common kernels.

/*! The Epetra_JadOperator class takes an existing Epetra_RowMatrix ojbect, analyzes it and 
    builds a jagged diagonal equivalent of it. Once constructed, it is also possible to 
    update the values of the matrix with values from another Epetra_RowMatrix that has
    the identical structure.
    
*/    

class Epetra_JadOperator: public Epetra_CompObject, public Epetra_Object, public virtual Epetra_Operator  {
      
 public:

  //@{ \name Constructors/Destructor.
  //! Epetra_JadOperator constuctor.
  /* The constructor for this class requires a fully constructed instance of an Epetra_RowMatrix
     object.
     \param Matrix (In) An existing Epetra_RowMatrix.
     \param UseFloats (In) Optional argument, by default is false.  If set to true, 
            matrix values will be stored as floats instead of doubles.
     \param UseShorts (In) Optional argument, by default is false.  If set to true and if the local indices are small enough, 
            matrix indices will be stored as unsigned shorts instead of ints.
	    \warning UseFloats and UseShorts are not supported in this version of the code.
     \pre Matrix must have Matrix.Filled()==true.
  */
  Epetra_JadOperator(const Epetra_RowMatrix & Matrix, bool UseFloats = false, bool UseShorts = false);

  //! Epetra_JadOperator Destructor
  virtual ~Epetra_JadOperator();
  //@}
  
  //@{ \name Post-construction modifications.
  //! Update values using a matrix with identical structure.
  /* Updates the values only using a matrix that has exactly the same structure as
     the matrix used to construct this Epetra_JadOperator object.  Once the constructor
     is called, the Matrix argument is no longer needed.
     \param Matrix (In) An existing Epetra_RowMatrix with \e identical structure to 
            the matrix used to create this Epetra_JadOperator.
     \param CheckStructure (In) Optional argument, by default is false.  If set to true, 
            the method will check to see if the structure of Matrix is compatible with
	    the structure of matrix used to create this Epetra_JadOperator.  Performing
	    this check has signficant overhead, so it should only be turned on for debugging.
     \pre Matrix must have Matrix.Filled()==true.
  */ 
  int UpdateValues(const Epetra_RowMatrix & Matrix, bool CheckStructure = false);
  //@}
  
  //@{ \name Additional methods required to support the Epetra_Operator interface.
    
    //! If set true, transpose of this operator will be applied.
    /*! This flag allows the transpose of the given operator to be used implicitly.  Setting this flag
        affects only the Apply() and ApplyInverse() methods.  If the implementation of this interface 
	does not support transpose use, this method should return a value of -1.
      
    \param In
	   UseTranspose -If true, multiply by the transpose of operator, otherwise just use operator.

    \return Always returns 0.
  */
  int SetUseTranspose(bool UseTranspose) {return(UseTranspose_ = UseTranspose);};

  //! Returns the infinity norm of the global matrix.
  /* Returns the quantity \f$ \| A \|_\infty\f$ such that
     \f[\| A \|_\infty = \max_{1\lei\lem} \sum_{j=1}^n |a_{ij}| \f].

     \warning This method is supported if and only if the Epetra_RowMatrix Object that was used to create this supports this method.
  */ 
  double NormInf() const {return(NormInf_);};
  
  //! Returns a character string describing the operator
  const char* Label() const {return(Epetra_Object::Label());}

  //! Returns the Epetra_Comm communicator associated with this matrix.
  const Epetra_Comm & Comm() const {return(*Comm_);};
  
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

    \return Integer error code = -1.
    \warning This method is NOT supported.
  */
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {return(-1);}

  //! Returns true because this class can compute an Inf-norm.
  bool HasNormInf() const {return(HasNormInf_);};
  
  //! Returns the current UseTranspose setting.
  bool UseTranspose() const {return(UseTranspose_);};
  
  //! Returns the Epetra_Map object associated with the domain of this matrix operator.
  const Epetra_Map & OperatorDomainMap() const {return(OperatorDomainMap_);};
  
  //! Returns the Epetra_Map object associated with the range of this matrix operator.
  const Epetra_Map & OperatorRangeMap() const {return(OperatorRangeMap_);};

  //@}
  
  //@{ \name Additional accessor methods.

  //! Returns the Epetra_Import object that contains the import operations for distributed operations.
  const Epetra_Import* Importer() const {return(Importer_);}
  
  //! Returns the Epetra_Export object that contains the export operations for distributed operations.
  const Epetra_Export* Exporter() const {return(Exporter_);}
  
  //! Returns the number of nonzero entries in the global matrix.
  /*
    Note that if maps are defined such that some nonzeros appear on
    multiple processors, then those nonzeros will be counted multiple
    times.
  */
  inline int NumGlobalNonzeros() const {return(NumGlobalNonzeros_);}

  //@}

  //@{ \name I/O Methods.

  //! Print method
  virtual void Print(ostream& os) const;
  //@}

 protected:

  void UpdateImportVector(int NumVectors) const;
  void UpdateExportVector(int NumVectors) const;
  void GeneralMV(bool TransA, double * x, double * y) const;
  void GeneralMM(bool TransA, double ** X, int LDX, double ** Y, int LDY, int NumVectors) const;
  void GeneralMM3RHS(bool TransA, double ** X, int LDX, double ** Y, int LDY, int NumVectors) const;
  void GeneralMM2RHS(bool TransA, double * x, int ldx, double * y, int ldy) const;
  double NormInf_;
  Epetra_Comm * Comm_;
  Epetra_Map OperatorDomainMap_;
  Epetra_Map OperatorRangeMap_;
  
  int Allocate(const Epetra_RowMatrix & Matrix, bool UseFloats);
  int NumMyRows_;
  int NumMyCols_;
  int NumMyNonzeros_;
  int NumGlobalNonzeros_;
  Epetra_SerialDenseVector Values_;
  float * FloatValues_;
  Epetra_IntSerialDenseVector Indices_;
  unsigned short * ShortIndices_;
  Epetra_IntSerialDenseVector IndexOffset_;
  Epetra_IntSerialDenseVector RowPerm_;

  bool UseTranspose_;
  bool HasNormInf_;
  bool UsingFloats_;
  bool UsingShorts_;
  int NumJaggedDiagonals_;
    

  mutable Epetra_MultiVector * ImportVector_;
  mutable Epetra_MultiVector * ExportVector_;
  Epetra_Import * Importer_;
  Epetra_Export * Exporter_;

};
#endif /* EPETRA_JADOPERATOR_H */
