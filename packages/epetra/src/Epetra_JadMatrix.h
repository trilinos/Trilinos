
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

#ifndef EPETRA_JADMATRIX_H
#define EPETRA_JADMATRIX_H

#include "Epetra_RowMatrix.h"
#include "Epetra_Object.h"
#include "Epetra_CompObject.h"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_IntSerialDenseVector.h"


class Epetra_Vector;
class Epetra_MultiVector;
class Epetra_Import;
class Epetra_Export;

//! Epetra_JadMatrix: A class for constructing matrix objects optimized for common kernels.

/*! The Epetra_JadMatrix class takes an existing Epetra_RowMatrix ojbect, analyzes it and 
    builds a jagged diagonal equivalent of it. Once constructed, it is also possible to 
    update the values of the matrix with values from another Epetra_RowMatrix that has
    the identical structure.
    
*/    

class Epetra_JadMatrix: public Epetra_CompObject, public Epetra_Object, public virtual Epetra_RowMatrix  {
      
 public:

  //@{ \name Constructors/Destructor.
  //! Epetra_JadMatrix constuctor.
  /* The constructor for this class requires a fully constructed instance of an Epetra_RowMatrix
     object.
     \param Matrix (In) An existing Epetra_RowMatrix.
     \pre Matrix must have Matrix.Filled()==true.
  */
  Epetra_JadMatrix(const Epetra_RowMatrix & Matrix);

  //! Epetra_JadMatrix Destructor
  virtual ~Epetra_JadMatrix();
  //@}
  
  //@{ \name Post-construction modifications.
  //! Update values using a matrix with identical structure.
  /* Updates the values only using a matrix that has exactly the same structure as
     the matrix used to construct this Epetra_JadMatrix object.  Once the constructor
     is called, the Matrix argument is no longer needed.
     \param Matrix (In) An existing Epetra_RowMatrix with \e identical structure to 
            the matrix used to create this Epetra_JadMatrix.
     \param CheckStructure (In) Optional argument, by default is false.  If set to true, 
            the method will check to see if the structure of Matrix is compatible with
	    the structure of matrix used to create this Epetra_JadMatrix.  Performing
	    this check has signficant overhead, so it should only be turned on for debugging.
     \pre Matrix must have Matrix.Filled()==true.
  */ 
  int UpdateValues(const Epetra_RowMatrix & Matrix, bool CheckStructure = false);
  //@}
  
  //@{ \name Extraction methods.

    //! Returns a copy of the specified local row in user-provided arrays.
    /*! 
    \param In
           MyRow - Local row to extract.
    \param In
	   Length - Length of Values and Indices.
    \param Out
	   NumEntries - Number of nonzero entries extracted.
    \param Out
	   Values - Extracted values for this row.
    \param Out
	   Indices - Extracted global column indices for the corresponding values.
	  
    \return Integer error code, set to 0 if successful, set to -1 if MyRow not valid, -2 if Length is too short (NumEntries will have required length).
  */
    int ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const;

    //! Returns a reference to the ith entry in the matrix, along with its row and column index.
    /*! 
    \param In
           CurEntry - Local entry to extract.
    \param Out
	   Value - Extracted reference to current values.
    \param Out
	   RowIndex - Row index for current entry.
    \param Out
	   ColIndex - Column index for current entry.
	  
    \return Integer error code, set to 0 if successful, set to -1 if CurEntry not valid.
  */
    inline int ExtractMyEntryView(int CurEntry, double *Value, int & RowIndex, int & ColIndex) { 
      if (CurEntry>=NumMyNonzeros_) EPETRA_CHK_ERR(-1); 
      Value = &(Values_[0])+CurEntry;
      ColIndex = Indices_[CurEntry];
      for (int j=0; j<NumJaggedDiagonals_; j++) if (CurEntry<IndexOffset_[j+1]) RowIndex = InvRowPerm_[IndexOffset_[j]-CurEntry];
      return(0);
    }

    //! Returns a const reference to the ith entry in the matrix, along with its row and column index.
    /*! 
    \param In
           CurEntry - Local entry to extract.
    \param Out
	   Value - Extracted reference to current values.
    \param Out
	   RowIndex - Row index for current entry.
    \param Out
	   ColIndex - Column index for current entry.
	  
    \return Integer error code, set to 0 if successful, set to -1 if CurEntry not valid.
  */
    inline int ExtractMyEntryView(int CurEntry, double const * Value, int & RowIndex, int & ColIndex) const { 
      if (CurEntry>=NumMyNonzeros_) EPETRA_CHK_ERR(-1); 
      Value = &Values_[0]+CurEntry;
      ColIndex = Indices_[CurEntry];
      for (int j=0; j<NumJaggedDiagonals_; j++) if (CurEntry<IndexOffset_[j+1]) RowIndex = InvRowPerm_[IndexOffset_[j]-CurEntry];
      return(0);
    }

    //! Returns a copy of the main diagonal in a user-provided vector.
    /*! 
    \param Out
	   Diagonal - Extracted main diagonal.

    \return Integer error code, set to 0 if successful.
  */
    int ExtractDiagonalCopy(Epetra_Vector & Diagonal) const;
    //@}

  //@{ \name Computational methods.

    //! Returns the result of a Epetra_JadMatrix multiplied by a Epetra_MultiVector X in Y.
    /*! 
    \param In
	   TransA -If true, multiply by the transpose of matrix, otherwise just use matrix.
    \param In
	   X - A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
	   Y -A Epetra_MultiVector of dimension NumVectorscontaining result.

    \return Integer error code, set to 0 if successful.
  */
    int Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    //! Returns the result of a Epetra_JadMatrix solve with a Epetra_MultiVector X in Y (not implemented).
    /*! 
    \param In
	   Upper -If true, solve Ux = y, otherwise solve Lx = y.
    \param In
	   Trans -If true, solve transpose problem.
    \param In
	   UnitDiagonal -If true, assume diagonal is unit (whether it's stored or not).
    \param In
	   X - A Epetra_MultiVector of dimension NumVectors to solve for.
    \param Out
	   Y -A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
  */
    int Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {return(-1);}

    //! Computes the sum of absolute values of the rows of the Epetra_JadMatrix, results returned in x.
    /*! The vector x will return such that x[i] will contain the inverse of sum of the absolute values of the 
        \e this matrix will be scaled such that A(i,j) = x(i)*A(i,j) where i denotes the global row number of A
        and j denotes the global column number of A.  Using the resulting vector from this function as input to LeftScale()
	will make the infinity norm of the resulting matrix exactly 1.
    \param Out
	   x -A Epetra_Vector containing the row sums of the \e this matrix. 
	   \warning It is assumed that the distribution of x is the same as the rows of \e this.

    \return Integer error code, set to 0 if successful.
  */
    int InvRowSums(Epetra_Vector& x) const;

    //! Scales the Epetra_JadMatrix on the left with a Epetra_Vector x.
    /*! The \e this matrix will be scaled such that A(i,j) = x(i)*A(i,j) where i denotes the row number of A
        and j denotes the column number of A.
    \param In
	   x -A Epetra_Vector to solve for.

    \return Integer error code, set to 0 if successful.
  */
    int LeftScale(const Epetra_Vector& x);

    //! Computes the sum of absolute values of the columns of the Epetra_JadMatrix, results returned in x.
    /*! The vector x will return such that x[j] will contain the inverse of sum of the absolute values of the 
        \e this matrix will be sca such that A(i,j) = x(j)*A(i,j) where i denotes the global row number of A
        and j denotes the global column number of A.  Using the resulting vector from this function as input to 
	RighttScale() will make the one norm of the resulting matrix exactly 1.
    \param Out
	   x -A Epetra_Vector containing the column sums of the \e this matrix. 
	   \warning It is assumed that the distribution of x is the same as the rows of \e this.

    \return Integer error code, set to 0 if successful.
  */
    int InvColSums(Epetra_Vector& x) const;

    //! Scales the Epetra_JadMatrix on the right with a Epetra_Vector x.
    /*! The \e this matrix will be scaled such that A(i,j) = x(j)*A(i,j) where i denotes the global row number of A
        and j denotes the global column number of A.
    \param In
	   x -The Epetra_Vector used for scaling \e this.

    \return Integer error code, set to 0 if successful.
  */
    int RightScale(const Epetra_Vector& x);
  //@}

  //@{ \name Matrix Properties Query Methods.


    //! If FillComplete() has been called, this query returns true, otherwise it returns false, presently always returns true.
    bool Filled() const {return(true);}

    //! If matrix is lower triangular, this query returns true, otherwise it returns false.
    bool LowerTriangular() const {return(LowerTriangular_);}

    //! If matrix is upper triangular, this query returns true, otherwise it returns false.
    bool UpperTriangular() const {return(UpperTriangular_);}

  //@}
  
  //@{ \name Atribute access functions

    //! Returns the infinity norm of the global matrix.
    /* Returns the quantity \f$ \| A \|_\infty\f$ such that
       \f[\| A \|_\infty = \max_{1\lei\lem} \sum_{j=1}^n |a_{ij}| \f].

     \warning This method is supported if and only if the Epetra_RowMatrix Object that was used to create this supports this method.

    */ 
    double NormInf() const;

    //! Returns the one norm of the global matrix.
    /* Returns the quantity \f$ \| A \|_1\f$ such that
       \f[\| A \|_1= \max_{1\lej\len} \sum_{i=1}^m |a_{ij}| \f].

     \warning This method is supported if and only if the Epetra_RowMatrix Object that was used to create this supports this method.

    */ 
    double NormOne() const;

    //! Returns the number of nonzero entries in the global matrix.
    int NumGlobalNonzeros() const {return(NumGlobalNonzeros_);}

    //! Returns the number of global matrix rows.
    int NumGlobalRows() const {return(OperatorRangeMap().NumGlobalPoints());}

    //! Returns the number of global matrix columns.
    int NumGlobalCols() const {return(OperatorDomainMap().NumGlobalPoints());}

    //! Returns the number of global nonzero diagonal entries.
    int NumGlobalDiagonals() const{return(OperatorDomainMap().NumGlobalPoints());}
    
    //! Returns the number of nonzero entries in the calling processor's portion of the matrix.
    int NumMyNonzeros() const {return(NumMyNonzeros_);}

    //! Returns the number of matrix rows owned by the calling processor.
    int NumMyRows() const {return(OperatorRangeMap().NumMyPoints());}

    //! Returns the number of matrix columns owned by the calling processor.
    int NumMyCols() const {return(RowMatrixColMap().NumMyPoints());}

    //! Returns the number of local nonzero diagonal entries.
    int NumMyDiagonals() const {return(OperatorRangeMap().NumMyPoints());}

    //! Returns the maximum number of nonzero entries across all rows on this processor.
    int MaxNumEntries() const {return(NumJaggedDiagonals_);}

    //! Return the current number of values stored for the specified local row.
    /*! Similar to NumMyEntries() except NumEntries is returned as an argument
      and error checking is done on the input value MyRow.
      \param MyRow - (In) Local row.
      \param NumEntries - (Out) Number of nonzero values.
      
      \return Integer error code, set to 0 if successful, set to -1 if MyRow not valid.
      \pre None.
      \post Unchanged.
    */
     int NumMyRowEntries(int MyRow, int& NumEntries) const;

    //! Returns the Epetra_Map object associated with the domain of this operator.
    const Epetra_Map & OperatorDomainMap() const {return(OperatorDomainMap_);}

    //! Returns the Epetra_Map object associated with the range of this operator (same as domain).
    const Epetra_Map & OperatorRangeMap() const  {return(OperatorRangeMap_);}

    //! Implement the Epetra_SrcDistObjec::Map() function.
    const Epetra_BlockMap& Map() const {return(RowMatrixRowMap());}

    //! Returns the Row Map object needed for implementing Epetra_RowMatrix.
    const Epetra_Map & RowMatrixRowMap() const {return(RowMatrixRowMap_);}

    //! Returns the Column Map object needed for implementing Epetra_RowMatrix.
    const Epetra_Map & RowMatrixColMap() const {return(RowMatrixColMap_);}

    //! Returns the Epetra_Import object that contains the import operations for distributed operations.
    const Epetra_Import * RowMatrixImporter() const {return(Importer_);}

    //! Returns a pointer to the Epetra_Comm communicator associated with this matrix.
    const Epetra_Comm & Comm() const {return(*Comm_);}
  //@}
  
  
  //@{ \name I/O Methods.

  //! Print method
  virtual void Print(ostream & os) const;
  //@}

  //@{ \name Additional methods required to support the Epetra_RowMatrix interface.
    
    //! If set true, transpose of this operator will be applied.
    /*! This flag allows the transpose of the given operator to be used implicitly.  Setting this flag
        affects only the Apply() and ApplyInverse() methods.  If the implementation of this interface 
	does not support transpose use, this method should return a value of -1.
      
    \param In
	   UseTranspose -If true, multiply by the transpose of operator, otherwise just use operator.

    \return Always returns 0.
  */
  int SetUseTranspose(bool UseTranspose) {return(UseTranspose_ = UseTranspose);}

  //! Returns a character string describing the operator
  const char* Label() const {return(Epetra_Object::Label());}
  
  //! Returns the result of a Epetra_RowMatrix applied to a Epetra_MultiVector X in Y.
  /*! 
    \param X (In) - A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Y (Out) - A Epetra_MultiVector of dimension NumVectors containing result.
    
    \return Integer error code, set to 0 if successful.
  */
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
    return(Epetra_JadMatrix::Multiply(Epetra_JadMatrix::UseTranspose(), X, Y));}

    //! Returns the result of a Epetra_RowMatrix inverse applied to an Epetra_MultiVector X in Y.
    /*! 

    \param X (In) - A Epetra_MultiVector of dimension NumVectors to solve for.
    \param Y (Out) - A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code = -1.
    \warning This method is NOT supported.
  */
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {return(-1);}

  //! Returns true because this class can compute an Inf-norm.
  bool HasNormInf() const {return(true);}
  
  //! Returns the current UseTranspose setting.
  bool UseTranspose() const {return(UseTranspose_);}

  //@}
  
  //@{ \name Additional accessor methods.

  //! Returns the Epetra_Import object that contains the import operations for distributed operations, returns zero if none.
    /*! If RowMatrixColMap!=OperatorDomainMap, then this method returns a pointer to an Epetra_Import object that imports objects
        from an OperatorDomainMap layout to a RowMatrixColMap layout.  This operation is needed for sparse matrix-vector
	multiplication, y = Ax, to gather x elements for local multiplication operations.

	If RowMatrixColMap==OperatorDomainMap, then the pointer will be returned as 0.

    \return Raw pointer to importer.  This importer will be valid as long as the Epetra_RowMatrix object is valid.
  */
  const Epetra_Import* Importer() const {return(Importer_);}
  
  //! Returns the Epetra_Export object that contains the export operations for distributed operations, returns zero if none.
    /*! If RowMatrixRowMap!=OperatorRangeMap, then this method returns a pointer to an Epetra_Export object that exports objects
        from an RowMatrixRowMap layout to a OperatorRangeMap layout.  This operation is needed for sparse matrix-vector
	multiplication, y = Ax, to scatter-add y elements generated during local multiplication operations.

	If RowMatrixRowMap==OperatorRangeMap, then the pointer will be returned as 0.  For a typical Epetra_RowMatrix object,
	this pointer will be zero since it is often the case that RowMatrixRowMap==OperatorRangeMap.

    \return Raw pointer to exporter.  This exporter will be valid as long as the Epetra_RowMatrix object is valid.
  */
  const Epetra_Export* Exporter() const {return(Exporter_);}

  //@}

 protected:

  void UpdateImportVector(int NumVectors) const;
  void UpdateExportVector(int NumVectors) const;
  void GeneralMV(bool TransA, double * x, double * y) const;
  void GeneralMM(bool TransA, double ** X, int LDX, double ** Y, int LDY, int NumVectors) const;
  void GeneralMM3RHS(bool TransA, double ** X, int LDX, double ** Y, int LDY, int NumVectors) const;
  void GeneralMM2RHS(bool TransA, double * x, int ldx, double * y, int ldy) const;
  Epetra_Comm * Comm_;
  Epetra_Map OperatorDomainMap_;
  Epetra_Map OperatorRangeMap_;
  Epetra_Map RowMatrixRowMap_;
  Epetra_Map RowMatrixColMap_;
  
  int Allocate(const Epetra_RowMatrix & Matrix);
  int NumMyRows_;
  int NumMyCols_;
  int NumMyNonzeros_;
  int NumGlobalNonzeros_;
  Epetra_SerialDenseVector Values_;
  Epetra_IntSerialDenseVector Indices_;
  Epetra_IntSerialDenseVector IndexOffset_;
  Epetra_IntSerialDenseVector Profile_;
  Epetra_IntSerialDenseVector RowPerm_;
  Epetra_IntSerialDenseVector InvRowPerm_;

  bool UseTranspose_;
  bool HasNormInf_;
  bool LowerTriangular_;
  bool UpperTriangular_;
  int NumJaggedDiagonals_;
    

  mutable Epetra_MultiVector * ImportVector_;
  mutable Epetra_MultiVector * ExportVector_;
  Epetra_Import * Importer_;
  Epetra_Export * Exporter_;

};
#endif /* EPETRA_JADMATRIX_H */
