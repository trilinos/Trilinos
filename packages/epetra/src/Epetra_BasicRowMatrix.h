
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

#ifndef EPETRA_BASICROWMATRIX_H
#define EPETRA_BASICROWMATRIX_H

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

//! Epetra_BasicRowMatrix: A class for constructing matrix objects optimized for common kernels.

/*! The Epetra_BasicRowMatrix class takes an existing Epetra_RowMatrix ojbect, analyzes it and 
    builds a jagged diagonal equivalent of it. Once constructed, it is also possible to 
    update the values of the matrix with values from another Epetra_RowMatrix that has
    the identical structure.
    
*/    

class Epetra_BasicRowMatrix: public Epetra_CompObject, public Epetra_Object, public virtual Epetra_RowMatrix  {
      
 public:

  //@{ \name Constructors/Destructor.
  //! Epetra_BasicRowMatrix constuctor.
  /* This constructor takes a row and column map.  On each processor these maps describe the global rows and columns, resp, 
     that the processor will care about.  Note that these maps do not have to be one-to-one.  In other words, a row ID can appear
     on more than one processor, as can a column ID.
     \param RowMap (In) An Epetra_Map containing on each processor a list of GIDs of rows that the processor cares about.
     \param ColMap (In) An Epetra_Map containing on each processor a list of GIDs of columns that the processor cares about.

     In this constructor, the domain and range maps are assumed to be the same as the row map.  Note that this requires that 
     the global matrix be square.  If the matrix is not square, or the domain vectors or range vectors do not have the same layout
     as the rows, then the second constructor should be called.
  */
  Epetra_BasicRowMatrix(const Epetra_Map & RowMap, const Epetra_Map & ColMap);

  //! Epetra_BasicRowMatrix constuctor.
  /* This constructor takes a row, column, domain and range map.  On each processor these maps describe the global rows, columns, domain
     and range, resp, that the processor will care about.  The domain and range maps must be one-to-one, but note that the row and column
     maps do not have to be one-to-one.  In other words, a row ID can appear
     on more than one processor, as can a column ID.
     \param RowMap (In) An Epetra_Map containing on each processor a list of GIDs of rows that the processor cares about.
     \param ColMap (In) An Epetra_Map containing on each processor a list of GIDs of columns that the processor cares about.
     \param DomainMap (In) An Epetra_Map describing the distribution of domain vectors and multivectors.
     \param RangeMap (In) An Epetra_Map describing the distribution of range vectors and multivectors.

  */
  Epetra_BasicRowMatrix(const Epetra_Map & RowMap, const Epetra_Map & ColMap, 
			const Epetra_Map & DomainMap, const Epetra_Map & RangeMap);

  //! Epetra_BasicRowMatrix Destructor
  virtual ~Epetra_BasicRowMatrix();
  //@}
  
  //@{ \name Post-construction modifications.
  //! Update the constants associated with the structure of the matrix.
  /* Several constants are pre-computed to save excess computations.  However, if the structure of the
     problem changes, specifically if the nonzero count in any given row changes, then this function should be called
     to update these constants.
  */ 
  virtual int ComputeStructureConstants();
  //! Update the constants associated with the values of the matrix.
  /* Several numeric constants are pre-computed to save excess computations.  However, if the values of the
     problem change, then this function should be called to update these constants.
  */ 
  virtual int ComputeNumericConstants();
  //@}
  
  //@{ \name User-required implementation methods.

    //! Returns a copy of the specified local row in user-provided arrays.
    /*! 
    \param MyRow (In) - Local row to extract.
    \param Length (In) - Length of Values and Indices.
    \param NumEntries (Out) - Number of nonzero entries extracted.
    \param Values (Out) - Extracted values for this row.
    \param Indices (Out) - Extracted global column indices for the corresponding values.
	  
    \return Integer error code, set to 0 if successful, set to -1 if MyRow not valid, -2 if Length is too short (NumEntries will have required length).
  */
  virtual int ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const = 0;

    //! Returns a reference to the ith entry in the matrix, along with its row and column index.
    /*! 
    \param CurEntry (In) - Local entry to extract.
    \param Value (Out) - Extracted reference to current values.
    \param RowIndex (Out) - Row index for current entry.
    \param ColIndex (Out) - Column index for current entry.
	  
    \return Integer error code, set to 0 if successful, set to -1 if CurEntry not valid.
  */
    virtual int ExtractMyEntryView(int CurEntry, double *Value, int & RowIndex, int & ColIndex) = 0;

    //! Returns a const reference to the ith entry in the matrix, along with its row and column index.
    /*! 
    \param CurEntry (In) - Local entry to extract.
    \param Value (Out) - Extracted reference to current values.
    \param RowIndex (Out) - Row index for current entry.
    \param ColIndex (Out) - Column index for current entry.
	  
    \return Integer error code, set to 0 if successful, set to -1 if CurEntry not valid.
  */
    virtual int ExtractMyEntryView(int CurEntry, double const * Value, int & RowIndex, int & ColIndex) const = 0;

    //! Return the current number of values stored for the specified local row.
    /*! Similar to NumMyEntries() except NumEntries is returned as an argument
      and error checking is done on the input value MyRow.
      \param MyRow (In) - Local row.
      \param NumEntries (Out) - Number of nonzero values.
      
      \return Integer error code, set to 0 if successful, set to -1 if MyRow not valid.
    */
    virtual int NumMyRowEntries(int MyRow, int & NumEntries) const = 0;
    //@}

  //@{ \name Computational methods.

    //! Returns the result of a Epetra_BasicRowMatrix multiplied by a Epetra_MultiVector X in Y.
    /*! 
    \param TransA (In) - If true, multiply by the transpose of matrix, otherwise just use matrix.
    \param X (Out) - An Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Y (Out) - An Epetra_MultiVector of dimension NumVectorscontaining result.

    \return Integer error code, set to 0 if successful.
  */
    virtual int Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    //! Returns the result of a Epetra_BasicRowMatrix solve with a Epetra_MultiVector X in Y (not implemented).
    /*! 
    \param Upper (In) - If true, solve Ux = y, otherwise solve Lx = y.
    \param Trans (In) - If true, solve transpose problem.
    \param UnitDiagonal (In) - If true, assume diagonal is unit (whether it's stored or not).
    \param X (In) - An Epetra_MultiVector of dimension NumVectors to solve for.
    \param Y (Out) - An Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
  */
    virtual int Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {return(-1);}

    //! Returns a copy of the main diagonal in a user-provided vector.
    /*! 
    \param Diagonal (Out) - Extracted main diagonal.

    \return Integer error code, set to 0 if successful.
  */
    virtual int ExtractDiagonalCopy(Epetra_Vector & Diagonal) const;

    //! Computes the sum of absolute values of the rows of the Epetra_BasicRowMatrix, results returned in x.
    /*! The vector x will return such that x[i] will contain the inverse of sum of the absolute values of the 
        \e this matrix will be scaled such that A(i,j) = x(i)*A(i,j) where i denotes the global row number of A
        and j denotes the global column number of A.  Using the resulting vector from this function as input to LeftScale()
	will make the infinity norm of the resulting matrix exactly 1.
    \param x (Out) - An Epetra_Vector containing the row sums of the \e this matrix. 
	   \warning It is assumed that the distribution of x is the same as the rows of \e this.

    \return Integer error code, set to 0 if successful.
  */
    virtual int InvRowSums(Epetra_Vector& x) const;

    //! Scales the Epetra_BasicRowMatrix on the left with a Epetra_Vector x.
    /*! The \e this matrix will be scaled such that A(i,j) = x(i)*A(i,j) where i denotes the row number of A
        and j denotes the column number of A.
    \param x (In) - An Epetra_Vector to solve for.

    \return Integer error code, set to 0 if successful.
  */
    virtual int LeftScale(const Epetra_Vector& x);

    //! Computes the sum of absolute values of the columns of the Epetra_BasicRowMatrix, results returned in x.
    /*! The vector x will return such that x[j] will contain the inverse of sum of the absolute values of the 
        \e this matrix will be sca such that A(i,j) = x(j)*A(i,j) where i denotes the global row number of A
        and j denotes the global column number of A.  Using the resulting vector from this function as input to 
	RighttScale() will make the one norm of the resulting matrix exactly 1.
    \param x (Out) - An Epetra_Vector containing the column sums of the \e this matrix. 
	   \warning It is assumed that the distribution of x is the same as the rows of \e this.

    \return Integer error code, set to 0 if successful.
  */
    virtual int InvColSums(Epetra_Vector& x) const;

    //! Scales the Epetra_BasicRowMatrix on the right with a Epetra_Vector x.
    /*! The \e this matrix will be scaled such that A(i,j) = x(j)*A(i,j) where i denotes the global row number of A
        and j denotes the global column number of A.
    \param x (In) - The Epetra_Vector used for scaling \e this.

    \return Integer error code, set to 0 if successful.
  */
    virtual int RightScale(const Epetra_Vector& x);
  //@}

  //@{ \name Matrix Properties Query Methods.


    //! If FillComplete() has been called, this query returns true, otherwise it returns false, presently always returns true.
    virtual bool Filled() const {return(true);}

    //! If matrix is lower triangular, this query returns true, otherwise it returns false.
    bool LowerTriangular() const {return(LowerTriangular_);}

    //! If matrix is upper triangular, this query returns true, otherwise it returns false.
    virtual bool UpperTriangular() const {return(UpperTriangular_);}

  //@}
  
  //@{ \name Atribute access functions

    //! Returns the infinity norm of the global matrix.
    /* Returns the quantity \f$ \| A \|_\infty\f$ such that
       \f[\| A \|_\infty = \max_{1\lei\lem} \sum_{j=1}^n |a_{ij}| \f].

     \warning This method is supported if and only if the Epetra_RowMatrix Object that was used to create this supports this method.

    */ 
    virtual double NormInf() const{return(NormInf_);}

    //! Returns the one norm of the global matrix.
    /* Returns the quantity \f$ \| A \|_1\f$ such that
       \f[\| A \|_1= \max_{1\lej\len} \sum_{i=1}^m |a_{ij}| \f].

     \warning This method is supported if and only if the Epetra_RowMatrix Object that was used to create this supports this method.

    */ 
    virtual double NormOne() const{return(NormOne_);}

    //! Returns the number of nonzero entries in the global matrix.
    virtual int NumGlobalNonzeros() const{return(NumGlobalNonzeros_);}

    //! Returns the number of global matrix rows.
    virtual int NumGlobalRows() const {return(OperatorRangeMap().NumGlobalPoints());}

    //! Returns the number of global matrix columns.
    virtual int NumGlobalCols() const {return(OperatorDomainMap().NumGlobalPoints());}

    //! Returns the number of global nonzero diagonal entries.
    virtual int NumGlobalDiagonals() const{return(OperatorDomainMap().NumGlobalPoints());}
    
    //! Returns the number of nonzero entries in the calling processor's portion of the matrix.
    virtual int NumMyNonzeros() const{return(NumMyNonzeros_);}

    //! Returns the number of matrix rows owned by the calling processor.
    virtual int NumMyRows() const {return(OperatorRangeMap().NumMyPoints());}

    //! Returns the number of matrix columns owned by the calling processor.
    virtual int NumMyCols() const {return(RowMatrixColMap().NumMyPoints());}

    //! Returns the number of local nonzero diagonal entries.
    virtual int NumMyDiagonals() const {return(OperatorRangeMap().NumMyPoints());}

    //! Returns the maximum number of nonzero entries across all rows on this processor.
    virtual int MaxNumEntries() const{return(MaxNumEntries_);}

    //! Returns the Epetra_Map object associated with the domain of this operator.
    virtual const Epetra_Map & OperatorDomainMap() const {return(OperatorDomainMap_);}

    //! Returns the Epetra_Map object associated with the range of this operator (same as domain).
    virtual const Epetra_Map & OperatorRangeMap() const  {return(OperatorRangeMap_);}

    //! Implement the Epetra_SrcDistObjec::Map() function.
    virtual const Epetra_BlockMap& Map() const {return(RowMatrixRowMap());}

    //! Returns the Row Map object needed for implementing Epetra_RowMatrix.
    virtual const Epetra_Map & RowMatrixRowMap() const {return(RowMatrixRowMap_);}

    //! Returns the Column Map object needed for implementing Epetra_RowMatrix.
    virtual const Epetra_Map & RowMatrixColMap() const {return(RowMatrixColMap_);}

    //! Returns the Epetra_Import object that contains the import operations for distributed operations.
    virtual const Epetra_Import * RowMatrixImporter() const {return(Importer_);}

    //! Returns a pointer to the Epetra_Comm communicator associated with this matrix.
    virtual const Epetra_Comm & Comm() const {return(*Comm_);}
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
      
    \param UseTranspose (In) - If true, multiply by the transpose of operator, otherwise just use operator.

    \return Always returns 0.
  */
  virtual int SetUseTranspose(bool UseTranspose) {return(UseTranspose_ = UseTranspose);}

  //! Returns a character string describing the operator
  virtual const char* Label() const {return(Epetra_Object::Label());}
  
  //! Returns the result of a Epetra_RowMatrix applied to a Epetra_MultiVector X in Y.
  /*! 
    \param X (In) - A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Y (Out) - A Epetra_MultiVector of dimension NumVectors containing result.
    
    \return Integer error code, set to 0 if successful.
  */
  virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
    return(Epetra_BasicRowMatrix::Multiply(Epetra_BasicRowMatrix::UseTranspose(), X, Y));}

    //! Returns the result of a Epetra_RowMatrix inverse applied to an Epetra_MultiVector X in Y.
    /*! 

    \param X (In) - A Epetra_MultiVector of dimension NumVectors to solve for.
    \param Y (Out) - A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code = -1.
    \warning This method is NOT supported.
  */
  virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {return(-1);}

  //! Returns true because this class can compute an Inf-norm.
  bool HasNormInf() const {return(true);}
  
  //! Returns the current UseTranspose setting.
  virtual bool UseTranspose() const {return(UseTranspose_);}

  //@}
  
  //@{ \name Additional accessor methods.

  //! Returns the Epetra_Import object that contains the import operations for distributed operations, returns zero if none.
    /*! If RowMatrixColMap!=OperatorDomainMap, then this method returns a pointer to an Epetra_Import object that imports objects
        from an OperatorDomainMap layout to a RowMatrixColMap layout.  This operation is needed for sparse matrix-vector
	multiplication, y = Ax, to gather x elements for local multiplication operations.

	If RowMatrixColMap==OperatorDomainMap, then the pointer will be returned as 0.

    \return Raw pointer to importer.  This importer will be valid as long as the Epetra_RowMatrix object is valid.
  */
  virtual const Epetra_Import* Importer() const {return(Importer_);}
  
  //! Returns the Epetra_Export object that contains the export operations for distributed operations, returns zero if none.
    /*! If RowMatrixRowMap!=OperatorRangeMap, then this method returns a pointer to an Epetra_Export object that exports objects
        from an RowMatrixRowMap layout to a OperatorRangeMap layout.  This operation is needed for sparse matrix-vector
	multiplication, y = Ax, to scatter-add y elements generated during local multiplication operations.

	If RowMatrixRowMap==OperatorRangeMap, then the pointer will be returned as 0.  For a typical Epetra_RowMatrix object,
	this pointer will be zero since it is often the case that RowMatrixRowMap==OperatorRangeMap.

    \return Raw pointer to exporter.  This exporter will be valid as long as the Epetra_RowMatrix object is valid.
  */
  virtual const Epetra_Export* Exporter() const {return(Exporter_);}

  //@}

 protected:

  void UpdateImportVector(int NumVectors) const;
  void UpdateExportVector(int NumVectors) const;
  Epetra_Comm * Comm_;
  Epetra_Map OperatorDomainMap_;
  Epetra_Map OperatorRangeMap_;
  Epetra_Map RowMatrixRowMap_;
  Epetra_Map RowMatrixColMap_;
  
  int NumMyNonzeros_;
  int NumGlobalNonzeros_;
  int MaxNumEntries_;
  double NormInf_;
  double NormOne_;
  int NumMyRows_;
  int NumMyCols_;

  bool UseTranspose_;
  bool HasNormInf_;
  bool LowerTriangular_;
  bool UpperTriangular_;
    

  mutable Epetra_MultiVector * ImportVector_;
  mutable Epetra_MultiVector * ExportVector_;
  Epetra_Import * Importer_;
  Epetra_Export * Exporter_;

};
#endif /* EPETRA_BASICROWMATRIX_H */
