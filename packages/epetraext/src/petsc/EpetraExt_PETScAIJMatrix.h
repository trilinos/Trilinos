//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER

/*#############################################################################
# CVS File Information
#    Current revision: $Revision$
#    Last modified:    $Date$
#    Modified by:      $Author$
#############################################################################*/

#ifndef _EPETRAEXT_PETSCAIJMATRIX_H_
#define _EPETRAEXT_PETSCAIJMATRIX_H_

#include "Epetra_Object.h"
#include "Epetra_CompObject.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Map.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
extern "C" {
//Petsc headers.
//Note: Petsc internally hard-codes paths to headers, relative to the PETSC home
//      directory.  This means that --with-incdirs must contain the full path(s)
//      to the header below plus the PETSc home directory.
#include "src/mat/impls/aij/mpi/mpiaij.h"
}

class Epetra_Import;
class Epetra_Export;
class Epetra_Vector;
class Epetra_MultiVector;

//! Epetra_PETScAIJMatrix: A class for constructing and using real-valued sparse compressed row matrices.

/*! The Epetra_PETScAIJMatrix is a wrapper class for PETSc sequential or parallel AIJ matrices.  It is
    derived from the Epetra_RowMatrix class, and so provides PETSc users access to Trilinos preconditioners.
    This class is lightweight, i.e., there are no deep copies of matrix data.  Whenever possible, class
    methods utilize callbacks to native PETSc functions.  Currently, only sequential and parallel point AIJ
    PETSc matrix types are supported.
*/    

class Epetra_PETScAIJMatrix: public Epetra_Object, public Epetra_CompObject, public virtual Epetra_RowMatrix  {
      
 public:

   //! @name Constructors/Destructor
  //@{ 
  //! Epetra_PETScAIJMatrix constructor.
  /*! Creates a Epetra_PETScAIJMatrix object by encapsulating an existing PETSc matrix.
    
    \param In
           Amat - A completely constructed PETSc SEQAIJ or MPIAIJ matrix.
  */
  Epetra_PETScAIJMatrix(Mat Amat);

  //! Epetra_PETScAIJMatrix Destructor
  virtual ~Epetra_PETScAIJMatrix();
  //@}
  
  //! @name Extraction methods
  //@{ 

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
	   Indices - Extracted local column indices for the corresponding values.
	  
    \return Integer error code, set to 0 if successful.
  */
    int ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const;

    //! Returns a copy of the main diagonal in a user-provided vector.
    /*! 
    \param Out
	   Diagonal - Extracted main diagonal.

    \return Integer error code, set to 0 if successful.
  */
    int ExtractDiagonalCopy(Epetra_Vector & Diagonal) const;
    //@}

    //! @name Computational methods
  //@{ 

    //! Returns the result of a Epetra_PETScAIJMatrix multiplied by a Epetra_MultiVector X in Y.
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

    //! Returns the result of a Epetra_PETScAIJMatrix multiplied by a Epetra_MultiVector X in Y.
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
    int Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    //! Computes the sum of absolute values of the rows of the Epetra_PETScAIJMatrix, results returned in x.
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

    //! Scales the Epetra_PETScAIJMatrix on the left with a Epetra_Vector x.
    /*! The \e this matrix will be scaled such that A(i,j) = x(i)*A(i,j) where i denotes the row number of A
        and j denotes the column number of A.
    \param In
	   x -A Epetra_Vector to solve for.

    \return Integer error code, set to 0 if successful.
  */
    int LeftScale(const Epetra_Vector& x);

    //! Computes the sum of absolute values of the columns of the Epetra_PETScAIJMatrix, results returned in x.
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

    //! Scales the Epetra_PETScAIJMatrix on the right with a Epetra_Vector x.
    /*! The \e this matrix will be scaled such that A(i,j) = x(j)*A(i,j) where i denotes the global row number of A
        and j denotes the global column number of A.
    \param In
	   x -The Epetra_Vector used for scaling \e this.

    \return Integer error code, set to 0 if successful.
  */
    int RightScale(const Epetra_Vector& x);
  //@}

  //! @name Matrix Properties Query Methods
  //@{ 


    //! If FillComplete() has been called, this query returns true, otherwise it returns false.
    bool Filled() const {return(true);};

    //! If matrix is lower triangular, this query returns true, otherwise it returns false.
    bool LowerTriangular() const {return(false);};

    //! If matrix is upper triangular, this query returns true, otherwise it returns false.
    bool UpperTriangular() const {return(false);};

  //@}
  
  //! @name Attribute access functions
  //@{ 

    //! Returns a pointer to the PETSc matrix used to create this object.
    Mat Amat() const {return(Amat_);};

    //! Returns the infinity norm of the global matrix.
    /* Returns the quantity \f$ \| A \|_\infty\f$ such that
       \f[\| A \|_\infty = \max_{1\lei\lem} \sum_{j=1}^n |a_{ij}| \f].
    */ 
    double NormInf() const;

    //! Returns the one norm of the global matrix.
    /* Returns the quantity \f$ \| A \|_1\f$ such that
       \f[\| A \|_1= \max_{1\lej\len} \sum_{i=1}^m |a_{ij}| \f].
    */ 
    double NormOne() const;

    //! Returns the number of nonzero entries in the global matrix.
    int NumGlobalNonzeros() const {return(NumGlobalNonzeros_);};

    //! Returns the number of global matrix rows.
    int NumGlobalRows() const {return(OperatorRangeMap().NumGlobalPoints());};

    //! Returns the number of global matrix columns.
    int NumGlobalCols() const {return(OperatorDomainMap().NumGlobalPoints());};

    //! Returns the number of global nonzero diagonal entries.
    int NumGlobalDiagonals() const{return(OperatorDomainMap().NumGlobalPoints());};
    
    //! Returns the number of nonzero entries in the calling processor's portion of the matrix.
    int NumMyNonzeros() const {return(NumMyNonzeros_);};

    //! Returns the number of matrix rows owned by the calling processor.
    int NumMyRows() const {return(OperatorRangeMap().NumMyPoints());};

    //! Returns the number of matrix columns owned by the calling processor.
    int NumMyCols() const {return(RowMatrixColMap().NumMyPoints());};

    //! Returns the number of local nonzero diagonal entries.
    int NumMyDiagonals() const {return(OperatorRangeMap().NumMyPoints());};

    //! Returns the Epetra_Map object associated with the domain of this operator.
    const Epetra_Map & OperatorDomainMap() const {return(*DomainMap_);};

    //! Returns the Epetra_Map object associated with the range of this operator (same as domain).
    const Epetra_Map & OperatorRangeMap() const  {return(*DomainMap_);};

    //! Implement the Epetra_SrcDistObjec::Map() function.
    const Epetra_BlockMap& Map() const {return(RowMatrixRowMap());}

    //! Returns the Row Map object needed for implementing Epetra_RowMatrix.
    const Epetra_Map & RowMatrixRowMap() const {return(OperatorRangeMap());};

    //! Returns the Column Map object needed for implementing Epetra_RowMatrix.
    const Epetra_Map & RowMatrixColMap() const {return(*ColMap_);};

    //! Returns the Epetra_Import object that contains the import operations for distributed operations.
    virtual const Epetra_Import * RowMatrixImporter() const {return(Importer_);};

    //! Returns a pointer to the Epetra_Comm communicator associated with this matrix.
    const Epetra_Comm & Comm() const {return(*Comm_);};
  //@}
  
  
  //! @name I/O Methods
  //@{ 

  //! Print method
  virtual void Print(ostream & os) const;
  //@}

  //! @name Additional methods required to support the Epetra_Operator interface
  //@{ 

    //! Returns a character string describing the operator
    const char * Label() const {return(Epetra_Object::Label());};
    
    //! If set true, transpose of this operator will be applied.
    /*! This flag allows the transpose of the given operator to be used implicitly.  Setting this flag
        affects only the Apply() and ApplyInverse() methods.  If the implementation of this interface 
	does not support transpose use, this method should return a value of -1.
      
    \param In
	   UseTranspose -If true, multiply by the transpose of operator, otherwise just use operator.

    \return Always returns 0.
  */
  int SetUseTranspose(bool UseTranspose)
    {(void)UseTranspose; return(-1);}

    //! Returns the result of a Epetra_Operator applied to a Epetra_MultiVector X in Y.
    /*! 
    \param In
	   X - A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
	   Y -A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
  */
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
    return(Epetra_PETScAIJMatrix::Multiply(Epetra_PETScAIJMatrix::UseTranspose(), X, Y));};

    //! Returns the result of a Epetra_Operator inverse applied to an Epetra_MultiVector X in Y.
    /*! In this implementation, we use several existing attributes to determine how virtual
        method ApplyInverse() should call the concrete method Solve().  We pass in the UpperTriangular(), 
	the Epetra_PETScAIJMatrix::UseTranspose(), and NoDiagonal() methods. The most notable warning is that
	if a matrix has no diagonal values we assume that there is an implicit unit diagonal that should
	be accounted for when doing a triangular solve.

    \param In
	   X - A Epetra_MultiVector of dimension NumVectors to solve for.
    \param Out
	   Y -A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
  */
  int ApplyInverse(const Epetra_MultiVector& X,
                   Epetra_MultiVector& Y) const
    {(void)X; (void)Y; return(-1);}

    //! Returns true because this class can compute an Inf-norm.
    virtual bool HasNormInf() const {return(true);}

    //! Returns the current UseTranspose setting.
    virtual bool UseTranspose() const {return(false);}

  //@}
  //! @name Additional methods required to implement RowMatrix interface
  //@{ 

    //! Return the current number of values stored for the specified local row.
    /*! Similar to NumMyEntries() except NumEntries is returned as an argument
        and error checking is done on the input value MyRow.
    \param In
           MyRow - Local row.
    \param Out
	   NumEntries - Number of nonzero values.
	  
    \return Integer error code, set to 0 if successful.
  */
    int NumMyRowEntries(int MyRow, int & NumEntries) const;

    //! Returns the maximum of NumMyRowEntries() over all rows.
    int MaxNumEntries() const;
  //@}

 private:

    int GetRow(int Row) const;
    Mat Amat_;                //general PETSc matrix type
    mutable double * Values_;
    mutable int * Indices_;
    mutable int MaxNumEntries_;
    
#ifdef HAVE_MPI
    Epetra_MpiComm * Comm_;
#else
    Epetra_SerialComm * Comm_;
#endif  
    Epetra_Map * DomainMap_;
    Epetra_Map * ColMap_;
    Epetra_Import * Importer_;
    mutable Epetra_MultiVector * ImportVector_;
 
    double NumGlobalNonzeros_;
    int NumMyNonzeros_;
    int NumMyRows_;
    int NumGlobalRows_;
    int NumMyCols_;
    int PetscRowStart_;
    int PetscRowEnd_;
    enum petscMatrixType {PETSC_SEQ_AIJ, PETSC_MPI_AIJ};
    const MatType MatType_;         //really const char*
    mutable double NormInf_;
    mutable double NormOne_;

    
 //! Copy constructor (not accessible to users).
  //FIXME we need a copy ctor
  //Epetra_PETScAIJMatrix(const Epetra_PETScAIJMatrix & Matrix) {(void)Matrix;}
};
#endif /* _EPETRAEXT_PETSCAIJMATRIX_H_ */
