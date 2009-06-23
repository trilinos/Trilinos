//@HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2009) Sandia Corporation
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
//@HEADER

#ifndef EPETRAEXT_HYPREIJMATRIX_H_
#define EPETRAEXT_HYPREIJMATRIX_H_

// Trilinos source files
#include "Epetra_Object.h"
#include "Epetra_CompObject.h"
#include "Epetra_BasicRowMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

//Hypre source files
#include "_hypre_parcsr_mv.h"
#include "HYPRE_parcsr_mv.h"
#include "HYPRE_IJ_mv.h"
#include "_hypre_IJ_mv.h"
#include "HYPRE.h"

class Epetra_Vector;
class Epetra_MultiVector;
class Epetra_Import;

//! EpetraExt_HypreIJMatrix: A class for constructing and using real-valued sparse compressed row matrices.

/*! The EpetraExt_HypreIJMatrix is a wrapper class for Hypre parallel IJ matrices.  It is
    derived from the Epetra_BasicRowMatrix class, and so provides Hypre users access to Trilinos solvers.
    This class is lightweight, i.e., there are no deep copies of matrix data.  Whenever possible, class
    methods utilize callbacks to native Hypre functions.  
*/    

class EpetraExt_HypreIJMatrix: public Epetra_BasicRowMatrix  {
      
 public:

   //! @name Constructors/Destructor
  //@{ 
  //! Epetra_HypreIJMatrix constructor.
  /*! Creates a Epetra_HypreIJMatrix object by encapsulating an existing Hypre matrix.
    
    \param matrix (In) - A completely constructed Hypre IJ matrix.
  */
  EpetraExt_HypreIJMatrix(HYPRE_IJMatrix matrix);

  //! EpetraExt_HypreIJMatrix Destructor
  virtual ~EpetraExt_HypreIJMatrix();
  //@}
  
  //! @name Extraction methods
  //@{ 

    //! Returns a copy of the specified local row in user-provided arrays.
    /*! 
    \param MyRow (In) - Local row to extract.
    \param Length (In) - Length of Values and Indices.
    \param NumEntries (Out) - Number of nonzero entries extracted.
    \param Values (Out) - Extracted values for this row.
    \param Indices (Out) - Extracted local column indices for the corresponding values.
	  
    \return Integer error code, set to 0 if successful.
  */
    int ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const;
    
   //! Returns a reference to the ith entry in the matrix, along with its row and column index.
    /*! 
    \param CurEntry (In) - Index of local entry (from 0 to NumMyNonzeros()-1) to extract. 
    \param Value (Out) - Extracted reference to current values. 
    \param RowIndex (Out) - Row index for current entry. 
    \param ColIndex (Out) - Column index for current entry.
	  
    \return Integer error code, set to 0 if successful.
   */
   
    int ExtractMyEntryView(int CurEntry, double *&Value, int &RowIndex, int &ColIndex);
    
    //! Returns a const reference to the ith entry in the matrix, along with its row and column index.
    /*! 
    \param CurEntry (In) - Index of local entry (from 0 to NumMyNonzeros()-1) to extract. 
    \param Value (Out) - Extracted reference to current values.
    \param RowIndex (Out) - Row index for current entry. 
    \param ColIndex (Out) - Column index for current entry.
	  
    \return Integer error code, set to 0 if successful.
   */
   
    int ExtractMyEntryView(int CurEntry, const double *&Value, int &RowIndex, int &ColIndex) const;
    //@}

    //! @name Computational methods
  //@{ 

    //! Returns the result of a EpetraExt_HypreIJMatrix multiplied by a Epetra_MultiVector X in Y.
    /*! 
    \param TransA (In) -If true, multiply by the transpose of matrix, otherwise just use matrix.
    \param X (In) - A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Y (Out) -A Epetra_MultiVector of dimension NumVectorscontaining result.

    \return Integer error code, set to 0 if successful.
  */
    int Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    //! Returns the result of a EpetraExt_HypreIJMatrix multiplied by a Epetra_MultiVector X in Y.
    /*! 
    \param Upper (In) -If true, solve Ux = y, otherwise solve Lx = y.
    \param Trans (In) -If true, solve transpose problem.
    \param UnitDiagonal (In) -If true, assume diagonal is unit (whether it's stored or not).
    \param X (In) - A Epetra_MultiVector of dimension NumVectors to solve for.
    \param Y (Out) -A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
  */
    int Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    //! Scales the EpetraExt_HypreIJMatrix on the left with a Epetra_Vector x.
    /*! The \e this matrix will be scaled such that A(i,j) = x(i)*A(i,j) where i denotes the row number of A
        and j denotes the column number of A.
    \param X (In) -A Epetra_Vector to solve for.

    \return Integer error code, set to 0 if successful.
  */
    int LeftScale(const Epetra_Vector& X);


    //! Scales the EpetraExt_HypreIJMatrix on the right with a Epetra_Vector x.
    /*! The \e this matrix will be scaled such that A(i,j) = x(j)*A(i,j) where i denotes the global row number of A
        and j denotes the global column number of A.
    \param X (In) -The Epetra_Vector used for scaling \e this.

    \return Integer error code, set to 0 if successful.
  */
    int RightScale(const Epetra_Vector& X);
  //@}

  //! @name Additional methods required to support the Epetra_Operator interface
  //@{ 

    //! Returns the result of a Epetra_Operator applied to a Epetra_MultiVector X in Y.
    /*! 
    \param X (In) - A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Y (Out) -A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
  */
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
    return(EpetraExt_HypreIJMatrix::Multiply(EpetraExt_HypreIJMatrix::UseTranspose(), X, Y));};

    //! Returns the result of a Epetra_Operator inverse applied to an Epetra_MultiVector X in Y.
    /*! In this implementation, we use several existing attributes to determine how virtual
        method ApplyInverse() should call the concrete method Solve().  We pass in the UpperTriangular(), 
	the EpetraExt_HypreIJMatrix::UseTranspose(), and NoDiagonal() methods. The most notable warning is that
	if a matrix has no diagonal values we assume that there is an implicit unit diagonal that should
	be accounted for when doing a triangular solve.

    \param X (In) - A Epetra_MultiVector of dimension NumVectors to solve for.
    \param Y (Out) -A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
  */
  int ApplyInverse(const Epetra_MultiVector& X,
                   Epetra_MultiVector& Y) const
    {(void)X; (void)Y; return(-1);}

    //! Returns the current UseTranspose setting.
    virtual bool UseTranspose() const {return(false);}

  //@}
  //! @name Additional methods required to implement RowMatrix interface
  //@{ 

    //! Return the current number of values stored for the specified local row.
    /*! Similar to NumMyEntries() except NumEntries is returned as an argument
        and error checking is done on the input value MyRow.
    \param MyRow (In) - Local row.
    \param NumEntries (Out) - Number of nonzero values.
	  
    \return Integer error code, set to 0 if successful.
  */
    int NumMyRowEntries(int MyRow, int & NumEntries) const;

  //@}

 protected:

   int InitializeDefaults();
   
    mutable HYPRE_IJMatrix Matrix_;
    mutable HYPRE_ParCSRMatrix ParMatrix_;
    mutable HYPRE_IJVector X_hypre;
    mutable HYPRE_IJVector Y_hypre;
    mutable HYPRE_ParVector par_x;
    mutable HYPRE_ParVector par_y;
    mutable hypre_ParVector *x_vec;
    mutable hypre_ParVector *y_vec;
    mutable hypre_Vector *x_local;
    mutable hypre_Vector *y_local;
     
    int NumMyRows_;
    int NumGlobalRows_;
    int NumGlobalCols_;
    int MyRowStart_;
    int MyRowEnd_;
    mutable int MatType_;
};
#endif /* EPETRAEXT_HYPREIJMATRIX_H_ */
