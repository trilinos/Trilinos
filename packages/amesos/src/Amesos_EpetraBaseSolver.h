/*
// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
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
*/

#ifndef _AMESOS_EPETRA_BASESOLVER_H_
#define _AMESOS_EPETRA_BASESOLVER_H_

class Epetra_Import;
class Epetra_CrsMatrix;
class Epetra_RowMatrix;
class Epetra_CrsMatrix;
class Epetra_VbrMatrix;
class Epetra_MultiVector;
#include "Epetra_SerialDenseVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_Map.h"

#include "Amesos_ConfigDefs.h"
#include "Amesos_BaseSolver.h"
#include "Epetra_LinearProblem.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_Comm.h"
#endif


#define AMESOS_ROW_MATRIX 0
#define AMESOS_CRS_MATRIX 1
#define AMESOS_VBR_MATRIX 2

// numbers here as required by MUMPS
#define AMESOS_UNSYM 0   
#define AMESOS_SPD   1
#define AMESOS_SYM   2

//! Amesos_EpetraBaseSolver: A generic implementation of Amesos_BaseSolver for Epetra Matrices.
/*! Amesos_EpetraBaseSolver is a concrete implementation of
  Amesos_BaseSolver for Epetra_RowMatrix, Epetra_Crsmatrix, and
  Epetra_VbrMatrix. The class furnishes a \c GetRow function, which
  operates in \c Copy mode for \c Epetra_RowMatrix, and in \c View mode
  (hence faster) for \c Epetra_CrsMatrix and \c Epetra_VbrMatrix.

  \author Marzio Sala, SNL 9214
*/

class Amesos_EpetraBaseSolver : public Amesos_BaseSolver {
  
public:
  
  Amesos_EpetraBaseSolver(const Epetra_LinearProblem & Problem);
  
  ~Amesos_EpetraBaseSolver();

  //! Sets the interface.  
  int SetInterface(Epetra_RowMatrix * Mat);

  //! Gets a given row of Epetra_RowMatrix.
  /*! Returns the pointers to the nonzero columns and values for
    specified row. If the linear system matrix is an Epetra_VbrMatrix,
    then the specified row refers to a given block row, and this
    function returns all the nonzero elements in the block
    row. (However, in this case, \c RowIndices and \c ColIndices refer
    to non-block rows and columns.)

    \param \in BlockRow: number of (local) row (or block row for \c Epetra_RowMatrix) to obtain

    \param \out NumIndices: number of nonzero elements in specified row

    \param \our RowIndices:  pointer to an integer vector, containing all the row
    indices of nonzero elements in the specified block row. For VBR
    matrices, output row indices refers to the NON-BLOCK map.

    \param \our ColIndices:  pointer to an integer vector, containing all the column
    indices of nonzero elements in the specified block row

    \param \out Values: pointer to a double vector, containing the values.
  */
  int GetRow(int BlockRow, int & NumIndices,
	     int * & RowIndices, 
	     int * & ColIndices, double * & Values);

  //! Gets the matrix type (SPD, symmetric, or general).
  inline int MatrixType() const
  {
    return MatrixType_;
  }

  //! Returns the number of  rows in the calling process.
  inline int NumMyRows() const
  {
    return NumMyRows_;
  }

  //! Returns the number of block rows in the calling process.
  inline int NumMyBlockRows() const 
  {
    return NumMyBlockRows_;
  }

  //! Returns the number of global rows.
  inline int NumGlobalRows() const
  {
    return NumGlobalRows_;
  }

  //! Returns the number of global block rows.
  inline int NumMyNonzeros() const 
  {
    return NumMyNonzeros_;
  }

  //! Returns the number of global nonzero elements.
  inline int NumGlobalNonzeros() const
  {
    return NumGlobalNonzeros_;
  }

  //! Returns the maximum number of nonzero entries in a row.
  inline int MaxNumEntries() const 
  {
    return MaxNumEntries_;
  }

  //! Returns the local numbering for local element \c i.
  inline int MyGlobalElements(int i) const 
  {
    return MyGlobalElements_[i];
  }

  //! Returns the number of points for each element in BlockMap.
  inline int NumPDEEqns() const 
  {
    return NumPDEEqns_;
  }
  
  inline int * GetRowIndices() const 
  {
    return RowIndices_->Values();
  }
  
  inline int * GetColIndices() const 
  {
    return ColIndices_->Values();
  }

  inline double * GetValues() const 
  {
    return Values_->Values();
  }

  //! Returns a pointer to the linear system matrix, as Epetra_RowMatrix
  inline Epetra_RowMatrix * RowA() const
  {
    return RowA_;
  }

  //! Returns a pointer to the linear system matrix, as Epetra_CrsMatrix (or 0 if cast fails)
  inline Epetra_CrsMatrix * CrsA() const
  {
    return CrsA_;
  }

  //! Returns a pointer to the linear system matrix, as Epetra_VbrMatrix (or 0 if cast fails)
  inline Epetra_VbrMatrix * VbrA() const
  {
    return VbrA_;
  }

  //! Returns the index base.
  int IndexBase() const;
  
  //! Gets a pointer to the Epetra_LinearProblem.
  const Epetra_LinearProblem * GetProblem() const { return(&Problem_); };

  //! Returns a pointer to the Epetra_Comm communicator associated with this matrix.
  const Epetra_Comm & Comm() const {return(GetProblem()->GetOperator()->Comm());};

  //! Returns true if the linear problem is defined only on the calling process.
  inline bool IsLocal() const
  {
    return( IsLocal_);
  }

  //! If true, ignores off-processor contributions.
  inline int SetIsLocal(const bool flag) 
  {
    IsLocal_ = flag;
    return 0;
  }
  
  //! Set the matrix property (unsymmetric, SPD, general symmetric).
  /*! Set the matrix property as follows:
     -# 0 : general unsymmetric matrix;
     -# 1 : SPD;
     -# 2 : general symmetric matrix.
  */
  int SetMatrixProperty(const int property) 
  {
    switch( property ) {
    case AMESOS_UNSYM:
    case AMESOS_SPD:
    case AMESOS_SYM:
      MatrixProperty_ = property;
      return 0;
      break;
    default:
      return -1;
    }
  }

  //! Returns the matrix property.
  int MatrixProperty() const
  {
    return(MatrixProperty_);
  }

  inline Epetra_RowMatrix * GetMatrix() const 
  {
    return( Problem_.GetMatrix() );
  }

  //! Returns a pointer to  LHS
  inline Epetra_MultiVector * GetLHS() const
  {
    return( Problem_.GetLHS() );
  }

  //! Returns a pointer to  RHS
  inline Epetra_MultiVector * GetRHS() const
  {
    return( Problem_.GetRHS() );
  }

  // do nothing now, could be useful in the future
  inline int UpdateLHS() 
  {
    return 0;
  }

  bool MatrixShapeOK() const;  

protected:

private:

  int MaxNumEntries_;
  int MatrixType_;
  int NumMyRows_;
  int NumMyBlockRows_;
  int NumGlobalRows_;
  int NumMyNonzeros_;
  int NumGlobalNonzeros_;
  int * MyGlobalElements_;
  
  Epetra_RowMatrix * RowA_; 
  Epetra_CrsMatrix * CrsA_;               // cast RowMatrix to Crs (if possible)
  Epetra_VbrMatrix * VbrA_;               // cast RowMatrix to Vbr (if possible)

  Epetra_IntSerialDenseVector * RowIndices_;
  Epetra_IntSerialDenseVector * ColIndices_;
  Epetra_SerialDenseVector    * Values_;

  Epetra_SerialDenseMatrix ** Entries_;

  int BlockIndices_;
  int NumPDEEqns_;
  
  bool IsSetInterfaceOK_;

  const Epetra_LinearProblem & Problem_;

  int MatrixProperty_;

  bool IsLocal_;            //  1 if Problem_->GetOperator() is stored entirely on process 0
                           //  Note:  Local Problems do not require redistribution of
                           //  the matrix A or vectors X and B.  
}; /* Amesos_EpetraBaseSolver */


#endif
