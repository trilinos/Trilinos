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

#ifndef _AMESOS_EPETRA_INTERFACE_H_
#define _AMESOS_EPETRA_INTERFACE_H_

class Epetra_Import;
class Epetra_CrsMatrix;
class Epetra_RowMatrix;
class Epetra_CrsMatrix;
class Epetra_VbrMatrix;
class Epetra_MultiVector;
#include "Epetra_SerialDenseVector.h"
class Epetra_IntSerialDenseVector;
class Epetra_SerialDenseMatrix;

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


class Amesos_EpetraInterface : public Amesos_BaseSolver {
  
public:
  
  Amesos_EpetraInterface(const Epetra_LinearProblem * Problem);
  
  ~Amesos_EpetraInterface();

  
  int SetOperator(Epetra_RowMatrix * Mat);
  int GetRow(int BlockRow, int & NumIndices,
	     int * & RowIndices, 
	     int * & ColIndices, double * & Values);

  inline int MatrixType() const
  {
    return MatrixType_;
  }
  
  inline int NumMyRows() const
  {
    return NumMyRows_;
  }

  inline int NumMyBlockRows() const 
  {
    return NumMyBlockRows_;
  }
  
  inline int NumGlobalRows() const
  {
    return NumGlobalRows_;
  }
  
  inline int NumMyNonzeros() const 
  {
    return NumMyNonzeros_;
  }
  
  inline int NumGlobalNonzeros() const
  {
    return NumGlobalNonzeros_;
  }
  
  inline int MaxNumEntries() const 
  {
    return MaxNumEntries_;
  }
  
  inline int MyGlobalElements(int i) const 
  {
    return MyGlobalElements_[i];
  }

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

  inline Epetra_RowMatrix * RowA() const
  {
    return RowA_;
  }
  
  inline Epetra_CrsMatrix * CrsA() const
  {
    return CrsA_;
  }

  inline Epetra_VbrMatrix * VbrA() const
  {
    return VbrA_;
  }

  //! Get a pointer to the Problem.
  const Epetra_LinearProblem *GetProblem() const { return(Problem_); };

  //! Returns a pointer to the Epetra_Comm communicator associated with this matrix.
  const Epetra_Comm & Comm() const {return(GetProblem()->GetOperator()->Comm());};

  int AddToSymFactTime(double);
  int AddToNumFactTime(double);
  int AddToSolFactTime(double);
  
  inline int GetSymFactTime() const
  {
    return SymFactTime_;
  }
  
  inline int GetNumFactTime() const
  {
    return SymFactTime_;
  }

  inline int GetSolFactTime() const 
  {
    return SymFactTime_;
  }

  //! Set the matrix property (unsymmetric, SPD, general symmetric).
  /*! Set the matrix property as follows:
     - 0 : general unsymmetric matrix;
     - 1 : SPD;
     - 2 : general symmetric matrix.
  */
  int SetMatrixProperty(int property) 
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

  //! Retrive the matrix property.
  int MatrixProperty() const
  {
    return(MatrixProperty_);
  }

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
  Epetra_CrsMatrix * CrsA_;               // MS // cast RowMatrix to Crs (if possible)
  Epetra_VbrMatrix * VbrA_;               // MS // cast RowMatrix to Vbr (if possible)

  Epetra_IntSerialDenseVector * RowIndices_;
  Epetra_IntSerialDenseVector * ColIndices_;
  Epetra_SerialDenseVector    * Values_;

  Epetra_SerialDenseMatrix ** Entries_;

  int BlockIndices_;
  int NumPDEEqns_;
  
  bool IsSetOperatorOK_;

  const Epetra_LinearProblem * Problem_;

  double SymFactTime_;                    // MS // keep trace of timing
  double NumFactTime_;
  double SolFactTime_;
  
  int MatrixProperty_;

}; /* Amesos_EpetraInterface */


#endif
