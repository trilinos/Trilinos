
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

#ifndef _EPETRA_MATRIXFREE_H_
#define _EPETRA_MATRIXFREE_H_

#include "Epetra_RowMatrix.h"	// base class

#include "Epetra_Import.h"
#include "Epetra_Vector.h"

// Forward Declarations
class Epetra_Comm;
class Epetra_Map;

namespace NOX {
namespace Epetra {
class Interface;
}
}

namespace NOX {

namespace Epetra {

//! Concrete implementation for creating an Epetra_RowMatrix in Matrix-Free mode.
/*! Inherits the pure virtual classes Epetra_RowMatrix.
 * NOTE: This will change when Mike adds the 
 * pure virtual "Operator" class to Trilinos.  This will be a subset 
 * of the Epetra_RowMatrix class that contains only the essentials 
 * for needed for linear solves and matrix operations.
 */ 
class MatrixFree : public Epetra_RowMatrix {
      
 public:

  //! Constructor 
  MatrixFree(Interface& i, const Epetra_Vector& x);

  //! Pure virtual destructor
  virtual ~MatrixFree();

  virtual int SetUseTranspose(bool UseTranspose);

  virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  virtual bool UseTranspose() const;

  virtual bool HasNormInf() const;

  virtual const Epetra_BlockMap & DomainMap() const;

  virtual const Epetra_BlockMap & RangeMap() const;

  virtual bool Filled() const;

  virtual int NumMyRowEntries(int MyRow, int & NumEntries) const;
  
  virtual int ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const;
  
  virtual int ExtractDiagonalCopy(Epetra_Vector & Diagonal) const;

  virtual int Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  virtual int Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X,  Epetra_MultiVector& Y) const;

  virtual int InvRowSums(Epetra_Vector& x) const;
  
  virtual int LeftScale(const Epetra_Vector& x);
  
  virtual int InvColSums(Epetra_Vector& x) const;
  
  virtual int RightScale(const Epetra_Vector& x);
  
  virtual double NormInf() const;

  virtual double NormOne() const;
  
  virtual int NumGlobalNonzeros() const;
  
  virtual int NumGlobalRows() const;
  
  virtual int NumGlobalCols() const;
  
  virtual int NumGlobalDiagonals() const;
  
  virtual int NumMyNonzeros() const;
  
  virtual int NumMyRows() const;
  
  virtual int NumMyCols() const;
  
  virtual int NumMyDiagonals() const;
  
  virtual bool LowerTriangular() const;

  virtual bool UpperTriangular() const;

  virtual const Epetra_Comm & Comm() const;

  virtual const Epetra_BlockMap & BlockRowMap() const;

  virtual const Epetra_BlockMap & BlockImportMap() const;
  
  virtual const Epetra_Import * Importer() const;

  //! Compute Jacobian given the specified input vector, x.  Returns true if computation was successful.
  virtual bool computeJacobian(const Epetra_Vector& x, Epetra_RowMatrix& Jac);

protected:
  
  //! Dummy Importer
  Epetra_Import importer;

  //! User provided interface function
  Interface& interface;

  //! The current solution vector 
  Epetra_Vector currentX;

  //! Perturbed solution vector 
  mutable Epetra_Vector perturbX;

  //! Function evaluation at currentX 
  mutable Epetra_Vector fo;

  //! Function evaluation at perturbX 
  mutable Epetra_Vector fp;

};
} // namespace Epetra
} // namespace NOX

#endif /* _EPETRA_ROWMATRIX_H_ */
