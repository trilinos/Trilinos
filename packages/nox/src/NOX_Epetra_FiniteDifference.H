
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

#ifndef _EPETRA_FINITEDIFFERENCE_H_
#define _EPETRA_FINITEDIFFERENCE_H_

// Forward Declarations
class Epetra_Comm;
class Epetra_Map;
class Epetra_Import;
class Epetra_Vector;
class Epetra_CrsGraph;
class Epetra_CrsMatrix;
class Epetra_RowMatrix;

namespace NOX {
  class Utils;
  namespace Epetra {
    class Interface;
  }
}

namespace NOX {

namespace Epetra {

//! Concrete implementation for creating a Jacobian Epetra_RowMatrix using finite differencing.
/*! Inherits the pure virtual classes Epetra_RowMatrix. 
 * NOTE: This will change when Mike adds the 
 * pure virtual "Operator" class to Trilinos.  This will be a subset 
 * of the Epetra_RowMatrix class that contains only the essentials 
 * for needed for linear solves and matrix operations.
 */ 
class FiniteDifference : public Epetra_RowMatrix {
      
 public:

  //! Constructor 
  FiniteDifference(Interface& i, const Epetra_Vector& x);

  //! Pure virtual destructor
  virtual ~FiniteDifference();

  //! Returns a character string describing the name of the operator
  virtual char* Label () const;

  //! If set true, the transpose of this operator will be applied
  virtual int SetUseTranspose(bool UseTranspose);

  //! Return the result on an Epetra_Operator applied to an Epetra_MultiVector X in Y.
  virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! Return the result on an Epetra_Operator inverse applied to an Epetra_MultiVector X in Y.
  virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! Returns the current use transpose setting
  virtual bool UseTranspose() const;

  //! Returns true if the this object can provide an approximate Inf-norm, false otherwise.
  virtual bool HasNormInf() const;

  //!Returns the Epetra_BlockMap object associated with the domain of this matrix operator.
  virtual const Epetra_BlockMap & DomainMap() const;

  //!Returns the Epetra_BlockMap object associated with the range of this matrix operator.
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

  //! Compute Jacobian given the specified input vector, x. Returns true if computation was successful.
  virtual bool computeJacobian(const Epetra_Vector& x, Epetra_RowMatrix& Jac);

protected:

  //! Construct the Epetra_RowMatrix for the Jacobian
  Epetra_CrsMatrix* createJacobian(Interface& i, const Epetra_Vector& x);

protected:
  
  //! Map that contains element distribution data
  const Epetra_BlockMap& map;

  //! Pointer to the Jacobian graph
  Epetra_CrsGraph* graph;

  //! Pointer to the Jacobian
  Epetra_CrsMatrix* jacobian;

  //! User provided interface function
  Interface& interface;

  //! Perturbed solution vector - a work array that needs to be mutable. 
  mutable Epetra_Vector x_perturb;

  //! Function evaluation at currentX - a work array that needs to be mutable.
  mutable Epetra_Vector fo;

  //! Function evaluation at perturbX - a work array that needs to be mutable.
  mutable Epetra_Vector fp;

  //! Column vector of the jacobian - a work array that needs to be mutable.
  mutable Epetra_Vector Jc;

};
}  // namespace Epetra
}  // namespace NOX

#endif /* _EPETRA_FINITEDIFFERENCE_H_ */
