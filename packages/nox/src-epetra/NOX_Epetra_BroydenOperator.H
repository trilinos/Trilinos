//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef NOX_EPETRA_BROYDENOPERATOR_H
#define NOX_EPETRA_BROYDENOPERATOR_H

#include "NOX.H"
#include "NOX_Common.H"
#include "NOX_Epetra.H"
#include "Teuchos_ParameterList.hpp"
#include "NOX_Abstract_PrePostOperator.H"

#include "Epetra_Vector.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"

namespace NOX {

namespace Epetra {

/*! \brief A concrete implementation of a Broyden-type operator for NOX.

  This operator is intended to allow cheap updates to an existing Jacobian
  or preconditioning matrix that would otherwise be difficult or impossible
  to obtain by other means.  It computes updates using secant approximations
  emobdied in truncated %Broyden updates that preserve matrix sparsity.

  This class derives from NOX::Abstract::PrePostOperator in order to
  perform a Broyden-type update on an existing matrix that it holds but does
  not own.  This update is performed after each nonlinear iteration within
  method runPostIterate(...) according to the recursive formula:

  \f[ \tilde{B}_{k+1} = \tilde{B}_k + \frac{({y_k -
    \tilde{B}_k s_k})s_k^T}{s^T s} \f]

  where \f[ y_k = F_{k+1} - F_k \f] and \f[ s_k = x_{k+1} - x_k \f]

  The tilde on the matrices \f$ B \f$ indicates that the updates
  are constrained so that the nonzero structure of the original matrix
  passed into the constructor is preserved.  Inasmuch as unconstrained
  %Broyden updates produce dense matrices, these constrained updates lead
  to a loss of Broyden-matrix properties, e.g.
  
  \f[ \tilde{B}_{k+1} s_k \ne \tilde{B}_k + s_k \f]

  \f[ \tilde{B}_{k+1} q \ne \tilde{B}_k q \quad \forall q : s_k^T q = 0 \f]

  One could recover these properties by passing into the constructor
  a dense %Epetra_CrsMatrix, though the cost of typical use of this
  matrix, e.g.  applying ILU to it, would be significant.  Additionally,
  "better" values obtained from another Jacobian or preconditioning
  matrix can be used to replace corresponding values in the updated
  %Broyden matrix by passing the Jacobian or preconditioning matrix and
  its associated interface to the constructor.  The structure of the
  Jacobain or preconditioning matrix typically represents a subset of
  the %Broyden matrix, e.g. a block diagonal matrix.

 */
class BroydenOperator : public Epetra_RowMatrix,
                        public NOX::Abstract::PrePostOperator,
                        public NOX::Epetra::Interface::Jacobian,
                        public NOX::Epetra::Interface::Preconditioner 
{

 public:

  //! Constructor taking an initial matrix to be updated
  BroydenOperator(Teuchos::ParameterList & nlParams, 
		  Epetra_Vector & solnVec,
                  const Teuchos::RefCountPtr<Epetra_CrsMatrix>& broydMat0, 
		  bool verbose = false);

  //! Constructor taking an initial matrix to be updated along with a Jacobian matrix whose values replace those of the %Broyden update.
  BroydenOperator(
       Teuchos::ParameterList & nlParams, 
       Epetra_Vector & solnVec,
       const Teuchos::RefCountPtr<Epetra_CrsMatrix>& broydMat0, 
       const Teuchos::RefCountPtr<NOX::Epetra::Interface::Jacobian>& jacInt,
       const Teuchos::RefCountPtr<Epetra_CrsMatrix>& jacMatrix, 
       bool verbose = false);

  //! Constructor taking an initial matrix to be updated along with a preconditioning matrix whose values replace those of the %Broyden update.
  BroydenOperator(
   Teuchos::ParameterList & nlParams, 
   Epetra_Vector & solnVec,
   const Teuchos::RefCountPtr<Epetra_CrsMatrix>& broydMat0,
   const Teuchos::RefCountPtr<NOX::Epetra::Interface::Preconditioner>& precInt,
   const Teuchos::RefCountPtr<Epetra_CrsMatrix>& precMatrix, 
   bool verbose = false);

  //! Copy Constructor
  BroydenOperator(const BroydenOperator &);

  //! Destructor
  virtual ~BroydenOperator();

  //! Returns a character string describing the name of the operator
  virtual const char* Label () const;

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
  virtual const Epetra_Map & OperatorDomainMap() const;

  //!Returns the Epetra_BlockMap object associated with the range of this matrix operator.
  virtual const Epetra_Map & OperatorRangeMap() const;

  //! See Epetra_RowMatrix documentation.
  virtual bool Filled() const;

  //! See Epetra_RowMatrix documentation.
  virtual int NumMyRowEntries(int MyRow, int & NumEntries) const;

  //! See Epetra_RowMatrix documentation.
  virtual int MaxNumEntries() const;

  //! See Epetra_RowMatrix documentation.
  virtual int ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const;

  //! See Epetra_RowMatrix documentation.
  virtual int ExtractDiagonalCopy(Epetra_Vector & Diagonal) const;

  //! See Epetra_RowMatrix documentation.
  virtual int Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! See Epetra_RowMatrix documentation.
  virtual int Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X,  Epetra_MultiVector& Y) const;

  //! See Epetra_RowMatrix documentation.
  virtual int InvRowSums(Epetra_Vector& x) const;

  //! See Epetra_RowMatrix documentation.
  virtual int LeftScale(const Epetra_Vector& x);

  //! See Epetra_RowMatrix documentation.
  virtual int InvColSums(Epetra_Vector& x) const;

  //! See Epetra_RowMatrix documentation.
  virtual int RightScale(const Epetra_Vector& x);

  //! See Epetra_RowMatrix documentation.
  virtual double NormInf() const;

  //! See Epetra_RowMatrix documentation.
  virtual double NormOne() const;

  //! See Epetra_RowMatrix documentation.
  virtual int NumGlobalNonzeros() const;

  //! See Epetra_RowMatrix documentation.
  virtual int NumGlobalRows() const;

  //! See Epetra_RowMatrix documentation.
  virtual int NumGlobalCols() const;

  //! See Epetra_RowMatrix documentation.
  virtual int NumGlobalDiagonals() const;

  //! See Epetra_RowMatrix documentation.
  virtual int NumMyNonzeros() const;

  //! See Epetra_RowMatrix documentation.
  virtual int NumMyRows() const;

  //! See Epetra_RowMatrix documentation.
  virtual int NumMyCols() const;

  //! See Epetra_RowMatrix documentation.
  virtual int NumMyDiagonals() const;

  //! See Epetra_RowMatrix documentation.
  virtual bool LowerTriangular() const;

  //! See Epetra_RowMatrix documentation.
  virtual bool UpperTriangular() const;

  //! See Epetra_RowMatrix documentation.
  virtual const Epetra_Comm & Comm() const;

  //! See Epetra_RowMatrix documentation.
  virtual const Epetra_Map & RowMatrixRowMap() const;

  //! See Epetra_RowMatrix documentation.
  virtual const Epetra_Map & RowMatrixColMap() const;

  //! See Epetra_RowMatrix documentation.
  virtual const Epetra_Import * RowMatrixImporter() const;

  //! See Epetra_SrcDistObj documentation.
  virtual const Epetra_BlockMap& Map() const;

  //! Set the current step vector, \f[ y_k = x_{k+1} - x_k \f]
  void setStepVector ( Epetra_Vector & vec );

  //! Set the current step vector, \f[ y_k = x_{k+1} - x_k \f]
  void setStepVector ( NOX::Epetra::Vector & vec );

  //! Set the current yield vector, \f[ y_k = F_{k+1} - F_k \f]
  void setYieldVector( Epetra_Vector & vec );

  //! Set the current yield vector, \f[ y_k = F_{k+1} - F_k \f]
  void setYieldVector( NOX::Epetra::Vector & vec );

  //! Compute the sparse Broyden update
  bool computeSparseBroydenUpdate();

  //! Return a reference to the %Broyden matrix.  The matrix is not owned but is obtained from the client at construction.
  Epetra_CrsMatrix & getBroydenMatrix()
  { return *crsMatrix; };

 protected:

  //! Initialize operator and data
  virtual bool initialize( Teuchos::ParameterList & nlParams, const Epetra_Vector & x );

  //! Allow for fill of Jacobian matrix whose values will replace corresponding entries in the %Broyden matrix.
  virtual bool computeJacobian( const Epetra_Vector & x, Epetra_Operator& Jac);

  //! Allow for fill of preconditioning matrix whose values will replace corresponding entries in the %Broyden matrix.
  virtual bool computePreconditioner( const Epetra_Vector &, 
				      Epetra_Operator& Prec,
                                      Teuchos::ParameterList * params = 0 );

  //! Old Update method
  virtual void old_runPostIterate( const NOX::Solver::Generic & solver);

  //! Set a flag to skip the first call to computeSparseBroydenUpdate as no valid step vector or residual vector is available 
  virtual void runPreSolve( const NOX::Solver::Generic & solver);

  //! Update the %Broyden matrix using changes in residuals the solution vector from the most recent nnlinear iteration.
  virtual void runPostIterate( const NOX::Solver::Generic & solver);

  //! Replace values in %Broyden matrix with either Jacobian or preconditioning matrix entries.
  void replaceBroydenMatrixValues( const Epetra_CrsMatrix & mat);

 protected:

  bool verbose;

  bool firstIteration;

  Teuchos::RefCountPtr<NOX::Epetra::Vector> stepVec;

  Teuchos::RefCountPtr<NOX::Epetra::Vector> yieldVec;

  Teuchos::RefCountPtr<NOX::Epetra::Vector> workVec;

  Teuchos::RefCountPtr<NOX::Epetra::Vector> updateVectorPtr;
  NOX::Epetra::Vector & updateVector;

  Teuchos::RefCountPtr<NOX::Epetra::Vector> broydenVecPtr;
  Teuchos::RefCountPtr<NOX::Epetra::Vector> residualVecPtr;

  Teuchos::RefCountPtr<NOX::Epetra::Vector> tempVecPtr;

  Teuchos::RefCountPtr<Epetra_CrsMatrix> crsMatrix;

  Teuchos::RefCountPtr<NOX::Epetra::Interface::Jacobian> jacIntPtr;
  Teuchos::RefCountPtr<Epetra_CrsMatrix> jacMatrixPtr;

  Teuchos::RefCountPtr<NOX::Epetra::Interface::Preconditioner> precIntPtr;
  Teuchos::RefCountPtr<Epetra_CrsMatrix> precMatrixPtr;

  //! label for the Epetra_RowMatrix
  string label;
};
} // namespace Epetra
} // namespace NOX

#endif /* _NOX_EPETRANEW_BROYDENOPERATOR */
