//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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

#ifndef NOX_EPETRA_FINITEDIFFERENCE_H
#define NOX_EPETRA_FINITEDIFFERENCE_H

#include "Epetra_RowMatrix.h"                       // base class
#include "NOX_Epetra_Interface_Jacobian.H"          // base class
#include "NOX_Epetra_Interface_Preconditioner.H"    // base class

#include "NOX_Epetra_Interface_Required.H" // for enum FillType
#include "NOX_Common.H"                    // for <string>
#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"

// Forward Declarations
class Epetra_Comm;
class Epetra_Map;
class Epetra_Import;
class Epetra_Vector;
class Epetra_CrsGraph;
class Epetra_CrsMatrix;

namespace NOX {
  namespace Abstract {
    class Group;
  }
  namespace Epetra {
    class Vector;
  }
}

namespace NOX {

namespace Epetra {

  /*! \brief Concrete implementation for creating an Epetra_RowMatrix Jacobian via finite differencing of the residual.

The Jacobian entries are calculated via 1st order finite differencing.  This requires \f$ N  + 1 \f$ calls to computeF() where \f$ N \f$ is the number of unknowns in the problem.

  \f[ J_{ij} = \frac{\partial F_i}{\partial x_j} = \frac{F_i(x+\delta\mathbf{e}_j) - F_i(x)}{\delta}  \f]

where \f$J\f$ is the Jacobian, \f$F\f$ is the function evaluation, \f$x\f$ is the solution vector, and \f$\delta\f$ is a small perturbation to the \f$x_j\f$ entry.

The perturbation, \f$ \delta \f$, is calculated based on one of the following equations:

\f[ \delta = \alpha * | x_j | + \beta \f]
\f[ \delta = \alpha * | x_j | + \beta_j \f]

where \f$ \alpha \f$ is a scalar value (defaults to 1.0e-4) and \f$ \beta \f$ can be either a scalar or a vector (defaults to a scalar value of 1.0e-6).  The choice is defined by the type of constructor used.  All parameters are supplied in the constructor.
In addition to the forward difference derivative approximation, backward or centered differences can be used via the setDifferenceMethod function.  Note that centered difference provides second order spatial accuracy but at the cost of twice as many function evaluations.

  Since this inherits from the Epetra_RowMatrix class, it can be used as the preconditioning matrix for AztecOO preconditioners.  This method is very inefficient when computing the Jacobian and is not recommended for large-scale systems but only for debugging purposes.
  */
class FiniteDifference : public Epetra_RowMatrix,
			 public NOX::Epetra::Interface::Jacobian,
			 public NOX::Epetra::Interface::Preconditioner {

 public:

  //! Define types for use of the perturbation parameter \f$ \delta\f$.
  enum DifferenceType {Forward, Backward, Centered};

  //! Constructor with scalar beta.
  FiniteDifference(Teuchos::ParameterList& printingParams,
                   const Teuchos::RCP<NOX::Epetra::Interface::Required>& i, 
		   const NOX::Epetra::Vector& initialGuess,
		   double beta = 1.0e-6, 
		   double alpha = 1.0e-4);

  //! Constructor with vector beta.
  FiniteDifference(Teuchos::ParameterList& printingParams,
                   const Teuchos::RCP<NOX::Epetra::Interface::Required>& i, 
		   const NOX::Epetra::Vector& initialGuess,
		   const Teuchos::RCP<const Epetra_Vector>& beta, 
		   double alpha = 1.0e-4);

  //! Constructor that takes a pre-constructed Epetra_CrsGraph so it does not have to determine the non-zero entries in the matrix.
  FiniteDifference(Teuchos::ParameterList& printingParams,
                   const Teuchos::RCP<NOX::Epetra::Interface::Required>& i, 
		   const NOX::Epetra::Vector& initialGuess,
		   const Teuchos::RCP<Epetra_CrsGraph>& g, 
		   double beta = 1.0e-6,
		   double alpha = 1.0e-4);

  //! Constructor with output control that takes a pre-constructed Epetra_CrsGraph so it does not have to determine the non-zero entries in the matrix.
  FiniteDifference(Teuchos::ParameterList& printingParams,
                   const Teuchos::RCP<NOX::Epetra::Interface::Required>& i, 
		   const NOX::Epetra::Vector& initialGuess,
		   const Teuchos::RCP<Epetra_CrsGraph>& g, 
		   const Teuchos::RCP<const Epetra_Vector>& beta,
		   double alpha = 1.0e-4);

  //! Pure virtual destructor
  virtual ~FiniteDifference();

  //! Returns a character std::string describing the name of the operator
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
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
   virtual int NumGlobalNonzeros() const;
#endif
   virtual long long NumGlobalNonzeros64() const;

  //! See Epetra_RowMatrix documentation.
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
   virtual int NumGlobalRows() const;
#endif
   virtual long long NumGlobalRows64() const;

  //! See Epetra_RowMatrix documentation.
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
   virtual int NumGlobalCols() const;
#endif
   virtual long long NumGlobalCols64() const;

  //! See Epetra_RowMatrix documentation.
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
   virtual int NumGlobalDiagonals() const;
#endif
   virtual long long NumGlobalDiagonals64() const;

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

  //! Compute Jacobian given the specified input vector, x. Returns true if computation was successful.
  virtual bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac);

  //! Compute Jacobian given the specified input vector, x. Returns true if computation was successful.
  virtual bool computeJacobian(const Epetra_Vector& x);

  //! Compute an Epetra_RowMatrix to be used by Aztec preconditioners given the specified input vector, x. Returns true if computation was successful.
  virtual bool computePreconditioner(const Epetra_Vector& x,
				     Epetra_Operator& Prec,
				     Teuchos::ParameterList* precParams = 0);

  //! Set the type of perturbation method used (default is Forward)
  virtual void setDifferenceMethod( DifferenceType type );

  //! An accessor method for the underlying Epetra_CrsMatrix
  virtual Epetra_CrsMatrix& getUnderlyingMatrix() const;

  //! Output the underlying matrix
  virtual void Print(ostream&) const;

  //! Register a NOX::Abstract::Group derived object and use the computeF() method of that group for the perturbation instead of the NOX::Epetra::Interface::Required::computeF() method.  This is required for LOCA to get the operators correct during homotopy.
  void setGroupForComputeF(NOX::Abstract::Group& group);

protected:

  //! Constructs an Epetra_CrsGraph and Epetra_RowMatrix for the Jacobian.  This is only called if the user does not supply an Epetra_CrsGraph.
  Teuchos::RCP<Epetra_CrsMatrix> 
  createGraphAndJacobian(Interface::Required& i, const Epetra_Vector& x);

  bool computeF(const Epetra_Vector& input, Epetra_Vector& result,
		NOX::Epetra::Interface::Required::FillType);

protected:

  //! Printing Utilities object
  const NOX::Utils utils;

  //! Pointer to the Jacobian graph.
  Teuchos::RCP<Epetra_CrsGraph> graph;

  //! Pointer to the Jacobian.
  Teuchos::RCP<Epetra_CrsMatrix> jacobian;

  //! User provided interface function.
  Teuchos::RCP<NOX::Epetra::Interface::Required> interface;

  //! Perturbed solution vector - a work array that needs to be mutable.
  mutable Epetra_Vector x_perturb;

  //! Function evaluation at currentX - a work array that needs to be mutable.
  mutable Epetra_Vector fo;

  //! Function evaluation at perturbX - a work array that needs to be mutable.
  mutable Epetra_Vector fp;

  //! Optional pointer to function evaluation at -perturbX - needed only for centered finite differencing
  Teuchos::RCP<Epetra_Vector> fmPtr;

  //! Column vector of the jacobian - a work array that needs to be mutable.
  mutable Epetra_Vector Jc;

  //! Constant for the perturbation calculation.
  double alpha;

  //! Constant for the perturbation calculation.
  double beta;

  //! Vector for the perturbation calculation.
  Teuchos::RCP<const Epetra_Vector> betaVector;

  //! Define types for the \f$ \beta \f$ parameter during the computation of the perturbation parameter \f$ \delta\f$.
  enum BetaType {Scalar, Vector};

  //! Flag that sets whether \f$ \beta \f$ is a scalar or a vector.
  BetaType betaType;

  //! Define types for use of the perturbation parameter \f$ \delta\f$.
  DifferenceType diffType;

  //! label for the Epetra_RowMatrix
  std::string label;

  //! Flag to enables the use of a group instead of the interface for the computeF() calls in the directional difference calculation.
  bool useGroupForComputeF;

  //! Pointer to the group for possible use in computeF() calls.
  Teuchos::RCP<NOX::Abstract::Group> groupPtr;

};
}  // namespace Epetra
}  // namespace NOX

#endif /* NOX_EPETRA_FINITEDIFFERENCE_H */
