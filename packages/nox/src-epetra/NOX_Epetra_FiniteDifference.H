
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

#include "Epetra_RowMatrix.h"    // base class

#include "NOX_Common.H"          // for <string>

// Forward Declarations
class Epetra_Comm;
class Epetra_Map;
class Epetra_Import;
class Epetra_Vector;
class Epetra_CrsGraph;
class Epetra_CrsMatrix;

namespace NOX {
  class Utils;
  namespace Epetra {
    class Interface;
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

  Since this inherits from the Epetra_RowMatrix class, it can be used as the preconditioning matrix for AztecOO preconditioners.  This method is very inefficient when computing the Jacobian and is not recommended for large-scale systems but only for debugging purposes.
  */
class FiniteDifference : public Epetra_RowMatrix {
      
 public:

  //! Constructor 
  FiniteDifference(Interface& i, const Epetra_Vector& initialGuess, double beta = 1.0e-6, double alpha = 1.0e-4);

  //! Constructor
  FiniteDifference(Interface& i, const Epetra_Vector& initialGuess, const Epetra_Vector& beta, double alpha = 1.0e-4);

  //! Constructor that takes a pre-constructed Epetra_CrsGraph so it does not have to determine the non-zero entries in the matrix.
  FiniteDifference(Interface& i, const Epetra_Vector& initialGuess, const Epetra_CrsGraph& g, double beta = 1.0e-6, double alpha = 1.0e-4);

  //! Constructor that takes a pre-constructed Epetra_CrsGraph so it does not have to determine the non-zero entries in the matrix.
  FiniteDifference(Interface& i, const Epetra_Vector& initialGuess, const Epetra_CrsGraph& g, const Epetra_Vector& beta, double alpha = 1.0e-4);

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
  virtual const Epetra_Map & OperatorDomainMap() const;

  //!Returns the Epetra_BlockMap object associated with the range of this matrix operator.
  virtual const Epetra_Map & OperatorRangeMap() const;

  //! See Epetra_RowMatrix documentation.
  virtual bool Filled() const;

  //! See Epetra_RowMatrix documentation.
  virtual int NumMyRowEntries(int MyRow, int & NumEntries) const;
  
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

  //! Compute Jacobian given the specified input vector, x. Returns true if computation was successful.
  virtual bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac);

  //! Compute an Epetra_RowMatrix to be used by Aztec preconditioners given the specified input vector, x. Returns true if computation was successful.
  virtual bool computePreconditioner(const Epetra_Vector& x, Epetra_RowMatrix& M);

protected:

  //! Constructs an Epetra_CrsGraph and Epetra_RowMatrix for the Jacobian.  This is only called if the user does not supply an Epetra_CrsGraph.
  Epetra_CrsMatrix* createGraphAndJacobian(Interface& i, const Epetra_Vector& x);

protected:

  //! Map that contains element distribution data
  const Epetra_BlockMap& map;

  //! Pointer to the Jacobian graph.
  Epetra_CrsGraph* graph;

  //! Pointer to the Jacobian.
  Epetra_CrsMatrix* jacobian;

  //! User provided interface function.
  Interface& interface;

  //! Perturbed solution vector - a work array that needs to be mutable. 
  mutable Epetra_Vector x_perturb;

  //! Function evaluation at currentX - a work array that needs to be mutable.
  mutable Epetra_Vector fo;

  //! Function evaluation at perturbX - a work array that needs to be mutable.
  mutable Epetra_Vector fp;

  //! Column vector of the jacobian - a work array that needs to be mutable.
  mutable Epetra_Vector Jc;

  //! Constant for the perturbation calculation.
  double alpha;

  //! Constant for the perturbation calculation.
  double beta;

  //! Vector for the perturbation calculation.
  const Epetra_Vector* betaVector;

  //! Define types for the \f$ \beta \f$ parameter during the computation of the perturbation parameter \f$ \delta\f$. 
  enum BetaType {Scalar, Vector};

  //! Flag that sets whether \f$ \beta \f$ is a scalar or a vector.
  BetaType betaType;

  //! label for the Epetra_RowMatrix
  string label;

};
}  // namespace Epetra
}  // namespace NOX

#endif /* _EPETRA_FINITEDIFFERENCE_H_ */
