// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#ifndef NOX_EPETRA_GROUP_H
#define NOX_EPETRA_GROUP_H

#include "NOX_Abstract_Group.H"	// base class
#include "NOX_Epetra_Vector.H"	// class data element
#include "NOX_Common.H" // class data element (string)

// Forward declares
namespace NOX {
namespace Epetra {
class Interface;
class SharedOperator;
}
namespace Parameter {
class List;
}
}
class Epetra_Vector;
class Epetra_Operator;
class Epetra_RowMatrix;
class AztecOO;

namespace NOX {
namespace Epetra {

//! Concrete implementation of NOX::Abstract::Group for %Trilinos/%Epetra.
/*! This group is set up to use the linear algebra services provided through the Trilinos/Epetra package with AztecOO for the linear %solver.
 
<B>Supplying a Jacobian %Operator</B><br>
A linear solve using this class requires at a minimum that the user supply a Jacobian object derived from the pure virtual Epetra_Operator class.  This does not have to be an explicit matrix since Newton-Krylov methods only require the action of the Jacobian on a vector (Jy).  The user has five options for providing the Epetra_Operator Jacobian.  This operator is passed in to through the group constructor.  Depending on the type of object, the user may have to implement additional functions in the NOX::Epetra::Interface object (the user must always implement the computeF() function) to actually fill/evaluate the object at the current solution.  Any additional functions are listed below:
  <ul>
  <li> An Epetra_Operator derived object - For this object the user must implement the computeJacobian() method in the NOX::Epetra::Interface.
  <li> An Epetra_RowMatrix derived object - For this object the user must implement the computeJacobian() method in the NOX::Epetra::Interface.
  <li> A NOX::Epetra::Operator derived object - For this object, the user is not required to implement any additional methods in the NOX::Epetra::Interface.
  <li> A NOX::Epetra::FiniteDifference object - Estimates the entire Jacobian by finite differencing.  For this object, the user is not required to implement any additional methods in the NOX::Epetra::Interface.
  <li> A NOX::Epetra::MatrixFree object - Estimates the action of the Jacobian on a vector via a directional derivative.  This requires minimal memory since no Jacobian is stored.  The Jacobian is stored in an Epetra_CrsMatirx.  This is not optimied and can be very time consuming since a ::Epetra::NOXInterface::computeF() call is made for each unknown. This was meant primarily as a validation tool to ensure Jacobian fills are correct.For this object, the user is not required to implement any additional methods in the NOX::Epetra::Interface.
  </ul>

<B>"Linear %Solver" sublist parameters</B><br>
A NOX::ParameterList is supplied in the constructor and during calls to the NOX::Epetra::Group::computeNewton() and NOX::Epetra::Group::applyJacobianInverse().  This parameter list is the "Linear %Solver" sublist.  The following parameters can be set in the linear Solver sublist for NOX::Epetra::Groups
<ul>

<li> "Aztec %Solver" - The following options are valid:<br>
  <ul> 
  <li> "GMRES" - Restarted generalized minimal residual (default).
  <li> "CG" - Conjugate gradient.
  <li> "CGS" - Conjugate gradient squared.
  <li> "TFQMR" - Transpose-free quasi-minimal reasidual.
  <li> "BiCGStab" - Bi-conjugate gradient with stabilization.
  <li> "LU" - Sparse direct solve (single processor only).
  </ul><br>

<li> "Size of Krylov Subspace" - When using restarted GMRES this sets the size of the Krylov subspace (defaults to 300). <br>
 
<li> "Orthogonalization" - The orthogonalization routine use for the Gram-Schmidt orthogonalization procedure in Aztec.The following options are valid:<br>
  <ul> 
  <li> "Classical" - (default).
  <li> "Modified" 
  </ul><br>

<li> "Scaling" - The scaling of the linear system before it is solved.  We can NOT use the Aztec scaling methods since AztecOO uses the matrix-free interface to Aztec.  Therefore we implement our own scaling and apply it to the Epetra_LinearProblem object.  Output from the Aztec solve will always say no scaling since we scale things before Aztec gets the problem.  The following options are valid:<br>
  <ul>
  <li> "None" - No scaling is used (default).
  <li> "Row Sum"
  </ul><br>

<li> "Convergence Criteria" - Algorithm used to calculate the residual that is used for determining the convergence of the linear solver.  See the Aztec 2.1 manual for more information.  The following options are valid:
  <ul>
  <li> "r0" - (default)
  <li> "rhs"
  <li> "norm"
  <li> "no scaling"
  <li> "sol"
  </ul><br>
 
<li> "Tolerance" - tolerance used by AztecOO to determine if a linear solve has converged.

<li> "Max Iterations" - maximum number of linear iterations in the linear solve.  Defaults to 

<li> "Output Frequency" - number of linear solve iterations between output of the linear solve residual. Takes an integer, or one of the AztecOO flags: AZ_none, AZ_last, or AZ_all as a value. Defaults to AZ_last.

<li> "Preconditioning" - Preconditioning of the linear system also has multiple options that may depend on the Jacobian. A preconditioner object is an Epetra_Operator derived class just like th Jacobian.  If the AztecOO preconditioners are to be used, they require some knowledge of the matrix and thus at a minimum the preconditioner object must also be an Epetra_RowMatrix. The valid options are listed below:<br>
  <ul> <br>
  <li> "None" - no preconditioning.  This is the default setting.

  <li> "AztecOO: Jacobian Matrix" - An AztecOO preconditioner will be used with the Jacobian matrix being used as the preconditioning matrix.  This requires that the Epetra_Operator supplied for the Jacobian MUST be derived from an Epetra_RowMatrix (NOX will test for this and throw an error if it is not).  AztecOO preconditioners need to know the graph pattern of the matrix.  This means a NOX::Epetra::MatrixFree Jacobian object can NOT be used with this option for preconditioning.

  <li> "AztecOO: User RowMatrix" - An AztecOO preconditioner will be used where the preconditioning matrix is a separate oject than the Jacobian.  The Group constructor with separate entries for the Jacobian and preconditioner objects must be used. Since AztecOO preconditioners require matrix information, the preconditioning object (M) must be derived from an Epetra_RowMatix (NOX will test for this and throw an error if it is not). 

  <li> "User Supplied %Preconditioner" - The user supplies an Epetra_Operator derived class or a NOX::Epetra::Operator derived class.  A vector is preconditioned through calls to the ApplyInverse() function of the Epetra_Operator class.  If an Epetra_Operator is passed in, the user must implement the function computePreconditioner() in the NOX::Epetra::Interface derived class.  If a NOX::Epetra::Operator object is passed in, the user does not need to implement the NOX::Epetra::Interface::computePreconditioner() function since it is implemented as part of the NOX::Epetra::Preconditioner class. 
  </ul><br>

<li> "Aztec %Preconditioner" - If an AztecOO preconditioner is used, the type of preconditioner can be set with the flag "Aztec Preconditioner" in the "Linear %Solver" sublist (defaults to ilu). Currently supported preconditioners and their corresponding parameters that can be set in the "Linear %Solver" sublist are shown below (See the Aztec 2.1 manual for more information):<br>
  <ul><br>
  <li> "ilu" with additional parameters:
  <ul>
  <li> "Overlap" - defaults to 0.
  <li> "Graph Fill" - defaults to 0.
  </ul> <br>

  <li> "ilut" with additional parameters:
  <ul>
  <li> "Overlap" - defaults to 0.
  <li> "Fill Factor" - defaults to 1.0.
  <li> "Drop Tolerance" - defaults to 1.0e-12
  </ul> <br>

  <li> "Jacobi" - k step Jacobi where k is set by the "Steps" flag: 
  <ul>
  <li> "Steps" - defaults to 3.
  </ul><br>

  <li> "Symmetric Gauss-Siedel" - Non-overlapping domain decomposition k step symmetric Gauss-Siedel where k is set by the "Steps" flag: 
  <ul>
  <li> "Steps" - defaults to 3.
  </ul><br>

  <li> "Polynomial" - Neumann polynomial with order set by the parameter:
  <ul>
  <li> "Polynomial Order" - defaults to 3.
  </ul><br>

  <li> "Least-squares Polynomial" - Least-squares polynomial with order set by the parameter:
  <ul>
  <li> "Polynomial Order" - defaults to 3.
  </ul><br>

  </ul>

<li> "RCM Reordering" - Enables RCM reordering in conjunction with domain decomp incomplete factorization preconditioning.  The following options are valid:<br>
  <ul> 
  <li> "Disabled" - (default).
  <li> "Enabled" 
  </ul><br>

</ul>
<B>"Output" sublist</B><br>
The parameter list passed into the group during a computeNewton() or ApplyJacobianInverse() will have an "Output" sublist created that contains the following parameters:
<ul>
<li> "Number of Linear Iterations" - Number of iterations used by the linear solver <br>
<li> "True Unscaled Residual"<br>
</ul><br>

 */

class Group : public Abstract::Group {

protected:
  //! List of types of operators that can be used for the Jacobian and/or Preconditioner.
  /*! Current options are:
    <ul> 
    <li> None - No operator was supplied (used only for preconditioner).
    <li> EpetraOperator - the user implements an Epetra_Operator derived object.
    <li> EpetraRowMatrix - the user implements an Epetra_RowMatrix derived object.
    <li> NoxOperator - the user implements a NOX::Epetra::Operator derived object.
    <li> NoxFiniteDifferenceRowMatrix - the operator is a NOX::Epetra::FiniteDifference object that is derived from an Epetra_CrsMatrix. The compute() function evaluates the Jacobian using finite differencing.
    <li> NoxMatrixFreeOperator - The operator is a NOX::Epetra::MatrixFree object that derives from the Epetra_Operator class.  The compute() function evaluates the Jacobian.  
    </ul>
  */ 
  enum OperatorType {None, 
		     EpetraOperator, 
		     EpetraRowMatrix,
		     NoxOperator, 
		     NoxFiniteDifferenceRowMatrix, 
		     NoxMatrixFreeOperator
  }; 

public:
  //! Constructor with Jacobian Operator only.  
  /*! Either there is no preconditioning or the Jacobian will be used
   * for preconditioning.  An Epetra_Operator must be supplied for 
   * the Jacobian even in Matrix-Free mode. linearSolverParams is 
   * the "Linear Solver" sublist of parameter list.
   */
  Group(const Parameter::List& linearSolverParams, Interface& i, 
 	Epetra_Vector& x, Epetra_Operator& J);

  //! Constructor with Jacobian Operator only.  
  /*! Either there is no preconditioning or the Jacobian will be used
   * for preconditioning.  An Epetra_Operator must be supplied for 
   * the Jacobian even in Matrix-Free mode. linearSolverParams is 
   * the "Linear Solver" sublist of parameter list.
   */
  Group(const Parameter::List& linearSolverParams, Interface& i, 
 	Vector& x, Epetra_Operator& J);

  //! Constructor with a separate Jacobian (J) and %Preconditioner (M).  
  //! linearSolverParams is the "Linear %Solver" sublist of parameter list.
  Group(const Parameter::List& linearSolverParams, Interface& i, 
 	Epetra_Vector& x, Epetra_Operator& J, Epetra_Operator& M);

  //! Constructor with a separate Jacobian (J) and %Preconditioner (M).  
  //! linearSolverParams is the "Linear %Solver" sublist of parameter list.
  Group(const Parameter::List& linearSolverParams, Interface& i, 
 	Vector& x, Epetra_Operator& J, Epetra_Operator& M);

  /*! \brief Copy constructor. If type is DeepCopy, takes ownership of
    valid shared Jacobian and shared preconditioning matrix. */
  Group(const Group& source, CopyType type = DeepCopy);

  //! Destructor.
  virtual ~Group();

  virtual Abstract::Group& operator=(const Group& source);
  //! See above.
  virtual Abstract::Group& operator=(const Abstract::Group& source);

  /** @name "Compute" functions. */
  //@{

  virtual bool setX(const Vector& y);

  //! See above
  virtual bool setX(const Abstract::Vector& y);

  virtual bool computeX(const Group& grp, 
			const Vector& d, 
			double step);
  //! See above.
  virtual bool computeX(const Abstract::Group& grp,
			const Abstract::Vector& d, 
			double step);

  virtual bool computeF();

  virtual bool computeJacobian();

  virtual bool computeGradient();

  virtual bool computeNewton(Parameter::List& params);

  //@}

  /** @name Jacobian operations.
   *
   * Operations using the Jacobian matrix. These may not be defined in
   * matrix-free scenarios. */

  //@{
  
  virtual bool applyJacobian(const Vector& input, Vector& result) const;
  //! See above
  virtual bool applyJacobian(const Abstract::Vector& input, Abstract::Vector& result) const;

  virtual bool applyJacobianTranspose(const Vector& input, Vector& result) const;
  //! See above
  virtual bool applyJacobianTranspose(const Abstract::Vector& input, Abstract::Vector& result) const;
  
  virtual bool applyJacobianInverse(Parameter::List &params, const Vector &input, Vector &result) const;

  //! See above
  virtual bool applyJacobianInverse(Parameter::List &params, const Abstract::Vector &input, Abstract::Vector &result) const;

  virtual bool applyJacobianDiagonalInverse(const Vector& input, Vector& result) const;
  //! See above
  virtual bool applyJacobianDiagonalInverse(const Abstract::Vector& input, Abstract::Vector& result) const;
  
  //@}

  /** @name "Is" functions
   *
   * Checks to see if various objects have been computed. Returns true
   * if the corresponding "compute" function has been called since the
   * last update to the solution vector (via instantiation or
   * computeX). */
  //@{

  virtual bool isF() const;
  virtual bool isJacobian() const;
  virtual bool isGradient() const;
  virtual bool isNewton() const;
  virtual bool isNormNewtonSolveResidual() const;

  //@}

  /** @name "Get" functions 
   *
   * Note that these function do not check whether or not the vectors
   * are valid. Must use the "Is" functions for that purpose. */
  //@{

  virtual const Abstract::Vector& getX() const;

  virtual const Abstract::Vector& getF() const;
  
  virtual double getNormF() const;

  virtual const Abstract::Vector& getGradient() const;

  virtual const Abstract::Vector& getNewton() const;

  double getNormNewtonSolveResidual() const;

  //@}

  virtual Abstract::Group* clone(CopyType type = DeepCopy) const;

  //! Return the SharedJacobian.
  virtual SharedOperator& getSharedJacobian();

  //! Return the SharedPreconditioner.
  virtual SharedOperator& getSharedPreconditioner();

  //! Return the userInterface.
  virtual Interface& getUserInterface();

protected:

  //! resets the isValid flags to false
  virtual void resetIsValid();

  /** @name "Linear Solver" functions 
   */
  //@{
  //! Set any required Aztec options.  NOTE: This does not handle all aztec options.  They will be added as needed. 
  virtual void setAztecOptions(const Parameter::List& params, AztecOO& aztec) const;
  
  //! Checks to make sure that the supplied operators can be used for the requested preconditioning. 
  virtual bool checkOperatorConsistency();
  
  //! Computes the preconditioner based on the current solution vector and sets the correct preconditioning option (Operator or RowMatrix) in the AztecOO object.
  virtual bool computePreconditioner(AztecOO& aztec) const;
  
  //! Computes the 2-norm of the residual of the linear model used in the Newton solve computation, ||Js+f||.
  virtual bool computeNormNewtonSolveResidual();
  //@}

  //! Returns the type of operator that is passed into the group constructors.
  /*! Uses dynamic casting to identify the underlying object type. */
  virtual OperatorType getOperatorType(const Epetra_Operator& o);

protected:

  /** @name Vectors */
  //@{
  //! Solution vector pointer.
  Vector* xVectorPtr;
  //! Solution vector.
  Vector& xVector;
  //! Right-hand-side vector (function evaluation).
  Vector* RHSVectorPtr;
  //! Right-hand-side vector pointer (function evaluation).
  Vector& RHSVector;
  //! Gradient vector pointer(steepest descent vector).
  Vector* gradVectorPtr;
  //! Gradient vector (steepest descent vector).
  Vector& gradVector;
  //! Newton direction vector pointer.
  Vector* NewtonVectorPtr;
  //! Newton direction vector.
  Vector& NewtonVector;
  //@}

  //! An extra temporary vector, only allocated if needed.
  mutable Epetra_Vector* tmpVectorPtr;

  /** @name IsValid flags 
   *  
   * True if the current solution is up-to-date with respect to the
   * currect xVector. */
  //@{
  bool isValidRHS;
  bool isValidJacobian;
  bool isValidGrad;
  bool isValidNewton;
  bool isValidNormNewtonSolveResidual;
  //@}
  
  //! 2-Norm of RHS
  double normRHS;

  //! 2-Norm of the Newton solve residual: ||Js+f||
  double normNewtonSolveResidual;

  /** @name Shared Operators */
  //@{
  //! Pointer to shared Jacobian matrix 
  SharedOperator* sharedJacobianPtr;

  //! Reference to shared Jacobian matrix 
  SharedOperator& sharedJacobian;

  //! Pointer to shared Preconditioning matrix.  We do not create a corresponding reference as we did with the sharedJacobian since this object is optionally supplied. 
  SharedOperator* sharedPreconditionerPtr;
  //@}

  //! Reference to shared Jacobian matrix 
  SharedOperator& sharedPreconditioner;

  //! Reference to the user supplied interface functions
  Interface& userInterface;

  //! Flag that tells the code how the Jacobian operator is implemented.
  /*! This flag is set in the constructor by calling the function getJacobianType().
  */ 
  OperatorType jacobianOperatorType;
  
  //! Flag that tells the code how the Preconditioning operator is implemented.
  /*! This flag is set in the constructor by calling the function 
    getPrecType(). 
  */ 
  OperatorType preconditionerOperatorType;
  
  //! Determines how the preconditioning is handled. 
  /*! This variable is set in the constructors from the "Preconditioning" parameter in the "Linear Solver" sublist.
   <li> Key in "Linear Solver" sublist: "Preconditioning"
   <li> Valid Options:
   <ul>
   <li> "None" - no precondtioning is used.
   <li> "AztecOO: Jacobian Matrix" - AztecOO preconditioner is used.  This requires an explicit Jacobian Matrix.  Therefore the Jacobian matrix must be an Epetra_RowMatrix or a NOX_FiniteDifference object.
   <li> "AztecOO: User RowMatrix" - AztecOO preconditioner is used. This requires an explicit Epetra_RowMatrix that is NOT the same as the Jacobian.  This was specifically written for "Matrix-Free" mode when there is no explicit Matrix for the Jacobian.  
   <li> "User Supplied Preconditioner" - user supplied routine computes preconditioner and applies it to vector.  The user must implement the routine as an Epetra_Operator class such that the Apply() function preconditions the incoming vector. 
   </ul>
  */ 
  string preconditioner;
  
};

} // namespace Epetra
} // namespace NOX


#endif
