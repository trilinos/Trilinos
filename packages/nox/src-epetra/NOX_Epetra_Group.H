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
#include "NOX_Utils.H"          // class data element
#include "NOX_Common.H"         // class data element (string)

// Forward declares
namespace NOX {
  namespace Epetra {
    class Interface;
    class SharedOperator;
    class Scaling;
  }
  namespace Parameter {
    class List;
  }
}
class Epetra_Vector;
class Epetra_Operator;
class Epetra_RowMatrix;
class AztecOO;
class Ifpack_IlukGraph;
class Ifpack_CrsRiluk;

namespace NOX {
namespace Epetra {

//! Concrete implementation of NOX::Abstract::Group for Trilinos/Epetra.
/*! 

This group is set up to use the linear algebra services provided
through the Trilinos/Epetra package with AztecOO for the linear
%solver.
 
<B>Supplying a Jacobian %Operator</B>

A linear solve using this class requires at a minimum that the user
supply a Jacobian object derived from the pure virtual Epetra_Operator
class.  This does not have to be an explicit matrix since
Newton-Krylov methods only require the action of the Jacobian on a
vector (Jy).  The user has five options for providing the
Epetra_Operator Jacobian.  This operator is passed in to through the
group constructor.  Depending on the type of object, the user may have
to implement additional functions in the NOX::Epetra::Interface object
(the user must always implement the computeF() function) to actually
fill/evaluate the object at the current solution.  Any additional
functions are listed below:

  <ul>
  
  <li> An Epetra_Operator derived object - For this object the user
  must implement the computeJacobian() method in the
  NOX::Epetra::Interface.

  <li> An Epetra_RowMatrix derived object - For this object the user
  must implement the computeJacobian() method in the
  NOX::Epetra::Interface.

  <li> A NOX::Epetra::Operator derived object - For this object, the
  user is not required to implement any additional methods in the
  NOX::Epetra::Interface.

  <li> A NOX::Epetra::FiniteDifference object - Estimates the entire
  Jacobian by finite differencing.  For this object, the user is not
  required to implement any additional methods in the
  NOX::Epetra::Interface.

  <li> A NOX::Epetra::MatrixFree object - Estimates the action of the
  Jacobian on a vector via a directional derivative.  This requires
  minimal memory since no Jacobian is stored.  The Jacobian is stored
  in an Epetra_CrsMatirx.  This is not optimied and can be very time
  consuming since a ::Epetra::NOXInterface::computeF() call is made
  for each unknown. This was meant primarily as a validation tool to
  ensure Jacobian fills are correct.For this object, the user is not
  required to implement any additional methods in the
  NOX::Epetra::Interface.

  </ul>

<B>"Linear %Solver" sublist parameters</B>

A NOX::ParameterList is supplied in the constructor and during calls
to the NOX::Epetra::Group::computeNewton() and
NOX::Epetra::Group::applyJacobianInverse().  This parameter list is
the "Linear %Solver" sublist.  The following parameters can be set in
the linear Solver sublist for NOX::Epetra::Groups <ul>

<li> "Aztec %Solver" - The following options are valid:

  <ul> 
  <li> "GMRES" - Restarted generalized minimal residual (default).
  <li> "CG" - Conjugate gradient.
  <li> "CGS" - Conjugate gradient squared.
  <li> "TFQMR" - Transpose-free quasi-minimal reasidual.
  <li> "BiCGStab" - Bi-conjugate gradient with stabilization.
  <li> "LU" - Sparse direct solve (single processor only).
  </ul>

<li> "Size of Krylov Subspace" - When using restarted GMRES this sets
the size of the Krylov subspace (defaults to 300).
 
<li> "Orthogonalization" - The orthogonalization routine use for the
Gram-Schmidt orthogonalization procedure in Aztec.  The following
options are valid:

  <ul> 
  <li> "Classical" - (default).
  <li> "Modified" 
  </ul>

<li> "Convergence Criteria" - Algorithm used to calculate the residual
that is used for determining the convergence of the linear solver.
See the Aztec 2.1 manual for more information.  The following options
are valid:

  <ul>
  <li> "r0" - (default)
  <li> "rhs"
  <li> "norm"
  <li> "no scaling"
  <li> "sol"
  </ul>
 
<li> "Tolerance" - tolerance used by AztecOO to determine if a linear solve has converged.

<li> "Ill-Conditioning Threshold" - If the upper hessenberg matrix
during GMRES generates a condition number greater than this parameter
value, aztec will exit the linear solve returning the it's current
solution.  The default is 1.0e11.
 

<li> "Preconditioner Iterations" - Number of iterations an
AztecOO_Operator should take when solving the preconditioner.  This is
only used if an AztecOO preconditioner is used and the solver makes a
call to NOX::Epetra::Group::applyRightPreconditioning().

<li> "Max Iterations" - maximum number of linear iterations in the linear solve.  Defaults to 

<li> "Output Frequency" - number of linear solve iterations between
output of the linear solve residual. Takes an integer, or one of the
AztecOO flags: AZ_none, AZ_last, or AZ_all as a value. Defaults to
AZ_last.

<li> "Preconditioning" - Preconditioning of the linear system also has
multiple options that may depend on the Jacobian. A preconditioner
object is an Epetra_Operator derived class just like th Jacobian.  If
the AztecOO preconditioners are to be used, they require some
knowledge of the matrix and thus at a minimum the preconditioner
object must also be an Epetra_RowMatrix. The valid options are listed
below:

  <ul> 
  <li> "None" - no preconditioning.  This is the default setting.

  <li> "AztecOO: Jacobian Matrix" - An AztecOO preconditioner will be
  used with the Jacobian matrix being used as the preconditioning
  matrix.  This requires that the Epetra_Operator supplied for the
  Jacobian MUST be derived from an Epetra_RowMatrix (NOX will test for
  this and throw an error if it is not).  AztecOO preconditioners need
  to know the graph pattern of the matrix.  This means a
  NOX::Epetra::MatrixFree Jacobian object can NOT be used with this
  option for preconditioning.

  <li> "AztecOO: User RowMatrix" - An AztecOO preconditioner will be
  used where the preconditioning matrix is a separate object than the
  Jacobian.  The Group constructor with separate entries for the
  Jacobian and preconditioner objects must be used. Since AztecOO
  preconditioners require matrix information, the preconditioning
  object (M) must be derived from an Epetra_RowMatix (NOX will test
  for this and throw an error if it is not).

  <li> "IFPACK: Jacobian Matrix" - An IFPACK preconditioner will be
  used where the Jacobian matrix is used as the preconditioning
  matrix.  Since IFPACK preconditioners require matrix information,
  the Jacobian must be derived from an Epetra_RowMatix (NOX will test
  for this and throw an error if it is not).

  <li> "IFPACK: User RowMatrix" - An IFPACK preconditioner will be
  used where the preconditioning matrix is a separate object than the
  Jacobian.  The Group constructor with separate entries for the
  Jacobian and preconditioner objects must be used. Since IFPACK
  preconditioners require matrix information, the preconditioning
  object (M) must be derived from an Epetra_RowMatix (NOX will test
  for this and throw an error if it is not).

  <li> "User Supplied %Preconditioner" - The user supplies an
  Epetra_Operator derived class or a NOX::Epetra::Operator derived
  class.  A vector is preconditioned through calls to the
  ApplyInverse() function of the Epetra_Operator class.  If an
  Epetra_Operator is passed in, the user must implement the function
  computePreconditioner() in the NOX::Epetra::Interface derived class.
  If a NOX::Epetra::Operator object is passed in, the user does not
  need to implement the
  NOX::Epetra::Interface::computePreconditioner() function since the
  same call is implemented as part of the NOX::Epetra::Preconditioner
  class.  


</ul>


<li> "Aztec %Preconditioner" - If an AztecOO preconditioner is used,
the type of preconditioner can be set with the flag "Aztec
Preconditioner" in the "Linear %Solver" sublist (defaults to
ilu). Currently supported preconditioners and their corresponding
parameters that can be set in the "Linear %Solver" sublist are shown
below (See the Aztec 2.1 manual for more information):

  <ul>
  <li> "ilu" with additional parameters:
  <ul>
  <li> "Overlap" - defaults to 0.
  <li> "Graph Fill" - defaults to 0.
  </ul> 

  <li> "ilut" with additional parameters:
  <ul>
  <li> "Overlap" - defaults to 0.
  <li> "Fill Factor" - defaults to 1.
  <li> "Drop Tolerance" - defaults to 1.0e-12
  </ul> 

  <li> "Jacobi" - k step Jacobi where k is set by the "Steps" flag: 
  <ul>
  <li> "Steps" - defaults to 3.
  </ul>

  <li> "Symmetric Gauss-Siedel" - Non-overlapping domain decomposition
  k step symmetric Gauss-Siedel where k is set by the "Steps" flag:

  <ul>
  <li> "Steps" - defaults to 3.
  </ul>

  <li> "Polynomial" - Neumann polynomial with order set by the parameter:
  <ul>
  <li> "Polynomial Order" - defaults to 3.
  </ul>

  <li> "Least-squares Polynomial" - Least-squares polynomial with order set by the parameter:
  <ul>
  <li> "Polynomial Order" - defaults to 3.
  </ul>

  </ul>

<li> "RCM Reordering" - Enables RCM reordering in conjunction with
domain decomp incomplete factorization preconditioning.  The following
options are valid:

  <ul> 
  <li> "Disabled" - (default).
  <li> "Enabled" 
  </ul>

<li> "Use Adaptive Linear Solve" - Enables the use of AztecOO's AdaptiveIterate() method instead of calling the Iterate() method.  This causes the preconditioning matrix to be modified to make the linear solves easier. AztecOO will attempt to solve the linear system multiple times now and if the solves are failing it will modify the preconditioner and try again.  Boolean value, defaults to false.  NOTE: This only works for internal Aztec preconditioners!  The "Preconditioning" parameter must be set to "AztecOO: Jacobian Matrix" or "AztecOO: User RowMatrix".

<li> "Max Adaptive Solve Iterations" - Maximum number of attempts that the linear solver will make when trying to solve a linear system. Defaults to 5.

</ul>

<B>"Output" sublist</B>

The parameter list passed into the group during a computeNewton() or
ApplyJacobianInverse() will have an "Output" sublist created that
contains the following parameters:

<ul>

<li> "Acheived Tolerance" - Actual tolerance achieved by the linear
solver computed via the convergence test requested.

<li> "Number of %Linear Iterations" - Number of iterations used by the
linear solver in the last call to applyJacobianInverse

<li> "Total Number of %Linear Iterations" - Total number of linear
solve iterations performed by groups that have used this input list

</ul>

 */

class Group : public virtual Abstract::Group {

protected:
  /*! 
    \brief 
    Types that define how the construction of the preconditioner
    should be handled during a call to createPreconditioner().

    Preconditioners are used in two methods: (1) during linear solves
    in applyJacobianInverse() and (2) for the method
    applyRightPreconditioner().  For the applyJacobianInverse() the
    preconditioner is created and destroyed "implicitly" in the linear
    solver.  The correct operator need only be set in the AztecOO
    solver and the solver will take care of creation and destruction
    internally.  The opposite is true for
    applyRightPreconditioner(). Here we must force an "explicit"
    construction and retention of the preconditioner since there will
    be multple calls to applyRightPreconditioner() that reuses the
    same preconditioner.  This flag type is used to specify the
    construction type.
   */

  enum PrecConstructionType {

    ExplicitConstruction, 

    ImplicitConstruction

  };

  //! List of types of operators that can be used for the Jacobian and/or Preconditioner.
  enum OperatorType {

    //! No operator was supplied (used only for preconditioner).
    None, 

    //! the user implements an Epetra_Operator derived object.
    EpetraOperator, 

    //! the user implements an Epetra_RowMatrix derived object.
    EpetraRowMatrix,

    //! the user implements a NOX::Epetra::Operator derived object.
    NoxOperator, 

    /*! 
      \brief
      the operator is a NOX::Epetra::FiniteDifference object that is
      derived from an Epetra_CrsMatrix. The compute() function
      evaluates the Jacobian using finite differencing.
    */
    NoxFiniteDifferenceRowMatrix, 

    /*!
      \brief

      The operator is a NOX::Epetra::MatrixFree object that derives
      from the Epetra_Operator class.  The compute() function
      evaluates the Jacobian.
    */
    NoxMatrixFreeOperator
  }; 

public:
  //! Constructor with Jacobian Operator only.  
  /*! Either there is no preconditioning or the Jacobian will be used
    for preconditioning.  An Epetra_Operator must be supplied for 
    the Jacobian even in Matrix-Free mode. linearSolverParams is 
    the "Linear Solver" sublist of parameter list.
   */
  Group(NOX::Parameter::List& printingParams, 
	NOX::Parameter::List& linearSolverParams, NOX::Epetra::Interface& i, 
 	Epetra_Vector& x, Epetra_Operator& J);

  //! Constructor with Jacobian Operator only.  
  /*! Either there is no preconditioning or the Jacobian will be used
    for preconditioning.  An Epetra_Operator must be supplied for 
    the Jacobian even in Matrix-Free mode. linearSolverParams is 
    the "Linear Solver" sublist of parameter list.
   */
  Group(NOX::Parameter::List& printingParams, 
	NOX::Parameter::List& linearSolverParams, NOX::Epetra::Interface& i, 
 	NOX::Epetra::Vector& x, Epetra_Operator& J);

  //! Constructor with a separate Jacobian (J) and %Preconditioner (M).  
  //! linearSolverParams is the "Linear %Solver" sublist of parameter list.
  Group(NOX::Parameter::List& printingParams, 
	NOX::Parameter::List& linearSolverParams, NOX::Epetra::Interface& i, 
 	Epetra_Vector& x, Epetra_Operator& J, Epetra_Operator& M);

  //! Constructor with a separate Jacobian (J) and %Preconditioner (M).  
  //! linearSolverParams is the "Linear %Solver" sublist of parameter list.
  Group(NOX::Parameter::List& printingParams, 
	NOX::Parameter::List& linearSolverParams, NOX::Epetra::Interface& i, 
 	NOX::Epetra::Vector& x, Epetra_Operator& J, Epetra_Operator& M);

  /*! \brief Copy constructor. If type is DeepCopy, takes ownership of
    valid shared Jacobian and shared preconditioning matrix. */
  Group(const NOX::Epetra::Group& source, NOX::CopyType type = NOX::DeepCopy);

  //! Destructor.
  virtual ~Group();

  virtual NOX::Abstract::Group& operator=(const NOX::Abstract::Group& source);
  virtual NOX::Abstract::Group& operator=(const NOX::Epetra::Group& source);

  /** @name "Compute" functions. */
  //@{

  virtual void setX(const NOX::Epetra::Vector& y);
  virtual void setX(const NOX::Abstract::Vector& y);

  virtual void computeX(const Group& grp, 
			const NOX::Epetra::Vector& d, 
			double step);
  virtual void computeX(const NOX::Abstract::Group& grp,
			const NOX::Abstract::Vector& d, 
			double step);

  virtual NOX::Abstract::Group::ReturnType computeF();

  virtual NOX::Abstract::Group::ReturnType computeJacobian();

  virtual NOX::Abstract::Group::ReturnType computeGradient();

  virtual NOX::Abstract::Group::ReturnType computeNewton(NOX::Parameter::List& params);

  //@}

  /** @name Jacobian operations.
   *
   * Operations using the Jacobian matrix. These may not be defined in
   * matrix-free scenarios. */

  //@{
  
  virtual NOX::Abstract::Group::ReturnType 
  applyJacobian(const NOX::Epetra::Vector& input, NOX::Epetra::Vector& result) const;
  virtual NOX::Abstract::Group::ReturnType 
  applyJacobian(const NOX::Abstract::Vector& input, NOX::Abstract::Vector& result) const;

  virtual NOX::Abstract::Group::ReturnType 
  applyJacobianTranspose(const NOX::Epetra::Vector& input, NOX::Epetra::Vector& result) const;
  virtual NOX::Abstract::Group::ReturnType 
  applyJacobianTranspose(const NOX::Abstract::Vector& input, NOX::Abstract::Vector& result) const;
  
  virtual NOX::Abstract::Group::ReturnType 
  applyJacobianInverse(NOX::Parameter::List &params, const NOX::Epetra::Vector &input, NOX::Epetra::Vector &result) 
    const;
  virtual NOX::Abstract::Group::ReturnType 
  applyJacobianInverse(NOX::Parameter::List &params, const NOX::Abstract::Vector &input, NOX::Abstract::Vector &result) 
    const;

  virtual NOX::Abstract::Group::ReturnType 
  applyRightPreconditioning(NOX::Parameter::List& params, const NOX::Epetra::Vector& input, 
			    NOX::Epetra::Vector& result) const;

  virtual NOX::Abstract::Group::ReturnType 
  applyRightPreconditioning(NOX::Parameter::List& params, const NOX::Abstract::Vector& input, 
			    NOX::Abstract::Vector& result) const;

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

  /*! \brief Returns true if the value of the Norm of the linear model
  for a full Newton step ||Js + f|| is valid with respect to the
  current solution vector. 
  */

  virtual bool isNormNewtonSolveResidual() const;

  /*! \brief Returns true if an explicitly constructed preconditioner
  exists (i.e. one that is computed and saved for further use in
  multiple calls to applyRightPreconditioner).  
  */
  virtual bool isPreconditioner() const;

  //@}

  /** @name "Get" functions 
   *
   * Note that these function do not check whether or not the vectors
   * are valid. Must use the "Is" functions for that purpose. */
  //@{

  virtual const NOX::Abstract::Vector& getX() const;

  virtual const NOX::Abstract::Vector& getF() const;
  
  virtual double getNormF() const;

  virtual const NOX::Abstract::Vector& getGradient() const;

  virtual const NOX::Abstract::Vector& getNewton() const;

  /*! 
    \brief
    Returns the 2-norm of the residual of the linear model used in the
    Newton solve computation, ||Js+f||.  This does not account for
    line search adjustments to the step length!
  */
  virtual NOX::Abstract::Group::ReturnType 
  getNormLastLinearSolveResidual(double & residual) const;

  //@}

  virtual NOX::Abstract::Group* clone(CopyType type = DeepCopy) const;

  //! Return the SharedJacobian.
  virtual NOX::Epetra::SharedOperator& getSharedJacobian();

  //! Return the SharedPreconditioner.
  virtual NOX::Epetra::SharedOperator& getSharedPreconditioner();

  //! Return the userInterface.
  virtual NOX::Epetra::Interface& getUserInterface();

  //! Sets the diagonal scaling vector(s) used in scaling the linear system.  See NOX::Epetra::Scaling for details on how to specify scaling of the linear system.
  virtual void setLinearSolveScaling(NOX::Epetra::Scaling& scalingObject);

protected:

  //! resets the isValid flags to false
  virtual void resetIsValid();

  /** @name "Linear Solver" functions 
   */
  //@{

  //! Set any required Aztec options.  NOTE: This does not handle all aztec options.  They will be added as needed. 
  virtual void setAztecOptions(const NOX::Parameter::List& params, AztecOO& aztec) const;
  
  //! Checks to make sure that the supplied operators are valid for the requested preconditioning options. 
  virtual bool checkOperatorConsistency();
  
  /*! 
    \brief
    Computes a precoditioning matrix based on the current solution
    vector.  It then sets the correct preconditioning object (Operator
    or Epetra_RowMatrix) in the AztecOO solver object.  The first
    argument tells the method whether or not to enforce explicit
    construction of the preconditioner so that it is saved for reuse.
  */
  virtual bool createPreconditioner(PrecConstructionType c, 
				    NOX::Parameter::List& p) const;
  
  /*! 
    \brief 
    Allocates the objects required for using ifpack preconditioners
    (NOX::Epetra::Group::ifpackGraph and
    NOX::Epetra::Group::ifpackPreconditioner).  This is called from
    NOX::Epetra::Group::computePreconditioner().
  */

  virtual bool createIfpackPreconditioner(NOX::Parameter::List& p) const;

  /*! 
    \brief 
    Deletes all objects associated with the chosen preconditioner.
    This is called during linear solves and when the solution vector
    changes to reset the preconditioner.
  */
  virtual bool destroyPreconditioner() const;

  /*! 
    \brief

    Deletes the AztecOO solver object.  This is called when the
    solution vector for the group is changed.  The preconditioner is
    no longer valid so the solver and preconditioner are destroyed by
    a call to this method.
  */
  virtual bool destroyAztecSolver() const;

  /*! 
    \brief 
    Computes the 2-norm of the residual of the linear model used in
    the Newton solve computation, ||Js+f||.
  */
  virtual bool computeNormNewtonSolveResidual();
  //@}

  //! Returns the type of operator that is passed into the group constructors.
  /*! Uses dynamic casting to identify the underlying object type. */
  virtual OperatorType getOperatorType(const Epetra_Operator& o);

protected:

  //! Printing Utilities object
  const NOX::Utils utils;

  /** @name Vectors */
  //@{
  //! Solution vector pointer.
  NOX::Epetra::Vector* xVectorPtr;
  //! Solution vector.
  NOX::Epetra::Vector& xVector;
  //! Right-hand-side vector (function evaluation).
  NOX::Epetra::Vector* RHSVectorPtr;
  //! Right-hand-side vector pointer (function evaluation).
  NOX::Epetra::Vector& RHSVector;
  //! Gradient vector pointer(steepest descent vector).
  NOX::Epetra::Vector* gradVectorPtr;
  //! Gradient vector (steepest descent vector).
  NOX::Epetra::Vector& gradVector;
  //! Newton direction vector pointer.
  NOX::Epetra::Vector* NewtonVectorPtr;
  //! Newton direction vector.
  NOX::Epetra::Vector& NewtonVector;
  //! An extra temporary vector, only allocated if needed.
  mutable Epetra_Vector* tmpVectorPtr;
  //@}

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
  mutable bool isValidPreconditioner;
  //@}
  
  //! 2-Norm of RHS
  double normRHS;

  //! 2-Norm of the Newton solve residual: ||Js+f||
  double normNewtonSolveResidual;

  /** @name Shared Operators */
  //@{
  //! Pointer to shared Jacobian matrix 
  NOX::Epetra::SharedOperator* sharedJacobianPtr;

  //! Reference to shared Jacobian matrix 
  NOX::Epetra::SharedOperator& sharedJacobian;

  //! Pointer to shared Preconditioning matrix.  
  NOX::Epetra::SharedOperator* sharedPreconditionerPtr;

  /*! 
    \brief 
    Reference to shared Preconditioning matrix.  If a Group
    constructor is used that only supplies a Jacobian (i.e. no
    separate preconditioning matrix is supplied) then this will point
    to the Jacobian! 
  */

  NOX::Epetra::SharedOperator& sharedPreconditioner;
  //@}

  //! Reference to the user supplied interface functions
  NOX::Epetra::Interface& userInterface;

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

   <li> "AztecOO: Jacobian Matrix" - AztecOO preconditioner is used.
   This requires an explicit Jacobian Matrix.  Therefore the Jacobian
   matrix must be an Epetra_RowMatrix or a NOX_FiniteDifference
   object.

   <li> "AztecOO: User RowMatrix" - AztecOO preconditioner is
   used. This requires an explicit Epetra_RowMatrix that is NOT the
   same as the Jacobian.  This was specifically written for
   "Matrix-Free" mode when there is no explicit Matrix for the
   Jacobian.

   <li> "IFPACK: Jacobian Matrix" - IFPACK preconditioner is used.
   This requires an explicit Jacobian Matrix.  Therefore the Jacobian
   matrix must be an Epetra_RowMatrix or a NOX_FiniteDifference
   object.

   <li> "IFPACK: User RowMatrix" - IFPACK preconditioner is used. This
   requires an explicit Epetra_RowMatrix that is NOT the same as the
   Jacobian.  This was specifically written for "Matrix-Free" mode
   when there is no explicit Matrix for the Jacobian.

   <li> "User Supplied Preconditioner" - user supplied routine that
   computes the preconditioner and applies it to a vector.  The user
   must implement the routine as an Epetra_Operator class such that
   the ApplyInverse() function preconditions the incoming vector.

   </ul>
  */ 
  string preconditioner;
  
  //! Aztec solver object used for preconditioning.

  /*! 

    If calling applyRightPreconditioning() and using an AztecOO
    solver, the AztecOO object must be saved.  Previously we created
    and deleted this object during linear solves.  The solver object
    must be stored so that we don't have to recompute the
    preconditioner for each call to applyRightPreconditioning().  We
    only need to recompute the preconditioner when the solution vector
    and it's corresponding preconditioning matrix is updated.

   */
  mutable AztecOO* aztecSolver;

  /*! 
    \brief 
    If using an IFPACK preconditioner, we must store the IFpack graph.
    This is mutable since the applyRightPreconditioner() is a const
    member.
  */
  mutable Ifpack_IlukGraph* ifpackGraph;

  //! If using an IFPACK preconditioner, we must store the IFpack preconditioner.
  mutable Ifpack_CrsRiluk* ifpackPreconditioner;

  //! Scaling object supplied by the user
  NOX::Epetra::Scaling* scaling;

};

} // namespace Epetra
} // namespace NOX


#endif
