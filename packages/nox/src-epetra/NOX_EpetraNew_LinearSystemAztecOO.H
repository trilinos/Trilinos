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

#ifndef NOX_EPETRANEW_LINEARSYSTEMAZTECOO_H
#define NOX_EPETRANEW_LINEARSYSTEMAZTECOO_H

#include "NOX_EpetraNew_LinearSystem.H"	// base class
#include "NOX_Epetra_Vector.H"	// class data element
#include "NOX_Utils.H"          // class data element
#include "NOX_Common.H"         // class data element (string)

// Forward declares
namespace NOX {
  namespace EpetraNew {
    class Scaling;
    namespace Interface {
      class Required;
      class Jacobian;
      class Preconditioner;
    }
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
//! Improved version of the Epetra support class.
namespace EpetraNew {

//! Concrete implementation of NOX::Epetra::LinearSolver for AztecOO.
/*! 

This group is set up to use the linear algebra services provided
through the Trilinos/Epetra package with AztecOO for the linear
%solver.
 
<B>Supplying a Jacobian %Operator</B>


if using an internally ctor prec - must be a rowmatrix.....


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

<li> "Convergence Test" - Algorithm used to calculate the residual
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

<li> "Jacobian Operator" - When a constructor does not require a Jacobian operator, the linear system withh create a default operator using:
  <ul>
  <li> "Matrix-Free" (default)
  <li> "Finite Difference"
  </ul>


<li> "Preconditioner" - Preconditioning of the linear system also has
multiple options that may depend on the types of Jacobian and 
preconditioner operators supplied.The valid options are listed
below:

  <ul> 
  <li> "None" - no preconditioning.  This is the default setting.

  <li> "AztecOO" - An AztecOO preconditioner will be
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

<li> "Zero Initial Guess" - boolean.  Zero out the initial guess for linear solves performed through applyJacobianInverse calls (i.e. zero out the result vector before the linear solve).  Defaults to false.

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

class LinearSystemAztecOO : public virtual NOX::EpetraNew::LinearSystem {

protected:

  //! List of types of epetra objects that can be used for the Jacobian and/or Preconditioner.
  enum OperatorType {
    //! An Epetra_Operator derived object.
    EpetraOperator, 
    //! An Epetra_RowMatrix derived object.
    EpetraRowMatrix,
    //! An Epetra_VbrMatrix object.
    EpetraVbrMatrix,
    //! An Epetra_CrsMatrix object.
    EpetraCrsMatrix
  };

public:
  //! Constructor with no Operators.  
  /*! Jacobian Operator will be constructed internally based on the
    parameter "Jacobian Operator".  Defaults to using a 
    NOX::EpetraNew::MatrixFree object.
   */
  LinearSystemAztecOO(NOX::Parameter::List& printingParams, 
		      NOX::Parameter::List& linearSolverParams,
		      NOX::EpetraNew::Interface::Required& iReq, 
		      Epetra_Vector& cloneVector,
		      NOX::EpetraNew::Scaling* scalingObject = 0);

  //! Constructor with a user supplied Jacobian Operator only.  
  /*! Either there is no preconditioning or the Jacobian will be used
    with an internal preconditioner.  linearSolverParams is 
    the "Linear Solver" sublist of parameter NOX list.
   */
  LinearSystemAztecOO(NOX::Parameter::List& printingParams, 
		      NOX::Parameter::List& linearSolverParams, 
		      NOX::EpetraNew::Interface::Required& iReq, 
		      NOX::EpetraNew::Interface::Jacobian& iJac, 
		      Epetra_Operator& J,
		      Epetra_Vector& cloneVector,
		      NOX::EpetraNew::Scaling* scalingObject = 0);

  //! Constructor with a user supplied Preconditioner Operator only.  
  /*! Jacobian operator will be constructed internally based on the
    parameter "Jacobian Operator".  Defaults to using a 
    NOX::EpetraNew::MatrixFree object.
   */
  LinearSystemAztecOO(NOX::Parameter::List& printingParams, 
		      NOX::Parameter::List& linearSolverParams, 
		      NOX::EpetraNew::Interface::Required& i, 
		      NOX::EpetraNew::Interface::Preconditioner& iPrec, 
		      Epetra_Operator& M,
		      Epetra_Vector& cloneVector,
		      NOX::EpetraNew::Scaling* scalingObject = 0);

  //! Constructor with user supplied separate objects for the
  //! Jacobian (J) and %Preconditioner (M).  
  //! linearSolverParams is the "Linear %Solver" sublist of parameter list.
  LinearSystemAztecOO(NOX::Parameter::List& printingParams, 
		      NOX::Parameter::List& linearSolverParams, 
		      NOX::EpetraNew::Interface::Jacobian& iJac, 
		      Epetra_Operator& J,  
		      NOX::EpetraNew::Interface::Preconditioner& iPrec, 
		      Epetra_Operator& M,
		      Epetra_Vector& cloneVector,
		      NOX::EpetraNew::Scaling* scalingObject = 0);

  //! Destructor.
  virtual ~LinearSystemAztecOO();

  virtual bool applyJacobian(const NOX::Epetra::Vector& input, 
			     NOX::Epetra::Vector& result) const;

  virtual bool applyJacobianTranspose(const NOX::Epetra::Vector& input, 
				      NOX::Epetra::Vector& result) const;
  
  virtual bool applyJacobianInverse(NOX::Parameter::List &params, 
				    const NOX::Epetra::Vector &input, 
				    NOX::Epetra::Vector &result) const;

  virtual bool applyRightPreconditioning(bool useTranspose,
					 NOX::Parameter::List& params, 
					 const NOX::Epetra::Vector& input, 
					 NOX::Epetra::Vector& result) const;

  //virtual void getNormLastLinearSolveResidual(double & residual) const;

  /*! 
    \brief
    Computes a precoditioning matrix based on the current solution
    vector.  It then sets the correct preconditioning object (Operator
    or Epetra_RowMatrix) in the AztecOO solver object.  The first
    argument tells the method whether or not to enforce explicit
    construction of the preconditioner so that it is saved for reuse.
  */
  virtual bool createPreconditioner(Epetra_Vector& x, 
				    NOX::Parameter::List& p,
				    bool recomputeGraph) const;
  
  /*! 
    \brief 
    Deletes all objects associated with the chosen preconditioner.
    This is called during linear solves and when the solution vector
    changes to reset the preconditioner.
  */
  virtual bool destroyPreconditioner() const;

  //! Sets the diagonal scaling vector(s) used in scaling the linear system.  See NOX::EpetraNew::Scaling for details on how to specify scaling of the linear system.  
  void resetScaling(NOX::EpetraNew::Scaling& scalingObject);

  //! Compute the Jacobian 
  virtual bool computeJacobian(Epetra_Vector& x);

  //! NOX::Interface::Jacobian accessor
  virtual const NOX::EpetraNew::Interface::Jacobian& 
  getJacobianInterface() const;

  //! NOX::Interface::Preconditioiner accessor
  virtual const NOX::EpetraNew::Interface::Preconditioner& 
  getPrecInterface() const;

  //! Jacobian Epetra_Operator accessor
  virtual const Epetra_Operator& getJacobianOperator() const;

  //! Preconditioner Epetra_Operator accessor (only the base matrix if using an internal preconditioner - aztecoo or ifpack). 
  virtual const Epetra_Operator& getPrecOperator() const;

protected:

  //! Reset the linear solver parameters
  virtual void reset(NOX::Parameter::List& linearSolverParams);

  //! Set any required Aztec options.
  virtual void setAztecOptions(const NOX::Parameter::List& lsParams, 
			       AztecOO& aztec) const;
  
  //! Creates an internally owned Epetra_Operator for the Jacobian. 
  virtual bool createJacobianOperator(NOX::Parameter::List& lsParams,
				   NOX::EpetraNew::Interface::Required& iReq, 
				   const Epetra_Vector& cloneVector); 

  //! Creates an internally owned Epetra_Operator for the Preconditioner. 
  virtual bool createPrecOperator(NOX::Parameter::List& lsParams,
				  NOX::EpetraNew::Interface::Required& iReq, 
				  const Epetra_Vector& cloneVector); 

  //! Checks to make sure that the supplied operators are valid for the requested preconditioning options set in the parameter list. 
  virtual bool checkPreconditionerValidity();
  
  /*! 
    \brief 
    Allocates the objects required for using ifpack preconditioners
    (NOX::EpetraNew::Group::ifpackGraph and
    NOX::EpetraNew::Group::ifpackPreconditioner).  This is called from
    NOX::EpetraNew::Group::computePreconditioner().
  */
  virtual bool createIfpackPreconditioner(NOX::Parameter::List& p) const;

  /*! 
    \brief

    Deletes the AztecOO solver object.  This is called when the
    solution vector for the group is changed.  The preconditioner is
    no longer valid so the solver and preconditioner are destroyed by
    a call to this method.
  */
  //virtual bool destroyAztecSolver() const;

  //! Returns the type of operator that is passed into the group constructors.
  /*! Uses dynamic casting to identify the underlying object type. */
  virtual OperatorType getOperatorType(const Epetra_Operator& o);

  virtual void throwError(const string& functionName, 
			  const string& errorMsg) const;

protected:

  //! Printing Utilities object
  const NOX::Utils utils;

  //! Reference to the user supplied Jacobian interface functions
  NOX::EpetraNew::Interface::Jacobian* jacInterfacePtr;

  //! Type of operator for the Jacobian.
  OperatorType jacType;
  
  //! Poniter to the Jacobian operator.
  mutable Epetra_Operator* jacPtr;

  //! True if this object created and owns the jacobian operator.
  bool ownsJacOperator;

  //! Reference to the user supplied preconditioner interface functions
  NOX::EpetraNew::Interface::Preconditioner* precInterfacePtr;

  //! Type of operator for the preconditioner. 
  OperatorType precType;

  //! Pointer to the preconditioner operator.
  mutable Epetra_Operator* precPtr;

  //! True if this object created and owns the preconditioner operator.
  bool ownsPrecOperator;

  //! Source of the RowMatrix if using an AztecOO native preconditioner 
  enum PreconditionerMatrixSourceType {UseJacobian, 
				       SeparateMatrix};
  PreconditionerMatrixSourceType precMatrixSource;

  enum PreconditionerType {None_, AztecOO_, Ifpack_, UserDefined_};
  PreconditionerType precAlgorithm;
  
  //! Aztec solver object.
  mutable AztecOO* aztecSolverPtr;
  mutable AztecOO& aztecSolver;

  /*! 
    \brief Stores the ifpack graph. 
    Mutable since the applyRightPreconditioner() is a const method.
  */
  mutable Ifpack_IlukGraph* ifpackGraphPtr;

  /*! \brief Stores the ifpack pareconditioner.
    Mutable since the applyRightPreconditioner() is a const method.
  */
  mutable Ifpack_CrsRiluk* ifpackPreconditionerPtr;

  //! Scaling object supplied by the user
  NOX::EpetraNew::Scaling* scaling;

  //! An extra temporary vector, only allocated if needed.
  mutable Epetra_Vector* tmpVectorPtr;

  mutable double conditionNumberEstimate;

  //! True if the preconditioner has been computed 
  mutable bool isPrecConstructed;

  //! If set to true, solver information is printed to the "Output" sublist of the "Linear Solver" list.
  bool outputSolveDetails;

  //! Zero out the initial guess for linear solves performed through applyJacobianInverse calls (i.e. zero out the result vector before the linear solve).
  bool zeroInitialGuess;

};

} // namespace Epetra
} // namespace NOX


#endif
