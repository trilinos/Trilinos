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
#include "Epetra_Time.h"        // class data element

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

/*! 

\brief Concrete implementation of NOX::Epetra::LinearSolver for AztecOO.

This solver provides the linear algebra services provided
through the AztecOO parallel iterative linear %solver.
 
The NOX::EpetraNew::LinearSystemAztecOO object provides a flexible
and efficient way to interface an Epetra based application code
to the Aztec linear solver.  This class handles construction 
of both the preconditioners and AztecOO solver.  All options are
determined through parameter lists and the basic constructors. 

<B>Constructing a %Linear System</B>

There are four different constructors that can be used.  The difference
between constructors is based on whether the user supplies a 
Jacobian, a preconditioner, neither or both.  

If a Jacobian is not supplied then this object can create an internally constructed Jacobian based on a Finite Difference or Matrif-Free object.  The user can specify which type of object to use by setting the parameter "Jacbian Operator" in the parameter list.  The choices are "Matrix-Free" or "Finite Difference".

The user can supply their own preconditioner as an Epetra_Operator, or they can supply their own matrix (an Epetra_RowMatrix derived object) that can be used by one of the internal preconditioner libraries (currently aztecoo or ifpack).  If they supply their own preconditioner the object must implement the Epetra_Operator::ApplyInverse method.  This is the method called during the linear solve to introduce preconditoning into aztecoo.  If the user supplies a matrix to be used with an internal preconditioner, it must be derived from the Epetra_RowMatrix class and must implement all functionality in the Epetra_RowMatrix.  If a Preconditioner is not supplied, then this object can create an internal preconditioner matrix by finite differencing or it can use the Jacobian operator if the Jacobian derives from the Epetra_RowMatrix class.   The user can specify which type of object to use by setting the parameter "Preconditioner Operator" in the parameter list.  The choices are "Use Jacobian" or "Finite Difference".

The Jacobian and preconditioner each require an interface to update the state of the operator with respect to the solution vector and any other parameters.  There are three interfaces that can be implemented, NOX::EpetraNew::Interface::Required,   NOX::EpetraNew::Interface::Jacobian, and NOX::EpetraNew::Interface::Preconditioner.  

NOX::EpetraNew::Interface::Required supplies the computeF() function so codes can tell NOX what the nonlinear equations are.  This is the minimum requirement to run nox through the epetra interface.  LinearSolverAztecOO requires this in some constructors so that if a Jacobian or preconditoner is not supplied, it will use computeF from the Required interface to estimate the Jacobian or preconditioner via finite differences or directional derivatives.

NOX::EpetraNew::Interface::Jacobian is used for updating a user supplied Jacobian opertor with respect to the solution vector and any other parameters.  It is required only in constructors in which a user supplies a Jacobian operator.

NOX::EpetraNew::Interface::Preconditioner is used for updating a user supplied preconditioner opertor/matrix with respect to the solution vector and any other parameters.  It is required only in constructors in which a user supplies a preconditioner operator.

<B>"Linear %Solver" sublist parameters</B>

A NOX::ParameterList called linearSolverParams is required in the various constructors and during some method calls such as applyJacobianInverse() and applyRightPreconditioning().  Typically, this list is the "Linear %Solver" sublist found in the nox parameter list.  The following parameters can be set in the linear solver sublist and are vaild for the NOX::EpetraNew::LinearSolverAztecOO object:<ul>

<li> "Aztec %Solver" - Determine the iterative technique used in the solve.  The following options are valid:

  <ul> 
  <li> "GMRES" - Restarted generalized minimal residual (default).
  <li> "CG" - Conjugate gradient.
  <li> "CGS" - Conjugate gradient squared.
  <li> "TFQMR" - Transpose-free quasi-minimal reasidual.
  <li> "BiCGStab" - Bi-conjugate gradient with stabilization.
  <li> "LU" - Sparse direct solve (single processor only).
  </ul>

<li> "Size of Krylov Subspace" - When using restarted GMRES this sets
the maximum size of the Krylov subspace (defaults to 300).
 
<li> "Orthogonalization" - The orthogonalization routine used for the
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
 
<li> "Tolerance" - Tolerance used by AztecOO to determine if an iterative linear solve has converged. 

<li> "Ill-Conditioning Threshold" - If the upper hessenberg matrix
during GMRES generates a condition number greater than this parameter
value, aztec will exit the linear solve returning the it's current
solution.  The default is 1.0e11.

<li> "Preconditioner Iterations" - Number of iterations an
AztecOO_Operator should take when solving the preconditioner.  This is
only used if an AztecOO preconditioner is used and the solver makes a
call to NOX::Epetra::Group::applyRightPreconditioning().  This is NOT a
recomended approach.

<li> "Max Iterations" - maximum number of iterations in the 
linear solve.  Default is 400.

<li> "Zero Initial Guess" - boolean.  Zero out the initial guess for linear solves performed through applyJacobianInverse calls (i.e. zero out the result vector before the linear solve).  Defaults to false.

<li> "Output Frequency" - number of linear solve iterations between
output of the linear solve residual. Takes an integer, or one of the
AztecOO flags: AZ_none, AZ_last, or AZ_all as a value. Defaults to
AZ_last.

<li> "Jacobian Operator" - When a constructor does not require a Jacobian operator, the linear system will create a default operator using:

  <ul>
  <li> "Matrix-Free" (default)
  <li> "Finite Difference"
  </ul>


<li> "Preconditioner" - Sets the choice of the preconditioner to use during linear solves. The validity of the choice of preconditioner will depend on the types of operators that are available for the Jacobian and preconditioner.  NOTE: This flag will override any constructor details.  For example, if you supply a preconditioner operator in the constructor, it will not be used if this flag is set to "None".  If you supply an Epetra_Operator for the preconditioner but the "Preconditioner" flag is set to "AztecOO" (this requires an Epetra_RowMatrix for the preconditioner operator), this object will exit with a failure.  The valid options and any requirements on the operator type are listed below:

  <ul> 
  <li> "None" - No preconditioning.  (default)
  <li> "AztecOO" - AztecOO internal preconditioner.  This requires a preconditioner operator that derives from the Epetra_RowMatrix class.  
  <li> "Ifpack" -  Ifpack internal preconditioner.  This requires a preconditioner object that derives from the Epetra_RowMatrix class or it can use a Jacobian if the Jacobian derives from an Epetra_RowMatrix.  
  <li> "User Defined" - The user supplies an Epetra_Operator derived 
  class.  Users must implement at a minimum the ApplyInverse() function of the Epetra_Operator class since preconditioning of vectors is accomplished through calls to this method.  
  </ul>

<li> "Jacobian Operator" - If a constructor is used that does not supply a Jacobian operator, nox will create an internal Jacobian operator.  This flag is ONLY valid in such cases.  This will determine which Operator is used:
  <ul>
  <li> "Matrix-Free" - Create a NOX::EpetraNew::MatrixFree object.
  <li> "Finite Difference" - Create a NOX::EpetraNew::FiniteDifference object.
  </ul>

<li> "Preconditioner Operator" - If a constructor is used that does not supply a preconditioner operator, nox will create an internal preconditioner operator.  This flag is ONLY valid in such cases.  This will determine which Operator is used:
  <ul>
  <li> "Use Jacobian" - Use the Jacobian Operator (it must be an Epetra_RowMatrix derived object).
  <li> "Finite Difference" - Create a NOX::EpetraNew::FiniteDifference object.
  </ul>

<li> "Aztec %Preconditioner" - If the "Preconditioner" flag is set to "AztecOO" then the specific AztecOO preconditioner is specified with this flag.
Currently supported preconditioners and their corresponding
parameters that can be set are shown
below (See the Aztec 2.1 manual for more information):

  <ul>
  <li> "ilu" - ilu preconditioning.  This choice allows the following  additional parameters to be specified:
  <ul>
  <li> "Overlap" - defaults to 0.0
  <li> "Graph Fill" - defaults to 0.0
  </ul> 

  <li> "ilut" - ilut preconditioning. This choice allows the following  additional parameters to be specified:
  <ul>
  <li> "Overlap" - defaults to 0.0
  <li> "Fill Factor" - defaults to 1.0
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

<li> "Reuse Preconditioner" - (boolean) If set to true, the group/linear solver will not refactor the preconditioner even if the group's solution vector changes.  This flag tells the group an dlinear system to turn off control of preconditioner recalculation.  This is a dangerous condition but can really speed up the computations if the user knows what they are doing.  Toggling this flag is left to the user (ideally it should be done through a status test).  Defaults to false.

<li> "RCM Reordering" - Enables RCM reordering in conjunction with
domain decomp incomplete factorization preconditioning.  The following
options are valid:

  <ul> 
  <li> "Disabled" - (default).
  <li> "Enabled" 
  </ul>

<li> "Use Adaptive %Linear %Solve" - Enables the use of AztecOO's AdaptiveIterate() method instead of calling the Iterate() method.  This causes the preconditioning matrix to be modified to make the linear solves easier. AztecOO will attempt to solve the linear system multiple times now and if the solves are failing it will modify the preconditioner and try again.  Boolean value, defaults to false.  NOTE: This only works for internal Aztec preconditioners!  The "Preconditioning" parameter must be set to "AztecOO: Jacobian Matrix" or "AztecOO: User RowMatrix".  (NOTE: This parameter is currently NOT supported)

<li> "Max Adaptive %Solve Iterations" - (integer) Maximum number of attempts that the linear solver will make when trying to solve a linear system. Defaults to 5.  (NOTE: This parameter is currently NOT supported)

<li> "Output %Solver Details" - (boolean) Write the output sublist below to the parameter list after each linear solve. default is true.

</ul>

<B>"Output" sublist</B>

The parameter list passed in during calls to ApplyJacobianInverse() 
will have an "Output" sublist created that contains the following 
parameters if the flag "Output %Solver Details" is set to true:

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
  /*! Either there is no preconditioning or the preconditioner will be 
    used/created internally.  The Jacobian (if derived from an 
    Epetra_RowMatrix class can be used with an internal preconditioner. 
    See the parameter key "Preconditioner Operator" for more details.
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
    parameter "Jacobian Operator" in the parameter list.  See the 
    parameter key "Jacobian Operator" for more details.  Defaults 
    to using a NOX::EpetraNew::MatrixFree object.
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
  
  virtual bool applyJacobianInverse(NOX::Parameter::List& linearSolverParams, 
				    const NOX::Epetra::Vector& input, 
				    NOX::Epetra::Vector& result);

  virtual bool applyRightPreconditioning(bool useTranspose,
				   NOX::Parameter::List& linearSolverParams, 
				   const NOX::Epetra::Vector& input, 
				   NOX::Epetra::Vector& result) const;

  virtual bool createPreconditioner(Epetra_Vector& x, 
				    NOX::Parameter::List& linearSolverParams,
				    bool recomputeGraph) const;
  
  /*! 
    \brief 
    Deletes all objects associated with the chosen preconditioner.
    This is called during linear solves and when the solution vector
    changes to reset the preconditioner.
  */
  virtual bool destroyPreconditioner() const;

  //! Reset the linear solver parameters
  virtual void reset(NOX::Parameter::List& linearSolverParams);

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

  //virtual void getNormLastLinearSolveResidual(double & residual) const;

  //! Returns the total time (sec.) spent in createPreconditioner(). 
  double getTimeCreatePreconditioner() const;

  //! Returns the total time (sec.) spent in applyJacobianInverse(). 
  double getTimeApplyJacobianInverse() const;

protected:

  //! Parse the parameter list and set the corresponding options in the AztecOO solver onject.
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

  // Allow reuse of Preconditioner
  virtual bool checkPreconditionerReuse();

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

  //! Counter for number of times called since reset or construction
  int precQueryCounter;
  
  //! Parameter to determine whether or not to recompute Preconditioner
  int maxAgeOfPrec;

  //! Epetra_Time object
  Epetra_Time timer;

  //! Total time spent in createPreconditioner (sec.).
  mutable double timeCreatePreconditioner;

  //! Total time spent in applyJacobianInverse (sec.).
  mutable double timeApplyJacbianInverse;

};

} // namespace Epetra
} // namespace NOX


#endif
