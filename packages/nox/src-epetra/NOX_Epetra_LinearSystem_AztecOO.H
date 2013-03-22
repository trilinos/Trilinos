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

#ifndef NOX_EPETRA_LINEARSYSTEMAZTECOO_H
#define NOX_EPETRA_LINEARSYSTEMAZTECOO_H

#include "NOX_Epetra_LinearSystem.H"	// base class
#include "NOX_Epetra_Vector.H"	// class data element
#include "NOX_Utils.H"          // class data element
#include "NOX_Common.H"         // class data element (std::string)
#include "Teuchos_ParameterList.hpp"
#include "Epetra_Time.h"        // class data element
#ifdef HAVE_NOX_ML_EPETRA
#include "Epetra_Map.h" // This appears to be needed to precede the ml include
#include "ml_epetra_preconditioner.h"
#endif
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

// Forward declares
namespace Teuchos {
  class ParameterList;
}
namespace NOX {
  namespace Epetra {
    class Scaling;
    namespace Interface {
      class Required;
      class Jacobian;
      class Preconditioner;
    }
  }
}
class Epetra_Vector;
class Epetra_Operator;
class Epetra_RowMatrix;
class AztecOO;
class Ifpack_IlukGraph;
class Ifpack_CrsRiluk;
class Ifpack_Preconditioner;

namespace NOX {
//! Improved version of the Epetra support class.
namespace Epetra {

/*! 

\brief Concrete implementation of NOX::Epetra::LinearSolver for AztecOO.

This solver provides the linear algebra services provided
through the AztecOO parallel iterative linear solver.
 
The NOX::Epetra::LinearSystemAztecOO object provides a flexible
and efficient way to interface an Epetra based application code
to the Aztec linear solver.  This class handles construction 
of both the preconditioners and AztecOO solver.  All options are
determined through parameter lists and the basic constructors. 

<B>Constructing a Linear System</B>

There are four different constructors that can be used.  The difference
between constructors is based on whether the user supplies a 
Jacobian, a preconditioner, neither or both.  

If a Jacobian is not supplied then this object can create an
internally constructed Jacobian based on a Finite Difference or
Matrif-Free object.  The user can specify which type of object to use
by setting the parameter "Jacobian Operator" in the parameter list.
The choices are "Matrix-Free" or "Finite Difference".

The user can supply their own preconditioner as an Epetra_Operator, or they can supply their own matrix (an Epetra_RowMatrix derived object) that can be used by one of the internal preconditioner libraries (currently aztecoo or ifpack).  If they supply their own preconditioner the object must implement the Epetra_Operator::ApplyInverse method.  This is the method called during the linear solve to introduce preconditoning into aztecoo.  If the user supplies a matrix to be used with an internal preconditioner, it must be derived from the Epetra_RowMatrix class and must implement all functionality in the Epetra_RowMatrix.  If a Preconditioner is not supplied, then this object can create an internal preconditioner matrix by finite differencing or it can use the Jacobian operator if the Jacobian derives from the Epetra_RowMatrix class.   The user can specify which type of object to use by setting the parameter "Preconditioner Operator" in the parameter list.  The choices are "Use Jacobian" or "Finite Difference".

The Jacobian and preconditioner each require an interface to update the state of the operator with respect to the solution vector and any other parameters.  There are three interfaces that can be implemented, NOX::Epetra::Interface::Required,   NOX::Epetra::Interface::Jacobian, and NOX::Epetra::Interface::Preconditioner.  

NOX::Epetra::Interface::Required supplies the computeF() function so codes can tell NOX what the nonlinear equations are.  This is the minimum requirement to run nox through the epetra interface.  LinearSolverAztecOO requires this in some constructors so that if a Jacobian or preconditoner is not supplied, it will use computeF from the Required interface to estimate the Jacobian or preconditioner via finite differences or directional derivatives.

NOX::Epetra::Interface::Jacobian is used for updating a user supplied Jacobian opertor with respect to the solution vector and any other parameters.  It is required only in constructors in which a user supplies a Jacobian operator.

NOX::Epetra::Interface::Preconditioner is used for updating a user supplied preconditioner opertor/matrix with respect to the solution vector and any other parameters.  It is required only in constructors in which a user supplies a preconditioner operator.

<B>"Linear Solver" sublist parameters</B>

A Teuchos::ParameterList called linearSolverParams is required in the various constructors and during some method calls such as applyJacobianInverse() and applyRightPreconditioning().  Typically, this list is the "Linear Solver" sublist found in the nox parameter list.  The following parameters can be set in the linear solver sublist and are vaild for the NOX::Epetra::LinearSolverAztecOO object:<ul>

<li> "Aztec Solver" - Determine the iterative technique used in the solve.  The following options are valid:

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

<li> "Throw Error on Prec Failure" - boolean.  If set to true, an exception will be thrown if the preconditioner fails to initialize or recompute/refactor.  If set to false, a warning will br printed if the NOX::Utils::Warning is enabled in the printing utilities (NOX::Utils). Defaults to true.

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
  <li> "Ifpack" -  Ifpack internal preconditioner.  This requires a preconditioner object that derives from the Epetra_RowMatrix class or it can use a Jacobian if the Jacobian derives from an Epetra_RowMatrix. This option is deprecated.  Please use "New Ifpack". 
  <li> "New Ifpack" -  Ifpack internal preconditioner.  This requires a preconditioner object that derives from the Epetra_RowMatrix class or it can use a Jacobian if the Jacobian derives from an Epetra_RowMatrix.  
  <li> "User Defined" - The user supplies an Epetra_Operator derived 
  class.  Users must implement at a minimum the ApplyInverse() function of the Epetra_Operator class since preconditioning of vectors is accomplished through calls to this method.  
  </ul>

<li> "Jacobian Operator" - If a constructor is used that does not supply a Jacobian operator, nox will create an internal Jacobian operator.  This flag is ONLY valid in such cases.  This will determine which Operator is used:
  <ul>
  <li> "Matrix-Free" - Create a NOX::Epetra::MatrixFree object.
  <li> "Finite Difference" - Create a NOX::Epetra::FiniteDifference object.
  </ul>

<li> "Preconditioner Operator" - If a constructor is used that does not supply a preconditioner operator, nox will create an internal preconditioner operator.  This flag is ONLY valid in such cases.  This will determine which Operator is used:
  <ul>
  <li> "Use Jacobian" - Use the Jacobian Operator (it must be an Epetra_RowMatrix derived object).
  <li> "Finite Difference" - Create a NOX::Epetra::FiniteDifference object.
  </ul>

<li> "Aztec Preconditioner" - If the "Preconditioner" flag is set to "AztecOO" then the specific AztecOO preconditioner is specified with this flag.
Currently supported preconditioners and their corresponding
parameters that can be set are shown
below (See the Aztec 2.1 manual for more information):

  <ul>
  <li> "ilu" - ilu preconditioning.  This choice allows the following  additional parameters to be specified:
  <ul>
  <li> "Overlap" - defaults to 0
  <li> "Graph Fill" - defaults to 0
  </ul> 

  <li> "ilut" - ilut preconditioning. This choice allows the following  additional parameters to be specified:
  <ul>
  <li> "Overlap" - defaults to 0
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

<li> "Ifpack" - If the "Preconditioner" flag is set to "New Ifpack" then any of the options supported by the Ifpack Create factory can be specified using a Teuchos::ParameterList containing the Ifpack options and then setting this as a parameter named "Ifpack" in the "Linear Solver" sublist.

<li> "ML" - If the "Preconditioner" flag is set to "ML" then any of the options supported by the ML factory can be specified using a Teuchos::ParameterList containing the ML options and then setting this as a parameter named "ML" in the "Linear Solver" sublist.

<li> "Preconditioner Reuse Policy" - (std::string) Allows the user to set how and when the preconditioner should be computed.  This flag supports native Aztec, Ifpack and ML preconditioners.  There are three options:
    <ul>
    <li>"Rebuild" - The "Rebuild" option always completely destroys and then rebuilds the preconditioner each time a linear solve is requested.
    <li>"Reuse" - The group/linear solver will not recompute the preconditioner even if the group's solution vector changes.  It just blindly reuses what has been constructed.  This turns off control of preconditioner recalculation.  This is a dangerous condition but can really speed up the computations if the user knows what they are doing.  We don't recommend users trying this.
    <li>"Recompute" - Recomputes the preconditioner, but will try to efficiently reuse any objects that don't need to be destroyed.  How efficient the "Recompute" option is depends on the type of preconditioner.  For example if we are using ILU from the Ifpack library, we would like to not destroy and reallocate the graph each solve.  With this option, we tell Ifpack to reuse the graph from last time - e.g the sparisty pattern has not changed between applications of the preconditioner. 
    </ul>

<li> "Max Age Of Prec" - (int) If the "Preconditioner Reuse Policy" is set to "Reuse", this integer tells the linear system how many times to reuse the preconditioner before rebuilding it.  Defaults to 1.

<li> "RCM Reordering" - Enables RCM reordering in conjunction with
domain decomp incomplete factorization preconditioning.  The following
options are valid:

  <ul> 
  <li> "Disabled" - (default).
  <li> "Enabled" 
  </ul>

<li> "Use Adaptive Linear Solve" - Enables the use of AztecOO's AdaptiveIterate() method instead of calling the Iterate() method.  This causes the preconditioning matrix to be modified to make the linear solves easier. AztecOO will attempt to solve the linear system multiple times now and if the solves are failing it will modify the preconditioner and try again.  Boolean value, defaults to false.  NOTE: This only works for internal Aztec preconditioners!  The "Preconditioning" parameter must be set to "AztecOO: Jacobian Matrix" or "AztecOO: User RowMatrix".  (NOTE: This parameter is currently NOT supported)

<li> "Max Adaptive Solve Iterations" - (integer) Maximum number of attempts that the linear solver will make when trying to solve a linear system. Defaults to 5.  (NOTE: This parameter is currently NOT supported)

<li> "Compute Scaling Manually" - (boolean) The linear system can be scaled if a NOX::Epetra::Scaling object is supplied to LinearSystemAztecOO.  When to compute the scaling can be handled either manually by the user, or this object can automatically compute the scaling prior to a linear solve.  By setting this flag to true, the user will call NOX::Epetra::Scaling::computeScaling() manually - on their own!  Setting this to false means the LinearSystemAztecOO object will call the computeScaling function right before it applies the scaling to the matrix in the applyJacobianInverse function.  Default is true (user will call compute scaling). 

<li> "Output Solver Details" - (boolean) Write the output sublist below to the parameter list after each linear solve. default is true.

<li> "Write Linear System" - (boolean) If set to true, the linear system (Epetra_Map, Jacobian, LHS and RHS) is printed to a set of files in matrix market format.  This option requires building nox with the flag --enable-nox-debug and building the EpetraExt library.

<li> "Write Linear System File Prefix" - (std::string) If writing of the linear system is enabled (see above parameter) users can change the name of the output file prefix.  The default is "NOX_LinSys".  This option requires building nox with the flag --enable-nox-debug and building the EpetraExt library.

</ul>

<B>"Output" sublist</B>

The parameter list passed in during calls to ApplyJacobianInverse() 
will have an "Output" sublist created that contains the following 
parameters if the flag "Output Solver Details" is set to true:

<ul>

<li> "Acheived Tolerance" - Actual tolerance achieved by the linear
solver computed via the convergence test requested.

<li> "Number of Linear Iterations" - Number of iterations used by the
linear solver in the last call to applyJacobianInverse

<li> "Total Number of Linear Iterations" - Total number of linear
solve iterations performed by groups that have used this input list

</ul>

 */

class LinearSystemAztecOO : public virtual NOX::Epetra::LinearSystem {

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
    NOX::Epetra::MatrixFree object.
   */
  LinearSystemAztecOO(
    Teuchos::ParameterList& printingParams, 
    Teuchos::ParameterList& linearSolverParams,
    const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq, 
    const NOX::Epetra::Vector& cloneVector,
    const Teuchos::RCP<NOX::Epetra::Scaling> scalingObject = 
    Teuchos::null);

  //! Constructor with a user supplied Jacobian Operator only.  
  /*! Either there is no preconditioning or the preconditioner will be 
    used/created internally.  The Jacobian (if derived from an 
    Epetra_RowMatrix class can be used with an internal preconditioner. 
    See the parameter key "Preconditioner Operator" for more details.
   */
  LinearSystemAztecOO(
    Teuchos::ParameterList& printingParams, 
    Teuchos::ParameterList& linearSolverParams, 
    const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq, 
    const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac, 
    const Teuchos::RCP<Epetra_Operator>& J,
    const NOX::Epetra::Vector& cloneVector,
    const Teuchos::RCP<NOX::Epetra::Scaling> scalingObject = 
    Teuchos::null);

  //! Constructor with a user supplied Preconditioner Operator only.  
  /*! Jacobian operator will be constructed internally based on the
    parameter "Jacobian Operator" in the parameter list.  See the 
    parameter key "Jacobian Operator" for more details.  Defaults 
    to using a NOX::Epetra::MatrixFree object.
   */
  LinearSystemAztecOO(
    Teuchos::ParameterList& printingParams, 
    Teuchos::ParameterList& linearSolverParams, 
    const Teuchos::RCP<NOX::Epetra::Interface::Required>& i, 
    const Teuchos::RCP<NOX::Epetra::Interface::Preconditioner>& iPrec,
    const Teuchos::RCP<Epetra_Operator>& M,
    const NOX::Epetra::Vector& cloneVector,
    const Teuchos::RCP<NOX::Epetra::Scaling> scalingObject = 
    Teuchos::null);
  
  //! Constructor with user supplied separate objects for the
  //! Jacobian (J) and Preconditioner (M).  
  //! linearSolverParams is the "Linear Solver" sublist of parameter list.
  LinearSystemAztecOO(
    Teuchos::ParameterList& printingParams, 
    Teuchos::ParameterList& linearSolverParams, 
    const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac, 
    const Teuchos::RCP<Epetra_Operator>& J,  
    const Teuchos::RCP<NOX::Epetra::Interface::Preconditioner>& iPrec, 
    const Teuchos::RCP<Epetra_Operator>& M,
    const NOX::Epetra::Vector& cloneVector,
    const Teuchos::RCP<NOX::Epetra::Scaling> scalingObject = 
    Teuchos::null);

  //! Destructor.
  virtual ~LinearSystemAztecOO();

  virtual bool applyJacobian(const NOX::Epetra::Vector& input, 
			     NOX::Epetra::Vector& result) const;

  virtual bool applyJacobianTranspose(const NOX::Epetra::Vector& input, 
				      NOX::Epetra::Vector& result) const;
  
  virtual bool applyJacobianInverse(Teuchos::ParameterList& linearSolverParams, 
				    const NOX::Epetra::Vector& input, 
				    NOX::Epetra::Vector& result);

  virtual bool applyRightPreconditioning(bool useTranspose,
				   Teuchos::ParameterList& linearSolverParams, 
				   const NOX::Epetra::Vector& input, 
				   NOX::Epetra::Vector& result) const;

  virtual bool createPreconditioner(const NOX::Epetra::Vector& x, 
				    Teuchos::ParameterList& linearSolverParams,
				    bool recomputeGraph) const;
  
  /*! 
    \brief 
    Deletes all objects associated with the chosen preconditioner.
    This is called during linear solves and when the solution vector
    changes to reset the preconditioner.
  */
  virtual bool destroyPreconditioner() const;

  /*! \brief Recalculates the preconditioner using an already allocated graph.
    
    Use this to compute a new preconditioner while using the same
    graph for the preconditioner.  This avoids deleting and
    reallocating the memory required for the preconditioner and
    results in a big speed-up for large-scale jobs.
  */
  virtual bool recomputePreconditioner(const NOX::Epetra::Vector& x, 
			     Teuchos::ParameterList& linearSolverParams) const;
  // Derived class
  virtual PreconditionerReusePolicyType 
  getPreconditionerPolicy(bool advanceReuseCounter=true);

  //! Reset the linear solver parameters.
  virtual void reset(Teuchos::ParameterList& linearSolverParams);

  //! Get the scaling object.
  virtual Teuchos::RCP<NOX::Epetra::Scaling> getScaling();

  //! Sets the diagonal scaling vector(s) used in scaling the linear system.  See NOX::Epetra::Scaling for details on how to specify scaling of the linear system.  
  void resetScaling(const Teuchos::RCP<NOX::Epetra::Scaling>& s);

  //! Compute the Jacobian 
  virtual bool computeJacobian(const NOX::Epetra::Vector& x);

  //! NOX::Interface::Jacobian accessor
  virtual Teuchos::RCP<const NOX::Epetra::Interface::Jacobian> 
  getJacobianInterface() const;

  //! NOX::Interface::Preconditioiner accessor
  virtual Teuchos::RCP<const NOX::Epetra::Interface::Preconditioner> 
  getPrecInterface() const;

  //! Indicates whether a preconditioner has been constructed
  virtual bool isPreconditionerConstructed() const;

  //! Indicates whether the linear system has a preconditioner
  virtual bool hasPreconditioner() const;

  //! Jacobian Epetra_Operator accessor
  virtual Teuchos::RCP<const Epetra_Operator> getJacobianOperator() const;

  //! Jacobian Epetra_Operator accessor
  virtual Teuchos::RCP<Epetra_Operator> getJacobianOperator();

  //! Preconditioner Epetra_Operator accessor (only the base matrix if using an internal preconditioner - aztecoo or ifpack). 
  virtual Teuchos::RCP<const Epetra_Operator> getPrecOperator() const;

  //! Return preconditioner operator generated and stored in AztecOO
  /*!
   * Note:  This should only be called if hasPreconditioner() returns true.
   */
  virtual Teuchos::RCP<const Epetra_Operator> 
  getGeneratedPrecOperator() const;

  //! Return preconditioner operator generated and stored in AztecOO
  virtual Teuchos::RCP<Epetra_Operator> getGeneratedPrecOperator();

  //virtual void getNormLastLinearSolveResidual(double & residual) const;

  //! Returns the total time (sec.) spent in createPreconditioner(). 
  double getTimeCreatePreconditioner() const;

  //! Returns the total time (sec.) spent in applyJacobianInverse(). 
  double getTimeApplyJacobianInverse() const;

  //! Set Jacobian operator for solve
  virtual void setJacobianOperatorForSolve(const Teuchos::RCP<const Epetra_Operator>& solveJacOp);

  //! Set preconditioner operator for solve
  /*!
   * Note:  This should only be called if hasPreconditioner() returns true.
   */
  virtual void setPrecOperatorForSolve(const Teuchos::RCP<const Epetra_Operator>& solvePrecOp);

protected:

  //! Parse the parameter list and set the corresponding options in the AztecOO solver onject.
  virtual void setAztecOptions(Teuchos::ParameterList& lsParams, 
			       AztecOO& aztec) const;
  
  //! Creates an internally owned Epetra_Operator for the Jacobian. 
  virtual bool createJacobianOperator(
       Teuchos::ParameterList& printParams,
       Teuchos::ParameterList& lsParams,
       const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq, 
       const NOX::Epetra::Vector& cloneVector); 

  //! Creates an internally owned Epetra_Operator for the Preconditioner. 
  virtual bool createPrecOperator(
       Teuchos::ParameterList& printParams,
       Teuchos::ParameterList& lsParams,
       const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq, 
       const NOX::Epetra::Vector& cloneVector); 

  //! Checks to make sure that the supplied operators are valid for the requested preconditioning options set in the parameter list. 
  virtual bool checkPreconditionerValidity();

  /*! 
    \brief 
    Allocates the objects required for using ifpack preconditioners
    (NOX::Epetra::Group::ifpackGraph and
    NOX::Epetra::Group::ifpackPreconditioner).  This is called from
    NOX::Epetra::Group::computePreconditioner().
  */
  virtual bool createIfpackPreconditioner(Teuchos::ParameterList& p) const;

  /*! 
    \brief 
    Allocates the objects required for using the new version of ifpack preconditioners via the Ifpack Create factory.
    This is called from NOX::Epetra::Group::computePreconditioner().
  */
  virtual bool createNewIfpackPreconditioner(Teuchos::ParameterList& p) const;

#ifdef HAVE_NOX_ML_EPETRA
  /*! 
    \brief 
    Allocates the objects required for using ML preconditioners
  */
  virtual bool createMLPreconditioner(Teuchos::ParameterList& p) const;
#endif

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

  /*! \brief Sets the epetra Jacobian operator in the AztecOO object.

      Makes a call to SetUserMatrix or SetUserOperator to set the Jacobian.
  */
  virtual void setAztecOOJacobian() const;

  /*! \brief Sets the epetra Preconditioner operator in the AztecOO object.

      Makes a call to SetUserOperator.  This must be done AFTER the
      Jacobian matrix/operators is set by setAztecOOJacobian(),
      otherwise the aztec object may ignore this operation.
  */
  virtual void setAztecOOPreconditioner() const;

  virtual void throwError(const std::string& functionName, 
			  const std::string& errorMsg) const;

  //! Prints a warning for ifpack preconditioner failures (error_code != 0).
  virtual void precError(int error_code, 
			 const std::string& nox_function,
			 const std::string& prec_type,
			 const std::string& prec_function) const;

protected:

  //! Printing Utilities object
  NOX::Utils utils;

  //! Reference to the user supplied Jacobian interface functions
  Teuchos::RCP<NOX::Epetra::Interface::Jacobian> jacInterfacePtr;

  //! Type of operator for the Jacobian.
  OperatorType jacType;
  
  //! Pointer to the Jacobian operator.
  mutable Teuchos::RCP<Epetra_Operator> jacPtr;

  //! Reference to the user supplied preconditioner interface functions
  Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> precInterfacePtr;

  //! Type of operator for the preconditioner. 
  OperatorType precType;

  //! Pointer to the preconditioner operator.
  mutable Teuchos::RCP<Epetra_Operator> precPtr;

  //! Source of the RowMatrix if using an AztecOO native preconditioner 
  enum PreconditionerMatrixSourceType {UseJacobian, 
				       SeparateMatrix};
  PreconditionerMatrixSourceType precMatrixSource;

  enum PreconditionerType {None_, AztecOO_, Ifpack_, NewIfpack_, ML_, UserDefined_};
  PreconditionerType precAlgorithm;
  
  //! Aztec solver object.
  mutable Teuchos::RCP<AztecOO> aztecSolverPtr;

  /*! 
    \brief Stores the ifpack graph. 
    Mutable since the applyRightPreconditioner() is a const method.
  */
  mutable Teuchos::RCP<Ifpack_IlukGraph> ifpackGraphPtr;

  /*! \brief Stores the ifpack preconditioner.
    Mutable since the applyRightPreconditioner() is a const method.
  */
  mutable Teuchos::RCP<Ifpack_CrsRiluk> 
  ifpackPreconditionerPtr;

  /*! \brief Stores the new ifpack preconditioner.
    Mutable since the applyRightPreconditioner() is a const method.
  */
  mutable Teuchos::RCP<Ifpack_Preconditioner> 
  newIfpackPreconditionerPtr;

#ifdef HAVE_NOX_ML_EPETRA
  /*! \brief Stores the ml preconditioner.
    Mutable since the applyRightPreconditioner() is a const method.
  */
  mutable Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> 
  MLPreconditionerPtr;
#endif

  //! Scaling object supplied by the user
  Teuchos::RCP<NOX::Epetra::Scaling> scaling;

  //! An extra temporary vector, only allocated if needed.
  mutable Teuchos::RCP<NOX::Epetra::Vector> tmpVectorPtr;

  mutable double conditionNumberEstimate;

  //! True if the preconditioner has been computed 
  mutable bool isPrecConstructed;

  //! If set to true, solver information is printed to the "Output" sublist of the "Linear Solver" list.
  bool outputSolveDetails;

  //! Zero out the initial guess for linear solves performed through applyJacobianInverse calls (i.e. zero out the result vector before the linear solve).
  bool zeroInitialGuess;

  //! Stores the parameter "Compute Scaling Manually".
  bool manualScaling;

  //! Policy for how to handle the preconditioner between nonlineaer iterations.
  PreconditionerReusePolicyType precReusePolicy;

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

  //! Preconditioner operator that will be used in solves.
  mutable Teuchos::RCP<Epetra_Operator> solvePrecOpPtr;

  //! If true, any preconditioner error will cause a throw instead of a warning.
  bool throwErrorOnPrecFailure;

// #ifdef HAVE_NOX_DEBUG
#ifdef HAVE_NOX_EPETRAEXT
  mutable int linearSolveCount;
#endif
// #endif

};

} // namespace Epetra
} // namespace NOX
#endif
