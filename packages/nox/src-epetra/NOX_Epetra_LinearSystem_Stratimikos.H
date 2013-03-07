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

#ifndef NOX_EPETRA_LINEARSYSTEMSTRATIMIKOS_H
#define NOX_EPETRA_LINEARSYSTEMSTRATIMIKOS_H

#include "NOX_Common.H" 

#ifdef HAVE_NOX_STRATIMIKOS

#include "NOX_Epetra_LinearSystem.H"	// base class
#include "NOX_Epetra_Vector.H"	// class data element
#include "NOX_Utils.H"          // class data element
#include "Teuchos_ParameterList.hpp"
#include "Epetra_Time.h"        // class data element
#ifdef HAVE_NOX_ML_EPETRA
#include "Epetra_Map.h" // This appears to be needed to precede the ml include
#endif
#include "Teuchos_RCP.hpp"

#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"

#include "Thyra_PreconditionerBase.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_LinearOpWithSolveFactoryBase.hpp"
#include "Thyra_PreconditionerFactoryBase.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"

#include "Teuchos_FancyOStream.hpp"

// Forward declares
namespace NOX {
  namespace Epetra {
    class Scaling;
    namespace Interface {
      class Jacobian;
      class Preconditioner;
    }
  }
}
class Epetra_Vector;
class Epetra_Operator;
class Epetra_RowMatrix;

namespace NOX {
//! Improved version of the Epetra support class.
namespace Epetra {

/*! 

\brief Concrete implementation of NOX::Epetra::LinearSolver for Stratimikos.

This solver class provides the linear algebra services 
through the Stratimikos linear solver strategies
package.
 
This class handles construction of both the preconditioners 
and solvers.  All options are determined through parameter 
lists and the basic constructors. 

<B>Constructing a Linear System</B>

There are two constructors for this linear system class.
The first is for supplying a Jacobian operator which, if
it is a Epetra_RowMatrix, can be used to generate a 
preconditioner. 

The second is for the cases when the user wants to supply
an operator other than the Jacobian for use in preconditioning.
This constructor allows for two use cases itself. The default
case is that the preconditioning operator is a matrix that is
an approximation of the Jacobian, which is used internally
to generate a preconditioner (such as by ILU). The second
use case is where the preconditioning operator supplied to
this constructor approximates the inverse of the Jacobian
matrix. This case requires that the bool option 
precIsAlreadyInverted is set as true (optional eighth argument).

As compared to the NOX::Epetra::LinearSystem::AztecOO class,
this implementation has dropped internal support for 
Matrix-Free and Finite Difference Jacobians. However, 
the Jacobian operator supplied to this class can be
constructed using these classes, as is done in the 
1DfemStratimikos examples.

<B>"Stratimikos Linear Solver" sublist parameters</B>

A Teuchos::ParameterList called "Stratimikos Linear Solver" is input to this class
to select all the options. This list itself contains only 2 sublists.

<li> "Stratimikos" sublist - This very extensive list is sent directly 
to the Stratimikos layer to allow selection of numerous linear solver
and preconditioning strategies. This list is documented here:
Stratimikos::DefaultLinearSolverBuilder .

<li> "NOX Stratimikos Options" sublist - This allows the NOX driver code
for Stratimikos to parse some options, such as for preconditioner reuse.

<ul>
<li> "Preconditioner Reuse Policy" - (std::string) Allows the user to set how and when the preconditioner should be computed.  This flag supports native Aztec, Ifpack and ML preconditioners.  There are three options:
    <ul>
    <li>"Rebuild" - The "Rebuild" option always completely destroys and then rebuilds the preconditioner each time a linear solve is requested.
    <li>"Reuse" - The group/linear solver will not recompute the preconditioner even if the group's solution vector changes.  It just blindly reuses what has been constructed.  This turns off control of preconditioner recalculation.  This is a dangerous condition but can really speed up the computations if the user knows what they are doing.  We don't recommend users trying this.
    <li>"Recompute" - Recomputes the preconditioner, but will try to efficiently reuse any objects that don't need to be destroyed.  How efficient the "Recompute" option is depends on the type of preconditioner.  For example if we are using ILU from the Ifpack library, we would like to not destroy and reallocate the graph each solve.  With this option, we tell Ifpack to reuse the graph from last time - e.g the sparisty pattern has not changed between applications of the preconditioner.
    </ul>

<li> "Max Age Of Prec" - (int) If the "Preconditioner Reuse Policy" is set to "Reuse", this integer tells the linear system how many times to reuse the preconditioner before rebuilding it.  Defaults to 1.

<li> "Zero Initial Guess" - boolean.  Zero out the initial guess for linear solves performed through applyJacobianInverse calls (i.e. zero out the result vector before the linear solve).  Defaults to false.


</ul>
 */

class LinearSystemStratimikos : public virtual NOX::Epetra::LinearSystem {

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

  //! Constructor with a user supplied Jacobian Operator.  
  /*! Either there is no preconditioning or the preconditioner will be 
    used/created internally.  The Jacobian (if derived from an 
    Epetra_RowMatrix class) can be used to internally to create a
    preconditioner. The stratimikosLinSolParam list needs to have 
    a sublist called Stratimikos which is sent directly to the 
    Stratimikos package. See Stratimikos for full documentation.
   */

  LinearSystemStratimikos(
    Teuchos::ParameterList& printingParams, 
    Teuchos::ParameterList& stratimikosLinSolParams, 
    const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac, 
    const Teuchos::RCP<Epetra_Operator>& J,
    const NOX::Epetra::Vector& cloneVector,
    const Teuchos::RCP<NOX::Epetra::Scaling> scalingObject = 
    Teuchos::null);

  //! Constructor with user supplied separate objects for the
  //! Jacobian (J) and Preconditioner (M). 
  /*! If precIsAlreadyInverted=false, then M is considered an
      approximate Jacobian. Preconditioner algorithms such as
      ILU are then used to generate an approximate M-inverse.
      Else, if precIsAlreadyInverted=true, then the M operator
      is already assumed to approximate J-inverse and it is
      just applied as-is. 
  
      The stratimikosLinSolParams list needs to have a 
      sublist called Stratimikos which is sent directly to the 
      Stratimikos package. See Stratimikos for full documentation.
  */
  LinearSystemStratimikos(
    Teuchos::ParameterList& printingParams, 
    Teuchos::ParameterList& stratimikosLinSolParams, 
    const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac, 
    const Teuchos::RCP<Epetra_Operator>& J,  
    const Teuchos::RCP<NOX::Epetra::Interface::Preconditioner>& iPrec, 
    const Teuchos::RCP<Epetra_Operator>& M,
    const NOX::Epetra::Vector& cloneVector,
    const bool& precIsAlreadyInverted = false,
    const Teuchos::RCP<NOX::Epetra::Scaling> scalingObject = 
    Teuchos::null);

  //! Destructor.
  virtual ~LinearSystemStratimikos();

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
  
  virtual bool destroyPreconditioner() const;

  /*! \brief recompute preconditioner not supported here. The
      current supported options are rebuild and reuse.
  */
  virtual bool recomputePreconditioner(const NOX::Epetra::Vector& x, 
			     Teuchos::ParameterList& linearSolverParams) const;
  // Derived class
  virtual PreconditionerReusePolicyType 
  getPreconditionerPolicy(bool advanceReuseCounter=true);

  //! Reset the NOX Stratimimos Options, such as preconditioner reuse.
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

  //! Statistics for number of times the linear solver has been called
  int getNumLinearSolves();

  //! Statistics for number of iterations taken in last linear solve 
  int getLinearItersLastSolve();
  
  //! Statistics for cumulative number of iterations in all linear solve
  int getLinearItersTotal();
  
  //! Statistics for the achieved tolerance of last linear solve 
  double getAchievedTol();

protected:

  //! Builds the linear op with solve factory from the Stratimikos sublist
  void initializeStratimikos(Teuchos::ParameterList& stratParams);

  //! Returns the type of operator that is passed into the group constructors.
  /*! Uses dynamic casting to identify the underlying object type. */
  virtual OperatorType getOperatorType(const Epetra_Operator& o);

  virtual void setStratimikosPreconditioner() const;

  /*! \brief Sets the epetra Preconditioner operator in the AztecOO object.

      Makes a call to SetUserOperator.  This must be done AFTER the
      Jacobian matrix/operators is set by setAztecOOJacobian(),
      otherwise the aztec object may ignore this operation.
  */

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

  //! Source of preconditioner: 
  //! SeparateMatrix is approximate Jacobian, to do ILU on, e.g.
  //! User_Defined_ is just an operator 
  enum PreconditionerMatrixSourceType
    {UseJacobian, SeparateMatrix, UserDefined_};
  PreconditionerMatrixSourceType precMatrixSource;

  //! Thyra linear op objects.
  Teuchos::RCP< ::Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory;
  Teuchos::RCP< ::Thyra::LinearOpWithSolveBase<double> > lows;
  Teuchos::RCP<const ::Thyra::LinearOpBase<double> > linearOp;
  mutable Teuchos::RCP< ::Thyra::PreconditionerBase<double> > precObj;

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

  //! If true, the linear solve tolerance comes from NOX (i.e. Inexact Newton)
  bool getLinearSolveToleranceFromNox;

  //! Iterative linear solver statistics
  int linearSolveCount;
  int linearSolveIters_last;
  int linearSolveIters_total;
  double linearSolveAchievedTol;

  //! Output stream for stratimikos
  Teuchos::RCP<Teuchos::FancyOStream> outputStream;

};

} // namespace Epetra
} // namespace NOX
#endif //HAVE_NOX_STRATIMIKOS
#endif
