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


// Forward declares
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

namespace NOX {
//! Improved version of the Epetra support class.
namespace Epetra {

/*! 

\brief Concrete implementation of NOX::Epetra::LinearSolver for Stratimikos.

This solver provides the linear algebra services provided
by Trilinos through the Stratimikos linear solver strategies
package.
 
This class handles construction of both the preconditioners 
and solvers.  All options are determined through parameter 
lists and the basic constructors. 

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
    Epetra_RowMatrix class can be used with an internal preconditioner. 
    The stratimikosParams object is sent directly to the Stratimikos
    package. See Stratimikos for full documentation.
   */

  LinearSystemStratimikos(
    Teuchos::ParameterList& printingParams, 
    Teuchos::ParameterList& stratimikosParams, 
    const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq, 
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
  
      The stratimikosParams object is sent directly to the Stratimikos
      package. See Stratimikos for full documentation.
  */
  LinearSystemStratimikos(
    Teuchos::ParameterList& printingParams, 
    Teuchos::ParameterList& linearSolverParams, 
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

  //! Creates an internally owned Epetra_Operator for the Jacobian. 
  virtual bool createJacobianOperator(
       Teuchos::ParameterList& printParams,
       Teuchos::ParameterList& lsParams,
       const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq, 
       const NOX::Epetra::Vector& cloneVector); 

  //! Builds the linear op with solve factory
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

  virtual void throwError(const string& functionName, 
			  const string& errorMsg) const;

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
  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory;
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<double> > lows;
  Teuchos::RCP<const Thyra::LinearOpBase<double> > linearOp;

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

#ifdef HAVE_NOX_DEBUG
#ifdef HAVE_NOX_EPETRAEXT
  mutable int linearSolveCount;
#endif
#endif

};

} // namespace Epetra
} // namespace NOX
#endif //HAVE_NOX_STRATIMIKOS
#endif
