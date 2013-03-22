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

#ifndef NOX_EPETRA_GROUP_H
#define NOX_EPETRA_GROUP_H

#include "NOX_Abstract_Group.H"	// base class
#include "NOX_Epetra_Vector.H"	// class data element
#include "NOX_Utils.H"          // class data element
#include "NOX_Common.H"         // class data element (std::string)
#include "NOX_Epetra_LinearSystem.H"  // class data element 
#include "NOX_SharedObjectTemplate.H"  // class data element 
#include "Teuchos_RCP.hpp"  // class data element

// Forward declares
namespace NOX {
  namespace Epetra {
    class Scaling;
    namespace Interface {
      class Required;
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
class AztecOOConditionNumber;
class Ifpack_IlukGraph;
class Ifpack_CrsRiluk;

namespace NOX {
namespace Epetra {

//! Concrete implementation of NOX::Abstract::Group for Trilinos/Epetra.
/*! This group is set up to use the linear algebra services provided
through the Trilinos/Epetra package with AztecOO for the linear
%solver.
 */
class Group : public virtual NOX::Abstract::Group {

public:
  //! Constructor with NO linear system (VERY LIMITED).
  /*! WARNING: If this constructor is used, then methods that require
    a Jacobian or preconditioning will not be available.  You will be
    limited to simple algorithms like nonlinear-CG with no
    preconditioning.
  */
  Group(Teuchos::ParameterList& printingParams, 
	const Teuchos::RCP<NOX::Epetra::Interface::Required>& i, 
	const NOX::Epetra::Vector& initialGuess);

  //! Standard Constructor.
  Group(Teuchos::ParameterList& printingParams, 
	const Teuchos::RCP<NOX::Epetra::Interface::Required>& i, 
	const NOX::Epetra::Vector& initialGuess, 
	const Teuchos::RCP<NOX::Epetra::LinearSystem>& linSys);

  /*! \brief Copy constructor. If type is DeepCopy, takes ownership of
    valid shared linear system. */
  Group(const NOX::Epetra::Group& source, 
	NOX::CopyType type = NOX::DeepCopy);

  //! Destructor.
  virtual ~Group();

  virtual NOX::Abstract::Group& operator=(const NOX::Abstract::Group& source);

  //! See operator=(const NOX::Abstract::Group&);
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

  virtual NOX::Abstract::Group::ReturnType computeNewton(Teuchos::ParameterList& params);

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

     /*!
     \brief Applies the inverse of the Jacobian matrix to the given
     input vector and puts the answer in result.
 
     Computes
     \f[ v = J^{-1} u, \f]
     where \f$J\f$ is the Jacobian, \f$u\f$ is the input vector, and \f$v\f$ is
the result vector.
 
     The "Tolerance" parameter specifies that the
     solution should be such that
     \f[
     \frac{\| J v - u \|_2}{\max \{ 1, \|u\|_2\} } < \mbox{Tolerance}
     \f]
 
     \return
     <ul>
     <li> NOX::Abstract::Group::NotDefined - Returned by default implementation
in NOX::Abstract::Group
     <li> NOX::Abstract::Group::BadDependency - If \f$J\f$ has not been computed     <li> NOX::Abstract::Group::NotConverged - If the linear solve fails to satisfy the "Tolerance"
          specified in \c params
     <li> NOX::Abstract::Group::Failed - If the computation fails
     <li> NOX::Abstract::Group::Ok - Otherwise
     </ul>
 
     The parameter "Tolerance" may be added/modified in the list of
     parameters - this is the ideal solution tolerance for an iterative
     linear solve.

     The parameter "Reuse Preconditioner" is a boolean that tells the group to turn off control of preconditioner recalculation.  This is a dangerous flag but can really speed the computations if the user knows what they are doing.  Toggling this flag is left to the user (ideally it should be done through a status test).  Defaults to false.
   */
  virtual NOX::Abstract::Group::ReturnType 
  applyJacobianInverse(Teuchos::ParameterList &params, const NOX::Epetra::Vector &input, NOX::Epetra::Vector &result) 
    const;
  virtual NOX::Abstract::Group::ReturnType 
  applyJacobianInverse(Teuchos::ParameterList &params, const NOX::Abstract::Vector &input, NOX::Abstract::Vector &result) 
    const;

  virtual NOX::Abstract::Group::ReturnType 
  applyRightPreconditioning(bool useTranspose,
			    Teuchos::ParameterList& params, 
			    const NOX::Epetra::Vector& input, 
			    NOX::Epetra::Vector& result) const;

  virtual NOX::Abstract::Group::ReturnType 
  applyRightPreconditioning(bool useTranspose, 
			    Teuchos::ParameterList& params, 
			    const NOX::Abstract::Vector& input, 
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

  /*! \brief Returns true if the condition number has been computed.  
  */
  virtual bool isConditionNumber() const;

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

  inline virtual Teuchos::RCP< const NOX::Abstract::Vector > getXPtr() const {return xVectorPtr;};

  inline virtual Teuchos::RCP< const NOX::Abstract::Vector > getFPtr() const {return RHSVectorPtr;};

  inline virtual Teuchos::RCP< const NOX::Abstract::Vector > getGradientPtr() const {return gradVectorPtr;};

  inline virtual Teuchos::RCP< const NOX::Abstract::Vector > getNewtonPtr() const {return NewtonVectorPtr;};

  /*! 
    \brief
    Returns the 2-norm of the residual of the linear model used in the
    Newton solve computation, ||Js+f||.  This does not account for
    line search adjustments to the step length!
  */
  virtual NOX::Abstract::Group::ReturnType 
  getNormLastLinearSolveResidual(double & residual) const;

  //@}

  virtual Teuchos::RCP<NOX::Abstract::Group> 
  clone(CopyType type = DeepCopy) const;

  //! Return the userInterface.
  virtual Teuchos::RCP<NOX::Epetra::Interface::Required> 
  getRequiredInterface();

  //! Return the Linear System.
  virtual Teuchos::RCP<const NOX::Epetra::LinearSystem> 
  getLinearSystem() const;

  //! Return the Linear System.
  virtual Teuchos::RCP<NOX::Epetra::LinearSystem> getLinearSystem();

  /* \brief Computes the condition number of the Jacobian matrix.

     Uses GMRES to estimate the condtion number.  The tolerance and
     maxIters are used to control the GMRES iterations.  Typically the
     solves do not have to be tight to get the estimate.
  */
  virtual NOX::Abstract::Group::ReturnType 
  computeJacobianConditionNumber(int maxIters, double tolerance,
				 int krylovSubspaceSize=100, 
				 bool printOutput=false);

  //! Returns the condition number of the Jacobian matrix.
  virtual double getJacobianConditionNumber() const;

  /*! 
    \brief 
    Sets option to disable linear resid computation. If disabled, 
    this saves on a MatVec per Newton but disallows inexact Newton methods
  */
  virtual void disableLinearResidualComputation(const bool disableChoice);
 
protected:

  //! resets the isValid flags to false
  virtual void resetIsValid();

  /*! 
    \brief 
    Computes the 2-norm of the residual of the linear model used in
    the Newton solve computation, ||Js+f||.
  */
  virtual bool computeNormNewtonSolveResidual();

protected:

  //! Printing Utilities object
  const NOX::Utils utils;

  /** @name Vectors */
  //@{
  //! Solution vector pointer.
  Teuchos::RCP<NOX::Epetra::Vector> xVectorPtr;
  //! Solution vector.
  NOX::Epetra::Vector& xVector;
  //! Right-hand-side vector (function evaluation).
  Teuchos::RCP<NOX::Epetra::Vector> RHSVectorPtr;
  //! Right-hand-side vector pointer (function evaluation).
  NOX::Epetra::Vector& RHSVector;
  //! Gradient vector pointer(steepest descent vector).
  Teuchos::RCP<NOX::Epetra::Vector> gradVectorPtr;
  //! Gradient vector (steepest descent vector).
  NOX::Epetra::Vector& gradVector;
  //! Newton direction vector pointer.
  Teuchos::RCP<NOX::Epetra::Vector> NewtonVectorPtr;
  //! Newton direction vector.
  NOX::Epetra::Vector& NewtonVector;
  //! An extra temporary vector, only allocated if needed.
  mutable Teuchos::RCP<Epetra_Vector> tmpVectorPtr;
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
  mutable bool isValidSolverJacOp;
  bool isValidConditionNumber;
  //@}
  
  //! 2-Norm of the Newton solve residual: ||Js+f||
  double normNewtonSolveResidual;

  //! condition number of Jacobian
  double conditionNumber;

  //! Pointer to the condition number object.
  Teuchos::RCP<AztecOOConditionNumber> azConditionNumberPtr;

  /** @name Shared Operators */
  //@{
  //! Pointer to shared Jacobian matrix 
  Teuchos::RCP< 
    NOX::SharedObject<NOX::Epetra::LinearSystem, NOX::Epetra::Group> 
    > sharedLinearSystemPtr;

  //! Reference to shared Jacobian matrix 
  NOX::SharedObject<NOX::Epetra::LinearSystem, NOX::Epetra::Group>& 
  sharedLinearSystem;

  //@}

  // Internal flag to disable linear resid computation. False unless set.
  bool linearResidCompDisabled;

  //! Reference to the user supplied interface functions
  Teuchos::RCP<NOX::Epetra::Interface::Required> userInterfacePtr;

};

} // namespace Epetra
} // namespace NOX


#endif
