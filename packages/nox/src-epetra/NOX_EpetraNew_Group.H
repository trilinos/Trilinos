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

#ifndef NOX_EPETRANEW_GROUP_H
#define NOX_EPETRANEW_GROUP_H

#include "NOX_Abstract_Group.H"	// base class
#include "NOX_Epetra_Vector.H"	// class data element
#include "NOX_Utils.H"          // class data element
#include "NOX_Common.H"         // class data element (string)
#include "NOX_EpetraNew_LinearSystem.H"  // class data element 
#include "NOX_SharedObjectTemplate.H"  // class data element 

// Forward declares
namespace NOX {
  namespace EpetraNew {
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
class Ifpack_IlukGraph;
class Ifpack_CrsRiluk;

namespace NOX {
namespace EpetraNew {

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
  Group(NOX::Parameter::List& printingParams, 
	NOX::EpetraNew::Interface::Required& i, 
	NOX::Epetra::Vector& initialGuess);

  //! Standard Constructor.
  Group(NOX::Parameter::List& printingParams, 
	NOX::EpetraNew::Interface::Required& i, 
	NOX::Epetra::Vector& initialGuess, 
	NOX::EpetraNew::LinearSystem& linSys);

  /*! \brief Copy constructor. If type is DeepCopy, takes ownership of
    valid shared lainear system. */
  Group(const NOX::EpetraNew::Group& source, 
	NOX::CopyType type = NOX::DeepCopy);

  //! Destructor.
  virtual ~Group();

  virtual NOX::Abstract::Group& operator=(const NOX::Abstract::Group& source);

  //! See operator=(const NOX::Abstract::Group&);
  virtual NOX::Abstract::Group& operator=(const NOX::EpetraNew::Group& source);

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
  applyJacobianInverse(NOX::Parameter::List &params, const NOX::Epetra::Vector &input, NOX::Epetra::Vector &result) 
    const;
  virtual NOX::Abstract::Group::ReturnType 
  applyJacobianInverse(NOX::Parameter::List &params, const NOX::Abstract::Vector &input, NOX::Abstract::Vector &result) 
    const;

  virtual NOX::Abstract::Group::ReturnType 
  applyRightPreconditioning(bool useTranspose,
			    NOX::Parameter::List& params, 
			    const NOX::Epetra::Vector& input, 
			    NOX::Epetra::Vector& result) const;

  virtual NOX::Abstract::Group::ReturnType 
  applyRightPreconditioning(bool useTranspose, 
			    NOX::Parameter::List& params, 
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

  //! Return the userInterface.
  virtual NOX::EpetraNew::Interface::Required& getRequiredInterface();

  virtual const NOX::EpetraNew::LinearSystem& getLinearSystem() const;

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
  mutable bool isValidSolverJacOp;
  //@}
  
  //! 2-Norm of RHS
  double normRHS;

  //! 2-Norm of the Newton solve residual: ||Js+f||
  double normNewtonSolveResidual;

  /** @name Shared Operators */
  //@{
  //! Pointer to shared Jacobian matrix 
  NOX::SharedObject<NOX::EpetraNew::LinearSystem, NOX::EpetraNew::Group>* sharedLinearSystemPtr;

  //! Reference to shared Jacobian matrix 
  NOX::SharedObject<NOX::EpetraNew::LinearSystem, NOX::EpetraNew::Group>& sharedLinearSystem;

  //@}

  //! Reference to the user supplied interface functions
  NOX::EpetraNew::Interface::Required& userInterface;

};

} // namespace Epetra
} // namespace NOX


#endif
