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

#ifndef LOCA_EPETRA_GROUP_H
#define LOCA_EPETRA_GROUP_H

#include "LOCA_Abstract_Group.H"    // base class
#include "NOX_Epetra_Group.H"       // base class

#include "LOCA_Epetra_Vector.H"	    // class data element
#include "LOCA_Parameter_Vector.H"  // class data element
#include "NOX_Common.H"             // class data element (string)

// Forward declares
namespace NOX {
  namespace Epetra {
    class SharedOperator;
  }
  namespace Parameter {
    class List;
  }
}
namespace LOCA {
  namespace Epetra {
    class Interface;
  }
}

namespace LOCA {
namespace Epetra {

//! Concrete implementation of LOCA::Abstract::Group for %Trilinos/%Epetra.
/*! This group is set up to use the linear algebra services provided through the Trilinos/Epetra package with AztecOO for the linear %solver. 
 */

class Group : public virtual LOCA::Abstract::Group, public virtual NOX::Epetra::Group {

public:
  //! Constructor with Jacobian Operator only.  
  /*! Either there is no preconditioning or the Jacobian will be used
   * for preconditioning.  An Epetra_Operator must be supplied for 
   * the Jacobian even in Matrix-Free mode. linearSolverParams is 
   * the "Linear Solver" sublist of parameter list.
   */
  Group(const NOX::Parameter::List& linearSolverParams, Interface& i, 
 	const ParameterVector& p, Vector& x, Epetra_Operator& J);

  //! Constructor with a separate Jacobian (J) and %Preconditioner (M).  
  //! linearSolverParams is the "Linear %Solver" sublist of parameter list.
  Group(const NOX::Parameter::List& linearSolverParams, Interface& i, 
 	const ParameterVector& p, Vector& x, Epetra_Operator& J, 
	Epetra_Operator& M);

  /*! \brief Copy constructor. If type is DeepCopy, takes ownership of
    valid shared Jacobian and shared preconditioning matrix. */
  Group(const Group& source, NOX::CopyType type = NOX::DeepCopy);

  //! Destructor.
  virtual ~Group();

  virtual Abstract::Group& operator=(const Group& source);
  //! See above.
  virtual Abstract::Group& operator=(const NOX::Abstract::Group& source);
  //! See above.
  virtual Abstract::Group& operator=(const LOCA::Abstract::Group& source);

  /** @name "Compute" functions. */
  //@{

  virtual bool setX(const Vector& y);

  //! See above
  virtual bool setX(const NOX::Abstract::Vector& y);

  virtual bool computeX(const Group& grp, 
			const Vector& d, 
			double step);
  //! See above.
  virtual bool computeX(const NOX::Abstract::Group& grp,
			const NOX::Abstract::Vector& d, 
			double step);

  //! Set the bifurcation/continuation parameter list in the group.
  virtual bool setParams(const ParameterVector& p);

  //! Set the bifurcation/continuation parameter list in the group.
  virtual bool computeParams(const ParameterVector& oldParams,
			     const ParameterVector& direction, 
			     double step);

  virtual bool computeF();

  virtual bool computeJacobian();

  virtual bool computeGradient();

  virtual bool computeNewton(NOX::Parameter::List& params);

  //@}

  /** @name Jacobian operations.
   *
   * Operations using the Jacobian matrix. These may not be defined in
   * matrix-free scenarios. */

  //@{
  
  virtual bool applyJacobian(const Vector& input, Vector& result) const;
  //! See above
  virtual bool applyJacobian(const NOX::Abstract::Vector& input, NOX::Abstract::Vector& result) const;

  virtual bool applyJacobianTranspose(const Vector& input, Vector& result) const;
  //! See above
  virtual bool applyJacobianTranspose(const NOX::Abstract::Vector& input, NOX::Abstract::Vector& result) const;
  
  virtual bool applyJacobianInverse(NOX::Parameter::List &params, const Vector &input, Vector &result) const;

  //! See above
  virtual bool applyJacobianInverse(NOX::Parameter::List &params, const NOX::Abstract::Vector &input, NOX::Abstract::Vector &result) const;

  virtual bool applyJacobianDiagonalInverse(const Vector& input, Vector& result) const;

  //! See above
  virtual bool applyJacobianDiagonalInverse(const NOX::Abstract::Vector& input, NOX::Abstract::Vector& result) const;
  
  virtual bool applyRightPreconditioning(NOX::Parameter::List& params, const Vector& input, Vector& result) const;

  //! See above
  virtual bool applyRightPreconditioning(NOX::Parameter::List& params, const NOX::Abstract::Vector& input, NOX::Abstract::Vector& result) const;

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

  //! Returns true if the value of the Norm of the linear model for a full Newton step ||Js + f|| is valid with respect to the current solution vector.
  virtual bool isNormNewtonSolveResidual() const;

  //! Returns true if an explicitly constructed preconditioner exists (i.e. one that is computed and saved for further use in multiple calls to applyRightPreconditioner).  
  virtual bool isPreconditioner() const;

  //@}

  /** @name "Get" functions 
   *
   * Note that these function do not check whether or not the vectors
   * are valid. Must use the "Is" functions for that purpose. */
  //@{

  virtual const Abstract::Vector& getX() const;

  virtual const ParameterVector& getParams() const;

  virtual const Abstract::Vector& getF() const;
  
  virtual double getNormF() const;

  virtual const Abstract::Vector& getGradient() const;

  virtual const Abstract::Vector& getNewton() const;

  //! Returns the 2-norm of the residual of the linear model used in the Newton solve computation, ||Js+f||.  This does not account for line search adjustments to the step length!
  double getNormNewtonSolveResidual() const;

  //@}

  virtual Abstract::Group* clone(NOX::CopyType type = NOX::DeepCopy) const;

  //! Return the SharedJacobian.
  virtual NOX::Epetra::SharedOperator& getSharedJacobian();

  //! Return the SharedPreconditioner.
  virtual NOX::Epetra::SharedOperator& getSharedPreconditioner();

  //! Return the userInterface.
  virtual NOX::Epetra::Interface& getUserInterface();

protected:

  //! resets the isValid flags to false
  //virtual void resetIsValid();

protected:

  //! Solution vector pointer.
  ParameterVector* pVectorPtr;

  /** @name Additional Shared Operators */
  //@{
  //! Pointer to shared mass matrix 
  //NOX::Epetra::SharedOperator* sharedMassMatrixPtr;

  //@}

  //! Reference to the user supplied interface functions
  Interface& userInterface;

};

} // namespace Epetra
} // namespace LOCA


#endif
