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

#ifndef LOCA_BIFURCATION_TPBORDGROUP_H
#define LOCA_BIFURCATION_TPBORDGROUP_H

#include "LOCA_Abstract_Group.H"
#include "LOCA_Bifurcation_TPBordVector.H"
#include "LOCA_DerivUtils.H"

namespace LOCA {
  class ParameterVector;
}

//! Nonlinear solvers package namespace
namespace LOCA { 

//! %NOX  bifurcation vector and group
namespace Bifurcation { 

/*! 
  \brief %LOCA's generic turning point group.  Implements the bordering
  algorithm to solve turning point equations.  Derived from the 
  LOCA::Abstract::Group class.

  \author Eric Phipps (SNL 9233)
*/

class TPBordGroup : public virtual LOCA::Abstract::Group {

public:
  
  //! Constructor to set the base group, and derivative objects.
  TPBordGroup(const Abstract::Group& g, 
	      const NOX::Abstract::Vector& lenVec, int paramId = 0,
	      const DerivUtils& d = DerivUtils());

  //! Copy constructor.
  TPBordGroup(const TPBordGroup& source, NOX::CopyType type = NOX::DeepCopy);

  //! Destructor.
  virtual ~TPBordGroup();
  
  /*! 
    \brief Copies the values of all vectors and any other data in
    the source group to this group.  (May invalidate shared data for
    source group.)
  */
  virtual LOCA::Abstract::Group& operator=(const LOCA::Abstract::Group& source);
  virtual NOX::Abstract::Group& operator=(const NOX::Abstract::Group& source);
  virtual TPBordGroup& operator=(const TPBordGroup& source);
     
  /*! 
    \brief Create a new %Group of the same derived type as this one by
    cloning this one, and return a pointer to the new group.  

    If type is \c DeepCopy, then we need to create an exact replica of
    "this". Otherwise, if type is \c ShapeCopy, we need only replicate
    the shape of "this". Returns \c NULL if clone is not supported.
  */
  virtual NOX::Abstract::Group* clone(NOX::CopyType type = NOX::DeepCopy) const;

  //! Set the parameter list in the group to p (pVector = p).  
  virtual void setParams(const ParameterVector& p);

  //! Return a const reference to the ParameterVector owned by the group. 
  virtual const ParameterVector& getParams() const;

  //! Set parameter indexed by paramID
  virtual void setParam(int paramID, double val);

  //! Return copy of parameter indexed by paramID
  virtual double getParam(int paramID) const;

  //! Set parameter indexed by paramID
  virtual void setParam(string paramID, double val);

  //! Return copy of parameter indexed by paramID
  virtual double getParam(string paramID) const;

  //@{ \name "Compute" functions.

  //! Set the solution vector, x, to y.
  virtual void setX(const NOX::Abstract::Vector& y);

  virtual void setX(const TPBordVector& y);

  //! Compute and return solution vector, x, where this.x = grp.x + step * d.
  /*! Let \f$x\f$ denote this group's solution and \f$\hat x\f$ denote
    the solution associated with group \c grp. Then we compute 
    \f[ x = \hat x + \mbox{step} \; d. \f]
  */
  virtual void computeX(const NOX::Abstract::Group& g, 
			const NOX::Abstract::Vector& d,
			double step);

  virtual void computeX(const TPBordGroup& g, const TPBordVector& d,
			double step);

  //! Compute and return F(x). 
  /*! In other words, return \f$F \equiv F(x)\f$. It's generally
    useful to also compute and stores norm of F(x) at this point for
    later access by the getNormF function. */
  virtual NOX::Abstract::Group::ReturnType computeF();

  //! Compute Jacobian.
  /*! In other words, compute \f$  J_{ij} = \frac{\partial F_i}{\partial x_j} (x)\f$.
    \note If not defined in derived object, returns \c false.
  */
  virtual NOX::Abstract::Group::ReturnType computeJacobian();

  //! Compute and return gradient.
  //! Returns \c false if either F or the Jacobian have not been computed.
  /*! In other words, compute \f$g = J^T F\f$. 
    \note If not defined in derived object, returns \c false.
  */
  virtual NOX::Abstract::Group::ReturnType computeGradient();

  //! Compute and return Newton direction, using parameters for nonlinear solve.
  //! Returns \c false if F and/or the Jacobian have not been computed.
  /*! Solve \f$Js = -F\f$ for \f$s\f$.  

  The parameters are from the "Linear %Solver" sublist of the
  "%Direction" sublist that is passed to solver during construction.

  The parameter "Tolerance" may be added/modified in the list of
  parameters - this is the ideal solution tolerance for an iterative
  linear solve.

  \note If not defined in derived object, returns \c false.
  */
  virtual NOX::Abstract::Group::ReturnType 
  computeNewton(NOX::Parameter::List& params);

  //@}

  /** @name Jacobian operations.
   
    Operations using the Jacobian matrix. These may not be defined in
    matrix-free scenarios.
  */

  //@{
  
  //! Applies Jacobian to the given input vector and puts the answer in the result.
  /*! 
    Computes 
    \f[ \mbox{result} = J \; \mbox{input}, \f]
    and returns \c true.

    Returns \c false if any errors are encountered.

    \note If not defined in derived object, returns \c false.
  */
  virtual NOX::Abstract::Group::ReturnType 
  applyJacobian(const NOX::Abstract::Vector& input, 
		NOX::Abstract::Vector& result) const;

  //! Applies Jacobian-Transpose to the given input vector and puts the answer in the result.

  /*!  
    Computes 
    \f[ \mbox{result} = J^T \; \mbox{input}, \f]
    and returns \c true.

    Returns \c false if any errors are encountered.

    \note If not defined in derived object, returns \c false.
  */
  virtual NOX::Abstract::Group::ReturnType 
  applyJacobianTranspose(const NOX::Abstract::Vector& input, 
			 NOX::Abstract::Vector& result) const;

  /*!
    \brief Applies the inverse of the Jacobian matrix to the given input vector.
  */
  /*! 
    Computes 
    \f[ \mbox{result} = J^{-1} \; \mbox{input}, \f]
    and returns \c true.

    Returns \c false if any errors are encountered.

  The parameters are from the "Linear Solver" sublist of the
  "Direction" sublist that is passed to solver during construction.

  The parameter "Tolerance" may be added/modified in the list of
  parameters - this is the ideal solution tolerance for an iterative
  linear solve.

    \note If not defined in derived object, returns \c false.
  */
  virtual NOX::Abstract::Group::ReturnType 
  applyJacobianInverse(NOX::Parameter::List& params, 
		       const NOX::Abstract::Vector& input, 
		       NOX::Abstract::Vector& result) const;

  //! Applies the inverse of the diagonal of the Jacobian to the given input vector.
  /*! 
    Only necessary if using NOX::Solver::NonlinearCG.
    \note If not defined in derived object, returns \c false.
  */
  virtual NOX::Abstract::Group::ReturnType 
  applyJacobianDiagonalInverse(const NOX::Abstract::Vector& input,
			       NOX::Abstract::Vector& result) const;

  //! Apply right preconditiong to the given input vector
  /*!
    Let \f$M\f$ be a right preconditioner for the Jacobian \f$J\f$; in
    other words, \f$M\f$ is a matrix such that
    \f[ JM \approx I. \f]
    Computes 
    \f[ \mbox{result} = M^{-1} \; \mbox{input}, \f]
    and returns \c true.

    Returns false if any errors are encountered.

    The parameters are from the "Linear %Solver" sublist of the
    "%Direction" sublist that is passed to solver during construction.
  */
  virtual NOX::Abstract::Group::ReturnType 
  applyRightPreconditioning(NOX::Parameter::List& params,
			    const NOX::Abstract::Vector& input, 
			    NOX::Abstract::Vector& result) const;
  //@}

  /** @name "Is" functions.
   
    Checks to see if various objects have been computed. Returns \c true
    if the corresponding "compute" function has been called since the
    last update to the solution vector (via instantiation or
    computeX).
  */

  //@{
  
  //! Return \c true if F is valid.
  virtual bool isF() const;
  //! Return \c true if the Jacobian is valid.
  /*! \note If not defined in derived object, returns \c false. */
  virtual bool isJacobian() const;
  //! Return \c true if the gradient is valid.
  /*! \note If not defined in derived object, returns \c false. */
  virtual bool isGradient() const;
  //! Return \c true if the Newton direction is valid.
  /*! \note If not defined in derived object, returns \c false. */
  virtual bool isNewton() const;
  //@}

  /** @name "Get" functions.
   
    Note that these function do not check whether or not the vectors
    are valid. Must use the "Is" functions for that purpose.
  */
  //@{ 

  //! Return solution vector.  
  virtual const NOX::Abstract::Vector& getX() const;

  //! Return F(x)
  virtual const NOX::Abstract::Vector& getF() const;

  //! Return 2-norm of F(x). 
  /*! In other words, \f[ \sqrt{\sum_{i=1}^n F_i^2} \f] */
  virtual double getNormF() const;

  //! Return gradient.
  virtual const NOX::Abstract::Vector& getGradient() const;

  //! Return Newton direction.
  virtual const NOX::Abstract::Vector& getNewton() const;

  /*!
    \brief Return the norm of the Newton solve residual.
  
    This is the error in the linear model for the particular Newton step.
    In other words, 

    \f[ ||Js - (-F)||_2 \f] 

    where \f$J\f$ is the Jacobian calculated by computeJacobian(),
    \f$s\f$ is the Newton step calculated by computeNewton(), and
    \f$F\f$ is the function value calulated by computeF().

    \note Returns -1 if there is an error in calculating the linear residual.
  */
  virtual double getNormNewtonSolveResidual() const;

  //@}

  virtual void print() const;

protected:

  //! Get bifurcation parameter
  double getBifParam() const;

  //! Set bifurcation parameter
  void setBifParam(double param);

protected:

  //! Stores the underlying group, which stores the solution vector,
  //! jacobian, and solution residual
  LOCA::Abstract::Group* grpPtr;   
  
  //! Stores the turning point solution vector
  TPBordVector tpXVec;

  //! Stores the turning point residual vector
  TPBordVector tpFVec;

  //! Stores the turning point Newton vector
  TPBordVector tpNewtonVec;

  //! Stores the constant length vector
  NOX::Abstract::Vector* lengthVecPtr;

  //! Stores the bifurcation parameter index
  int bifParamId;

  //! Stores the derivative of the solution residual w.r.t. bifparam
  NOX::Abstract::Vector* derivResidualParamPtr;

  //! Stores the derivative of the null vector residual w.r.t. bifparam
  NOX::Abstract::Vector* derivNullResidualParamPtr;

  //! Stores the directional derivative object used to compute derivatives
  //! of residual equations
  DerivUtils* derivPtr;

  //! Is residual vector valid
  bool isValidF;

  //! Is Jacobian matrix valid
  bool isValidJacobian;

  //! Is Newton vector valid
  bool isValidNewton;

  //! Is Tangent vector valid
  bool isValidTangent;
};
} // namespace Bifurcation
} // namespace LOCA

#endif
