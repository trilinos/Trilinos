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

#ifndef LOCA_CONTINUATION_ARCLENGTHGROUP_H
#define LOCA_CONTINUATION_ARCLENGTHGROUP_H

#include "LOCA_Continuation_Group.H"

//! Nonlinear solvers package namespace
namespace LOCA { 

//! %NOX  bifurcation vector and group
namespace Continuation { 

class ArcLengthGroup : public LOCA::Continuation::Group {

public:
  
  //! Constructor to set the base group, and derivative objects.
  ArcLengthGroup(const Abstract::Group& g, int paramID,
		 const NOX::Parameter::List& linSolverParams,
		 NOX::Parameter::List& params);
  ArcLengthGroup(const Abstract::Group& g, string paramID,
		 const NOX::Parameter::List& linSolverParams,
		 NOX::Parameter::List& params);

  //! Copy constructor.
  ArcLengthGroup(const ArcLengthGroup& source, NOX::CopyType type = NOX::DeepCopy);

  //! Destructor.
  virtual ~ArcLengthGroup();
  
  /*! 
    \brief Copies the values of all vectors and any other data in
    the source group to this group.  (May invalidate shared data for
    source group.)
  */
  virtual LOCA::Continuation::Group& 
  operator=(const LOCA::Continuation::Group& source);

  virtual NOX::Abstract::Group& operator=(const NOX::Abstract::Group& source);
  virtual ArcLengthGroup& operator=(const ArcLengthGroup& source);
     
  /*! 
    \brief Create a new %Group of the same derived type as this one by
    cloning this one, and return a pointer to the new group.  

    If type is \c DeepCopy, then we need to create an exact replica of
    "this". Otherwise, if type is \c ShapeCopy, we need only replicate
    the shape of "this". Returns \c NULL if clone is not supported.
  */
  virtual NOX::Abstract::Group* clone(NOX::CopyType type = NOX::DeepCopy) const;

  //! Set step size for continuation constraint equation
  virtual void setStepSize(double deltaS);

  //@{ \name "Compute" functions.

  //! Set the solution vector, x, to y.
  virtual void setX(const NOX::Abstract::Vector& y);
  virtual void setX(const LOCA::Continuation::Vector& y);

  //! Set the previous solution vector, x_old, to y.
  virtual void setPrevX(const NOX::Abstract::Vector& y);
  virtual void setPrevX(const LOCA::Continuation::Vector& y);

  //! Gets the previous solution vector
  virtual const LOCA::Continuation::Vector& getPrevX() const;

  virtual bool isPrevXVec() const;

  //! Compute and return solution vector, x, where this.x = grp.x + step * d.
  virtual void computeX(const NOX::Abstract::Group& g, 
			const NOX::Abstract::Vector& d,
			double step);
  virtual void computeX(const ArcLengthGroup& g, 
			const LOCA::Continuation::Vector& d,
			double step);
  
  //! Compute and return F(x). 
  virtual NOX::Abstract::Group::ReturnType computeF();

  //! Compute Jacobian.
  virtual NOX::Abstract::Group::ReturnType computeJacobian();

  //! Compute and return gradient.
  virtual NOX::Abstract::Group::ReturnType computeGradient();

  //! Compute and return Newton direction, using parameters for nonlinear solve
  virtual NOX::Abstract::Group::ReturnType 
  computeNewton(NOX::Parameter::List& params);

  //! Compute tangent direction
  virtual NOX::Abstract::Group::ReturnType
  computeTangent(NOX::Parameter::List& params);

  //@}

  /** @name Jacobian operations.
   
    Operations using the Jacobian matrix. These may not be defined in
    matrix-free scenarios.
  */

  //@{
  
  //! Applies Jacobian to the given input vector and puts the answer in the result.
  virtual NOX::Abstract::Group::ReturnType 
  applyJacobian(const NOX::Abstract::Vector& input, 
		NOX::Abstract::Vector& result) const;

  //! Applies Jacobian-Transpose to the given input vector and puts the answer in the result.
  virtual NOX::Abstract::Group::ReturnType 
  applyJacobianTranspose(const NOX::Abstract::Vector& input, 
			 NOX::Abstract::Vector& result) const;

  /*!
    \brief Applies the inverse of the Jacobian matrix to the given input vector.
  */
  virtual NOX::Abstract::Group::ReturnType 
  applyJacobianInverse(NOX::Parameter::List& params, 
		       const NOX::Abstract::Vector& input, 
		       NOX::Abstract::Vector& result) const;

  //! Apply right preconditiong to the given input vector
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

  virtual void setContinuationParameter(double val);

  virtual void print() const;

  virtual void recalculateScaleFactor(double dpds);

protected:

  void resetIsValid();

protected:
  
  //! Stores the arclength solution vector
  LOCA::Continuation::Vector xVec;

  //! Stores the arclength residual vector
  LOCA::Continuation::Vector fVec;

  //! Stores the arclength Newton vector
  LOCA::Continuation::Vector newtonVec;

  //! Stores the previous arclength solution vector
  LOCA::Continuation::Vector prevXVec;

  //! Stores the derivative of the solution residual w.r.t. arcparam
  NOX::Abstract::Vector* derivResidualParamPtr;

  //! Stores the parameter index
  double arclengthStep;

  //! Is residual vector valid
  bool isValidF;

  //! Is Jacobian matrix valid
  bool isValidJacobian;

  //! Is Newton vector valid
  bool isValidNewton;

  //! Is previous solution vector valid
  bool isValidPrevXVec;

  //! Goal value of dp/ds squared
  double gGoal;

  //! Minimum value for dp/ds for which rescaling is applied
  double gMax;

  //! Maximum value of scale factor
  double thetaMin;

};
} // namespace Bifurcation
} // namespace LOCA

#endif
