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

namespace LOCA { 

  namespace Continuation { 

    /*! 
      \brief Concrete implementation of a LOCA::Continuation::Group for 
      pseudo arc-length continuation.
    */
    /*!
      Pseudo arc-length continuation corresponds to a continuation equation
      \f$g(x,p,x_0,p_0,x^\ast,p^\ast,v,\Delta s)=0\f$ with \f$g\f$ given by
      \f$g(x,p,x_0,p_0,x^\ast,p^\ast,v,\Delta s)=(x-x^\ast)^Tv_x + (p-p^\ast) v_p - \Delta s\f$
      where \f$v_x\f$, \f$v_p\f$ are the solution and  parameter components 
      of the predictor direction \f$v\f$ respectively.
      This corresponds geometrically to constraining the nonlinear solver 
      steps used in calculating \f$F(x,p)=0\f$ to be orthogonal to the 
      predictor direction \f$v\f$.

      This group is a complete implementation of a NOX::Abstract::Group
      for the extended set of equations \f$F = 0\f$, \f$g = 0\f$ and defines
      all relevant methods for this set of equations (computeF, 
      computeJacobian, computeNewton, etc).
    */

    class ArcLengthGroup : public LOCA::Continuation::Group {

    public:
  
      //! Constructor with integer continuation parameter id
      ArcLengthGroup(const Abstract::Group& g, int paramID,
		     const NOX::Parameter::List& linSolverParams,
		     NOX::Parameter::List& params);

      //! Constructor with string continuation parameter id
      ArcLengthGroup(const Abstract::Group& g, string paramID,
		     const NOX::Parameter::List& linSolverParams,
		     NOX::Parameter::List& params);

      //! Copy constructor.
      ArcLengthGroup(const ArcLengthGroup& source, NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      virtual ~ArcLengthGroup();
  
      //! Assignment operator
      virtual LOCA::Continuation::Group& 
      operator=(const LOCA::Continuation::Group& source);

      //! Assignment operator
      virtual NOX::Abstract::Group& operator=(const NOX::Abstract::Group& source);

      //! Assignment operator
      virtual ArcLengthGroup& operator=(const ArcLengthGroup& source);
     
      //! Cloning function
      virtual NOX::Abstract::Group* clone(NOX::CopyType type = NOX::DeepCopy) const;

      //! Set step size for continuation constraint equation
      virtual void setStepSize(double deltaS);

      //@{ \name "Compute" functions.

      //! Set the solution vector to y.
      virtual void setX(const NOX::Abstract::Vector& y);

      //! Set the solution vector to y.
      virtual void setX(const LOCA::Continuation::Vector& y);

      //! Set the previous solution vector to y.
      virtual void setPrevX(const NOX::Abstract::Vector& y);

      //! Set the previous solution vector to y.
      virtual void setPrevX(const LOCA::Continuation::Vector& y);

      //! Gets the previous solution vector
      virtual const LOCA::Continuation::Vector& getPrevX() const;

      //! Returns true if the previous solution vector has been set
      virtual bool isPrevXVec() const;

      //! Compute and return solution vector, x, where this.x = grp.x + step * d.
      virtual void computeX(const NOX::Abstract::Group& g, 
			    const NOX::Abstract::Vector& d,
			    double step);

      //! Compute and return solution vector, x, where this.x = grp.x + step * d.
      virtual void computeX(const ArcLengthGroup& g, 
			    const LOCA::Continuation::Vector& d,
			    double step);
  
      //! Compute extended continuation equations 
      virtual NOX::Abstract::Group::ReturnType computeF();

      //! Compute extended continuation jacobian
      virtual NOX::Abstract::Group::ReturnType computeJacobian();

      //! Gradient is not defined for this system
      virtual NOX::Abstract::Group::ReturnType computeGradient();

      //! Compute Newton direction for extended continuation system
      virtual NOX::Abstract::Group::ReturnType 
      computeNewton(NOX::Parameter::List& params);

      //! Overloaded Compute tangent from LOCA::Continuation::Group
      virtual NOX::Abstract::Group::ReturnType
      computeTangent(NOX::Parameter::List& params);

      //@}

      /** @name Jacobian operations. */

      //@{
  
     //! Applies Jacobian for extended system
      virtual NOX::Abstract::Group::ReturnType 
      applyJacobian(const NOX::Abstract::Vector& input, 
		    NOX::Abstract::Vector& result) const;

      //! Jacobian transpose not defined for this system
      virtual NOX::Abstract::Group::ReturnType 
      applyJacobianTranspose(const NOX::Abstract::Vector& input, 
			     NOX::Abstract::Vector& result) const;

      //! Applies Jacobian inverse for extended system
      virtual NOX::Abstract::Group::ReturnType 
      applyJacobianInverse(NOX::Parameter::List& params, 
			   const NOX::Abstract::Vector& input, 
			   NOX::Abstract::Vector& result) const;

      //! Apply right preconditiong for the extended system
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
  
      //! Return \c true if extended residual is valid.
      virtual bool isF() const;

      //! Return \c true if the extended Jacobian is valid.
      virtual bool isJacobian() const;
      
      //! Always returns false
      virtual bool isGradient() const;
      
      //! Return \c true if the extended Newton direction is valid
      virtual bool isNewton() const;
      //@}

      /** @name "Get" functions.
   
      Note that these function do not check whether or not the vectors
      are valid. Must use the "Is" functions for that purpose.
      */
      //@{ 

      //! Return extended solution vector.  
      virtual const NOX::Abstract::Vector& getX() const;

      //! Return extended residual
      virtual const NOX::Abstract::Vector& getF() const;

      //! Return 2-norm of extended residual. 
      virtual double getNormF() const;

      //! Gradient is never valid
      virtual const NOX::Abstract::Vector& getGradient() const;

      //! Return extended Newton direction.
      virtual const NOX::Abstract::Vector& getNewton() const;

      //! Returns 2-norm of extended Newton solve residual
      virtual double getNormNewtonSolveResidual() const;

      //@}

      //! Sets the continuation parameter value
      virtual void setContinuationParameter(double val);

      //! Prints the group
      virtual void print() const;

      //! Calculates scale factor as described in the loca v1.0 manual
      virtual void recalculateScaleFactor(double dpds);

    protected:

      //! Resets all isValid flags
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
