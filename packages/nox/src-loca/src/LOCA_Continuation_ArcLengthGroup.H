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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#ifndef LOCA_CONTINUATION_ARCLENGTHGROUP_H
#define LOCA_CONTINUATION_ARCLENGTHGROUP_H

#include "LOCA_Continuation_ExtendedGroup.H"

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

    class ArcLengthGroup : public LOCA::Continuation::ExtendedGroup {

    public:
  
      //! Constructor with integer continuation parameter id
      ArcLengthGroup(LOCA::Continuation::AbstractGroup& g, int paramID,
		     NOX::Parameter::List& params);

      //! Constructor with string continuation parameter id
      ArcLengthGroup(LOCA::Continuation::AbstractGroup& g, string paramID,
		     NOX::Parameter::List& params);

      //! Copy constructor.
      ArcLengthGroup(const ArcLengthGroup& source, 
		     NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      virtual ~ArcLengthGroup();
  
      //! Assignment operator
      virtual LOCA::Continuation::ExtendedGroup& 
      operator=(const LOCA::Continuation::ExtendedGroup& source);

      //! Assignment operator
      virtual NOX::Abstract::Group& 
      operator=(const NOX::Abstract::Group& source);

      //! Assignment operator
      virtual ArcLengthGroup& 
      operator=(const ArcLengthGroup& source);
     
      //! Cloning function
      virtual NOX::Abstract::Group* 
      clone(NOX::CopyType type = NOX::DeepCopy) const;

      //! Set step size for continuation constraint equation
      virtual void setStepSize(double deltaS);

      //! Get step size for continuation constraint equation
      virtual double getStepSize() const;

      //@{ \name "Compute" functions.

      //! Set the solution vector to y.
      virtual void setX(const NOX::Abstract::Vector& y);

      //! Set the solution vector to y.
      virtual void setX(const LOCA::Continuation::ExtendedVector& y);

      //! Set the previous solution vector to y.
      virtual void setPrevX(const NOX::Abstract::Vector& y);

      //! Set the previous solution vector to y.
      virtual void setPrevX(const LOCA::Continuation::ExtendedVector& y);

      //! Gets the previous solution vector
      virtual const LOCA::Continuation::ExtendedVector& getPrevX() const;

      //! Returns true if the previous solution vector has been set
      virtual bool isPrevXVec() const;

      //! Scale predictor direction
      virtual void scalePredictor(LOCA::Continuation::ExtendedVector& v);

      //! Compute and return solution vector, x, where this.x = grp.x + step * d.
      virtual void computeX(const NOX::Abstract::Group& g, 
			    const NOX::Abstract::Vector& d,
			    double step);

      //! Compute and return solution vector, x, where this.x = grp.x + step * d.
      virtual void computeX(const ArcLengthGroup& g, 
			    const LOCA::Continuation::ExtendedVector& d,
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

      //! Calculates scale factor as described in the loca v1.0 manual
      virtual void recalculateScaleFactor(double dpds);

      //! Returns step size scale factor
      virtual double getStepSizeScaleFactor() const;

    protected:

      //! Resets all isValid flags
      void resetIsValid();

    protected:
  
      //! Stores the arclength solution vector
      LOCA::Continuation::ExtendedVector xVec;

      //! Stores the arclength residual vector
      LOCA::Continuation::ExtendedVector fVec;

      //! Stores the arclength Newton vector
      LOCA::Continuation::ExtendedVector newtonVec;

      //! Stores the arclength gradient vector
      LOCA::Continuation::ExtendedVector gradientVec;

      //! Stores the previous arclength solution vector
      LOCA::Continuation::ExtendedVector prevXVec;

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

      //! Is Gradient vector valid
      bool isValidGradient;

      //! Is previous solution vector valid
      bool isValidPrevXVec;

      //! Flag indicating whether to do arc-length scaling
      bool doArcLengthScaling;

      //! Goal value of dp/ds squared
      double gGoal;

      //! Minimum value for dp/ds for which rescaling is applied
      double gMax;

      //! Maximum value of scale factor
      double thetaMin;

      //! Scaling factor for step size
      double stepSizeScaleFactor;

      //! Flag indicating whether this is the first rescaling of predictor
      bool isFirstRescale;

    };
  } // namespace Bifurcation
} // namespace LOCA

#endif
