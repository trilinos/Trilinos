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

#ifndef LOCA_CONTINUATION_GROUP_H
#define LOCA_CONTINUATION_GROUP_H

#include "NOX_Parameter_List.H"
#include "NOX_Abstract_Group.H"        // base class
#include "LOCA_Abstract_Group.H"       // class data element
#include "LOCA_Continuation_Vector.H"  // class data element
#include "NOX_Common.H"                // for string 

namespace LOCA { 

  //! %Continuation vectors and groups
  namespace Continuation { 

    /*!
      \brief %Generic wrapper group for computing continuation curves.
    */
    /*!
      %Continuation is defined as computing some curve 
      \f$(x(s),p(s))\in\Re^{n+1}\f$ such that \f$F(x(s),p(s))=0\f$ for some
      parameterization \f$s\f$.  Given some point \f$(x_0,p_0)\f$ on the
      curve, another nearby point on the curve is calculated by first 
      computing a predictor direction \f$v\f$ and the approximate point
      \f$(x^\ast,p^\ast) = (x_0,p_0) + v\Delta s\f$ where \f$\Delta s\f$ is
      the step size.  Then the next point on the curve is computed by
      solving the extended set of equations
      \f[
          \begin{array}{cc}
             F(x,p) &= 0 \\
	     g(x,p,x_0,p_0,x^\ast,p^\ast,v,\Delta s) &= 0
	  \end{array}
      \f]
      for \f$(x,p)\f$.  The equation 
      \f$g(x,p,x_0,p_0,x^\ast,p^\ast,v,\Delta s)=0\f$ is called the 
      continuation equation
      and different choices of \f$g\f$ yield different continuation methods.

      This group, derived from the NOX::Abstract::Group, defines an additional
      interface that continuation groups should implement, as well as some
      generic algorithms for commom computations among continuation groups.
    */

    class Group : public virtual NOX::Abstract::Group {

    public:
  
      //! Constructor with integer continuation parameter id
      Group(const LOCA::Abstract::Group& grp, int paramID,
	    NOX::Parameter::List& linSolverParams,
	    NOX::Parameter::List& params);

      //! Constructor with string continuation parameter id
      Group(const LOCA::Abstract::Group& grp, string paramID,
	    NOX::Parameter::List& linSolverParams,
	    NOX::Parameter::List& params);

      //! Copy constructor
      Group(const Group& source, NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      virtual ~Group();

      //! Assignment operator
      virtual Group& operator=(const Group& source);

      //! Set the previous solution vector y
      virtual void setPrevX(const NOX::Abstract::Vector& y) = 0;

      //! Set the previous solution vector y
      virtual void setPrevX(const LOCA::Continuation::Vector& y) = 0;

      //! Gets the previous solution vector
      virtual const LOCA::Continuation::Vector& getPrevX() const = 0;

      //! Returns true if the previous solution vector has been set
      virtual bool isPrevXVec() const = 0;

      //! Set step size for continuation constraint equation
      virtual void setStepSize(double deltaS) = 0;

      //! Get step size for continuation constraint equation
      virtual double getStepSize() = 0;

      //! Compute tangent direction of continuation curve
      virtual NOX::Abstract::Group::ReturnType
      computeTangent();

      //! Compute secant direction of continuation curve
      virtual NOX::Abstract::Group::ReturnType
      computeSecant();

      //! Get predictor direction (tangent or secant)
      virtual const LOCA::Continuation::Vector& getPredictorDirection() const;

      //! Set predictor direction (tangent or secant)
      virtual void setPredictorDirection(const LOCA::Continuation::Vector&);

      //! Returns true if a predictor direction has been computed
      virtual bool isPredictorDirection() const;

      //! Sets the continuation parameter value
      virtual void setContinuationParameter(double val);

      //! Returns the continuation parameter value
      virtual double getContinuationParameter() const;

      //! Returns the base group (that defines \f$F\f$)
      virtual const LOCA::Abstract::Group& getGroup() const;

      //! Sets the scale factor \f$\theta\f$ 
      virtual void setScaleFactor(double theta);
      
      //! Returns the scale factor \f$\theta\f$
      virtual double getScaleFactor() const;

      //! Calculates scale factor as described in the loca v1.0 manual
      virtual void recalculateScaleFactor(double dpds) {}

      //! Returns step size scale factor
      virtual double getStepSizeScaleFactor()const;

      //! Prints the group
      virtual void print() const = 0;

      //! Computes a scaled dot product between two continuation vectors
      virtual double scaledDotProduct(const LOCA::Continuation::Vector& x, 
				      const LOCA::Continuation::Vector& y) const;

    protected:

      //! Resets all isValid flags to false
      void resetIsValid();

    protected:

      //! Pointer to base group that defines \f$F\f$
      LOCA::Abstract::Group* grpPtr;

      //! integer id of continuatin parameter
      int conParamID;

      //! Vector to hold predictor direction
      LOCA::Continuation::Vector predictorVec;

      //! True if predictor has been calculated
      bool isValidPredictor;

      //! Stores linear solver parameters for computing tangent vector
      NOX::Parameter::List& linearSolverParams;

      //! Stores (extended) scaling vector
      LOCA::Continuation::Vector scaleVec;

      /*! 
	Stores flag indicating whether a constant predictor was used 
	in place of the secant predictor (when no previous x-vec available)
      */
      bool usedConstantPredictor;
  
    };

    //! Computes a scaled dot product
    double scaledDotProduct(const NOX::Abstract::Vector& x, 
			    const NOX::Abstract::Vector& y,
			    const NOX::Abstract::Vector& s);

  } // namespace Continuation
} // namespace LOCA

#endif
