// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//                  LOCA Continuation Algorithm Package
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#ifndef LOCA_CONTINUATION_EXTENDEDGROUP_H
#define LOCA_CONTINUATION_EXTENDEDGROUP_H

#include "NOX_Parameter_List.H"
#include "LOCA_Extended_AbstractGroup.H"           // base class
#include "LOCA_Continuation_AbstractGroup.H"       // class data element
#include "LOCA_Continuation_ExtendedVector.H"      // class data element
#include "NOX_Common.H"                            // for string 

namespace LOCA { 

  namespace Continuation { 

    /*!
     * \brief Base class for all continuation groups. 
     */
    /*!
     * %Continuation is defined as computing some curve 
     * \f$(x(s),p(s))\in\Re^{n+1}\f$ such that \f$F(x(s),p(s))=0\f$ for some
     * parameterization \f$s\f$.  Given some point \f$(x_0,p_0)\f$ on the
     * curve, another nearby point on the curve is calculated by first 
     * computing a predictor direction \f$v\in\Re^{n+1}\f$ and the 
     * approximate point \f$(x^\ast,p^\ast) = (x_0,p_0) + v\Delta s\f$ 
     * where \f$\Delta s\f$ is the step size.  Then the next point on the 
     * curve is computed by solving the extended set of equations
     *  \f[
     *      \begin{array}{cc}
     *         F(x,p) &= 0 \\
     *         g(x,p,x_0,p_0,x^\ast,p^\ast,v,\Delta s) &= 0
     *      \end{array}
     *  \f]
     *  for \f$(x,p)\f$.  The equation 
     *  \f$g(x,p,x_0,p_0,x^\ast,p^\ast,v,\Delta s)=0\f$ is called the 
     *  continuation equation
     *  and different choices of \f$g\f$ yield different continuation methods.
     *
     * This group, derived from the LOCA::Extended::AbstractGroup, defines an 
     * additional interface that continuation groups should implement
     * so that all continuation groups can be manipulated by the 
     * LOCA::Stepper and LOCA::StepSize objects in a consistent fashion.
     * This group also defines some generic algorithms for commom 
     * computations among continuation groups and manages some of the data
     * needed by all continuation groups..  It stores a pointer to the 
     * underlying group (representing \f$F\f$),
     * a class of the type LOCA::Continuation::AbstractGroup.
     *
     * Ofter duing a continuation run, the continuation path must be
     * reparameterized due to changing relative magnitudes of the solution
     * component \f$x\f$ and parameter component \f$p\f$.  Since any 
     * reparameterization will impact step size control, this group defines
     * an interface for computing a "scale factor" for the continuation
     * equation \f$g = 0\f$ and a step size scale factor which can be
     * multiplied by the step size to have a consitent step size across
     * reparameterizations.  See LOCA::Continuation::ArcLengthGroup for an
     * example on how this can be used.  The default implementation here
     * is for no reparameterization (scale factors are identically 1).
     *
     * These parameters relavant to this group are:
     * <ul>
     * <li> "Initial Scale Factor" The initial scale factor for parameter 
     *      componenet scaling (defaults to 1.0).
     * </ul>
     */
    class ExtendedGroup : public virtual LOCA::Extended::AbstractGroup {

    public:
  
      //! Constructor with integer continuation parameter id
      /*!
       * The argument \em grp is the underlying group which implements \f$F\f$,
       * \em paramID is the index of the continuation parameter, 
       * \em params
       * is the parameter list for the continuation groups described above.
       *
       * The \em grp argument is a non-const reference and may be modified.
       */
      ExtendedGroup(LOCA::Continuation::AbstractGroup& grp, int paramID,
		    NOX::Parameter::List& params);

      //! Constructor with const grp, integer continuation parameter id
      /*!
       * The argument \em grp is the underlying group which implements \f$F\f$,
       * \em paramID is the index of the continuation parameter, 
       * \em params
       * is the parameter list for the continuation groups described above.
       */
      ExtendedGroup(const LOCA::Continuation::AbstractGroup& grp, int paramID,
		    NOX::Parameter::List& params);

      //! Constructor with string continuation parameter id
     /*!
       * The argument \em grp is the underlying group which implements \f$F\f$,
       * \em paramID is the index of the continuation parameter, 
       * \em params
       * is the parameter list for the continuation groups described above.
       *
       * The \em grp argument is a non-const reference and may be modified.
       */ 
      ExtendedGroup(LOCA::Continuation::AbstractGroup& grp, string paramID,
		    NOX::Parameter::List& params);

      //! Constructor with const grp, string continuation parameter id
      /*!
       * The argument \em grp is the underlying group which implements \f$F\f$,
       * \em paramID is the index of the continuation parameter, 
       * and \em params
       * is the parameter list for the continuation groups described above.
       */
      ExtendedGroup(const LOCA::Continuation::AbstractGroup& grp, 
		    string paramID,
		    NOX::Parameter::List& params);

      //! Copy constructor
      ExtendedGroup(const ExtendedGroup& source, 
		    NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      virtual ~ExtendedGroup();

       //! Assignment operator
      virtual NOX::Abstract::Group& 
      operator=(const NOX::Abstract::Group& source);

      //! Assignment operator
      virtual LOCA::Extended::AbstractGroup& 
      operator=(const LOCA::Extended::AbstractGroup& source);

      //! Assignment operator
      virtual ExtendedGroup& 
      operator=(const ExtendedGroup& source);

      //! Set the previous solution vector y
      virtual void setPrevX(const NOX::Abstract::Vector& y) = 0;

      //! Set the previous solution vector y
      virtual void setPrevX(const LOCA::Continuation::ExtendedVector& y) = 0;

      //! Gets the previous solution vector
      virtual const LOCA::Continuation::ExtendedVector& getPrevX() const = 0;

      //! Returns true if the previous solution vector has been set
      virtual bool isPrevXVec() const = 0;

      //! Set step size for continuation constraint equation
      virtual void setStepSize(double deltaS) = 0;

      //! Get step size for continuation constraint equation
      virtual double getStepSize() const = 0;

      //! Scale predictor direction
      virtual void scalePredictor(LOCA::Continuation::ExtendedVector& v) = 0;

      //! Get predictor direction (tangent or secant)
      virtual const LOCA::Continuation::ExtendedVector& 
      getPredictorDirection() const;

      //! Set predictor direction (tangent or secant)
      virtual void setPredictorDirection(
				   const LOCA::Continuation::ExtendedVector&);

      //! Returns true if a predictor direction has been computed
      virtual bool isPredictorDirection() const;

      //! Sets the continuation parameter value
      virtual void setContinuationParameter(double val);

      //! Returns the continuation parameter value
      virtual double getContinuationParameter() const;

      //! Get the continuation parameter id
      virtual int getContinuationParameterID() const;

      //! Sets the scale factor \f$\theta\f$ 
      virtual void setScaleFactor(double theta);
      
      //! Returns the scale factor \f$\theta\f$
      virtual double getScaleFactor() const;

      //! Calculates scale factor as described in the loca v1.0 manual
      virtual void recalculateScaleFactor(double dpds) {}

      //! Returns step size scale factor
      virtual double getStepSizeScaleFactor() const;

      //! Prints the group
      virtual void printSolution() const;

      //! Computes a scaled dot product between two continuation vectors
      virtual double computeScaledDotProduct(
					const NOX::Abstract::Vector& x, 
					const NOX::Abstract::Vector& y) const;

      //! Computes a scaled dot product between two continuation vectors
      virtual double computeScaledDotProduct(
			   const LOCA::Continuation::ExtendedVector& x, 
			   const LOCA::Continuation::ExtendedVector& y) const;

      //! Return underlying group
      virtual const LOCA::Continuation::AbstractGroup& 
      getUnderlyingGroup() const;
      
      //! Return underlying group
      virtual LOCA::Continuation::AbstractGroup& 
      getUnderlyingGroup();

    protected:

      //! Resets all isValid flags to false
      void resetIsValid();

    protected:

      //! Pointer to base group that defines \f$F\f$
      LOCA::Continuation::AbstractGroup* grpPtr;

      //! integer id of continuatin parameter
      int conParamID;

      //! Vector to hold predictor direction
      LOCA::Continuation::ExtendedVector predictorVec;

      //! Flag indicating whether we have our own copy of the underlying group
      bool ownsGroup;

      //! True if predictor has been calculated
      bool isValidPredictor;

      //! Stores scaling factor
      double theta;
    };

  } // namespace Continuation
} // namespace LOCA

#endif
