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

#ifndef LOCA_MULTICONTINUATION_EXTENDEDGROUP_H
#define LOCA_MULTICONTINUATION_EXTENDEDGROUP_H

#include "LOCA_Extended_AbstractGroup.H"                // base class
#include "LOCA_MultiContinuation_AbstractGroup.H"       // class data element
#include "LOCA_MultiContinuation_ExtendedVector.H"      // class data element
#include "LOCA_MultiContinuation_ExtendedMultiVector.H" // class data element
#include "LOCA_BorderedSystem_Manager.H"                // class data element
#include "LOCA_Predictor_Manager.H"                     // class data element

namespace LOCA { 

  namespace MultiContinuation { 

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
      ExtendedGroup(LOCA::MultiContinuation::AbstractGroup& grp, 
		    const vector<int>& paramIDs,
		    NOX::Parameter::List& params);

      //! Copy constructor
      ExtendedGroup(const ExtendedGroup& source, 
		    NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      virtual ~ExtendedGroup();

       //! Assignment operator
      virtual ExtendedGroup& 
      operator=(const ExtendedGroup& source);

      /*! 
       * @name Virtual interface for all continuation methods 
       */
      //@{

      //! Compute continuation constraint equations
      virtual NOX::Abstract::Group::ReturnType
      computeConstraints() = 0;

      //! Compute continuation constraint equations
      virtual NOX::Abstract::Group::ReturnType
      computeConstraintDerivatives() = 0;

      //! Return \c true if constraint residuals are valid
      virtual bool isConstraints() const = 0;

      //! Return \c true if constraint derivatives are valid
      virtual bool isConstraintDerivatives() const = 0;

      //! Return constraint residuals
      virtual const NOX::Abstract::MultiVector::DenseMatrix&
      getConstraints() const = 0;

      //! Return solution component of constraint derivatives
      /*!
       * May return NULL if constraint derivative is zero
       */
      virtual const NOX::Abstract::MultiVector*
      getConstraintDerivativesX() const = 0;

      //! Return parameter component of constraint derivatives
      /*!
       * May return NULL if constraint derivative is zero
       */
      virtual const NOX::Abstract::MultiVector::DenseMatrix*
      getConstraintDerivativesP() const = 0;

      /*! 
       * \brief Return \c true if solution component of constraint 
       * derivatives is zero
       */
      virtual bool isConstraintDerivativesXZero() const = 0;

      /*! 
       * \brief Return \c true if solution component of constraint 
       * derivatives is zero
       */
      virtual bool isConstraintDerivativesPZero() const = 0;

      //@}

      /*! 
       * @name Implementation of NOX::Abstract::Group virtual methods 
       */
      //@{

      //! Assignment operator
      virtual NOX::Abstract::Group& 
      operator=(const NOX::Abstract::Group& source);

      //! Set the solution vector to y.
      virtual void setX(const NOX::Abstract::Vector& y);

      /*! 
       * \brief Compute and return solution vector, x, where 
       * this.x = grp.x + step * d.
       */
      virtual void computeX(const NOX::Abstract::Group& g, 
			    const NOX::Abstract::Vector& d,
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

      //! Applies Jacobian for extended system
      virtual NOX::Abstract::Group::ReturnType 
      applyJacobianMultiVector(const NOX::Abstract::MultiVector& input, 
			       NOX::Abstract::MultiVector& result) const;

      //! Jacobian transpose not defined for this system
      virtual NOX::Abstract::Group::ReturnType 
      applyJacobianTransposeMultiVector(
				     const NOX::Abstract::MultiVector& input, 
				     NOX::Abstract::MultiVector& result) const;

      //! Applies Jacobian inverse for extended system
      virtual NOX::Abstract::Group::ReturnType 
      applyJacobianInverseMultiVector(
				    NOX::Parameter::List& params, 
				    const NOX::Abstract::MultiVector& input, 
				    NOX::Abstract::MultiVector& result) const;

      //! Return \c true if extended residual is valid.
      virtual bool isF() const;

      //! Return \c true if the extended Jacobian is valid.
      virtual bool isJacobian() const;
      
      //! Always returns false
      virtual bool isGradient() const;
      
      //! Return \c true if the extended Newton direction is valid
      virtual bool isNewton() const;

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

      /*! 
       * @name Implementation of LOCA::Extended::AbstractGroup virtual methods 
       */
      //@{

      //! Assignment operator
      virtual LOCA::Extended::AbstractGroup& 
      operator=(const LOCA::Extended::AbstractGroup& source);

      //! Return underlying group
      virtual const LOCA::Continuation::AbstractGroup& 
      getUnderlyingGroup() const;
      
      //! Return underlying group
      virtual LOCA::Continuation::AbstractGroup& 
      getUnderlyingGroup();

      //@}

      /*! 
       * @name Implementation of LOCA::Stepper::AbstractGroup virtual methods 
       */
      //@{

      //! Compute predictor directions
      virtual NOX::Abstract::Group::ReturnType 
      computePredictor();

      //! Get predictor direction \em i
      virtual const LOCA::MultiContinuation::ExtendedVector& 
      getPredictorDirection(int i) const;

      //! Get all predictor directions
      virtual const LOCA::MultiContinuation::ExtendedMultiVector& 
      getPredictorDirections() const;

      //! Set the previous solution vector y
      virtual void setPrevX(const LOCA::MultiContinuation::ExtendedVector& y);

      //! Gets the previous solution vector
      virtual const LOCA::MultiContinuation::ExtendedVector& getPrevX() const ;

      //! Set step size for continuation constraint equation \em i
      virtual void setStepSize(int i, double deltaS);

      //! Get step size for continuation constraint equation \em i
      virtual double getStepSize(int i) const;

      //! Sets the value for continuation parameter \em i
      virtual void setContinuationParameter(int i, double val);

      //! Returns the value for continuation parameter \em i
      virtual double getContinuationParameter(int i) const;

      //! Get the continuation parameter id for parameter \em i
      virtual int getContinuationParameterID(int i) const;

      //! Get the continuation parameter ids
      virtual const vector<int>& getContinuationParameterIDs() const;

      //! Returns step size scale factor for constraint equation \em i
      virtual double getStepSizeScaleFactor(int i) const;

      //! Prints the group
      virtual void printSolution() const;

      //@}

    protected:

      //! applyJacobianInverseMultiVector specialized to Newton solves
      virtual NOX::Abstract::Group::ReturnType 
      applyJacobianInverseNewton(NOX::Parameter::List& params);

      //! Resets all isValid flags to false
      virtual void resetIsValid();

      //! Sets up multivector views
      virtual void setupViews();

    protected:

      //! Pointer to base group that defines \f$F\f$
      LOCA::MultiContinuation::AbstractGroup* grpPtr;

      //! Number of parameters
      int numParams;

      //! Stores the extended solution vector and J^-1 df/dp
      LOCA::MultiContinuation::ExtendedMultiVector xMultiVec;

      //! Stores the extended residual vector and df/dp
      LOCA::MultiContinuation::ExtendedMultiVector fMultiVec;

      //! Stores the extended Newton vector
      LOCA::MultiContinuation::ExtendedMultiVector newtonMultiVec;

      //! Stores the extended gradient vector
      LOCA::MultiContinuation::ExtendedMultiVector gradientMultiVec;

      //! Vector to hold predictor direction
      LOCA::MultiContinuation::ExtendedMultiVector predictorMultiVec;

      //! Stores the previous extended solution vector
      LOCA::MultiContinuation::ExtendedMultiVector prevXMultiVec;

      //! Stores view of first column of xMultiVec
      LOCA::MultiContinuation::ExtendedVector *xVec;

      //! Stores view of first column of fMultiVec
      LOCA::MultiContinuation::ExtendedVector *fVec;
      
      //! Stores view of first column of fMultiVec as a multivec
      LOCA::MultiContinuation::ExtendedMultiVector* ffMultiVec;

      //! Stores view of df/dp columns of fMultiVec
      LOCA::MultiContinuation::ExtendedMultiVector* dfdpMultiVec;

      //! Stores view of first column of newtonMultiVec
      LOCA::MultiContinuation::ExtendedVector *newtonVec;

      //! Stores view of first column of gradientMultiVec
      LOCA::MultiContinuation::ExtendedVector *gradientVec;

      //! Stores view of first column of predictorMultiVec
      LOCA::MultiContinuation::ExtendedVector *predictorVec;

      //! Stores view of first column of prevXMultiVec
      LOCA::MultiContinuation::ExtendedVector *prevXVec;

      //! Stores manager for bordered solve methods
      LOCA::BorderedSystem::Manager borderedSolver;

      //! Stores manager for predictor computations
      LOCA::Predictor::Manager predictorManager;

      //! Stores indices for getting f part of fMultiVec
      vector<int> index_f;

      //! Stores indices for getting df/dp part of fMultiVec
      vector<int> index_dfdp;

      //! integer id of continuation parameters
      vector<int> conParamIDs;

      //! continuation step size
      vector<double> stepSize;

      //! step size scale factors
      vector<double> stepSizeScaleFactor;

      //! Flag indicating whether we have our own copy of the underlying group
      bool ownsGroup;

      //! Is residual vector valid
      bool isValidF;

      //! Is Jacobian matrix valid
      bool isValidJacobian;

      //! Is Newton vector valid
      bool isValidNewton;

      //! Is Gradient vector valid
      bool isValidGradient;

    };

  } // namespace MultiContinuation
} // namespace LOCA

#endif
