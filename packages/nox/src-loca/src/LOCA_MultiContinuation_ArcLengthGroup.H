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

#ifndef LOCA_MULTICONTINUATION_ARCLENGTHGROUP_H
#define LOCA_MULTICONTINUATION_ARCLENGTHGROUP_H

#include "LOCA_MultiContinuation_ExtendedGroup.H"       // base class

namespace LOCA { 

  namespace MultiContinuation { 

    /*! 
     * \brief Arclength continuation implementation of 
     * LOCA::MultiContinuation::ExtendedGroup
     */
    class ArcLengthGroup : 
      public virtual LOCA::MultiContinuation::ExtendedGroup {

    public:
  
      //! Constructor with integer continuation parameter ids
      /*!
       * The argument \em grp is the underlying group which implements \f$F\f$,
       * \em paramID is the index of the continuation parameter, 
       * \em params
       * is the parameter list for the continuation groups described above.
       *
       * The \em grp argument is a non-const reference and may be modified.
       */
      ArcLengthGroup(LOCA::MultiContinuation::AbstractGroup& grp, 
		     const vector<int>& paramIDs,
		     NOX::Parameter::List& params);

      //! Constructor for single parameter continuation
      /*!
       * The argument \em grp is the underlying group which implements \f$F\f$,
       * \em paramID is the name of the continuation parameter, 
       * \em params
       * is the parameter list for the continuation groups described above.
       *
       * The \em grp argument is a non-const reference and may be modified.
       */
      ArcLengthGroup(LOCA::MultiContinuation::AbstractGroup& grp, 
		     const string& paramID,
		     NOX::Parameter::List& params);

      //! Copy constructor
      ArcLengthGroup(const ArcLengthGroup& source, 
		     NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      virtual ~ArcLengthGroup();

      //! Assignment operator
      virtual ArcLengthGroup& 
      operator=(const ArcLengthGroup& source);

      /*! 
       * @name Implementation of LOCA::MultiContinuation::ExtendedGroup virtual methods
       */
      //@{

      //! Assignment operator
      virtual LOCA::MultiContinuation::ExtendedGroup& 
      operator=(const LOCA::MultiContinuation::ExtendedGroup& source);

      //! Compute continuation constraint equations
      virtual NOX::Abstract::Group::ReturnType
      computeConstraints();

      //! Compute continuation constraint equations
      virtual NOX::Abstract::Group::ReturnType
      computeConstraintDerivatives();

      //! Return \c true if constraint residuals are valid
      virtual bool isConstraints() const;

      //! Return \c true if constraint derivatives are valid
      virtual bool isConstraintDerivatives() const;

      //! Return constraint residuals
      virtual const NOX::Abstract::MultiVector::DenseMatrix&
      getConstraints() const;

      //! Return solution component of constraint derivatives
      virtual const NOX::Abstract::MultiVector*
      getConstraintDerivativesX() const;

      //! Return parameter component of constraint derivatives
      virtual const NOX::Abstract::MultiVector::DenseMatrix*
      getConstraintDerivativesP() const;

      /*! 
       * \brief Return \c true if solution component of constraint 
       * derivatives is zero
       */
      virtual bool isConstraintDerivativesXZero() const;

      /*! 
       * \brief Return \c true if solution component of constraint 
       * derivatives is zero
       */
      virtual bool isConstraintDerivativesPZero() const;

      //@}

      /*! 
       * @name Implementation of LOCA::Extended::AbstractGroup virtual methods 
       */
      //@{

      //! Assignment operator
      virtual LOCA::Extended::AbstractGroup& 
      operator=(const LOCA::Extended::AbstractGroup& source);

      //@}

      /*! 
       * @name Implementation of NOX::Abstract::Group virtual methods 
       */
      //@{

      //! Assignment operator
      virtual NOX::Abstract::Group& 
      operator=(const NOX::Abstract::Group& source);

      //! Clone function
      virtual NOX::Abstract::Group* 
      clone(NOX::CopyType type = NOX::DeepCopy) const;
      
      //@}

      /*! 
       * @name Implementation of LOCA::Stepper::AbstractGroup virtual methods 
       */
      //@{

      //! Notifies group that continuation step is completed
      virtual void notifyCompletedStep();

      //! Scales predictor
      virtual void scalePredictor();

      //! Computes a scaled dot product between two continuation vectors
      virtual double computeScaledDotProduct(
			       const NOX::Abstract::Vector& x, 
			       const NOX::Abstract::Vector& y) const;

      //! Calculates scale factors
      virtual void recalculateScaleFactor(double dpds, double thetaOld,
					  double& thetaNew);

      //@}

    protected:

      //! Resets all isValid flags to false
      virtual void resetIsValid();

    protected:

      //! Stores arclength constraint residuals
      NOX::Abstract::MultiVector::DenseMatrix constraints;

      //! Flag indicating whether constraints are valid
      bool isValidConstraints;

      //! Flag indicating whether constraint derivatives are valid
      bool isValidConstraintDerivatives;

      //! Stores scaling factor for each arclength equation
      vector<double> theta;

      //! Flag indicating whether to do arc-length scaling
      bool doArcLengthScaling;

      //! Goal value of dp/ds squared
      double gGoal;

      //! Minimum value for dp/ds for which rescaling is applied
      double gMax;

      //! Maximum value of scale factor
      double thetaMin;

      //! Flag indicating whether this is the first rescaling of predictor
      bool isFirstRescale;
    };

  } // namespace MultiContinuation
} // namespace LOCA

#endif
