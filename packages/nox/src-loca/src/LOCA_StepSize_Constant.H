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

#ifndef LOCA_STEPSIZE_CONSTANT_H
#define LOCA_STEPSIZE_CONSTANT_H

#include "LOCA_StepSize_Generic.H"  // base class

namespace LOCA {

  namespace StepSize {

    //! %Constant step size control strategy
    /*!
     * This class implements a roughly constant step size control strategy.
     * If the previous step was sucessful, the new step size is set equal
     * to the old, otherwise the step size is cut by a supplied factor.  
     * Once a sucessful
     * step is made, the step size is increased by a supplied factor 
     * until the initial step size is reached.
     *
     * This class also incorporates rescaling of the continuation parameter
     * when calculating a step size (common in arc-length continuation).
     * For the first continuation step, the step size is chosen so that
     * step size times the parameter component of the predictor is equal
     * to the initial step size.  From then on, the step size is multiplied
     * by the step size scale factor (see LOCA::Continuation::ExtendedGroup)
     * which incorporates rescaling of the continuation parameter.
     *
     * The parameters used by this class supplied in the constructor or reset
     * method are:
     * <ul>
     * <li> "Max Step Size" - Largest valid step size (Default 1.0e+12)
     * <li> "Min Step Size" - Smallest valid step size (Default 1.0e-12)
     * <li> "Initial Step Size" - Initial step size (Default 1.0)
     * <li> "Failed Step Reduction Factor" - Factor by which step size is 
     *      reduced after a failed step. (Default 0.5)
     * <li> "Successful Step Increase Factor" - Factor by which step size is
     *      increased after a successful step when the step size is smaller
     *      than the initial step size. (Default 1.26 = \f$2^{1/3}\f$)
     * </ul>
     */
    class Constant : public LOCA::StepSize::Generic {

    public:

      //! Constructor. 
      Constant(NOX::Parameter::List& params);

      //! Destructor
      virtual ~Constant();

      //! Reset parameters given above
      virtual NOX::Abstract::Group::ReturnType 
      reset(NOX::Parameter::List& params);
  
      //! Compute the step size as described above
      /*! 
       * Returns NOX::Abstract::Group::Failed if the computed step size is
       * smaller than the minimum step size
       */
      virtual NOX::Abstract::Group::ReturnType 
      compute(LOCA::Continuation::ExtendedGroup& curGroup,
	      const LOCA::Continuation::ExtendedVector& predictor,
	      const NOX::Solver::Generic& solver,
	      const LOCA::Abstract::Iterator::StepStatus& stepStatus,
	      const LOCA::Stepper& stepper,
	      double& stepSize);

      //! Compute the step size as described above
      /*! 
       * Returns NOX::Abstract::Group::Failed if the computed step size is
       * smaller than the minimum step size
       */
      virtual NOX::Abstract::Group::ReturnType 
      compute(LOCA::MultiContinuation::ExtendedGroup& curGroup,
	      const LOCA::MultiContinuation::ExtendedVector& predictor,
	      const NOX::Solver::Generic& solver,
	      const LOCA::Abstract::Iterator::StepStatus& stepStatus,
	      const LOCA::NewStepper& stepper,
	      double& stepSize);

      //! Returns previous step size
      virtual double getPrevStepSize() const;

      //! Returns initial step size
      virtual double getStartStepSize() const;

    protected:

      /*! 
       * \brief Clip the computed step size to the bounds given by the 
       * maximum and minimum step sizes
       */
      /*! 
       * Returns NOX::Abstract::Group::Failed if the computed step size is
       * smaller than the minimum step size
       */
      virtual NOX::Abstract::Group::ReturnType 
      clipStepSize(double& stepSize);

    protected:
      //! Maximum step size
      double maxStepSize;

      //! Minimum step size
      double minStepSize;

      //! Initial step size
      double startStepSize;

      //! Factor by which step size is reduced after a failed step
      double failedFactor;

      //! Factor by which step size is increased after a successful step
      double successFactor;

      //! Previous step size
      double prevStepSize;

      //! Flag indicating if this is the first step
      bool isFirstStep;
    };
  } // namespace StepSize
} // namespace LOCA

#endif
