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

#ifndef LOCA_STEPSIZE_GENERIC_H
#define LOCA_STEPSIZE_GENERIC_H

#include "NOX_Abstract_Group.H"  // for ReturnType, Parameter::List
#include "LOCA_Abstract_Iterator.H" // for StepStatus

// Forward declarations
namespace NOX {
  namespace Solver {
    class Generic;
  }
}
namespace LOCA {
  namespace Continuation {
    class ExtendedGroup;
    class ExtendedVector;
  }
  class Stepper;
  namespace MultiContinuation {
    class ExtendedGroup;
    class ExtendedVector;
  }
  class NewStepper;
}

namespace LOCA {

  //! Step size control strategies
  namespace StepSize {

    //! %Generic step size control interface.
    /*! 
     * %Generic interface for computing a step size.  Any step size control
     * strategy class should be derived from this class and implement the
     * pure virtual methods of this class.
     */
    class Generic {

    public:

      //! Constructor. 
      /*! 
	Constructors of derived objects should look like reset.
      */
      Generic() {};

      //! Destructor
      virtual ~Generic() {};

      //! Reset parameters
      virtual NOX::Abstract::Group::ReturnType 
      reset(NOX::Parameter::List& params) = 0;
  
      /*! 
	\brief Compute the step size.
      */
      /*!
       * Returns NOX::Abstract::Group::Failed if the computed step size
       * fails certain criteria (e.g., if it is too small)
       */
      virtual NOX::Abstract::Group::ReturnType 
      compute(LOCA::Continuation::ExtendedGroup& curGroup,
	      const LOCA::Continuation::ExtendedVector& predictor,
	      const NOX::Solver::Generic& solver,
	      const LOCA::Abstract::Iterator::StepStatus& stepStatus,
	      const LOCA::Stepper& stepper,
	      double& stepSize) = 0;

      /*! 
	\brief Compute the step size.
      */
      /*!
       * Returns NOX::Abstract::Group::Failed if the computed step size
       * fails certain criteria (e.g., if it is too small)
       */
      virtual NOX::Abstract::Group::ReturnType 
      compute(LOCA::MultiContinuation::ExtendedGroup& curGroup,
	      const LOCA::MultiContinuation::ExtendedVector& predictor,
	      const NOX::Solver::Generic& solver,
	      const LOCA::Abstract::Iterator::StepStatus& stepStatus,
	      const LOCA::NewStepper& stepper,
	      double& stepSize) = 0;

      //! Return the previous step size
      virtual double getPrevStepSize() const = 0;

      //! Return the initial step size
      virtual double getStartStepSize() const = 0;

    };
  } // namespace StepSize
} // namespace LOCA

#endif
