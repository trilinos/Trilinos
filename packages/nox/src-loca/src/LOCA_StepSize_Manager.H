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

#ifndef LOCA_STEPSIZE_MANAGER_H
#define LOCA_STEPSIZE_MANAGER_H

#include "LOCA_StepSize_Generic.H"  // base class

namespace LOCA {

  namespace StepSize {

    //! %Manager for step size control strategies
    /*! 
      The parameters passed to the constructor or reset should specify
      the "Method", as described below, as well as any additional
      parameters for that particular method.
 
      <ul>
      <li> "Method" - Name of the step size control method. Valid choices are
      <ul> 
      <li> "Constant" (LOCA::StepSize::Constant) [Default]
      <li> "Adaptive" (LOCA::StepSize::Adaptive) 
      </ul>
      </ul>
    */  
    class Manager : public LOCA::StepSize::Generic {

    public:

      //! Constructor. 
      Manager(
	 const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	 const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	 NOX::Parameter::List& params);

      //! Destructor
      virtual ~Manager();

      //! Reset parameters
      virtual NOX::Abstract::Group::ReturnType 
      reset(
	const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	NOX::Parameter::List& params);

      //! Computes step size based on step size control method
      virtual NOX::Abstract::Group::ReturnType 
      compute(LOCA::MultiContinuation::AbstractStrategy& curGroup,
	      const LOCA::MultiContinuation::ExtendedVector& predictor,
	      const NOX::Solver::Generic& solver,
	      const LOCA::Abstract::Iterator::StepStatus& stepStatus,
	      const LOCA::NewStepper& stepper,
	      double& stepSize);

      //! Returns step size control method
      const string& getMethod() const;

      //! Returns previous step size
      virtual double getPrevStepSize() const;

      //! Returns initial step size
      virtual double getStartStepSize() const;

    protected:

      //! Global data
      Teuchos::RefCountPtr<LOCA::GlobalData> globalData;

      //! Stores step size control method name
      string method;

      //! Stores pointer to current step size control object
      LOCA::StepSize::Generic *stepSizePtr;

    };
  } // namespace StepSize
} // namespace LOCA

#endif
