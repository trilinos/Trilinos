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

#ifndef LOCA_STEPPERARCLENGTH_H
#define LOCA_STEPPERARCLENGTH_H

#include "LOCA_Parameter_Vector.H"   // class data element
#include "NOX_Parameter_List.H"	     // class data element
#include "NOX_Solver_Manager.H"	     // class data element
#include "NOX_StatusTest_Generic.H"  // for enum StatusType
#include "NOX_Common.H"              // for <string> class
#include "LOCA_Abstract_DataOutput.H"
#include "LOCA_Stepper.H"   // class data element

// Forward Declarations
namespace LOCA {
  namespace Abstract {
    class Vector;
    class Group;
  }
  namespace Bifurcation {
    class ArcLengthVector;
    class ArcLengthGroup;
  }
}

namespace LOCA {

class StepperArclength : public LOCA::Stepper {

public:

  //! Constructor 
  StepperArclength(const LOCA::Stepper* stepper,
       LOCA::Bifurcation::ArcLengthGroup& alGroup);

  //! Destructor
  virtual ~StepperArclength();

  //! Solve for the entire continuation run
  virtual NOX::StatusTest::StatusType solve();

  //! Solve for one step of the continuation problem 
  //! (i.e. predictor, stepper, then call nonlinearSolve())
  virtual NOX::StatusTest::StatusType step();

protected:

  
  //! Compute a new step size based on results from the previous step
  virtual double computeStepSize(NOX::StatusTest::StatusType solverStatus);

private:
  //! Copy constructor: prevent from being automatically generated
  StepperArclength(const StepperArclength& s);


  //! Reset the Stepper for new run: prevent being used
  virtual bool reset(LOCA::Abstract::Group& initialGuess,
		     NOX::StatusTest::Generic& t,
		     const NOX::Parameter::List& p);

  //! Not used in derived class: prevent misuse
  virtual bool init();

protected:

  //! Pointer to original stepper used to create this
  //! augmented arclength stepper
  // const Stepper* stepperPtr;

  //! Pointer to current solution group, forced to be ArcLength Group
  LOCA::Bifurcation::ArcLengthGroup* curGroupPtr;

  //! Pointer to current solution group, forced to be ArcLength Group
  LOCA::Bifurcation::ArcLengthGroup* prevGroupPtr;

  //! Variable set to true for first step so that arclength
  //! step size can be calculated from original parameter step size
  //! after tangent calculation.
  bool isFirstArclengthStep;

  //! Ending value of original parameter in stepping routine
  const double originalParamFinalValue;
};
} // end namespace LOCA

#endif

