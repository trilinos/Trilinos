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

#ifndef LOCA_STEPPERMULTI_H
#define LOCA_STEPPERMULTI_H

#include "LOCA_Abstract_Iterator.H"    // for status
#include "NOX_Parameter_List.H"	       // class data element
#include "NOX_Solver_Manager.H"	       // class data element
#include "LOCA_Bifurcation_Manager.H"  // class data element
#include "NOX_StatusTest_Generic.H"    // for enum StatusType
#include "NOX_Common.H"                // for <string> class

#include "MFLOCA.H"

// Forward Declarations
namespace LOCA {
  namespace MultiContinuation {
    class ExtendedGroup;
    class AbstractGroup;
  }
}

namespace LOCA {

  /*! 
    \brief Implementation of LOCA::Abstract::Iterator for computing points
    along a continuation curve.
  */
  /*!
    The %Stepper class implements the pure virtual methods of the 
    LOCA::Abstract::Iterator for iteratively computing points along a 
    continuation curve.
  */

  class MultiStepper {

  public:

    //! Constructor 
    MultiStepper(LOCA::MultiContinuation::AbstractGroup& initialGuess,
		 NOX::StatusTest::Generic& t,
		 NOX::Parameter::List& p);
    
    //! Destructor
    virtual ~MultiStepper();

    //! Reset the Stepper to start a new continuation run
    virtual bool reset(LOCA::MultiContinuation::AbstractGroup& initialGuess,
		       NOX::StatusTest::Generic& t,
		       NOX::Parameter::List& p);

    //! Run continuation
    virtual LOCA::Abstract::Iterator::IteratorStatus run();

    //! Return the current solution group.
    virtual LOCA::MultiContinuation::AbstractGroup& getSolutionGroup();

    //! Return the output parameters from the stepper algorithm. 
    virtual const NOX::Parameter::List& getParameterList() const;

  protected:

    //! Copy Constructor: Not implemented, but put here to prevent compiler form generating one
    MultiStepper(const MultiStepper& s);

    //! Print to the screen the initialization information
    virtual void printInitializationInfo();

    //! Parses parameter list to get continuation parameter data
    virtual void getConParamData();

  protected:

    //! Bifurcation group manager
    LOCA::Bifurcation::Manager* bifGroupManagerPtr;

    //! Pointer to bifurcation group
    LOCA::MultiContinuation::AbstractGroup* bifGroupPtr;

    //! Pointer to the Group with the initial guess for the nonlinear solver.  This is used for resetting the nonlinear solver in between successful steps.
    LOCA::MultiContinuation::ExtendedGroup* curGroupPtr;

    //! Pointer to the status test.  
    NOX::StatusTest::Generic* statusTestPtr;

    //! Pointer to parameter list passed in constructor/reset method
    NOX::Parameter::List* paramListPtr;

    //! Solver (a.k.a corrector)
    NOX::Solver::Manager* solverPtr;

    //! Stores continuation parameter vector
    LOCA::ParameterVector paramVec;

    //! Stores vector of continuation parameter IDs
    vector<int> conParamIDVec;

    //! Stores list of continuation parameter data
    list<ParamData> conParamData;
  };
} // end namespace LOCA

#endif

