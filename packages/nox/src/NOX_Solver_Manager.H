// $Id$ 
// $Source$ 

// NOX: An Object-Oriented Nonlinear Solver Package
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#ifndef _NOX_SOLVER_MANAGER_H
#define _NOX_SOLVER_MANAGER_H

#include "NOX_Solver_Generic.H"	// base class
#include "NOX_Common.H"		// class data element (string)
namespace NOX {
namespace Solver{

/*! 

  \brief Concrete class to control the instatiation of the derived 
  NLS_Method class.  

  This is the class that takes control from the user code and hands it
  over to the nonlinear solver component.

  The parameters may specify the method:

   <ul>
   <li> "Nonlinear %Solver" - Name of the solver method. Valid choices are
   <ul> 
   <li> "%Newton" (NOX::Solver::LineSearch with step length of 1.0) [Default]
   <li> "Line Search" (NOX::Solver::LineSearch)
   <li> "Nonlinear CG" (NOX::Solver::NonlinearCG)
   <li> "Trust Region" (NOX::Solver::TrustRegion)
   </ul>
   </ul>
*/

class Manager : public Generic {

public:
  /*!
    \brief Generic interface. On input, the grp should contain the
    initial guess. On output it will contain the solution. 
  */
  Manager(Abstract::Group& grp, StatusTest::Generic &t, const Parameter::List& p);

  //! Destructor.
  ~Manager();

  bool reset(Abstract::Group& grp, StatusTest::Generic &t, const Parameter::List& p);

  StatusTest::StatusType getStatus();

  StatusTest::StatusType iterate();

  StatusTest::StatusType solve();

  const Abstract::Group& getSolutionGroup() const;
  
  const Abstract::Group& getPreviousSolutionGroup() const;

  int getNumIterations() const;

  const Parameter::List& getOutputParameters() const;

private:

  //! Name of method being used
  string method;

  //! Pointer to the solver object
  Generic* ptr;

};
}
}
#endif
