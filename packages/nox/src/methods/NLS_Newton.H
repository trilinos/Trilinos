// $Id$ 
// $Source$ 

// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#ifndef _NLS_NEWTON_H
#define _NLS_NEWTON_H

#include "NLS_Method.H"		// base class
#include "NLS_LineSearchManager.H" // line search
#include "NLS_StatusTest.H"	// status test for convergence

//! Implementation of plain Newton solver
class NLS_Newton : public NLS_Method {

public:

  //! Newton-method constructor. 
  NLS_Newton(NLS_Group& initialguess, 
	     NLS_Group& workspace, 
	     NLS_ParameterList& parameters);

  //! Destructor
  virtual ~NLS_Newton();

  virtual void resetParameters(NLS_ParameterList& p);
  virtual STATUS getStatus();
  virtual STATUS iterate();
  virtual STATUS solve();
  virtual NLS_Group& getSolutionGroup() const;
  virtual bool getProfileInfo(string& name, NLS_Parameter& p) const;

protected:
  
  NLS_Group& soln;		//!< Solution
  NLS_Group& oldsoln;		//!< Old solution
  NLS_ParameterList& params;    //!< Parameter List
  NLS_LineSearchManager linesearch; //!< Line search manager
  NLS_StatusTest statustest;	// !< Status tester
  double step;			//!< Step length for last iteration
  double defaultstep;		//!< Default step length

  // Statistics
  int niter;                    //!< Number of nonlinear iterations

};


#endif

