// $Id$ 
// $Source$ 

// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#ifndef NLS_METHOD_MANAGER_H
#define NLS_METHOD_MANAGER_H

#include "NLS_Method.H"		// Base class

//! Concrete class to control the instatiation of the derived 
//! NLS_Method class.  
/*! 
  This is the class that takes control from the user code and hands it over
  to the nonlinear solver component. 
*/

#include "NLS_Utilities.H"
//-----------------------------------------------------------------------------
class NLS_MethodManager : public NLS_Method {

public:
  //! Generic interface to any of the various methods. 
  NLS_MethodManager(NLS_Group& initial_guess, 
		    NLS_Group& solution_space, 
		    NLS_ParameterList& p);

  //! Destructor
  ~NLS_MethodManager();


  void resetParameters(NLS_ParameterList& p);
  bool isConverged();
  int iterate();
  int solve();
  NLS_Group& getSolutionGroup() const;
  bool getProfileInfo(string& name, NLS_Parameter& p) const;

private:

  //! Pointer to the solver object
  NLS_Method *ptr;
};

#endif
