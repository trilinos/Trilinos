// $Id$ 
// $Source$ 

// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#ifndef _NLS_METHOD_H
#define _NLS_METHOD_H

#include <string>		// argument to function
#include "NLS_ParameterList.H"	// argument to function
#include "NLS_Group.H"		// return value from function

//! Abstract nonlinear solver method interface.
class NLS_Method 
{

public:

  //! Status of solver.
  /*! Positive values indicate convergence, negative values indicate
    failure, and zero indicates neutral (i.e., neither converged nor
    failed). */
  enum STATUS {
    NotConverged = 0,
    Converged = 1,
    ConvAbsResidTol = 2,
    ConvRelResidTol = 3,
    ConvMaxResidTol = 4,
    ConvAbsUpdateTol = 5,
    ConvRelUpdateTol = 6,
    Failed = -1,
    MaxItersExceeded = -2,
    TooSlowConv = -3
  };

  //! Constructor 
  NLS_Method() {};

  //! Deconstructor
  virtual ~NLS_Method() {};

  //! Change parameters for the nonlinear method.
  virtual void resetParameters(NLS_ParameterList& p) = 0;

  //! Check current convergence/failure status. 
  virtual STATUS getStatus() = 0;

  //! Do one nonlinear method iteration and return status.
  virtual STATUS iterate() = 0;

  //! Solve the nonlinear problem and return final status.
  virtual STATUS solve() = 0;

  //! Return a reference to the current solution group.
  virtual NLS_Group& getSolutionGroup() const = 0;

  //! Fill in the parameter with the named profile information.
  virtual bool getProfileInfo(string& name, NLS_Parameter& p) const = 0;

};



#endif
