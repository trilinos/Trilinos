#ifndef _NLS_NEWTON_H
#define _NLS_NEWTON_H

#include "NLS_Method.H"

//! Implementation of plain Newton solver
class NLS_Newton : public NLS_Method {

public:

  //! Newton-method constructor
  NLS_Newton(NLS_AbstractFunction& function, NLS_ParameterList& params);

  void resetParameters(NLS_ParameterList& p);
  int iterate();
  int solve();
  bool isConverged();
  NLS_Group& getSolutionGroup() const;
  bool getProfileInfo(string& name, Parameter& p) const;

}
private:
  
  NLS_Function& f;		//!< Problem 
  NLS_Group& soln;		//!< Solution
  NLS_Group& oldsoln;		//!< Old solution
  double normrhs;		//!< Norm of RHS
  int niter;			//!< Number of nonlinear iterations
};


#endif

