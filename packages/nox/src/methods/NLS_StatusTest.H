// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#ifndef _NLS_STATUSTEST_H
#define _NLS_STATUSTEST_H
  
#include "NLS_ParameterList.H"
#include "NLS_Vector.H"
#include "NLS_Group.H"
#include "NLS_Method.H"

class NLS_StatusTest {

public:

  NLS_StatusTest();
  virtual ~NLS_StatusTest();
  virtual void setup(const NLS_ParameterList& p, double normrhs);
  virtual NLS_Method::STATUS operator()(const NLS_Group& oldsoln, const NLS_Group& soln, 
				      const NLS_Vector& updatedir, double step, int niter);
  
private:
  
  int maxiter;			//<! Max number of nonlinear solver iterations

  double absresidtol;		//<! Residual norm tolerance
  double relresidtol;		//<! Relative (R/R0) residual norm tolerance
  double maxresidtol;		//<! Residual max norm tolerance
  double ratioresidtol;		//<! Ratio of old norm RHS to normRHS lower bound

  double absupdatetol;		//<! Update norm tolerance
  double relupdate_epsilon_r;	//<! Relative update formula variable
  double relupdate_epsilon_a;	//<! Relative update formula variable
  
  double normrhszero;

};

#endif
