
// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#ifndef _NLS_FORCINGTERM_H
#define _NLS_FORCINGTERM_H
  
#include "NLS_ParameterList.H"
#include "NLS_Group.H"

//! Calculates the linear solver tolerance using an various adaptive
//! strategies (Walker et al. SIAM J. Scientific Computing 1996).
/*!
 * This routine calculates the residual reduction convergence criterion eta_k
 * (the forcing term). The algorithm determines the criterion based on how
 * close the linear model approximates the nonlinear function.  The algorithm
 * is taken from:
 *
 * S. C. Eisenstat and H. F. Walker, Choosing the forcing terms in an
 * inexact Newton methods, SIAM J. Scientific Computing, 17, (1996), pp 16-32
 * Authors: J. N. Shadid
 *          H. F. Walker
 */

class NLS_ForcingTerm {
public:

  //! Constructor
  NLS_ForcingTerm();

  //! Destructor
  ~NLS_ForcingTerm();

  //! Set the parameters in the utilities class.
  double getForcingTerm(NLS_ParameterList& p);

protected:
  double eta_k;
  double eta_km1;
  double eta_min;
  double eta_max;
  
  double alpha;
  double eta_k_alpha;
  double gamma;
  double residual_ratio;

  double residual_norm_k;
  double residual_norm_km1;
  double residual_norm_linear_model;
};

#endif
