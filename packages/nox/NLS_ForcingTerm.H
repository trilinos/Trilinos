
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

//! Functor to calculate the linear solver tolerance using an various
//! adaptive strategies.
/*!
 *
 * This routine calculates the residual reduction convergence
 * tolerance (the forcing term). The algorithm determines the
 * criterion based on how close the linear model approximates the
 * nonlinear function.  The algorithm is taken from:
 * S. C. Eisenstat and H. F. Walker, Choosing the forcing terms in an
 * inexact Newton methods, SIAM J. Scientific Computing, 17, (1996),
 * pp 16-32 
 *
 * Authors: J. N. Shadid, H. F. Walker */

class NLS_ForcingTerm {

public:

  //! Calculate the forcing term.
  /*! The normrhs is the norm of the current solution, the normoldrhs
    is the norm of the previous solution, and the normpredrhs is the
    norm of the predicted RHS. The last value does not need to be
    passed unless the "Forcing Term Method" is "Type 1". */
  double operator()(NLS_ParameterList& p, double normrhs, double normoldrhs, double normpredrhs = 0.0);

};

#endif
