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

#ifndef NOX_MORETHUENTE_H
#define NOX_MORETHUENTE_H

#include "NOX_LineSearch_Generic.H" // base class

namespace NOX {
namespace LineSearch {

//! More'-Thuente Line Search.
/*!

  This code is based on the More'-Thuente line search from the 1983
  MINPACK Project. More specifically, this code is based on Dianne
  O'Leary's 1991 Matlab-implementation of the More'-Thuente line
  search. The original comments are preserved in the descriptions of
  the individual subroutines. What follows is an updated summary.

  The function we are minimizing is given by

  f(x) = 0.5 * F(x)'* F(x)

  The purpose of the More'-Thuente line search is to find a step which
  satisfies a sufficient decrease condition and a curvature condition.
  At each stage the subroutine updates an interval of uncertainty with
  endpoints stx and sty. The interval of uncertainty is initially
  chosen so that it contains a minimizer of the modified function

  f(x+stp*s) - f(x) - ftol*stp*(gradf(x)'s).

  If a step is obtained for which the modified function has a
  nonpositive function value and nonnegative derivative, then the
  interval of uncertainty is chosen so that it contains a minimizer of
  f(x+stp*s).

  The algorithm is designed to find a step which satisfies the
  sufficient decrease condition

  f(x+stp*s) <= f(x) + ftol*stp*(gradf(x)'s),

  and the curvature condition

  abs(gradf(x+stp*s)'s)) <= gtol*abs(gradf(x)'s).

  If ftol is less than gtol and if, for example, the function is
  bounded below, then there is always a step which satisfies both
  conditions. If no step can be found which satisfies both conditions,
  then the algorithm usually stops when rounding errors prevent
  further progress. In this case stp only satisfies the sufficient
  decrease condition.

  This line search can be called via NOX::LineSearch::Manager.
 
  The following parameters can be specified for this line search:
 
  - "Sufficient Decrease" - The ftol in the sufficient decrease
    condition (defaults to 1.0e-4)

  - "Curvature Condition" - The gtol in the curvature condition
    (defaults to 0.9999)

  - "Interval Width" - The maximum width of the interval containing
    the minimum of the modified function (defaults to 1.0e-15)

  - "Maximum Step" - maximum allowable step length (defaults to 1.0e6)

  - "Minimum Step" - minimum allowable step length (defaults to 1.0e-12)

  - "Max Iters" - maximum number of right-hand-side and corresponding
    Jacobian evaluations (defaults to 20)

  - "Default Step" - starting step length (defaults to 1.0)

  - "Recovery Step" - step to take when the line search fails
    (defaults to value for "Default Step")


*/

class MoreThuente : public Generic {

public:

  //! Constructor
  MoreThuente(Parameter::List& params);

  //! Destructor
  ~MoreThuente();
  
  // derived
  bool reset(Parameter::List& params);

  // derived
  bool operator()(Abstract::Group& newgrp, double& step, 
		  const Abstract::Group& oldgrp, 
		  const Abstract::Vector& dir) ;


private:

  //! More'-Thuente's cvsrch subroutine.
  /*!
    This translation into C++ is based on a translation into Matlab by Dianne O'Leary.
    The return code corresponds to the "info" code described below.
    The original documentation follows.

\verbatim
%   Translation of minpack subroutine cvsrch
%   Dianne O'Leary   July 1991
%     **********
%
%     Subroutine cvsrch
%
%     The purpose of cvsrch is to find a step which satisfies 
%     a sufficient decrease condition and a curvature condition.
%     The user must provide a subroutine which calculates the
%     function and the gradient.
%
%     At each stage the subroutine updates an interval of
%     uncertainty with endpoints stx and sty. The interval of
%     uncertainty is initially chosen so that it contains a 
%     minimizer of the modified function
%
%          f(x+stp*s) - f(x) - ftol*stp*(gradf(x)'s).
%
%     If a step is obtained for which the modified function 
%     has a nonpositive function value and nonnegative derivative, 
%     then the interval of uncertainty is chosen so that it 
%     contains a minimizer of f(x+stp*s).
%
%     The algorithm is designed to find a step which satisfies 
%     the sufficient decrease condition 
%
%           f(x+stp*s) <= f(x) + ftol*stp*(gradf(x)'s),
%
%     and the curvature condition
%
%           abs(gradf(x+stp*s)'s)) <= gtol*abs(gradf(x)'s).
%
%     If ftol is less than gtol and if, for example, the function
%     is bounded below, then there is always a step which satisfies
%     both conditions. If no step can be found which satisfies both
%     conditions, then the algorithm usually stops when rounding
%     errors prevent further progress. In this case stp only 
%     satisfies the sufficient decrease condition.
%
%     The subroutine statement is
%
%        subroutine cvsrch(fcn,n,x,f,g,s,stp,ftol,gtol,xtol,
%                          stpmin,stpmax,maxfev,info,nfev,wa)
%     where
%
%	fcn is the name of the user-supplied subroutine which
%         calculates the function and the gradient.  fcn must 
%      	  be declared in an external statement in the user 
%         calling program, and should be written as follows.
%
%	  subroutine fcn(n,x,f,g)
%         integer n
%         f
%         x(n),g(n)
%	  ----------
%         Calculate the function at x and
%         return this value in the variable f.
%         Calculate the gradient at x and
%         return this vector in g.
%	  ----------
%	  return
%	  end
%
%       n is a positive integer input variable set to the number
%	  of variables.
%
%	x is an array of length n. On input it must contain the
%	  base point for the line search. On output it contains 
%         x + stp*s.
%
%	f is a variable. On input it must contain the value of f
%         at x. On output it contains the value of f at x + stp*s.
%
%	g is an array of length n. On input it must contain the
%         gradient of f at x. On output it contains the gradient
%         of f at x + stp*s.
%
%	s is an input array of length n which specifies the
%         search direction.
%
%	stp is a nonnegative variable. On input stp contains an
%         initial estimate of a satisfactory step. On output
%         stp contains the final estimate.
%
%       ftol and gtol are nonnegative input variables. Termination
%         occurs when the sufficient decrease condition and the
%         directional derivative condition are satisfied.
%
%	xtol is a nonnegative input variable. Termination occurs
%         when the relative width of the interval of uncertainty 
%	  is at most xtol.
%
%	stpmin and stpmax are nonnegative input variables which 
%	  specify lower and upper bounds for the step.
%
%	maxfev is a positive integer input variable. Termination
%         occurs when the number of calls to fcn is at least
%         maxfev by the end of an iteration.
%
%	info is an integer output variable set as follows:
%	  
%	  info = 0  Improper input parameters.
%
%	  info = 1  The sufficient decrease condition and the
%                   directional derivative condition hold.
%
%	  info = 2  Relative width of the interval of uncertainty
%		    is at most xtol.
%
%	  info = 3  Number of calls to fcn has reached maxfev.
%
%	  info = 4  The step is at the lower bound stpmin.
%
%	  info = 5  The step is at the upper bound stpmax.
%
%	  info = 6  Rounding errors prevent further progress.
%                   There may not be a step which satisfies the
%                   sufficient decrease and curvature conditions.
%                   Tolerances may be too small.
%
%       nfev is an integer output variable set to the number of
%         calls to fcn.
%
%	wa is a work array of length n.
%
%     Subprograms called
%
%	user-supplied......fcn
%
%	MINPACK-supplied...cstep
%
%	FORTRAN-supplied...abs,max,min
%	  
%     Argonne National Laboratory. MINPACK Project. June 1983
%     Jorge J. More', David J. Thuente
%
\endverbatim


   */
  int cvsrch(Abstract::Group& newgrp, double& step, 
	     const Abstract::Group& oldgrp, const Abstract::Vector& dir);


  //! More'-Thuente's cstep subroutine.
  /*!
    \verbatim
%   Translation of minpack subroutine cstep 
%   Dianne O'Leary   July 1991
%     **********
%
%     Subroutine cstep
%
%     The purpose of cstep is to compute a safeguarded step for
%     a linesearch and to update an interval of uncertainty for
%     a minimizer of the function.
%
%     The parameter stx contains the step with the least function
%     value. The parameter stp contains the current step. It is
%     assumed that the derivative at stx is negative in the
%     direction of the step. If brackt is set true then a
%     minimizer has been bracketed in an interval of uncertainty
%     with endpoints stx and sty.
%
%     The subroutine statement is
%
%       subroutine cstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,
%                        stpmin,stpmax,info)
% 
%     where
%
%       stx, fx, and dx are variables which specify the step,
%         the function, and the derivative at the best step obtained
%         so far. The derivative must be negative in the direction
%         of the step, that is, dx and stp-stx must have opposite 
%         signs. On output these parameters are updated appropriately.
%
%       sty, fy, and dy are variables which specify the step,
%         the function, and the derivative at the other endpoint of
%         the interval of uncertainty. On output these parameters are 
%         updated appropriately.
%
%       stp, fp, and dp are variables which specify the step,
%         the function, and the derivative at the current step.
%         If brackt is set true then on input stp must be
%         between stx and sty. On output stp is set to the new step.
%
%       brackt is a logical variable which specifies if a minimizer
%         has been bracketed. If the minimizer has not been bracketed
%         then on input brackt must be set false. If the minimizer
%         is bracketed then on output brackt is set true.
%
%       stpmin and stpmax are input variables which specify lower 
%         and upper bounds for the step.
%
%       info is an integer output variable set as follows:
%         If info = 1,2,3,4,5, then the step has been computed
%         according to one of the five cases below. Otherwise
%         info = 0, and this indicates improper input parameters.
%
%     Subprograms called
%
%       FORTRAN-supplied ... abs,max,min,sqrt
%                        ... dble
%
%     Argonne National Laboratory. MINPACK Project. June 1983
%     Jorge J. More', David J. Thuente
    \endverbatim
  */
  int cstep(double& stx, double& fx, double& dx,
	    double& sty, double& fy, double& dy,
	    double& stp, double& fp, double& dp,
	    bool& brackt, double stmin, double stmax);

  //! max
  double max(double a, double b);
  //! min
  double min(double a, double b);
  //! absmax - returns the max of the absolute value of the three inputs
  double absmax(double a, double b, double c);


  //! Computes d'*g where d is the given search direction (dir) and g
  //! is the gradient corresponding to the given group (grp). If the
  //! gradient has not been computed, then we compute F'*(J*d)
  //! instead. Here F is the right-hand-side and J is the Jacobian
  //! corresponding to the given group. 
  double dgcompute(const Abstract::Vector& dir, const Abstract::Group& grp);


private:

  //! Sufficient decrease parameter
  double ftol;
  //! Curvature condition parameter
  double gtol;
  //! Threshold on width of the uncertainty
  double xtol;
  //! Minimum step
  double stpmin;
  //! Maximum step
  double stpmax;
  //! Maximum number of function evaluations
  int maxfev;

  //! Default step length
  double defaultstep;
  //! Recovery step if linesearch fails
  double recoverystep;

  //! Temporary abstract vector pointer
  Abstract::Vector* tmpvecptr;


};
} // namespace LineSearch
} // namespace NOX
#endif
