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

#include "NOX_LineSearch_MoreThuente.H"	// class definition

#include "NOX_Common.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Parameter_List.H"
#include "NOX_Utils.H"

using namespace NOX;
using namespace NOX::LineSearch;

MoreThuente::MoreThuente(Parameter::List& params) 
{
  reset(params);
}

MoreThuente::~MoreThuente()
{
}

bool MoreThuente::reset(Parameter::List& params)
{ 
  ftol = params.getParameter("Sufficient Decrease", 1.0e-4);
  gtol = params.getParameter("Curvature Condition", 0.9999);
  xtol = params.getParameter("Interval Width", 1.0e-15);
  stpmin = params.getParameter("Minimum Step", 1.0e-12);
  stpmax = params.getParameter("Maximum Step", 1.0e+6);
  maxfev = params.getParameter("Max Iters", 20);
  defaultstep = params.getParameter("Default Step", 1.0);
  recoverystep = params.getParameter("Recovery Step", defaultstep);

  // Check the input parameters for errors.
  if ((ftol < 0.0) || (gtol < 0.0) || 
      (xtol < 0.0) || (stpmin < 0.0) || (stpmax < stpmin) || 
      (maxfev <= 0) || (defaultstep <= 0)) {
    cout << "NOX::LineSearch::MoreThuente::reset - Error in Input Parameter!" << endl;
    throw "NOX Error";
  }

  return true;
}


bool MoreThuente::compute(Abstract::Group& grp, double& step, 
			  const Abstract::Vector& dir,
			  const Solver::Generic& s) 
{
  const Abstract::Group& oldGrp = s.getPreviousSolutionGroup();
  int info = cvsrch(grp, step, oldGrp, dir);
  return (info == 1);
}

int MoreThuente::cvsrch(Abstract::Group& newgrp, double& stp, 
			const Abstract::Group& oldgrp, const Abstract::Vector& dir)
{
  if (Utils::doPrint(Utils::InnerIteration)) {
   cout << "\n" << Utils::fill(72) << "\n" << "-- More'-Thuente Line Search -- \n";
  }

  // Set default step
  stp = defaultstep;

  int info = 0;			// return code
  int infoc = 1;		// return code for subroutine cstep

  // Compute the initial gradient in the search direction and check
  // that s is a descent direction.

  double dginit = computeSlope(dir, oldgrp);

  if (dginit >= 0.0) {
    if (Utils::doPrint(Utils::Warning)) {
      cout << "NOX::LineSearch::MoreThuente::cvsrch - Non-descent direction (dginit = " << dginit << ")" << endl;
    }
    stp = recoverystep;
    newgrp.computeX(oldgrp, dir, stp);
    return 7;
  }

  // Initialize local variables.

  bool brackt = false;		// has the soln been bracketed?
  bool stage1 = true;		// are we in stage 1?
  int nfev = 0;			// number of function evaluations
  double finit = 0.5 * oldgrp.getNormF() * oldgrp.getNormF(); // initial function value
  double dgtest = ftol * dginit; // f for curvature condition
  double width = stpmax - stpmin; // interval width
  double width1 = 2 * width;	// ???

  // The variables stx, fx, dgx contain the values of the step,
  // function, and directional derivative at the best step.  The
  // variables sty, fy, dgy contain the value of the step, function,
  // and derivative at the other endpoint of the interval of
  // uncertainty.  The variables stp, f, dg contain the values of the
  // step, function, and derivative at the current step.

  double stx = 0.0;
  double fx = finit;
  double dgx = dginit;
  double sty = 0.0;
  double fy = finit;
  double dgy = dginit;

  // Start of iteration.

  double stmin, stmax;
  double fm, fxm, fym, dgm, dgxm, dgym;

  while (1) {

    // Set the minimum and maximum steps to correspond to the present
    // interval of uncertainty.

    if (brackt) {
      stmin = min(stx, sty);
      stmax = max(stx, sty);
    }
    else {
      stmin = stx;
      stmax = stp + 4 * (stp - stx);
    }
    
    // Force the step to be within the bounds stpmax and stpmin.
    stp = max(stp, stpmin);
    stp = min(stp, stpmax);

    // If an unusual termination is to occur then let stp be the
    // lowest point obtained so far.

    if ((brackt && ((stp <= stmin) || (stp >= stmax))) ||
	(nfev >= maxfev - 1) || (infoc == 0) ||
	(brackt && (stmax - stmin <= xtol * stmax))) {
      stp = stx;
    }

    // Evaluate the function and gradient at stp
    // and compute the directional derivative.

    newgrp.computeX(oldgrp, dir, stp);
    newgrp.computeF();
    double f = 0.5 * newgrp.getNormF() * newgrp.getNormF();
    newgrp.computeJacobian();
    nfev ++;
    string message = "";

    double dg = computeSlope(dir, newgrp);

    double ftest1 = finit + stp * dgtest;

    // Test for convergence.

    if ((brackt && ((stp <= stmin) || (stp >= stmax))) || (infoc == 0))	
      info = 6;			// Rounding errors

    if ((stp == stpmax) && (f <= ftest1) && (dg <= dgtest))
      info = 5;			// stp=stpmax

    if ((stp == stpmin) && ((f > ftest1) || (dg >= dgtest))) 
      info = 4;			// stp=stpmin

    if (nfev >= maxfev) 
      info = 3;			// max'd out on fevals

    if (brackt && (stmax-stmin <= xtol*stmax)) 
      info = 2;			// bracketed soln


    //cout << "f=" << f << " ftest1=" << ftest1 << " fabs(dg)=" << fabs(dg) 
    //	 << " gtol*(-dginit)=" << gtol*(-dginit) << endl;

    if ((f <= ftest1) && (fabs(dg) <= gtol*(-dginit))) 
      info = 1;			// Success!!!!

    if (info != 0) {		// Line search is done

      if (info != 1) {		// Line search failed 
	
	stp = recoverystep;
	newgrp.computeX(oldgrp, dir, stp);
	
	message = "(USING RECOVERY STEP!)";
	
	/*
	if (Utils::doPrint(Utils::Details))
	  message += "[Failure info flag = " + info + "]";
	*/
	    
      }
      else {			// Line search succeeded
	
	message = "(STEP ACCEPTED!)";
	
      }
      
      printStep(nfev, stp, finit, f, message);

      // Returning the line search flag
      return info;

    } // info != 0
    
    printStep(nfev, stp, finit, f, message);


    // In the first stage we seek a step for which the modified
    // function has a nonpositive value and nonnegative derivative.

    if (stage1 && (f <= ftest1) && (dg >= min(ftol, gtol) * dginit)) 
      stage1 = false;

    // A modified function is used to predict the step only if we have
    // not obtained a step for which the modified function has a
    // nonpositive function value and nonnegative derivative, and if a
    // lower function value has been obtained but the decrease is not
    // sufficient.

    if (stage1 && (f <= fx) && (f > ftest1)) {

      // Define the modified function and derivative values.

      fm = f - stp * dgtest;
      fxm = fx - stx * dgtest;
      fym = fy - sty * dgtest;
      dgm = dg - dgtest;
      dgxm = dgx - dgtest;
      dgym = dgy - dgtest;

      // Call cstep to update the interval of uncertainty 
      // and to compute the new step.

      infoc = cstep(stx,fxm,dgxm,sty,fym,dgym,stp,fm,dgm, 
		    brackt,stmin,stmax);

      // Reset the function and gradient values for f.

      fx = fxm + stx*dgtest;
      fy = fym + sty*dgtest;
      dgx = dgxm + dgtest;
      dgy = dgym + dgtest;

    }

    else {

      // Call cstep to update the interval of uncertainty 
      // and to compute the new step.

      infoc = cstep(stx,fx,dgx,sty,fy,dgy,stp,f,dg,
		    brackt,stmin,stmax);

    }

    // Force a sufficient decrease in the size of the
    // interval of uncertainty.

    if (brackt) {
      if (fabs(sty - stx) >= 0.66 * width1) 
	stp = stx + 0.5 * (sty - stx);
      width1 = width;
      width = fabs(sty-stx);
    }

  } // while-loop
}


int MoreThuente::cstep(double& stx, double& fx, double& dx,
		       double& sty, double& fy, double& dy,
		       double& stp, double& fp, double& dp,
		       bool& brackt, double stmin, double stmax)
{
  int info = 0;

  // Check the input parameters for errors.

  if ((brackt && ((stp <= min(stx, sty)) || (stp >= max(stx, sty)))) || 
      (dx * (stp - stx) >= 0.0) || (stmax < stmin))
    return info;

  // Determine if the derivatives have opposite sign.

  double sgnd = dp * (dx / fabs(dx));

  // First case. A higher function value.  The minimum is
  // bracketed. If the cubic step is closer to stx than the quadratic
  // step, the cubic step is taken, else the average of the cubic and
  // quadratic steps is taken.

  bool bound;
  double theta;
  double s;
  double gamma;
  double p,q,r;
  double stpc, stpq, stpf;

  if (fp > fx) {
    info = 1;
    bound = 1;
    theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
    s = absmax(theta, dx, dp);
    gamma = s * sqrt(((theta / s) * (theta / s)) - (dx / s) * (dp / s));
    if (stp < stx) 
      gamma = -gamma;

    p = (gamma - dx) + theta;
    q = ((gamma - dx) + gamma) + dp;
    r = p / q;
    stpc = stx + r * (stp - stx);
    stpq = stx + ((dx / ((fx - fp) / (stp - stx) + dx)) / 2) * (stp - stx);
    if (fabs(stpc - stx) < fabs(stpq - stx)) 
      stpf = stpc;
    else 
      stpf = stpc + (stpq - stpc) / 2;

    brackt = true;
  }

  // Second case. A lower function value and derivatives of opposite
  // sign. The minimum is bracketed. If the cubic step is closer to
  // stx than the quadratic (secant) step, the cubic step is taken,
  // else the quadratic step is taken.

  else if (sgnd < 0.0) {
    info = 2;
    bound = false;
    theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
    s = absmax(theta,dx,dp);
    gamma = s * sqrt(((theta/s) * (theta/s)) - (dx / s) * (dp / s));
    if (stp > stx) 
      gamma = -gamma;
    p = (gamma - dp) + theta;
    q = ((gamma - dp) + gamma) + dx;
    r = p / q;
    stpc = stp + r * (stx - stp);
    stpq = stp + (dp / (dp - dx)) * (stx - stp);
    if (fabs(stpc - stp) > fabs(stpq - stp))
      stpf = stpc;
    else
      stpf = stpq;
    brackt = true;
  }

  // Third case. A lower function value, derivatives of the same sign,
  // and the magnitude of the derivative decreases.  The cubic step is
  // only used if the cubic tends to infinity in the direction of the
  // step or if the minimum of the cubic is beyond stp. Otherwise the
  // cubic step is defined to be either stmin or stmax. The
  // quadratic (secant) step is also computed and if the minimum is
  // bracketed then the the step closest to stx is taken, else the
  // step farthest away is taken.

  else if (fabs(dp) < fabs(dx)) {
    info = 3;
    bound = true;
    theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
    s = absmax(theta, dx, dp);

    // The case gamma = 0 only arises if the cubic does not tend
    // to infinity in the direction of the step.

    gamma = s * sqrt(max(0,(theta / s) * (theta / s) - (dx / s) * (dp / s)));
    if (stp > stx) 
      gamma = -gamma;
      
    p = (gamma - dp) + theta;
    q = (gamma + (dx - dp)) + gamma;
    r = p / q;
    if ((r < 0.0) && (gamma != 0.0))
      stpc = stp + r * (stx - stp);
    else if (stp > stx)
      stpc = stmax;
    else
      stpc = stmin;
      
    stpq = stp + (dp/ (dp - dx)) * (stx - stp);
    if (brackt) {
      if (fabs(stp - stpc) < fabs(stp - stpq)) 
	stpf = stpc;
      else
	stpf = stpq;
    }
    else {
      if (fabs(stp - stpc) > fabs(stp - stpq)) 
	stpf = stpc;
      else
	stpf = stpq;
    }
  }

  // Fourth case. A lower function value, derivatives of the same
  // sign, and the magnitude of the derivative does not decrease. If
  // the minimum is not bracketed, the step is either stmin or
  // stmax, else the cubic step is taken.

  else {
    info = 4;
    bound = false;
    if (brackt) {
      theta = 3 * (fp - fy) / (sty - stp) + dy + dp;
      s = absmax(theta, dy, dp);
      gamma = s * sqrt(((theta/s)*(theta/s)) - (dy / s) * (dp / s));
      if (stp > sty) 
	gamma = -gamma;
      p = (gamma - dp) + theta;
      q = ((gamma - dp) + gamma) + dy;
      r = p / q;
      stpc = stp + r * (sty - stp);
      stpf = stpc;
    }
    else if (stp > stx)
      stpf = stmax;
    else
      stpf = stmin;
  }

  // Update the interval of uncertainty. This update does not depend
  // on the new step or the case analysis above.

  if (fp > fx) {
    sty = stp;
    fy = fp;
    dy = dp;
  }
  else {
    if (sgnd < 0.0) {
      sty = stx;
      fy = fx;
      dy = dx;
    }
    stx = stp;
    fx = fp;
    dx = dp;
  }

  // Compute the new step and safeguard it.

  stpf = min(stmax, stpf);
  stpf = max(stmin, stpf);
  stp = stpf;
  if (brackt && bound) {
    if (sty > stx) 
      stp = min(stx + 0.66 * (sty - stx), stp);
    else 
      stp = max(stx + 0.66 * (sty - stx), stp);
  }

  return info;

}

double MoreThuente::min(double a, double b)
{
  return (a < b ? a : b);
}

double MoreThuente::max(double a, double b)
{
  return (a > b ? a : b);
}


double MoreThuente::absmax(double a, double b, double c)
{
  a = fabs(a);
  b = fabs(b);
  c = fabs(c);

  if (a > b)
    return (a > c) ? a : c;
  else
    return (b > c) ? b : c;
}




