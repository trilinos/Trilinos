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
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include "NOX_LineSearch_MoreThuente.H"	// class definition

#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "Teuchos_ParameterList.hpp"
#include "NOX_Utils.H"
#include "NOX_MeritFunction_Generic.H"
#include "NOX_GlobalData.H"

NOX::LineSearch::MoreThuente::
MoreThuente(const Teuchos::RCP<NOX::GlobalData>& gd, 
	    Teuchos::ParameterList& params) :
  globalDataPtr(gd),
  print(gd->getUtils()),
  slope(gd),
  paramsPtr(0)
{
  reset(gd, params);
}

NOX::LineSearch::MoreThuente::~MoreThuente()
{
}

bool NOX::LineSearch::MoreThuente::
reset(const Teuchos::RCP<NOX::GlobalData>& gd,
      Teuchos::ParameterList& params)
{ 
  globalDataPtr = gd;
  meritFuncPtr = gd->getMeritFunction();
  print.reset(gd->getUtils());
  slope.reset(gd);

  paramsPtr = &params;
  Teuchos::ParameterList& p = params.sublist("More'-Thuente");
  ftol = p.get("Sufficient Decrease", 1.0e-4);
  gtol = p.get("Curvature Condition", 0.9999);
  xtol = p.get("Interval Width", 1.0e-15);
  stpmin = p.get("Minimum Step", 1.0e-12);
  stpmax = p.get("Maximum Step", 1.0e+6);
  maxfev = p.get("Max Iters", 20);
  defaultstep = p.get("Default Step", 1.0);
  recoverystep = p.get("Recovery Step", defaultstep);

  // Check the input parameters for errors.
  if ((ftol < 0.0) || 
      (gtol < 0.0) || 
      (xtol < 0.0) || 
      (stpmin < 0.0) || 
      (stpmax < stpmin) || 
      (maxfev <= 0) || 
      (defaultstep <= 0)) 
  {
    print.out() << "NOX::LineSearch::MoreThuente::reset - Error in Input Parameter!" << std::endl;
    throw "NOX Error";
  }

  counter.reset();

  std::string choice = p.get("Sufficient Decrease Condition", "Armijo-Goldstein");
  if (choice == "Ared/Pred") 
    suffDecrCond = AredPred;
  else if (choice == "Armijo-Goldstein") 
    suffDecrCond = ArmijoGoldstein;
  else {
    print.out() << "ERROR: NOX::LineSearch::MoreThuente::reset() - the choice of "
	 << "\"Sufficient Decrease Condition\" is invalid." << std::endl;
    throw "NOX Error";
  }

  choice = p.get("Recovery Step Type", "Constant");
  if (choice == "Constant")
    recoveryStepType = Constant;
  else if (choice == "Last Computed Step") {
    recoveryStepType = LastComputedStep;
  }
  else {
    print.out() << "NOX::LineSearch::MoreThuente::reset - Invalid "
	 << "\"Recovery Step Type\"" << std::endl;
    throw "NOX Error";
  }

  useOptimizedSlopeCalc = p.get("Optimize Slope Calculation", false);

  return true;
}


bool NOX::LineSearch::MoreThuente::compute(Abstract::Group& grp, double& step, 
			  const Abstract::Vector& dir,
			  const Solver::Generic& s) 
{
  counter.incrementNumLineSearches();
  const Abstract::Group& oldGrp = s.getPreviousSolutionGroup();
  int info = cvsrch(grp, step, oldGrp, dir, s);

  if (step != 1.0)
    counter.incrementNumNonTrivialLineSearches();    

  counter.setValues(*paramsPtr);

  return (info == 1);
}

int NOX::LineSearch::MoreThuente::
cvsrch(Abstract::Group& newgrp, double& stp, const Abstract::Group& oldgrp, 
       const Abstract::Vector& dir, const Solver::Generic& s)
{
  if (print.isPrintType(NOX::Utils::InnerIteration)) 
  {
   print.out() << "\n" << NOX::Utils::fill(72) << "\n" << "-- More'-Thuente Line Search -- \n";
  }

  // Set default step
  stp = defaultstep;

  int info = 0;			// return code
  int infoc = 1;		// return code for subroutine cstep

  // Compute the initial gradient in the search direction and check
  // that s is a descent direction.
  double dginit = 0.0;
  if (useOptimizedSlopeCalc)
    dginit = slope.computeSlopeWithOutJac(dir, oldgrp);
  else
    dginit = slope.computeSlope(dir, oldgrp);
  
  if (dginit >= 0.0) 
  {
    if (print.isPrintType(NOX::Utils::Warning)) 
    {
      print.out() << "NOX::LineSearch::MoreThuente::cvsrch - Non-descent direction (dginit = " << dginit << ")" << std::endl;
    }
    stp = recoverystep;
    newgrp.computeX(oldgrp, dir, stp);
    return 7;
  }

  // Initialize local variables.

  bool brackt = false;		// has the soln been bracketed?
  bool stage1 = true;		// are we in stage 1?
  int nfev = 0;			// number of function evaluations
  double dgtest = ftol * dginit; // f for curvature condition
  double width = stpmax - stpmin; // interval width
  double width1 = 2 * width;	// ???

  // initial function value
  double finit = meritFuncPtr->computef(oldgrp);

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

  // Get the linear solve tolerance for adjustable forcing term
  const Teuchos::ParameterList& p = s.getList();
  double eta_original = 0.0;
  double eta = 0.0;
  eta_original = const_cast<Teuchos::ParameterList&>(p).sublist("Direction").
    sublist("Newton").sublist("Linear Solver").
    get(std::string("Tolerance"), -1.0);
  eta = eta_original;
  

  // Start of iteration.

  double stmin, stmax;
  double fm, fxm, fym, dgm, dgxm, dgym;

  while (1) 
  {

    // Set the minimum and maximum steps to correspond to the present
    // interval of uncertainty.

    if (brackt) 
    {
      stmin = min(stx, sty);
      stmax = max(stx, sty);
    }
    else 
    {
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
	(brackt && (stmax - stmin <= xtol * stmax))) 
    {
      stp = stx;
    }

    // Evaluate the function and gradient at stp
    // and compute the directional derivative.

    newgrp.computeX(oldgrp, dir, stp);

    NOX::Abstract::Group::ReturnType rtype;
    rtype = newgrp.computeF();
    if (rtype != NOX::Abstract::Group::Ok) 
    {
      print.err() << "NOX::LineSearch::MoreThuente::cvrch - Unable to compute F" << std::endl;
      throw "NOX Error";
    }

    double f = meritFuncPtr->computef(newgrp);

    if (!useOptimizedSlopeCalc) {

      rtype = newgrp.computeJacobian();
      if (rtype != NOX::Abstract::Group::Ok) 
	{
	  print.err() << "NOX::LineSearch::MoreThuente::cvrch - Unable to compute Jacobian" << std::endl;
	  throw "NOX Error";
	}

      rtype = newgrp.computeGradient();
      if (rtype != NOX::Abstract::Group::Ok) 
	{
	  print.err() << "NOX::LineSearch::MoreThuente::cvrch - Unable to compute Gradient" << std::endl;
	  throw "NOX Error";
	}
    }

    nfev ++;
    std::string message = "";

    double dg = 0.0;
    if (useOptimizedSlopeCalc)
      dg = slope.computeSlopeWithOutJac(dir, newgrp);
    else 
      dg = slope.computeSlope(dir, newgrp);

    // Armijo-Goldstein sufficient decrease
    double ftest1 = finit + stp * dgtest;

    // Ared/Pred suffiecient decrease (could use a user defined norm)
    double ftest2 = 0.0;
    ftest2 = oldgrp.getNormF() * (1.0 - ftol * (1.0 - eta));

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


    //print.out() << "f=" << f << " ftest1=" << ftest1 << " fabs(dg)=" << fabs(dg) 
    //	 << " gtol*(-dginit)=" << gtol*(-dginit) << std::endl;

    // RPP sufficient decrease test can be different
    bool sufficientDecreaseTest = false;
    if (suffDecrCond == ArmijoGoldstein)
      sufficientDecreaseTest = (f <= ftest1);  // Armijo-Golstein
    else {
      double ap_normF = 0.0;
      ap_normF = newgrp.getNormF();

      sufficientDecreaseTest = (ap_normF <= ftest2); // Ared/Pred
    }

    if ((sufficientDecreaseTest) && (fabs(dg) <= gtol*(-dginit))) 
      info = 1;			// Success!!!!

    if (info != 0) 		// Line search is done
    {
      if (info != 1) 		// Line search failed 
      {
	// RPP add
	counter.incrementNumFailedLineSearches();
	
	if (recoveryStepType == Constant)
	  stp = recoverystep;

	newgrp.computeX(oldgrp, dir, stp);
	
	message = "(USING RECOVERY STEP!)";
	
	/*
	if (print.isPrintType(Utils::Details))
	  message += "[Failure info flag = " + info + "]";
	*/
	    
      }
      else 			// Line search succeeded
      {
	message = "(STEP ACCEPTED!)";
      }
      
      print.printStep(nfev, stp, finit, f, message);

      // Returning the line search flag
      return info;

    } // info != 0
    
    print.printStep(nfev, stp, finit, f, message);

    // RPP add
    counter.incrementNumIterations();

    // In the first stage we seek a step for which the modified
    // function has a nonpositive value and nonnegative derivative.

    if (stage1 && (f <= ftest1) && (dg >= min(ftol, gtol) * dginit)) 
      stage1 = false;

    // A modified function is used to predict the step only if we have
    // not obtained a step for which the modified function has a
    // nonpositive function value and nonnegative derivative, and if a
    // lower function value has been obtained but the decrease is not
    // sufficient.

    if (stage1 && (f <= fx) && (f > ftest1)) 
    {

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

    else 
    {

      // Call cstep to update the interval of uncertainty 
      // and to compute the new step.

      infoc = cstep(stx,fx,dgx,sty,fy,dgy,stp,f,dg,
		    brackt,stmin,stmax);

    }

    // Force a sufficient decrease in the size of the
    // interval of uncertainty.

    if (brackt) 
    {
      if (fabs(sty - stx) >= 0.66 * width1) 
	stp = stx + 0.5 * (sty - stx);
      width1 = width;
      width = fabs(sty-stx);
    }

  } // while-loop

}


int NOX::LineSearch::MoreThuente::cstep(double& stx, double& fx, double& dx,
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

  if (fp > fx) 
  {
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

  else if (sgnd < 0.0) 
  {
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

  else if (fabs(dp) < fabs(dx)) 
  {
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
    if (brackt) 
    {
      if (fabs(stp - stpc) < fabs(stp - stpq)) 
	stpf = stpc;
      else
	stpf = stpq;
    }
    else 
    {
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
    if (brackt) 
    {
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

  if (fp > fx) 
  {
    sty = stp;
    fy = fp;
    dy = dp;
  }
  else 
  {
    if (sgnd < 0.0) 
    {
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
  if (brackt && bound) 
  {
    if (sty > stx) 
      stp = min(stx + 0.66 * (sty - stx), stp);
    else 
      stp = max(stx + 0.66 * (sty - stx), stp);
  }

  return info;

}

double NOX::LineSearch::MoreThuente::min(double a, double b)
{
  return (a < b ? a : b);
}

double NOX::LineSearch::MoreThuente::max(double a, double b)
{
  return (a > b ? a : b);
}


double NOX::LineSearch::MoreThuente::absmax(double a, double b, double c)
{
  a = fabs(a);
  b = fabs(b);
  c = fabs(c);

  if (a > b)
    return (a > c) ? a : c;
  else
    return (b > c) ? b : c;
}


