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

#include "NOX_LineSearch_Polynomial.H"

#include "NOX_Common.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Parameter_List.H"
#include "NOX_Utils.H"

using namespace NOX;
using namespace NOX::LineSearch;

Polynomial::Polynomial(Parameter::List& params) 
{
  reset(params);
}

Polynomial::~Polynomial()
{

}

bool Polynomial::reset(Parameter::List& params)
{ 
  minstep = params.getParameter("Minimum Step", 1.0e-12);
  defaultstep = params.getParameter("Default Step", 1.0);
  recoverystep = params.getParameter("Recovery Step", defaultstep);
  maxiters = params.getParameter("Max Iters", 100);
  return true;
}

bool Polynomial::compute(Abstract::Group& newgrp, double& step, 
			 const Abstract::Vector& dir,
			 const Solver::Generic& s) 
{
  const Abstract::Group& oldgrp = s.getPreviousSolutionGroup();
  
  double oldf = 0.5*oldgrp.getNormF()*oldgrp.getNormF();  
                            // Redefined f(), RH

  // General computation of directional derivative used in curvature condition
  // Note that for Newton direction, oldfprime = -2.0*oldf
  Abstract::Vector* tmpvecptr = oldgrp.getX().clone(ShapeCopy);
  oldgrp.applyJacobian(dir,*tmpvecptr);
  double oldfprime = tmpvecptr->dot(oldgrp.getF());
  delete tmpvecptr;

  double newf, prevf;
  double tempStep, previousStep;
  double a,b,term1,term2,disc ;
  bool isfailed = false;
  bool firstPass = true ;

  step = defaultstep;
  newgrp.computeX(oldgrp, dir, step);
  newgrp.computeF();    
  newf = 0.5*newgrp.getNormF()*newgrp.getNormF();  
                            // Redefined f(), RH

  int niters = 1;

  if (Utils::doPrint(Utils::InnerIteration)) {
   cout << "\n" << Utils::fill(72) << "\n" << "-- Polynomial Line Search -- \n";
  }

  while (newf >= oldf+0.0001*step*oldfprime) {  //Armijo-Goldstein condition, RH

    if (Utils::doPrint(Utils::InnerIteration)) {
      cout << setw(3) << niters << ":";
      cout << " step = " << Utils::sci(step);
      cout << " oldf = " << Utils::sci(sqrt(2.*oldf));
      cout << " newf = " << Utils::sci(sqrt(2.*newf));
      cout << endl;
    }

    if( firstPass == true)
    {

     /*   First try quadratic  */

      tempStep = -oldfprime/(2.0*(newf - oldf - oldfprime)) ;
      firstPass = false;
     }

     else {

     /*   Do cubic as many times as needed */

       term1 = newf - oldf - step*oldfprime ;
       term2 = prevf - oldf - previousStep*oldfprime ;
       a = 1.0/(step-previousStep)*( term1/step/step -
                                    term2/previousStep/previousStep) ;
       b = 1.0/(step-previousStep)*( -term1*previousStep/step/step +
                                    term2*step/previousStep/previousStep) ;
       disc = b*b - 3.0*a*oldfprime ;
       if(disc < 0) {
         step = recoverystep;
         newgrp.computeX(oldgrp, dir, step);
         newgrp.computeF();    
         newf = 0.5*newgrp.getNormF()*newgrp.getNormF();  
         cout << Utils::fill(5,' ') << "step = " << Utils::sci(step);
         cout << Utils::fill(1,' ') << "oldf = " << Utils::sci(sqrt(2.*oldf));
         cout << Utils::fill(1,' ') << "newf = " << Utils::sci(sqrt(2.*newf));
         cout << " (USING RECOVERY STEP!)" << endl;
         cout << Utils::fill(72) << "\n" << endl;
         isfailed = true;
         return(!isfailed);
       }
       if( fabs(a) < 1.e-12)
          tempStep = -oldfprime/2.0/b ;
       else
          tempStep = (-b + sqrt(disc))/3.0/a ;

       if(tempStep > 0.5*step) tempStep = 0.5*step ;

     }
    
     previousStep = step ;
     prevf = newf ; 

     /*   Respect bounds, which are hard-coded here  */ 

    if(tempStep < 0.10*step) step *= 0.10;
    else step = tempStep ;


    if ((step < minstep) || (niters > maxiters))
    {
      step = recoverystep;
      newgrp.computeX(oldgrp, dir, step);
      newgrp.computeF();    
      newf = 0.5*newgrp.getNormF()*newgrp.getNormF(); 
      cout << Utils::fill(5,' ') << "step = " << Utils::sci(step);
      cout << Utils::fill(1,' ') << "oldf = " << Utils::sci(sqrt(2.*oldf));
      cout << Utils::fill(1,' ') << "newf = " << Utils::sci(sqrt(2.*newf));
      cout << " (USING RECOVERY STEP!)" << endl;
      cout << Utils::fill(72) << "\n" << endl;
      isfailed = true;
      return(!isfailed);
    }
    
    niters ++;
    newgrp.computeX(oldgrp, dir, step);
    newgrp.computeF();    
    newf = 0.5*newgrp.getNormF()*newgrp.getNormF();  
                            // Redefined f(), RH
  } 

  
  if (Utils::doPrint(Utils::InnerIteration)) {
      cout << setw(3) << niters << ":";
      cout << " step = " << Utils::sci(step);
      cout << " oldf = " << Utils::sci(sqrt(2.*oldf));
      cout << " newf = " << Utils::sci(sqrt(2.*newf));
      cout << " (STEP ACCEPTED!)" << endl;
      cout << Utils::fill(72) << "\n" << endl;

  }

  return (!isfailed);
}

