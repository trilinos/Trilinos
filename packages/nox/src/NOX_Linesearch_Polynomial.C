// $Id$ 
// $Source$ 

// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#include "NOX_Linesearch_Polynomial.H"

#include "NOX_Utils.H"		// for static doPrint function
#include <iomanip>		// for setw
#include <math.h>		// for abs, sqrt
#include <stdio.h>		// for getchar

using namespace NOX;
using namespace NOX::Linesearch;

Polynomial::Polynomial(const Parameter::List& params) 
{
  reset(params);
}

Polynomial::~Polynomial()
{

}

void Polynomial::reset(const Parameter::List& params)
{ 
  minstep = params.getParameter("Minimum Step", 1.0e-12);
  defaultstep = params.getParameter("Default Step", 1.0);
  recoverystep = params.getParameter("Recovery Step", defaultstep);
}

bool Polynomial::operator()(Abstract::Group& newgrp, double& step, 
			 const Abstract::Group& oldgrp, const Abstract::Vector& dir) 
{

  double oldf = 0.5*oldgrp.getNormRHS()*oldgrp.getNormRHS();  
                            // Redefined f(), RH
  double oldfprime = -2.0*oldf;  //  This holds for Newton direction, RH
//  Check that this indeed holds.....DONE !!  12-12-2001  RH
//  const NLS_Vector& chkvec = oldgrp.getGrad();
//  double chk = dir.dot(chkvec);
//  cout << "\n Computed Grad ....."<<oldfprime<<", "<<chk;
//  getchar();
  double newf, prevf;
  double tempStep, previousStep;
  double a,b,term1,term2,disc ;
  bool isfailed = false;
  bool firstPass = true ;

  step = defaultstep;
  newgrp.computeX(oldgrp, dir, step);
  newgrp.computeRHS();    
  newf = 0.5*newgrp.getNormRHS()*newgrp.getNormRHS();  
                            // Redefined f(), RH

  if (Utils::doPrint(1)) {
   cout << "\n" << Utils::fill(72) << "\n" << " -- Polynomial Line Search -- \n";
  }

  while (newf >= oldf+0.0001*step*oldfprime) {  //Armijo-Goldstein condition, RH

    if (Utils::doPrint(1)) {
      cout << Utils::fill(5,' ') << "step = " << Utils::sci(step);
      cout << Utils::fill(1,' ') << "oldf = " << Utils::sci(oldf);
      cout << Utils::fill(1,' ') << "newf = " << Utils::sci(newf);
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


    if (step < minstep)
    {
      step = recoverystep;
      cout << Utils::fill(5,' ') << "step = " << Utils::sci(step);
      cout << Utils::fill(1,' ') << "oldf = " << Utils::sci(oldf);
      cout << Utils::fill(1,' ') << "newf = " << Utils::sci(newf);
      cout << " (USING MINSTEP!)" << endl;
      cout << Utils::fill(72) << "\n" << endl;
      isfailed = true;
      return(!isfailed);
    }

    newgrp.computeX(oldgrp, dir, step);
    newgrp.computeRHS();    
    newf = 0.5*newgrp.getNormRHS()*newgrp.getNormRHS();  
                            // Redefined f(), RH
  } 

  
  if (Utils::doPrint(1)) {
      cout << Utils::fill(5,' ') << "step = " << Utils::sci(step);
      cout << Utils::fill(1,' ') << "oldf = " << Utils::sci(oldf);
      cout << Utils::fill(1,' ') << "newf = " << Utils::sci(newf);
      cout << " (STEP ACCEPTED!)" << endl;
      cout << Utils::fill(72) << "\n" << endl;

  }

  return (!isfailed);
}

