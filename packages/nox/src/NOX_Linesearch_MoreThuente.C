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

#include "NOX_Linesearch_MoreThuente.H"

#include "NOX_Utils.H"		// for static doPrint function
#include <iomanip>		// for setw
#include <math.h>		// for abs, sqrt
#include <stdio.h>		// for getchar

using namespace NOX;
using namespace NOX::Linesearch;

MoreThuente::MoreThuente(const Parameter::List& params) 
{
  reset(params);
}

MoreThuente::~MoreThuente()
{

}

void MoreThuente::reset(const Parameter::List& params)
{ 
  minstep = params.getParameter("Minimum Step", 1.0e-12);
  defaultstep = params.getParameter("Default Step", 1.0);
  recoverystep = params.getParameter("Recovery Step", defaultstep);
}

bool MoreThuente::operator()(Abstract::Group& newgrp, double& step, 
			 const Abstract::Group& oldgrp, const Abstract::Vector& dir) 
{

/*  First pass -->  copy alogrithm from Xyce/Xyce/src/NonlinearSolverPKG/src
                                           /N_NLS...........
    10-29-2001   RH
*/
   double xtrapf = 4.0 ;
   double stepmin=1.e-4; 
   double stepmax=1.e6 ;
   double ftol = 0.0001 ;
   double gtol = 0.500 ;
   double rtol = 1.e-10 ;
   int maxfev = 30 ;

   bool isfailed = false;
   step = defaultstep;

   int info ;
   double factor ;
   double dginit, dgtest ;
   double stx, sty ;
   double finit, fx, fy ;
   double fm, fxm, fym ;
   double dgm, dgxm, dgym ;
   double dg, dgx, dgy ;
   double width, width1 ;
   double value ;
   double gnorm, ftest1 ;
   int bracket ;
   int stage ; /*  1 - use auxillary function
                   0 - use objective function  */
   int i, nfev ;

   const Abstract::Vector & initialGradient = oldgrp.getGrad();
   dginit = dir.dot(initialGradient);

/* Check sign of dginit:  to do, RH  */

/* ...........  Initialization  ..........*/

   bracket = 0 ;
   stage = 1 ;
   finit = 0.5*oldgrp.getNormRHS()*oldgrp.getNormRHS();
                              /* initial value of objective function */
   dgtest = ftol*dginit ;  /* sufficient decrease condition  */
   width = stepmax - stepmin ;  /* width of search interval  */
   width1 = 2.0*width ;

   stx = 0.0 ;
   fx = finit ;
   dgx = dginit ;
   sty = 0.0 ;
   fy = finit ;
   dgy = dginit ;

//   printf("\n\nInitialization: finit, dginit, dgtest: %e, %e, %e\n\n",
//                               finit, dginit, dgtest);
//   getchar();

   nfev = 0 ; /* initialize # function evaluations  */
    
  if (Utils::doPrint(1)) {
   cout << "\n" << Utils::fill(72) << "\n" << " -- More'-Thuente Line Search -- \n";
  }

/*  Begin loop up to max # function evaluations  */

   for(i=0; i<maxfev; i++)
   {
      if(bracket)
      {
        stepmin = mymin(stx, sty);
        stepmax = mymax(stx, sty);
      }
      else
      {
         stepmin = stx ;
         stepmax = step + xtrapf*(step-stx) ;
      }

/* Constrain step to lie within bounds  */
      step = mymax(step, stepmin);
      step = mymin(step, stepmax);

/* Return best step so far if termination occurs  */
      if( (bracket&&(step <= stepmin || step >= stepmax)) ||
          (bracket&&(stepmax-stepmin <= rtol*stepmax))   ||
          (nfev >= maxfev-1) ) step = stx ;

 
/*  Find new Residual and objective function value  */

      newgrp.computeX(oldgrp, dir, step);
      newgrp.computeRHS();
      value = 0.5*newgrp.getNormRHS()*newgrp.getNormRHS();

      nfev++ ;
 
/*  Find gradient at new location  */

//      newgrp.computeRHS();
      newgrp.computeJacobian();
      newgrp.computeGrad();
      const Abstract::Vector & gradient = newgrp.getGrad();
/*
      factor = -1.0/value ;
      dscal_(&n, &factor, gradient, &iOne) ;
*/

//      gnorm = mydot(n, gradient, gradient) ;

      dg = dir.dot(gradient);
      ftest1 = finit + step*dgtest ;

//      printf("[i]\tstep\t\tdphi\t\tphi\t\tstepmin\t\tstepmax\n");
//      printf("[%d]\t%e\t%e\t%e\t%e\t%e",
//            i,step,dg,value,stepmin,stepmax);
//      getchar();
 
/*  Implement convergence checks, to do RH  */

      if((value <= ftest1) && (fabs(dg) <= gtol*-dginit))
      {
         if (Utils::doPrint(1)) {
            cout << Utils::fill(5,' ') << "step = " << Utils::sci(step);
            cout << Utils::fill(1,' ') << "oldf = " << Utils::sci(finit);
            cout << Utils::fill(1,' ') << "newf = " << Utils::sci(value);
            cout << " (STEP ACCEPTED!)" << endl;
            cout << Utils::fill(72) << "\n" << endl;
         }
         break;  /*  CONVERGENCE !!!  */
      }

     if (Utils::doPrint(1)) {
       cout << Utils::fill(5,' ') << "step = " << Utils::sci(step);
       cout << Utils::fill(1,' ') << "oldf = " << Utils::sci(finit);
       cout << Utils::fill(1,' ') << "newf = " << Utils::sci(value);
       cout << endl;
     }


/* Determine which objective function to use  */
      if( stage&&(value <= ftest1)&&(dg >= dginit*gtol)) stage=0;

     if(stage)
     {
//        printf("\nUsing auxillary function\n\n");
        fm = value - finit - step*dgtest ; /* Use auxillary function  */
        fxm = fx - finit - stx*dgtest ;    /* and derivatives  */
        fym = fy - finit - sty*dgtest ;
        dgm = dg - dgtest ;
        dgxm = dgx - dgtest ;
        dgym = dgy - dgtest ;

/*  Update interval and compute new step length using auxillary function */
        info = MTStep(&stx, &fxm, &dgxm, &sty, &fym, &dgym, &step, &fm, &dgm,
                      &stepmin, &stepmax, &bracket);

        fx = fxm + finit + stx*dgtest ;
        fy = fym + finit + sty*dgtest ;
        dgx = dgxm + dgtest ;
        dgy = dgym + dgtest ;
     }
     else
     {
//        printf("\nUsing original function\n\n");
/*  Update interval and compute new step length using original function */
        info = MTStep(&stx, &fx, &dgx, &sty, &fy, &dgy, &step, &value, &dg,
                      &stepmin, &stepmax, &bracket);
     }

//      printf("\nInfo : %d\tStep : %f\n\n",info,step);


/*  Put in checks based on returned info, to do RH  */

      if(bracket)
      {
         if(fabs(sty-stx) >= 0.66*width1) step = stx + 0.5*(sty-stx) ;
         width1 = width ;
         width = fabs(sty-stx) ;
      }
   }

   return(!isfailed);
}
/* --------------  End of lineSearchMT -----------------*/



int MoreThuente::MTStep(double *stx, double *fx, double *dx, 
           double *sty, double *fy,
           double *dy, double *stp, double *fp, double *dp, double *stepmin,
           double *stepmax, int *bracket)
{
/*  First pass -->  copy alogrithm from Xyce/Xyce/src/NonlinearSolverPKG/src
                                           /N_NLS...........
    10-29-2001   RH
*/

   double gamma1, p, q, r, s, sgnd, stpc, stpf, stpq, theta ;
   int bound, info ;

/*  Check input parameters for errors, to do RH  */

/*  Check if derivatives have opposite sign  */
   sgnd = *dp * (*dx/fabs(*dx)) ;

/*  For all of the following cases, the quadratic and cubic fits need to be
    checked to make sure they are correct, to do RH  */

/*
   printf("\n\nInside MTStep,  bracket %d\n\nstx, fx, dx :\t%e\t%e\t%e\n",
             *bracket,*stx,*fx,*dx);
   printf("stp, fp, dp :\t%e\t%e\t%e\n",*stp,*fp,*dp);
   printf("sty, fy, dy :\t%e\t%e\t%e\n",*sty,*fy,*dy);
   printf("stepmin, stepmax:\t%e\t%e\t%e\n",*stepmin,*stepmax);
   getchar();
*/
/*  Case 1:  */
   if(*fp > *fx)
   {
      info = 1 ;
      bound = 1 ;
      theta = 3.0*(*fx-*fp)/(*stp-*stx) + *dx + *dp ;
      s = mymax(fabs(theta),fabs(*dx)) ;
      s = mymax(s,fabs(*dp)) ;
      gamma1 = s * sqrt(pow(theta/s,2.0) - (*dx/s)*(*dp/s)) ;

      if(*stp < *stx) gamma1 = -gamma1 ;

      p = (gamma1 - *dx) + theta ;
      q = ((gamma1 - *dx) + gamma1) + *dp ;
      r = p / q ;
      stpc = *stx + r*(*stp-*stx) ;
      stpq = *stx + ((*dx/((*fx-*fp)/(*stp-*stx) + *dx)) * 0.5) * (*stp-*stx);

      if(fabs(stpc-*stx) < fabs(stpq-*stx)) stpf = stpc ;
      else stpf = stpc + 0.5*(stpq-stpc) ;

      *bracket = 1 ;
     
/*
//      printf("Case 1:\nCubic :\t%e\nQuad :\t%e\n\n",stpc,stpq);
//      getchar();
*/
   }

/*  Case 2:  */
   else if(sgnd < 0.0)
   {
      info = 2 ;
      bound = 0 ;
      theta = 3.0*(*fx-*fp)/(*stp-*stx) + *dx + *dp ;
      s = mymax(fabs(theta),fabs(*dx)) ;
      s = mymax(s,fabs(*dp)) ;
      gamma1 = s * sqrt(pow(theta/s,2.0) - (*dx/s)*(*dp/s)) ;

      if(*stp > *stx) gamma1 = -gamma1 ;

      p = (gamma1 - *dp) + theta ;
      q = ((gamma1 - *dp) + gamma1) + *dx ;
      r = p / q ;
      stpc = *stp + r*(*stx-*stp) ;
      stpq = *stp + (*dp/(*dp-*dx)) * (*stx-*stp);

      if(fabs(stpc-*stp) > fabs(stpq-*stp)) stpf = stpc ;
      else stpf = stpq ;

      *bracket = 1 ;
     
/*
//      printf("Case 2:\nCubic :\t%e\nQuad :\t%e\n\n",stpc,stpq);
//      getchar();
*/
   }

/*  Case 3:  */
   else if(fabs(*dp) < fabs(*dx))
   {
      info = 3 ;
      bound = 1 ;
      theta = 3.0*(*fx-*fp)/(*stp-*stx) + *dx + *dp ;
      s = mymax(fabs(theta),fabs(*dx)) ;
      s = mymax(s,fabs(*dp)) ;
/*  Comments....  */
      gamma1 = s * sqrt(mymax(0.0,pow(theta/s,2.0) - (*dx/s)*(*dp/s))) ;

      if(*stp > *stx) gamma1 = -gamma1 ;

      p = (gamma1 - *dp) + theta ;
      q = (gamma1 + (*dx-*dp)) + gamma1 ;
      r = p / q ;

      if(r < 0.0 && gamma1 != 0.0) stpc = *stp + r*(*stx-*stp) ;
      else if(*stp > *stx)         stpc = *stepmax ;
      else                         stpc = *stepmin ;

      stpq = *stp + (*dp/(*dp-*dx)) * (*stx-*stp);

     
/*
      printf("Case 3:\nCubic :\t%e\nQuad :\t%e\n\n",stpc,stpq);
*/

      if(*bracket)
      {
/*
         printf("fabs(stpc-*stp): %e ,\tfabs(stpq-*stp): %e\n", 
                         fabs(stpc-*stp),fabs(stpq-*stp));
*/
         if(fabs(stpc-*stp) < fabs(stpq-*stp)) stpf = stpc ;
         else stpf = stpq ;
      }
      else
      {
/*
         printf("fabs(stpc-*stp): %e ,\tfabs(stpq-*stp): %e\n", 
                         fabs(stpc-*stp),fabs(stpq-*stp));
*/
         if(fabs(stpc-*stp) > fabs(stpq-*stp)) stpf = stpc ;
         else stpf = stpq ;
      }
/*
      printf("stpf :\t%e\n\n",stpf);
      getchar();
*/
   }

/*  Case 4:  */
   else 
   {
      info = 4 ;
      bound = 0 ;

/*
      printf("Case 4:  bracket %d\n\n\n",*bracket);
*/
      if(*bracket)
      {
         theta = 3.0*(*fp-*fy)/(*sty-*stp) + *dy + *dp ;
         s = mymax(fabs(theta),fabs(*dy)) ;
         s = mymax(s,fabs(*dp)) ;
         gamma1 = s * sqrt(pow(theta/s,2.0) - (*dy/s)*(*dp/s)) ;

         if(*stp > *sty) gamma1 = -gamma1 ;

         p = (gamma1 - *dp) + theta ;
         q = ((gamma1 - *dp) + gamma1) + *dy ;
         r = p / q ;
         stpc = *stp + r*(*sty-*stp) ;
         stpf = stpc ;
     
/*
         printf("Case 4:  bracket %d\nCubic :\t%e\n\n",*bracket,stpc);
         getchar();
*/
      }

      else if(*stp > *stx)
      {
         stpf = *stepmax ;
      }
      else
      {
         stpf = *stepmin ;
      }
   }

/*  Update interval using ad hoc error trapping  */
   if(*fp==*fp)
   {
      if(*fp > *fx)
      {
         *sty = *stp ;
         *fy = *fp ;
         *dy = *dp ;
      }
      else
      {
         if(sgnd < 0.0)
         {
            *sty = *stx ;
            *fy = *fx ;
            *dy = *dx ;
         }
         *stx = *stp ;
         *fx = *fp ;
         *dx = *dp ;
      }
   }
   else
   {  /* *fp is NaN  */
      *sty = *stp ;
      *fy = *fp ;
      *dy = *dp ;
   }

/*  Compute new step and safeguard it  */
   if(*fp==*fp)
      *stp = 0.02231964 ;
   else
      stpf *= 1.0e-3;  /*  *fp is NaN  */

   stpf = mymin(*stepmax,stpf) ;
   stpf = mymax(*stepmin,stpf) ;
   *stp = stpf ;

   if(*bracket && bound)
   {
      if(*sty > *stx) *stp = mymin(*stx+0.66*(*sty-*stx),*stp) ;
      else            *stp = mymax(*stx+0.66*(*sty-*stx),*stp) ;
   }
         
/*
   printf("Final vals: stp, stx, sty\t%e, %e, %e\n",*stp,*stx,*sty);
   getchar();
*/
   return(info) ;
}
/* --------------  End of MTStep -----------------*/






double MoreThuente::mymin(double x, double y)
{
   if(x < y) return x ;
   return y ;
}
/* --------------  End of mymin -----------------*/




double MoreThuente::mymax(double x, double y)
{
   if(x < y) return y ;
   return x ;
}
/* --------------  End of mymax -----------------*/



