//
// @HEADER
// ***********************************************************************
// 
//                           Capo Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

/************************************************************ 
File:      vpl_integrator.cpp
Purpose:   VanderPol Integrator
Date:      6-30-05
Author:    Joseph Simonis

IMPORTANT NOTE ON LICENCE

This integrator has been modified from main.cpp for the VanderPol 
ODE found in the integration 
program written by Blake Ashby at Stanford.  His files included here: 
decsol.cpp, decsol.h, IntegratorT.cpp, IntegratorT.h, 
StiffIntegratorT.cpp and  StiffIntegratorT.h.  The files were 
downloaded from http://www.unige.ch/math/folks/hairer/software.html
It is his C++ version of RADAU5 details of which are found in
Solving Ordinary Differential Equations. Stiff and 
Differential-Algebraic Problems. 2nd edition. Springer Series in 
Comput. Math., vol. 14.  The following Licence must accompany
this code:

Copyright (c) 2004, Ernst Hairer

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are 
met:

- Redistributions of source code must retain the above copyright 
notice, this list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright 
notice, this list of conditions and the following disclaimer in the 
documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ?AS 
IS? AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR 
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**************************************************************/


/**** Include Files ****/
#include "vpl_integrator.hpp"
#include "StiffIntegratorT.h"


/**** Necessary Function Declerations ****/
void Function(double x, double *y, double *f, double eps);
void Jacobian(double x, double *y, double **J, double eps);
void Mass(double **M);

//-----------------------------------------------------------------
// Function      : VplIntegrator::VplIntegrator
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/30/05
//------------------------------------------------------------------
VplIntegrator::VplIntegrator()
{
  num_elements=2;
}
//-----------------------------------------------------------------
// Function      : VplIntegrator::VplIntegrator
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/30/05
//------------------------------------------------------------------
VplIntegrator::~VplIntegrator()
{

}
//-----------------------------------------------------------------
// Function      : VplIntegrator::Integrate
// Purpose       : The integration function
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/30/05
//------------------------------------------------------------------
bool VplIntegrator::Integrate(double *ynew, const double *y, const double t, 
			       const double param)
{
  int n = num_elements;
  double* z;
  z = new double[n];
  // initial values for y
  for (int i=0;i<2;i++)
    z[i]=y[i];
  // initial value for T
  double xbeg = 0.0;
  // final value for T
  const double xend = t;
  // interval of x for printing output
  double dx = 0.001;
  // The continuation parameter
  double eps = param;
  // rtoler and atoler are scalars
  int itoler = 0;
  // relative tolerance
  double *rtoler = new double(1.0e-10);
  // absolute tolerance
  double *atoler = new double(1.0e-10);
  // use SolutionOutput routine
  const int iout = 0;
  // initial step size
  double hinit = 0.0;
  // analytical Jacobian function provided
  const int ijac = 1;
  // number of non-zero rows below main diagonal of Jacobian
  int mljac = n;
  // number of non-zero rows above main diagonal of Jacobian
  int mujac = n;
  // Mass matrix routine is identity
  const int imas = 0;
  int mlmas = 0;
  int mumas = 0;
	
  // Use default values (see header files) for these parameters:
  double hmax(0.0);
  int nmax(0);
  double uround(0.0), safe(0.0), facl(0.0), facr(0.0);
  int nit(0);
  bool startn(false);
  int nind1(0), nind2(0), nind3(0), npred(0), m1(0), m2(0);
  bool hess(false);
  double fnewt(0.0), quot1(0.0), quot2(0.0), thet(0.0);

  StiffIntegratorT stiffT(n, z, xbeg, xend, dx, itoler, rtoler, atoler,
			  iout, hinit, hmax, nmax, uround, safe, facl, facr, ijac, mljac,
			  mujac, imas, mlmas, mumas, nit, startn, nind1, nind2, nind3, npred,
			  m1, m2, hess, fnewt, quot1, quot2, thet, eps);
  
  //cout << "\n\n*******Problem integrated with RADAU5*******\n\n";
  
  stiffT.Integrate(eps);
  
  // print statistics
  /*
  cout << "fcn = " << stiffT.NumFunction() <<
    " jac = " << stiffT.NumJacobian() <<
    " step = " << stiffT.NumStep() <<
    " accpt = " << stiffT.NumAccept() <<
    " rejct = " << stiffT.NumReject() <<
    " dec = " << stiffT.NumDecomp() <<
    " sol = " << stiffT.NumSol() << endl;
  */
  for (int i=0;i<2;i++)
    ynew[i]=z[i];

  delete [] z;
  delete rtoler;
  delete atoler;
  
  

  return true;

}

void Function(double x, double *y, double *f, double eps)
{
	double prod;

	f[0] = y[1];
	prod = 1.0 - y[0]*y[0];
	f[1] = (prod*y[1] - y[0])/eps;

	return;

} // fEval

void Jacobian(double x, double *y, double **J, double eps)
{

	J[0][0] = 0.0;
	J[0][1] = 1.0;
	J[1][0] = (-2.0*y[0]*y[1] - 1.0)/eps;
	J[1][1] = (1.0 - y[0]*y[0])/eps;

	return;

} // Jacobian

void Mass(double **M)
{

} // Mass
