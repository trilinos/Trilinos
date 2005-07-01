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
File:      rpmtest_integrator.cpp
Purpose:   A. Salinger's orginal test without psuedo-arclength
           continuation
Date:      7-01-05
Author:    Joseph Simonis
**************************************************************/


/**** Include Files ****/
#include "rpmtest_integrator.hpp"

//-----------------------------------------------------------------
// Function      : RpmtestIntegrator::RpmtestIntegrator
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 07/01/05
//------------------------------------------------------------------
RpmtestIntegrator::RpmtestIntegrator()
{
  num_elements=202;
}
//-----------------------------------------------------------------
// Function      : RpmtestIntegrator::RpmtestIntegrator
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 07/01/05
//------------------------------------------------------------------
RpmtestIntegrator::~RpmtestIntegrator()
{

}
//-----------------------------------------------------------------
// Function      : ReedIntegrator::Integrate
// Purpose       : The integration function
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/30/05
//------------------------------------------------------------------
bool RpmtestIntegrator::Integrate(double *ynew, const double *y, const double t, 
			       const double param)
{
  int n = num_elements;
  int m = n/2;
  double deltat = 0.002;
  double invdx_sq = (double) ((m-1)*(m-1)) / 400.0;
  
  // v is temporary for y at each explicit time step j
  double* v;
  v = new double[n];

  for (int j=0;j<n;j++)
    v[j]=y[j];


  // multiple explicit time steps per step
  for (int j=0; j<2000; j++) {

    for (int i=1; i<m-1; i++) 
      {
	
	ynew[i]   = v[i] + deltat * (
				     invdx_sq * (v[i-1] - 2*v[i] + v[i+1]) 
				     + v[i]* (1.0 - v[i]*v[i]) - v[i+m]  ); 
	ynew[m+i] = v[i+m] + deltat * (
				       4.0 * invdx_sq * (v[i-1+m] - 2*v[i+m] + v[i+1+m]) 
				       + param * (v[i] - 2.0 * v[i+m] +0.03)  );
      }
    
    ynew[0] = ynew[1];
    ynew[m-1] = ynew[m-2];
    ynew[m] = ynew[m+1];
    ynew[n-1] = ynew[n-2];
    
    // Increment to new explicit step
   for (int k=0;k<n;k++)
    v[k]=ynew[k];

   return true;
  }
};
