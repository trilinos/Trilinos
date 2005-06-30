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
File:      reed_integrator.cpp
Purpose:   Integrate the Clarinet Reed Problem
Date:      6-30-05
Author:    Joseph Simonis
**************************************************************/


/**** Include Files ****/
#include "reed_integrator.hpp"

//-----------------------------------------------------------------
// Function      : ReedIntegrator::ReedIntegrator
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/30/05
//------------------------------------------------------------------
ReedIntegrator::ReedIntegrator()
{
  num_elements=2;
}
//-----------------------------------------------------------------
// Function      : ReedIntegrator::ReedIntegrator
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/30/05
//------------------------------------------------------------------
ReedIntegrator::~ReedIntegrator()
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
bool ReedIntegrator::Integrate(double *ynew, const double *y, const double t, 
			       const double param)
{
  int n = num_elements;

  // Constants:
  double a=5.0;
  double b=4.0;
  double k=param;

  double* temp;
  temp = new double[n];

  double dt;
  dt=t/1000;

  // Set ynew=y;
  for (int i=0;i<n;i++) ynew[i]=y[i];



  // integrate over time
  for (int i=0; i<1000;i++)
    {
      temp[0]=dt*ynew[1];
      temp[1]=dt*(a*ynew[1]-b*(ynew[1]*ynew[1]*ynew[1])-k*ynew[0]);
      ynew[0]+=temp[0];
      ynew[1]+=temp[1];
    }
}
