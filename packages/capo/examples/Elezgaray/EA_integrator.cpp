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
File:      EA_integrator.cpp
Purpose:   The Elezgaray & Arneodo model problem.
Date:      8-01-05
Author:    Joseph Simonis
**************************************************************/


/**** Include Files ****/
#include "EA_integrator.hpp"

//-----------------------------------------------------------------
// Function      : EAIntegrator::EAIntegrator
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 08/01/05
//------------------------------------------------------------------
EAIntegrator::EAIntegrator(int n)
{
  numElements_=n;
}
//-----------------------------------------------------------------
// Function      : EAIntegrator::EAIntegrator
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 08/01/05
//------------------------------------------------------------------
EAIntegrator::~EAIntegrator()
{

}
//-----------------------------------------------------------------
// Function      : EAIntegrator::Integrate
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 08/01/05
//------------------------------------------------------------------
bool EAIntegrator::Integrate(double *u_T, const double *u,double t, double lambda)
{
  int n = numElements_;
  // Constants:
  double alpha=0.01;
  double epsilon=0.01;
  //  double C1=.008/(lambda*lambda);
  //  double C2=.004/(lambda*lambda);
  int k=(n+1)/2; 
  double h=1.0/(k+1);

  // Necessary storage
  double *a;
  double *b;
  double *c;
  double *r;
  double *temp;
    
  a=new double[k];
  b=new double[k];
  c=new double[k];
  r=new double[k];
  temp=new double[n];
 
  double dt;
  dt=t/10000;

  // Set u_T=u;
  for (int i=0;i<n;i++) 
    {
      u_T[i]=u[i];
    }

  // integrate over time
  for (int i=0;i<10000;i++)
    {
      for (int m=0;m<n;m++) temp[m]=u_T[m];

      // Use runga-kutta to integrate the reaction part.
      ea_rungakutta(dt,temp,n);

      for (int m=0;m<n;m++) u_T[m]=temp[m];

      // Now integrate the diffusion part:
      
      // First half....
      for (int m=0;m<k;m++)
	{
	  a[m]=(1+2*lambda*dt/(h*h));
	  r[m]=u_T[m];
	}

      for (int m=1;m<k;m++)
	{
	  b[m]=-(dt*lambda/(h*h));
	  c[m-1]=-(dt*lambda/(h*h));
	}
      b[0]=0;
      c[k-1]=0;
      r[0]+=(-2.0*dt*lambda/(h*h));
      r[k-1]+=(-2.0*dt*lambda/(h*h));

      thomas(b,a,c,r,k);

      for (int m=0;m<k;m++)
	{
	  u_T[m]=r[m];
	}

      // Second Half....
      for (int m=0;m<k;m++)
	{
	  a[m]=(1+2*lambda*dt/(h*h));
	  r[m]=u_T[m+k];
	}

      for (int m=1;m<k;m++)
	{
	  b[m]=-(dt*lambda/(h*h));
	  c[m-1]=-(dt*lambda/(h*h));
	}
      b[0]=0;
      c[k-1]=0;
      r[0]+=((-4.0)*dt*lambda/(h*h));
      r[k-1]+=((-4.0)*dt*lambda/(h*h));

      thomas(b,a,c,r,k);

 
      for (int m=0;m<k;m++)
	{
	  u_T[m+k]=r[m];
	}
      
    }

  delete [] a;
  delete [] b;
  delete [] c;
  delete [] r;
  delete [] temp;
  return true;
}


bool EAIntegrator::thomas(double* a, double* b, double* c, double* r, int n)
{
  c[0]=c[0]/b[0];
  r[0]=r[0]/b[0];
  for (int i=1;i<n;i++)
    {
      c[i]=c[i]/(b[i]-a[i]*c[i-1]);
      r[i]=(r[i]-a[i]*r[i-1])/(b[i]-a[i]*c[i-1]);
    }
  for (int i=n-2;i>-1;i--)
    {
      r[i]=r[i]-c[i]*r[i+1];
    }
  return true;
}


bool EAIntegrator::ea_rungakutta(double T,double* u, int n)
{

  double* kone;
  double* ktwo;
  double* kthree;
  double* kfour;

  kone=new double[n];
  ktwo=new double[n];
  kthree=new double[n];
  kfour=new double[n];

  double A=-2.0;
  double B=-4.0;

  int k=n/2;

  double h=T;

  // Calculate k1=h*f(u)
  for (int j=0; j<k; j++)
    {
      kone[j]=h*(100*(u[j+k]-(u[j]*u[j]+u[j]*u[j]*u[j])));
      kone[j+k]=h*(0.01-u[j]);
    }
  
  // Calculate k2=h*f(u+k1/2)
  for (int j=0; j<k; j++)
    {
      ktwo[j]=h*(100*( (u[j+k]+kone[j+k]/2.0)-((u[j]+kone[j]/2.0)*(u[j]+kone[j]/2.0)+(u[j]+kone[j]/2.0)*(u[j]+kone[j]/2.0)*(u[j]+kone[j]/2.0))));
      ktwo[j+k]=h*(0.01-u[j]-kone[j]/2.0);
    }
  
  // Calculate k3=h*f(u+k2/2)
  for (int j=0; j<k; j++)
    {
      kthree[j]=h*(100*( (u[j+k]+ktwo[j+k]/2.0)-((u[j]+ktwo[j]/2.0)*(u[j]+ktwo[j]/2.0)+(u[j]+ktwo[j]/2.0)*(u[j]+ktwo[j]/2.0)*(u[j]+ktwo[j]/2.0))));
      kthree[j+k]=h*(0.01-u[j]-ktwo[j]/2.0);
    }
  
  // Calculate k4=h*f(u+k3)
  for (int j=0; j<k; j++)
    {
      kfour[j]=h*(100*( (u[j+k]+kthree[j+k])-((u[j]+kthree[j])*(u[j]+kthree[j])+(u[j]+kthree[j])*(u[j]+kthree[j])*(u[j]+kthree[j]))));
      kfour[j+k]=h*(0.01-u[j]-kthree[j]);
    }
  
  for (int j=0;j<n;j++)
    {
      u[j]+=kone[j]/6+ktwo[j]/3+kthree[j]/3+kfour[j]/6;
    }
 
  delete [] kone;
  delete [] ktwo;
  delete [] kthree;
  delete [] kfour;
  
  return true;
} 



