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
File:      reed_interface.cpp
Purpose:   Create a class to wrap the clarinet reed integrator.
Date:      6-30-05
Author:    Joseph Simonis
**************************************************************/

/**** Include Files ****/
#include "reed_interface.hpp"
#include "reed_integrator.hpp"
#include "Epetra_Map.h"
#include "Thyra_VectorBase.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"


//-----------------------------------------------------------------
// Function      : ThyraIntegrator::ThyraIntegrator
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/30/05
//------------------------------------------------------------------
ThyraIntegrator::ThyraIntegrator(Teuchos::RefCountPtr<ReedIntegrator> problem)
{
  AppIntegrator = problem;
}

//-----------------------------------------------------------------
// Function      : ThyraIntegrator::Integrate
// Purpose       : Integrator which uses Thyra vectors for input 
//                 output.
// Special Notes : This function satisfies the virtual function
//                 declared in Capo_Integrate.hpp.
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/30/05
//------------------------------------------------------------------
bool ThyraIntegrator::Integrate(const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& y,const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& x,const double T, const double lambda)
{
  double *u;
  double *u_T;

  u = new double[2];
  u_T = new double[2];

  for (int i=0;i<2;i++)
    {
      u[i]=Thyra::get_ele(*y,i+1);
    }

  AppIntegrator->Integrate(u_T, u, T, lambda);

  for (int i=0;i<2;i++)
    {
      Thyra::set_ele(i+1,u_T[i],&*x);
    }
  delete [] u;
  delete [] u_T;

  return true;

}

