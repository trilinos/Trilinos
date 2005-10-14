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
File:      ExampleIntegrator.cpp
Purpose:   Simple example integrator using epetra vectors.
Date:      6-09-05
Author:    Joseph Simonis
**************************************************************/

/**** Include Files ****/
#include "Epetra_Map.h"
#include "Epetra_SerialComm.h"
#include "Epetra_Vector.h"
#include "ExampleIntegrator.hpp"

//-----------------------------------------------------------------
// Function      : ExampleIntegrator::ExampleIntegrator
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/09/05
//------------------------------------------------------------------
ExampleIntegrator::ExampleIntegrator(int numelements)
{  
  numElements_ = numelements;

}

//-----------------------------------------------------------------
// Function      : ExampleIntegrator::~ExampleIntegrator
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/09/05
//------------------------------------------------------------------
ExampleIntegrator::~ExampleIntegrator()
{
}

//-----------------------------------------------------------------
// Function      : ExampleIntegrator::Integrate
// Purpose       : The Integration function for the problem
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/09/05
//------------------------------------------------------------------

bool ExampleIntegrator::Integrate(Epetra_Vector *y, const Epetra_Vector &x, \
				   double t, double lambda)
{
  y->Scale(t*lambda,x);
  return true;
}
