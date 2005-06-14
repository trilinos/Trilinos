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
// Foundation, In(Bc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef Capo_EXAMPLE_INTEGRATOR_H
#define Capo_EXAMPLE_INTEGRATOR_H

/************************************************************ 
File:      ExampleIntegrator.hpp
Purpose:   Simple example integrator using epetra vectors.
Date:      6-09-05
Author:    Joseph Simonis
**************************************************************/


/**** Include Files ****/
#include "Epetra_Vector.h"
#include "Teuchos_RefCountPtr.hpp"


/**** ExampleIntegrator Class ****/
class ExampleIntegrator
{
public:
  // Constructor
  ExampleIntegrator(int numelements);
  
  // Destructor
  ~ExampleIntegrator();
  
  // Integration function
  bool Integrate(Epetra_Vector *y, const Epetra_Vector &x,double t, double lambda);
  
private:
  int numElements_;
  
};

#endif //Capo_EXAMPLE_INTEGRATOR_H
