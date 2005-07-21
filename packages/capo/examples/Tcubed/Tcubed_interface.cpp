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
File:      Tcubed_interface.cpp
Purpose:   Create a class to wrap the brusselator integroator.
Date:      7-15-05
Author:    Joseph Simonis
**************************************************************/

/**** Include Files ****/
#include "Tcubed_interface.hpp"
#include "Tcubed_FiniteElementProblem.H"
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
// Creation Date : 07/15/05
//------------------------------------------------------------------
ThyraIntegrator::ThyraIntegrator(Teuchos::RefCountPtr<Tcubed_FiniteElementProblem>& problem, Teuchos::RefCountPtr<const Epetra_Map>& emap)
{
  AppIntegrator = problem;
  Emap = emap;
}

//-----------------------------------------------------------------
// Function      : ThyraIntegrator::Integrate
// Purpose       : Integrator which uses Thyra vectors for input 
//                 output.
// Special Notes : This function satisfies the virtual function
//                 declared in Capo_Integrate.hpp.
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 07/15/05
//------------------------------------------------------------------
bool ThyraIntegrator::Integrate(const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& y,const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >& x,const double T, const double lambda)
{

  // Input y : Output x 
  Teuchos::RefCountPtr<const Epetra_Vector> epetra_y;
  Teuchos::RefCountPtr<Epetra_Vector> epetra_x;
  Teuchos::RefCountPtr<Epetra_Vector> epetra_temp;

  epetra_x = Thyra::get_Epetra_Vector(*Emap,x);
  epetra_y = Thyra::get_Epetra_Vector(*Emap,y);
  epetra_temp = Teuchos::rcp(new Epetra_Vector(*epetra_x));

  string temp_string;
  temp_string = "Right BC";
  AppIntegrator->setParameter(temp_string,lambda);

  epetra_x->Update(1.0,*epetra_y,0.0);
  epetra_temp->Update(1.0,*epetra_y,0.0);


  for (int i=0;i<100;i++)
    {
      AppIntegrator->evaluate(F_ONLY, &*epetra_x, &*epetra_temp, NULL);
      epetra_x->Update(1.0/1000.0,*epetra_temp,1.0);
    }


  /*
  cout << "Print the epetra vector after" << endl;
  for (int i=0;i<8;i++)
    cout <<(*epetra_x)[i] << " " << endl;
  */
  

  // In theory, when the objects epetra_y and epetra_x are destroyed,
  // the values in y and x are modified, therefore I don't actually
  // have to copy anything back...
  return true;

}

