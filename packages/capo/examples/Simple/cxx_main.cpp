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
File:      cxx_main.cpp
Purpose:   A simple example to call a Capo function on a 
           "epetra-integrator".
Date:      6-13-05
Author:    Joseph Simonis
**************************************************************/

/**** Include Files ****/
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Version.h"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "ExampleIntegrator.hpp"
#include "Example_Interface.hpp"
#include "Capo_Integrator.hpp"
#include "Capo_Parameter_List.hpp"
#include "Capo_Npgs.hpp"
#include "Capo_Stepper.hpp"
#include <string>

int main(int argc, char *argv[])
{

  typedef double Scalar; // Scalar type = double


#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  //Teuchos::RefCountPtr<const Epetra_Comm> epetra_comm = Teuchos::rcp( new Epetra_MpiComm(MPI_COMM_WORLD) );
  Teuchos::RefCountPtr<const Epetra_Comm> epetra_comm = Teuchos::rcp( new Epetra_MpiComm );
#else
  Teuchos::RefCountPtr<const Epetra_Comm> epetra_comm = Teuchos::rcp( new Epetra_SerialComm );
#endif

  int MyPID = epetra_comm->MyPID();

  int NumElements = 4;

  // Construct a Map with NumElements and index base of 0
  Teuchos::RefCountPtr<const Epetra_Map> epetra_map = Teuchos::rcp( new Epetra_Map(NumElements, 0, *epetra_comm) );

  // Construct a Thyra vector space
  Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<Scalar> > thyra_vs = Thyra::create_VectorSpace(epetra_map);

  // Create x and xn vectors
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > x  = Thyra::createMember(thyra_vs);
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > xn = Thyra::createMember(thyra_vs);

  if (MyPID == 0) 
  {
    cout << "Constructed a map with " << NumElements << " elements." << endl
        << "Constructed an Thyra vector space." << endl
        << "Created vectors x and xn." << endl;
  }



  Thyra::set_ele(1,1.0,&*xn);
  Thyra::set_ele(2,2.0,&*xn);
  Thyra::set_ele(3,1.0,&*xn);
  Thyra::set_ele(4,2.0,&*xn);

  Thyra::set_ele(1,0.0,&*x);
  Thyra::set_ele(2,0.0,&*x);
  Thyra::set_ele(3,0.0,&*x);
  Thyra::set_ele(4,0.0,&*x);

  if (MyPID == 0)
    {
      cout << "xn(0) = " << Thyra::get_ele(*xn,1) << endl;
      cout << "xn(1) = " << Thyra::get_ele(*xn,2) << endl;
      cout << "xn(2) = " << Thyra::get_ele(*xn,3) << endl;
      cout << "xn(3) = " << Thyra::get_ele(*xn,4) << endl;
    }

  double t0=3.2;
  double lambda0=1.2;

  // Create an Example Integrator Object...
  Teuchos::RefCountPtr<ExampleIntegrator> application;
  application = Teuchos::rcp(new ExampleIntegrator(4));
  cerr << "Created an Integrator" << endl;

  // Create two epetra vectors...
  Epetra_Vector xn_E = *(Thyra::get_Epetra_Vector(*&*epetra_map,xn));
  Epetra_Vector x_E = *(Thyra::get_Epetra_Vector(*&*epetra_map,x));
  cerr << "Created two Epetra Vectors out of the Thyra Vectors. " << endl;

  // Call the application integrator.
  (*application).Integrate(&x_E,xn_E,t0,lambda0);

  // Call the application integrator a different way...
  (*application).Integrate(&*(Thyra::get_Epetra_Vector(*&*epetra_map,x)),*(Thyra::get_Epetra_Vector(*&*epetra_map,xn)),t0,lambda0);

  if (MyPID == 0)
    {
      cout << "After Integration the vector returns:" << endl;
      cout << "x(0) = " << Thyra::get_ele(*x,1) << endl;
      cout << "x(1) = " << Thyra::get_ele(*x,2) << endl;
      cout << "x(2) = " << Thyra::get_ele(*x,3) << endl;
      cout << "x(3) = " << Thyra::get_ele(*x,4) << endl;
    }

  

  // Now Test the interface...

  // Create an ThyraIntegrate object...
  Teuchos::RefCountPtr<ThyraIntegrator> interfacetest;
  interfacetest = Teuchos::rcp(new ThyraIntegrator(epetra_map,application));
  cerr << "Created an exIntegrate object" << endl;

  /*
  t0=3.5;
  // Now try integrating through the interface...
  (*interfacetest).Integrate(x,xn,t0,lambda0);

  // Create object NP of the CAPO namespace.
  Teuchos::RefCountPtr<CAPO::NP> v;
  v = Teuchos::rcp(new CAPO::NP(xn,t0,lambda0,interfacetest));

  (*v).Print_Version();
  //Now call the twice_integrate function of NP.
  (*v).Twice_Integrate(x);

  if (MyPID == 0)
    {
      cout << "After Integration the vector returns:" << endl;
      cout << "x(0) = " << Thyra::get_ele(*x,1) << endl;
      cout << "x(1) = " << Thyra::get_ele(*x,2) << endl;
    }
  */

  // Build a parameter list

  Teuchos::RefCountPtr<CAPO::Parameter_List> PL;
  PL = Teuchos::rcp(new CAPO::Parameter_List);
  cout << "Successfully Created a Parameter List" << endl;

  string pm = "FloquetTolerence";
  (*PL).set_param(pm,0.4);
  cout << "Set the floquet tolerence to .4" << endl;

  pm = "mplus2tol";
  (*PL).set_param(pm,.2);
  cout << "Tried and failed to set the Mplus2tol" << endl;

  cout << "************************************" << endl;

  cout << "I've created an Integrator which uses Thyra vectors " << endl;
  cout << "and a parameter list.  Next I want to create a solver object." << endl;
  Teuchos::RefCountPtr<CAPO::Npgs> MySolver;
  MySolver = Teuchos::rcp(new CAPO::Npgs(PL,interfacetest,x,t0,lambda0));

  cout << "The solver has been created.  Now time to test the solver " << endl;
  cout << "functionality." << endl;

  MySolver->Initialize();
  MySolver->Predictor();
  cout << MySolver->Get_Tfinal() << endl;
  cout << MySolver->Get_lambdafinal() << endl;

}
