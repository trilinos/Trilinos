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
File:      EA_main.cpp
Purpose:   Solve for the period doubling of the E-A model as defined 
           in K.Engelborghs et. al.
Date:      8-1-05
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
#include "EA_integrator.hpp"
#include "EA_interface.hpp"
#include "Capo_Integrator.hpp"
#include "Capo_Parameter_List.hpp"
#include "Capo_NppD.hpp"
#include "Capo_Stepper.hpp"
#include <string>

int main(int argc, char *argv[])
{

  typedef double Scalar; // Scalar type = double


#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Teuchos::RefCountPtr<const Epetra_Comm> epetra_comm = Teuchos::rcp( new Epetra_MpiComm(MPI_COMM_WORLD) );
#else
  Teuchos::RefCountPtr<const Epetra_Comm> epetra_comm = Teuchos::rcp( new Epetra_SerialComm );
#endif

  int MyPID = epetra_comm->MyPID();
  // Get the number of elements from the command line
  if (argc!=2) { 
    cout << "Usage: " << argv[0] << " number_of_elements" << endl;
    exit(1);
  }
  int NumGlobalElements = atoi(argv[1]) + 1;


  // Construct a Map with NumElements and index base of 0
  Teuchos::RefCountPtr<const Epetra_Map> epetra_map = Teuchos::rcp( new Epetra_Map(NumGlobalElements, 0, *epetra_comm) );

  // Construct a Thyra vector space
  Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<Scalar> > thyra_vs = Thyra::create_VectorSpace(epetra_map);

  // Create x vector
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > x  = Thyra::createMember(thyra_vs);

  double h = 1.0/(NumGlobalElements/2+1);
  for (int i=1;i<(NumGlobalElements/2)+1;i++)
    {
      Thyra::set_ele(i,-2.0+exp((1-i*h)*(i*h)),&*x);
      Thyra::set_ele(NumGlobalElements/2+i,-4.0+(1-i*h)*(i*h),&*x);
    }
  // Set Stepper Parameters
  double Param0 = .0318;
  double t0 = 1.59;


  if (MyPID == 0) 
  {
    cout << "Setting up for Capo Run: Set Initial Values." << endl;
  }


  // Create Elezgaray Integrator
  Teuchos::RefCountPtr<EAIntegrator> application;
  application = Teuchos::rcp(new EAIntegrator(NumGlobalElements));
  cout << "Created the Elezgaray Arneodo Integrator." << endl;

  // Build the Integrator Interface
  Teuchos::RefCountPtr<ThyraIntegrator> Integrator_Interface;
  Integrator_Interface = Teuchos::rcp(new ThyraIntegrator(application,NumGlobalElements));
  cout << "Created the Integrator Interface." << endl;

  
  //To get an initial solution, we integrate x for awhile, then use it.
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > y  = Thyra::createMember(thyra_vs);
  Thyra::assign(&*y,0.0);
  Integrator_Interface->Integrate(x,y,5.5*t0,Param0);
  Thyra::assign(&*x,*y);
  Integrator_Interface->Integrate(x,y,1.0*t0,Param0);

  Thyra::Vp_StV(&*y,-1.0,*x);
  for (int j=0;j<NumGlobalElements;j++)
    cout << Thyra::get_ele(*y,j+1) << endl;
  

  /*
  // Build a parameter list
  Teuchos::RefCountPtr<CAPO::Parameter_List> PL;
  PL = Teuchos::rcp(new CAPO::Parameter_List);

  string temp_string;
  temp_string = "MaxInnerIts";
  PL->set_param(temp_string,20);
  temp_string ="FloquetTolerence";
  PL->set_param(temp_string,0.5); // Not properly implemented yet.
  temp_string ="tol";
  PL->set_param(temp_string,1.0e-6);
  temp_string ="SubspaceIterations";
  PL->set_param(temp_string, 3);
  temp_string ="NumberXtraVecsSubspace";
  PL->set_param(temp_string,4);
  temp_string ="printproc";
  PL->set_param(temp_string,1);
  temp_string ="lambda_stepsize";
  PL->set_param(temp_string,-.01);
  temp_string ="MaxOuterIts";
  PL->set_param(temp_string,325);
  temp_string ="EnablePeriodicity";
  PL->set_param(temp_string, true);
  temp_string ="EnableArclength";
  PL->set_param(temp_string, true);

  cout << "Successfully Created a Parameter List" << endl;


  // Build Solver
  Teuchos::RefCountPtr<CAPO::NppD> MySolver;
  MySolver = Teuchos::rcp(new CAPO::NppD(PL,Integrator_Interface,x,Param0,t0));
  cout << "Built the Solver." << endl;

  // Build a Stepper
  Teuchos::RefCountPtr<CAPO::Stepper> MyStepper;
  MyStepper = Teuchos::rcp(new CAPO::Stepper(PL,MySolver));
   cout << "Built the Problem Stepper. " << endl;

   // Solve
   MyStepper->Run();
  */
}

