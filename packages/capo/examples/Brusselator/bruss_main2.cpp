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
File:      bruss_main2.cpp
Purpose:   Solve for Branch 2 of the Brusselator as defined 
           in Lust et. al.
Date:      6-29-05
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
#include "bruss_integrator.hpp"
#include "bruss_interface.hpp"
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
  Teuchos::RefCountPtr<const Epetra_Comm> epetra_comm = Teuchos::rcp( new Epetra_MpiComm(MPI_COMM_WORLD) );
#else
  Teuchos::RefCountPtr<const Epetra_Comm> epetra_comm = Teuchos::rcp( new Epetra_SerialComm );
#endif

  int MyPID = epetra_comm->MyPID();

  int NumElements = 62;

  // Construct a Map with NumElements and index base of 0
  Teuchos::RefCountPtr<const Epetra_Map> epetra_map = Teuchos::rcp( new Epetra_Map(NumElements, 0, *epetra_comm) );

  // Construct a Thyra vector space
  Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<Scalar> > thyra_vs = Thyra::create_VectorSpace(epetra_map);

  // Create x vector
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > x  = Thyra::createMember(thyra_vs);

Thyra::set_ele(1, 2.039300000000000, &*x);
Thyra::set_ele(32, 2.7272000000000000, &*x);
Thyra::set_ele(2, 2.0393000000000000, &*x);
Thyra::set_ele(33, 2.7272000000000000, &*x);
Thyra::set_ele(3, 2.0393000000000000, &*x);
Thyra::set_ele(34, 2.7272000000000000, &*x);
Thyra::set_ele(4, 2.0743000000000000, &*x);
Thyra::set_ele(35, 2.7300000000000000, &*x);
Thyra::set_ele(5, 2.0743000000000000, &*x);
Thyra::set_ele(36, 2.7300000000000000, &*x);
Thyra::set_ele(6, 2.0988000000000000, &*x);
Thyra::set_ele(37, 2.7330000000000000, &*x);
Thyra::set_ele(7, 2.0988000000000000, &*x);
Thyra::set_ele(38, 2.7330000000000000, &*x);
Thyra::set_ele(8, 2.1080000000000000, &*x);
Thyra::set_ele(39, 2.7342000000000000, &*x);
Thyra::set_ele(9, 2.1080000000000000, &*x);
Thyra::set_ele(40, 2.7342000000000000, &*x);
Thyra::set_ele(10, 2.098800000000000, &*x);
Thyra::set_ele(41, 2.7320000000000000, &*x);
Thyra::set_ele(11, 2.0988000000000000, &*x);
Thyra::set_ele(42, 2.7320000000000000, &*x);
Thyra::set_ele(12, 2.0755000000000000, &*x);
Thyra::set_ele(43, 2.7281000000000000, &*x);
Thyra::set_ele(13, 2.0755000000000000, &*x);
Thyra::set_ele(44, 2.7281000000000000, &*x);
Thyra::set_ele(14, 2.0392000000000000, &*x);
Thyra::set_ele(45, 2.7257000000000000, &*x);
Thyra::set_ele(15, 2.0392000000000000, &*x);
Thyra::set_ele(46, 2.7257000000000000, &*x);
Thyra::set_ele(16, 1.9970000000000000, &*x);
Thyra::set_ele(47, 2.7270000000000000, &*x);
Thyra::set_ele(17, 1.997000000000000, &*x);
Thyra::set_ele(48, 2.7270000000000000, &*x);
Thyra::set_ele(18, 1.9561000000000000, &*x);
Thyra::set_ele(49, 2.7316000000000000, &*x);
Thyra::set_ele(19, 1.9561000000000000, &*x);
Thyra::set_ele(50, 2.7316000000000000, &*x);
Thyra::set_ele(20, 1.9226000000000000, &*x);
Thyra::set_ele(51, 2.7374000000000000, &*x);
Thyra::set_ele(21, 1.9226000000000000, &*x);
Thyra::set_ele(52, 2.7374000000000000, &*x);
Thyra::set_ele(22, 1.9008000000000000, &*x);
Thyra::set_ele(53, 2.7420000000000000, &*x);
Thyra::set_ele(23, 1.9008000000000000, &*x);
Thyra::set_ele(54, 2.7420000000000000, &*x);
Thyra::set_ele(24, 1.8930000000000000, &*x);
Thyra::set_ele(55, 2.7440000000000000, &*x);
Thyra::set_ele(25, 1.8930000000000000, &*x);
Thyra::set_ele(56, 2.7440000000000000, &*x);
Thyra::set_ele(26, 1.9000000000000000, &*x);
Thyra::set_ele(57, 2.7428000000000000, &*x);
Thyra::set_ele(27, 1.9000000000000000, &*x);
Thyra::set_ele(58, 2.7428000000000000, &*x);
Thyra::set_ele(28, 1.9216000000000000, &*x);
Thyra::set_ele(59, 2.7390000000000000, &*x);
Thyra::set_ele(29, 1.9216000000000000, &*x);
Thyra::set_ele(60, 2.7390000000000000, &*x);
Thyra::set_ele(30, 1.9563000000000000, &*x);
Thyra::set_ele(61, 2.732900000000000, &*x);
Thyra::set_ele(31, 1.9563000000000000, &*x);
Thyra::set_ele(62, 2.7329000000000000, &*x);

  // Set Stepper Parameters
  //double Param0 =1.0361;
 double Param0 =1.04;
 //double t0 = 2.9529;
  double t0 = 2.948;


  if (MyPID == 0) 
  {
    cout << "Setting up for Capo Run: Set Initial Values." << endl;
  }


  // Create Brusselator Integrator
  Teuchos::RefCountPtr<BrussIntegrator> application;
  application = Teuchos::rcp(new BrussIntegrator);
  cout << "Created the Brusselator Integrator." << endl;

  // Build the Integrator Interface
  Teuchos::RefCountPtr<ThyraIntegrator> Integrator_Interface;
  Integrator_Interface = Teuchos::rcp(new ThyraIntegrator(application));
  cout << "Created the Integrator Interface." << endl;


  // Build a parameter list
  Teuchos::RefCountPtr<CAPO::Parameter_List> PL;
  PL = Teuchos::rcp(new CAPO::Parameter_List);

  string temp_string;
  temp_string = "MaxInnerIts";
  PL->set_param(temp_string,20);
  temp_string ="FloquetTolerence";
  PL->set_param(temp_string,0.5); // Not properly implemented yet.
  temp_string ="tol";
  PL->set_param(temp_string,1.0e-5);
  temp_string ="SubspaceIterations";
  PL->set_param(temp_string, 1);
  temp_string ="NumberXtraVecsSubspace";
  PL->set_param(temp_string,4);
  temp_string ="printproc";
  PL->set_param(temp_string,1);
  temp_string ="lambda_stepsize";
  PL->set_param(temp_string,.005);
  temp_string ="lambda_extend_tol";
  PL->set_param(temp_string,4);
  temp_string ="MaxOuterIts";
  PL->set_param(temp_string,190);
  temp_string ="EnablePeriodicity";
  PL->set_param(temp_string, true);
  temp_string ="EnableArclength";
  PL->set_param(temp_string, false);

  cout << "Successfully Created a Parameter List" << endl;


  // Build Solver
  Teuchos::RefCountPtr<CAPO::Npgs> MySolver;
  MySolver = Teuchos::rcp(new CAPO::Npgs(PL,Integrator_Interface,x,Param0,t0));
  cout << "Built the Solver." << endl;

  // Build a Stepper
  Teuchos::RefCountPtr<CAPO::Stepper> MyStepper;
  MyStepper = Teuchos::rcp(new CAPO::Stepper(PL,MySolver));
   cout << "Built the Problem Stepper. " << endl;

   // Solve
   MyStepper->Run();

}
