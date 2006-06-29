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
File:      rpmtest_main.cpp
Purpose:   Solve A. Salinger's integration problem 
Date:      7-01-05
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
#include "rpmtest_integrator.hpp"
#include "rpmtest_interface.hpp"
#include "Capo_Integrator.hpp"
#include "Capo_Parameter_List.hpp"
#include "Capo_Rpm.hpp"
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

  int NumElements = 202;

  // Construct a Map with NumElements and index base of 0
  Teuchos::RefCountPtr<const Epetra_Map> epetra_map = Teuchos::rcp( new Epetra_Map(NumElements, 0, *epetra_comm) );

  // Construct a Thyra vector space
  Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<Scalar> > thyra_vs = Thyra::create_VectorSpace(epetra_map);

  // Create x vector
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > x  = Thyra::createMember(thyra_vs);

  double xp049[]={
    0.861701, 0.225987,
    0.861701, 0.225987,
    0.861866, 0.225772,
    0.862195, 0.22534,
    0.862689, 0.224693,
    0.863344, 0.223828,
    0.864158, 0.222744,
    0.86513, 0.221441,
    0.866254, 0.219916,
    0.867526, 0.218167,
    0.86894, 0.216193,
    0.870488, 0.213989,
    0.872161, 0.211554,
    0.873947, 0.208885,
    0.87583, 0.205977,
    0.877793, 0.202827,
    0.879811, 0.19943,
    0.881856, 0.195784,
    0.883889, 0.191882,
    0.885864, 0.187721,
    0.887721, 0.183295,
    0.889384, 0.178599,
    0.890756, 0.173627,
    0.891713, 0.168374,
    0.892099, 0.162834,
    0.891713, 0.157003,
    0.890301, 0.150873,
    0.887539, 0.14444,
    0.883019, 0.137699,
    0.876227, 0.130646,
    0.866522, 0.123277,
    0.853112, 0.115589,
    0.835039, 0.107582,
    0.811158, 0.0992559,
    0.780151, 0.0906152,
    0.740556, 0.0816663,
    0.690852, 0.0724198,
    0.6296, 0.0628911,
    0.555663, 0.0531008,
    0.468486, 0.0430755,
    0.368406, 0.0328482,
    0.256903, 0.0224579,
    0.136699, 0.0119489,
    0.0116068, 0.00136999,
    -0.113895, -0.00922799,
    -0.235269, -0.0197939,
    -0.348545, -0.0302787,
    -0.450782, -0.040637,
    -0.540278, -0.050829,
    -0.616502, -0.0608208,
    -0.679872, -0.0705848,
    -0.731439, -0.0800995,
    -0.772605, -0.089349,
    -0.804888, -0.0983222,
    -0.829766, -0.107012,
    -0.848586, -0.115415,
    -0.862521, -0.12353,
    -0.872564, -0.131358,
    -0.879532, -0.138901,
    -0.88409, -0.146165,
    -0.886771, -0.153153,
    -0.888001, -0.159872,
    -0.888115, -0.166327,
    -0.887378, -0.172524,
    -0.885996, -0.178471,
    -0.884134, -0.184172,
    -0.881918, -0.189636,
    -0.879448, -0.194868,
    -0.876803, -0.199875,
    -0.874043, -0.204663,
    -0.871218, -0.209237,
    -0.868364, -0.213605,
    -0.865512, -0.217771,
    -0.862685, -0.221741,
    -0.859901, -0.22552,
    -0.857175, -0.229114,
    -0.85452, -0.232527,
    -0.851944, -0.235763,
    -0.849455, -0.238828,
    -0.847058, -0.241726,
    -0.84476, -0.24446,
    -0.842563, -0.247034,
    -0.84047, -0.249452,
    -0.838485, -0.251717,
    -0.83661, -0.253833,
    -0.834845, -0.255802,
    -0.833194, -0.257628,
    -0.831657, -0.259312,
    -0.830235, -0.260858,
    -0.828928, -0.262267,
    -0.827739, -0.263542,
    -0.826667, -0.264684,
    -0.825713, -0.265695,
    -0.824877, -0.266577,
    -0.82416, -0.26733,
    -0.823562, -0.267956,
    -0.823083, -0.268456,
    -0.822724, -0.26883,
    -0.822484, -0.269079,
    -0.822364, -0.269204,
    -0.822364, -0.269204
  };

  for (int i=0; i<NumElements/2; i++) {
    Thyra::set_ele(i+1, xp049[2*i], &*x);
    Thyra::set_ele(i+1+NumElements/2, xp049[2*i+1], &*x);
  }


  // Set Stepper Parameters
  double Param0 = 0.0485;
  double t0 = 1.0; //dummy variable for this problem.


  if (MyPID == 0) 
  {
    cout << "Setting up for Capo Run: Set Initial Values." << endl;
  }


  // Create Integrator
  Teuchos::RefCountPtr<RpmtestIntegrator> application;
  application = Teuchos::rcp(new RpmtestIntegrator);
  cout << "Created the Rpm Test Integrator." << endl;

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
  PL->set_param(temp_string,1.0e-6);
  temp_string ="SubspaceIterations";
  PL->set_param(temp_string, 1);
  temp_string ="NumberXtraVecsSubspace";
  PL->set_param(temp_string,4);
  temp_string ="printproc";
  PL->set_param(temp_string,1);
  temp_string ="lambda_stepsize";
  PL->set_param(temp_string,-.001);
  temp_string ="MaxOuterIts";
  PL->set_param(temp_string,40);
  temp_string ="EnablePeriodicity";
  PL->set_param(temp_string, false);
  temp_string ="EnableArclength";
  PL->set_param(temp_string, true);


  cout << "Successfully Created a Parameter List" << endl;


  // Build Solver
  Teuchos::RefCountPtr<CAPO::Rpm> MySolver;
  MySolver = Teuchos::rcp(new CAPO::Rpm(PL,Integrator_Interface,x,Param0,t0));
  cout << "Built the Solver." << endl;

  // Build a Stepper
  Teuchos::RefCountPtr<CAPO::Stepper> MyStepper;
  MyStepper = Teuchos::rcp(new CAPO::Stepper(PL,MySolver));
   cout << "Built the Problem Stepper. " << endl;

   // Solve
   MyStepper->Run();

}

