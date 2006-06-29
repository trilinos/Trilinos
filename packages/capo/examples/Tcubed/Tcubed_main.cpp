//@HEADER
// ************************************************************************
// 
//                  LOCA Continuation Algorithm Package
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER
                                                                                
// 1D Finite Element Test Problem
/* Solves continuation problem (Parameter c="Right BC")
 *
 * d2u 
 * --- + a * u**3 = 0
 * dx2
 *
 * subject to @ x=0, u=b
 * subject to @ x=1, u=c
 */


// Trilinos Objects
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_LinearProblem.h"
//#include "AztecOO.h"

// User's application specific files 
//#include "Tcubed_integrator.hpp" // Problem Integrator
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Tcubed_interface.hpp" // Interface file to Capo Integrator
#include "Tcubed_FiniteElementProblem.H"              
#include "Capo_Integrator.hpp"
#include "Capo_Parameter_List.hpp"
#include "Capo_Rpm.hpp"
#include "Capo_Stepper.hpp"
using namespace std;

int main(int argc, char *argv[])
{
  typedef double Scalar; // Scalar type = double

  int ierr = 0, i;

  // Initialize MPI
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  // Create a communicator for Epetra objects
#ifdef HAVE_MPI
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif


  // Get the process ID and the total number of processors
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();
  
  // Get the number of elements from the command line
  if (argc!=2) { 
    cout << "Usage: " << argv[0] << " number_of_elements" << endl;
    exit(1);
  }
  int NumGlobalElements = atoi(argv[1]) + 1;



  // The number of unknowns must be at least equal to the 
  // number of processors.
  if (NumGlobalElements < NumProc) {
    cout << "numGlobalBlocks = " << NumGlobalElements 
	 << " cannot be < number of processors = " << NumProc << endl;
    exit(1);
  }
  Teuchos::RefCountPtr<const Epetra_Map> emap = Teuchos::rcp( new Epetra_Map(NumGlobalElements, 0, Comm) );


  // Create the FiniteElementProblem class.  This creates all required
  // Epetra objects for the problem and allows calls to the 
  // function (RHS) and Jacobian evaluation routines.
  Teuchos::RefCountPtr<Tcubed_FiniteElementProblem> Problem = Teuchos::rcp(new Tcubed_FiniteElementProblem(NumGlobalElements, Comm));

  string temp_string;
  temp_string = "Nonlinear Factor";
  Problem->setParameter(temp_string,1.0);
  temp_string = "Left BC";
  Problem->setParameter(temp_string,0.0);
  temp_string = "Right BC";
  Problem->setParameter(temp_string,0.1);

  // Get the vector from the Problem
  Epetra_Vector& soln = Problem->getSolution();
  soln.PutScalar(1.0);


  // Create x vector
  Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<Scalar> > xs  = Thyra::create_VectorSpace(emap);
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > x  = Thyra::createMember(xs);
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > removeme  = Thyra::createMember(xs);


  Thyra::assign(&*x,1.0);
  Thyra::assign(&*removeme,1.0);
  // Set Stepper Parameters
  double Param0 = 0.1;
  double dummy_time = 1.0; // This integrator doesn't have a time component.

  // Build a parameter list
  Teuchos::RefCountPtr<CAPO::Parameter_List> PL;
  PL = Teuchos::rcp(new CAPO::Parameter_List);

  temp_string = "MaxInnerIts";
  PL->set_param(temp_string,15);
  temp_string ="FloquetTolerence";
  PL->set_param(temp_string,0.5); // Not properly implemented yet.
  temp_string ="tol";
  PL->set_param(temp_string,1.0e-7);
  temp_string ="SubspaceIterations";
  PL->set_param(temp_string, 2);
  temp_string ="NumberXtraVecsSubspace";
  PL->set_param(temp_string,4);
  temp_string ="printproc";
  PL->set_param(temp_string,1);
  temp_string ="lambda_stepsize";
  PL->set_param(temp_string,.01);
  temp_string ="MaxOuterIts";
  PL->set_param(temp_string,130);
  temp_string ="EnablePeriodicity";
  PL->set_param(temp_string, false);
  temp_string ="EnableArclength";
  PL->set_param(temp_string, true);

  cout << "Successfully Created a Parameter List" << endl;

  // Build the Integrator Interface
  Teuchos::RefCountPtr<ThyraIntegrator> Integrator_Interface;
  Integrator_Interface = Teuchos::rcp(new ThyraIntegrator(Problem,emap));
  cerr << "Created the Integrator Interface." << endl;

  // Build Solver
  Teuchos::RefCountPtr<CAPO::Rpm> MySolver;
  MySolver = Teuchos::rcp(new CAPO::Rpm(PL,Integrator_Interface,x,Param0,dummy_time));
  cerr << "Built the Solver." << endl;

  // Build a Stepper
  Teuchos::RefCountPtr<CAPO::Stepper> MyStepper;
  MyStepper = Teuchos::rcp(new CAPO::Stepper(PL,MySolver));
   cerr << "Built the Problem Stepper. " << endl;

   //for (int i = 0; i<10 ; i++)
   //cout << soln[i] << " " <<endl;
   //cout << "***************" << endl;  
  //temp_string = "Right BC";
  //Problem->setParameter(temp_string,Param0);
  //Problem->evaluate(F_ONLY, &soln, &soln, NULL);
  //soln.Update(1.0,soln,1.0);
  //Integrator_Interface->Integrate(x,removeme,1.0,0.1);
  //for (int i = 0; i<10 ; i++)
  //cout << Thyra::get_ele(*removeme,i+1);

  //    cout << soln[i] << " " <<endl;
  


   // Solve
   MyStepper->Run();

}
