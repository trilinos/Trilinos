//@HEADER
// ************************************************************************
// 
//                          Rythmos Package 
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
// ************************************************************************
//@HEADER

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Version.h"

// Includes for Rythmos:
#include "Rythmos_ConfigDefs.h"
#include "ExampleApplicationRythmosInterface.hpp"
#include "Stepper_ForwardEuler.hpp"
#include "Stepper_ExplicitRK.hpp"

// Includes for Thyra:
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"

#include <string>

// Includes for Teuchos:
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

int main(int argc, char *argv[])
{

  double lambda = -0.5;  // ODE coefficient
  int numElements = 1; // number of elements in vector
  double x0 = 10.0; // ODE initial condition
  int N = 10;  // number of steps to take
  string method = "FE";  // other choice is method="ERK4"


  // Parse the command-line options:
  Teuchos::CommandLineProcessor  clp(false); // Don't throw exceptions
  clp.setOption( "lambda", &lambda, "ODE coefficient of decay or growth." );
  clp.setOption( "x0", &x0, "ODE initial condition." );
  clp.setOption( "method", &method, "Integration method:  FE, ERK4." );
  clp.setOption( "numsteps", &N, "Number of integration steps to take" );
  Teuchos::CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
  if( parse_return != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;


  // Set up the parameter list for the application:
  Teuchos::ParameterList params;
  params.set( "Lambda", lambda );
  params.set( "NumElements", numElements );
  params.set( "x0", x0 );
  
  // create interface to problem
  Teuchos::RefCountPtr<ExampleApplicationRythmosInterface> problem = Teuchos::rcp(new ExampleApplicationRythmosInterface(params));
  
  // Create Stepper object depending on command-line input
  Teuchos::RefCountPtr<Rythmos::Stepper<double> > stepper_ptr;
  if (method == "ERK4") // Explicit Runge-Kutta 4 stage
  {
    stepper_ptr = Teuchos::rcp(new Rythmos::ExplicitRK<double>(problem));
    method = "Explicit Runge-Kutta of order 4";
  } else // Forward Euler
  {
    stepper_ptr = Teuchos::rcp(new Rythmos::ForwardEuler<double>(problem));
    method = "Forward Euler";
  }
  Rythmos::Stepper<double> &stepper = *stepper_ptr;


  double t0 = 0.0;
  double t1 = 1.0;
  double dt = (t1-t0)/N;

  // Integrate forward with fixed step sizes:
  for (int i=1 ; i<=N ; ++i)
  {
    double dt_taken = stepper.TakeStep(dt);
    if (dt_taken != dt)
    {
      cerr << "Error, stepper took step of dt = " << dt_taken << " when asked to take step of dt = " << dt << endl;
      break;
    }
  }
  // Get solution out of stepper:
  Teuchos::RefCountPtr<const Thyra::VectorBase<double> > x_t = stepper.get_solution();
  // Convert Thyra::VectorBase to Epetra_Vector
  Teuchos::RefCountPtr<const Epetra_Vector> x = Thyra::get_Epetra_Vector(*(problem->get_Epetra_Map()),x_t);

  // These values should be passed by parameter list:
  // hard-coded values in ExampleApplicationRythmosInterface:
  // double numelements = 1;
  double x_initial = x0;

  // check exact answer
  double x_star = x_initial*exp(lambda*t1);

  // 06/03/05 tscoffe to get an Epetra_Map associated with an Epetra_Vector:
  // x.Map()
  // to get an Epetra_Comm associated with an Epetra_Vector:
  // x.Comm()
  
//  Teuchos::RefCountPtr<const Epetra_Comm> epetra_comm = (*problem).get_epetra_comm();
  //int MyPID = problem->get_Epetra_Map()->Comm()->MyPID();
  int MyPID = (*x).Comm().MyPID();
//  int MyPID = epetra_comm->MyPID();
  if (MyPID == 0)
  {
    cout << Rythmos::Rythmos_Version() << endl << endl;
    cout << "Integrating \\dot{x}=\\lambda x from t = " << t0 
         << " to t = " << t1 << endl
         << "with initial x_0 = " << x_initial 
         << ", \\Delta t = " << dt 
         << ", and \\lambda = " << lambda << endl
         << "using " << method << "." << endl;

    cout << "Computed: x(" << t1 << ") = " << (*x)[0] << endl;
    cout << "Exact:    x(" << t1 << ") = " << x_star << endl;
  }
  
  return 0;
}

