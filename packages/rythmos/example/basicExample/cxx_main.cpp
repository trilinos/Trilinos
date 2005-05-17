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

#include "Rythmos_ConfigDefs.h"

#include "ExampleApplication.hpp"

// Includes for Thyra:
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"

int main(int argc, char *argv[])
{

  typedef double Scalar; // Scalar type = double

  Teuchos::RefCountPtr<const Epetra_Comm> epetra_comm;
  Teuchos::RefCountPtr<const Epetra_Map> epetra_map;
  Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<Scalar> > epetra_vs;

#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  epetra_comm = Teuchos::rcp( new Epetra_MpiComm(MPI_COMM_WORLD) );
#else
  epetra_comm = Teuchos::rcp( new Epetra_SerialComm );
#endif

  int MyPID = epetra_comm->MyPID();

  if (MyPID == 0)
    cout << Rythmos::Rythmos_Version() << endl << endl;
	
  int NumElements = 1;

  // Construct a Map with NumElements and index base of 0
  epetra_map = Teuchos::rcp( new Epetra_Map(NumElements, 0, *epetra_comm) );

  // Construct a Thyra vector space
  epetra_vs = Thyra::create_MPIVectorSpaceBase(epetra_map);

  // Create x and xn vectors
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > x  = Thyra::createMember(epetra_vs);
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > xn = Thyra::createMember(epetra_vs);

  cout << "Integrating \\dot{x}=\\lambda x from t=0 to t=1" << endl
       << "with initial x_0 = 10, and \\Delta t=0.1" << endl
       << "using forward Euler." << endl;
  double lambda = -0.9;
  ExampleApplication problem(lambda);

  double t0 = 0.0;
  double t1 = 1.0;
  double dt = 0.1;
  double N = (t1-t0)/dt;
  double x_initial = 10.0; // initial condition
  Thyra::assign(&*xn,x_initial); // xn = x_initial
  cout << "x(0.0) = " << Thyra::get_ele(*xn,1) << endl;
  double t = t0;
  for (int i=1 ; i<N+1 ; ++i)
  {
    t = t0 + i*dt;
    Thyra::assign(&*x,*xn); // x = xn;
    problem.evalResidual( x, t ); 
    Thyra::Vp_StV(&*xn,dt,*x); // xn = xn + dt*x
    cout << "x(" << t << ") = " << Thyra::get_ele(*xn,1) << endl; 
  }
  cout << "       " << x_initial*exp(problem.getCoeff()*t) << " = Exact solution" << endl;
  return 0;
}

