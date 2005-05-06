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

#include "Rythmos_ConfigDefs.h"

#include "ExampleApplication.hpp"

// Includes for Thyra:
//#include "Thyra_VectorSpaceBase.hpp"
//#include "Thyra_VectorBase.hpp"

int main(int argc, char *argv[])
{

#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm (MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int MyPID = Comm.MyPID();

  if (MyPID == 0)
    cout << Rythmos::Rythmos_Version() << endl << endl;
	
  int NumElements = 1;

  // Construct a Map with NumElements and index base of 0
  Epetra_Map Map(NumElements, 0, Comm);

  // Construct a Thyra vector space
//  Thyra::Thyra_VectorSpace vectorspace(NumElements, 0, Comm);
  // Create two vectors in this vector space
//  Thyra::Thyra_Vector x = vectorspace.createMember();
//  Thyra::Thyra_Vector xn = vectorspace.createMember();

  // Create x and xn vectors
  Epetra_Vector x(Map);
  Epetra_Vector xn(Map);

  cout << "Integrating \\dot{x}=\\lambda x from t=0 to t=1" << endl
       << "with initial x_0 = 10, and \\Delta t=0.1" << endl
       << "using forward Euler." << endl;
  double lambda = -0.9;
  ExampleApplication problem(lambda);
  double t0 = 0;
  double t1 = 1;
  double dt = 0.1;
  double N = (t1-t0)/dt;
  double x_initial = 10; // initial condition
  xn.PutScalar(x_initial); 
//  Thyra::assign(xn,x_initial);
  cout << "x(0.0) = " << xn[0] << endl;
  double t = t0;
  for (int i=1 ; i<N+1 ; ++i)
  {
    t = t0+i*dt;
    x.Scale(1.0,xn); // x = xn;
//    Thyra::assign(x,xn);
    problem.evalResidual(x,t);
    xn.Update(dt,x,1.0); // xn = xn + dt*x
//    Thyra::Vp_StV(xn,dt,x);
    cout << "x(" << t << ") = " << xn[0] << endl;
  }
  cout << "       " << x_initial*exp(problem.getCoeff()*t) << " = Exact solution" << endl;


  return 0;
}




