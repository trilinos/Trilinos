//
// @HEADER
// ***********************************************************************
// 
//                           Rythmos Package
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


#include "ExampleApplication.hpp"


ExampleApplication::ExampleApplication(Teuchos::ParameterList &params)
{
  lambda_ = params.get( "Lambda", -0.5 );
  numElements_ = params.get( "NumElements", 1 );
  x0_ = params.get( "x0", 10.0 );
  // Serial only implementation here:
  // 05/26/05 tscoffe:  I haven't figured out how to get MPI_Init called with
  // argc and argv in such a way that MPI_COMM_WORLD is passed down here.
  epetra_comm_ = Teuchos::rcp( new Epetra_SerialComm  );
//  std::cout << "ExampleApplication::ExampleApplication(double lambda, int numElements)" << std::endl;
  
  // Construct a Map with NumElements and index base of 0
  epetra_map_ = Teuchos::rcp( new Epetra_Map(numElements_, 0, *epetra_comm_) );
//  std::cout << "Epetra_Map address = " << std::endl;
//  std::cout << epetra_map_.get() << std::endl;
  
}

int ExampleApplication::evalResidual(Epetra_Vector *y, const Epetra_Vector &x, double t)
{
  y->Scale(lambda_,x); // y = lambda*x
  return 0;
}


double ExampleApplication::getCoeff()
{
  return lambda_;
}

Teuchos::RefCountPtr<Epetra_Map> ExampleApplication::get_epetra_map()
{
  return(epetra_map_);
}

Teuchos::RefCountPtr<Epetra_Comm> ExampleApplication::get_epetra_comm()
{
  return(epetra_comm_);
}


Teuchos::RefCountPtr<Epetra_Vector> ExampleApplication::get_x0()
{
  Teuchos::RefCountPtr<Epetra_Vector> x0 = Teuchos::rcp(new Epetra_Vector(*epetra_map_));
  (*x0)[0] = x0_;
//  x0->Random();
  return(x0);
}

