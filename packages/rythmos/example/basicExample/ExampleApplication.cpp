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
  lambda_min_ = params.get( "Lambda_min", -0.9 );
  lambda_max_ = params.get( "Lambda_max", -0.01 );
  numElements_ = params.get( "NumElements", 1 );
  x0_ = params.get( "x0", 10.0 );
  // Serial only implementation here:
  // 05/26/05 tscoffe:  I haven't figured out how to get MPI_Init called with
  // argc and argv in such a way that MPI_COMM_WORLD is passed down here.
#ifdef HAVE_MPI
  int   argc = params.get( "main_argc" );
  char *argv = params.get( "main_argv" );
  MPI_Init(&argc,&argv);
  MPI_Comm mpiComm = MPI_COMM_WORLD;
  int procRank = 0;
  int numProc;
  MPI_Comm_size( mpiComm, &numProc );
  MPI_Comm_rank( mpiComm, &procRank );
  epetra_comm_ptr_ = Teuchos::rcp( new Epetra_Comm(mpiComm) );
#else
  epetra_comm_ptr_ = Teuchos::rcp( new Epetra_SerialComm  );
#endif // HAVE_MPI


  // Construct a Map with NumElements and index base of 0
  epetra_map_ptr_ = Teuchos::rcp( new Epetra_Map(numElements_, 0, *epetra_comm_ptr_) );

  lambda_ptr_ = Teuchos::rcp(new Epetra_Vector(*epetra_map_ptr_));
  Epetra_Vector &lambda = *lambda_ptr_;
  unsigned int seed = time(NULL); 
  seed *= seed;
  lambda.SetSeed(seed);
  lambda.Random(); // fill with random numbers in (-1,1)
  // Scale random numbers to (lambda_min_,lambda_max_)
  lambda.Scale( (lambda_min_ - lambda_max_)/2.0);
  double tmp = (lambda_max_ + lambda_min_)/2.0;
  for (int i=0 ; i<lambda.MyLength() ; ++i)
  {
    lambda[i] += tmp;
  }
  
}

int ExampleApplication::evalResidual(Epetra_Vector *y, const Epetra_Vector &x, double t)
{
  Epetra_Vector &lambda = *lambda_ptr_;
  int localNumElements = x.MyLength();
  for (int i=0 ; i<localNumElements ; ++i)
  {
    (*y)[i] = lambda[i]*x[i];
  }
  return 0;
}


Teuchos::RefCountPtr<const Epetra_Vector> ExampleApplication::get_coeff() const
{
  return(lambda_ptr_);
}

Teuchos::RefCountPtr<Epetra_Map> ExampleApplication::get_epetra_map() const
{
  return(epetra_map_ptr_);
}

Teuchos::RefCountPtr<Epetra_Comm> ExampleApplication::get_epetra_comm() const
{
  return(epetra_comm_ptr_);
}


Teuchos::RefCountPtr<Epetra_Vector> ExampleApplication::get_x0() const
{
  Teuchos::RefCountPtr<Epetra_Vector> x0_ptr = Teuchos::rcp(new Epetra_Vector(*epetra_map_ptr_));
  Epetra_Vector &x0 = *x0_ptr;
  for (int i=0 ; i<numElements_ ; ++i)
  {
    x0[i] = x0_;
  }
  return(x0_ptr);
}

