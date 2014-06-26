//@HEADER
// ***********************************************************************
//
//                     Rythmos Package
//                 Copyright (2006) Sandia Corporation
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#include "Epetra_SerialComm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Version.h"

#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"


int main(int argc, char *argv[])
{
  cout << Epetra_Version() << endl << endl;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  int NumElements = 1000;

  typedef double Scalar;

  Teuchos::RCP<const Epetra_Comm> epetra_comm;
  epetra_comm = Teuchos::rcp(new Epetra_SerialComm);

  // Construct a Map with NumElements and index base of 0
  Teuchos::RCP<const Epetra_Map> epetra_map;
  epetra_map = Teuchos::rcp(new Epetra_Map(NumElements, 0, *epetra_comm));
  // Construct a VectorSpace from the map
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > epetra_vs;
  epetra_vs = Thyra::create_VectorSpace(epetra_map);
  

  // Create x and b vectors
  Teuchos::RCP<Thyra::VectorBase<Scalar> > b = Thyra::createMember(epetra_vs);
  Teuchos::RCP<Thyra::VectorBase<Scalar> > x = Thyra::createMember(epetra_vs);

  Thyra::randomize(0.0,1.0,&*b); // b = random
  Thyra::V_StV(&*x,2.0,*b); // x = 2*b

  double bnorm, xnorm;
  xnorm = Thyra::norm_2(*x);
  bnorm = Thyra::norm_2(*b);

  cout << "2 norm of x = " << xnorm << endl
       << "2 norm of b = " << bnorm << endl;

//  Test out error messages from bad pointer arithmetic
//  cout << "Element 0 of x = " << (&*x)[0] << endl; // bad error message
//  cout << "Element 0 of x = " << (*x)[0] << endl; // good error message
//  for (int i=1;i<=NumElements;++i)
//    cout << "Element " << i << " of x = " << Thyra::get_ele(*x,i) << endl;

// Test out get_Epetra_Vector:
//  Teuchos::RCP<Epetra_Vector> x_epetra = Thyra::get_Epetra_Vector(*epetra_map,x);
//  for (int i=0;i<NumElements;++i)
//    cout << "Element " << i << " of x = " << (*x_epetra)[i] << endl;

  return 0;
}

