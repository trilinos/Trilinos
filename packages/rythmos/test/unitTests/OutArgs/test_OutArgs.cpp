//@HEADER
// ***********************************************************************
//
//                           Rythmos Package
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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

// Includes for Rythmos:
#include "EpetraExt_ModelEvaluator.h"

// Includes for Thyra:
#include "Thyra_VectorBase.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"

// Includes for Teuchos:
#include "Teuchos_RefCountPtr.hpp"

// Includes for Epetra:
#include "Epetra_SerialComm.h"
#include "Epetra_Map.h"

#include <iostream>

bool test1()
{
    bool test = false;
    Epetra_SerialComm epetra_comm;
    Teuchos::RefCountPtr<Epetra_Map> epetra_map = Teuchos::rcp(new Epetra_Map(1,0,epetra_comm));
    Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<double> > thyra_vs = Thyra::create_VectorSpace(epetra_map);
    Teuchos::RefCountPtr<Thyra::VectorBase<double> > in_vector = Thyra::createMember(thyra_vs);
    EpetraExt::OutArgs<double> outargs;
    outargs.request_F(in_vector);
    Teuchos::RefCountPtr<Thyra::VectorBase<double> > out_vector = outargs.get_F();
    if ( out_vector.get() == in_vector.get() )
      test = true; 
    return(test);
}

int main(int argc, char *argv[])
{
  bool success = false; 

  try { // catch exceptions

    success = test1(); 
    
   } // end try
    catch( const std::exception &excpt ) {
    std::cerr << "*** Caught a standard exception : " << excpt.what() << std::endl;
    success = false;
  }
  catch( ... ) {
    std::cerr << "*** Caught an unknown exception!\n";
    success = false;
  }

  if (success)
    std::cout << "Test successfull." << std::endl;
  else
    std::cout << "Test failed." << std::endl;

  return success ? 0 : 1;
} // end main() [Doxygen looks for this!]
