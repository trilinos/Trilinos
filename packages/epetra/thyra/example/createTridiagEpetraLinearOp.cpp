// @HEADER
// ***********************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright (2004) Sandia Corporation
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

#include "createTridiagEpetraLinearOp.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#ifdef HAVE_MPI
# include "Epetra_MpiComm.h"
#else
# include "Epetra_SerialComm.h"
#endif

Teuchos::RefCountPtr<Thyra::LinearOpBase<double> > createTridiagEpetraLinearOp(
  const int      globalDim
#ifdef HAVE_MPI
  ,MPI_Comm      mpiComm
#endif
  ,const double  diagScale
  ,const bool    verbose
  ,std::ostream  &out
  )
{
  using Teuchos::RefCountPtr;
  using Teuchos::rcp;

  //
  // (A) Create Epetra_Map
  //

#ifdef HAVE_MPI
  if(verbose) out << "\nCreating Epetra_MpiComm ...\n";
  Epetra_MpiComm epetra_comm(mpiComm); // Note, Epetra_MpiComm is just a handle class to Epetra_MpiCommData object!
#else
  if(verbose) out << "\nCreating Epetra_SerialComm ...\n";
  Epetra_SerialComm epetra_comm; // Note, Epetra_SerialComm is just a handle class to Epetra_SerialCommData object!
#endif
  // Create the Epetra_Map object giving it the handle to the Epetra_Comm object
  const Epetra_Map epetra_map(globalDim,0,epetra_comm); // Note, Epetra_Map is just a handle class to an Epetra_BlockMapData object!
  // Note that the above Epetra_Map object "copies" the Epetra_Comm object in some
  // since so that memory mangaement is performed safely.

  //
  // (B) Create the tridiagonal Epetra object
  //
  //       [  2  -1             ]
  //       [ -1   2  -1         ]
  //  A =  [      .   .   .     ]
  //       [          -1  2  -1 ]
  //       [             -1   2 ]
  //
  
  // (B.1) Allocate the Epetra_CrsMatrix object.
  RefCountPtr<Epetra_CrsMatrix> A_epetra = rcp(new Epetra_CrsMatrix(::Copy,epetra_map,3));
  // Note that Epetra_CrsMatrix is *not* a handle object so have to use RefCountPtr above.
  // Also note that the Epetra_Map object is copied in some sence by the Epetra_CrsMatrix
  // object so the memory managment is handled safely.

  // (B.2) Get the indexes of the rows on this processor
  const int numMyElements = epetra_map.NumMyElements();
  std::vector<int> myGlobalElements(numMyElements);
  epetra_map.MyGlobalElements(&myGlobalElements[0]);

  // (B.3) Fill the local matrix entries one row at a time.
  // Note, we set up Epetra_Map above to use zero-based indexing and that is what we must use below:
  const double offDiag = -1.0, diag = 2.0*diagScale;
  int numEntries; double values[3]; int indexes[3];
  for( int k = 0; k < numMyElements; ++k ) {
    const int rowIndex = myGlobalElements[k];
    if( rowIndex == 0 ) {                     // First row
      numEntries = 2;
      values[0]  = diag;             values[1]  = offDiag;
      indexes[0] = 0;                indexes[1] = 1; 
    }
    else if( rowIndex == globalDim - 1 ) {    // Last row
      numEntries = 2;
      values[0]  = offDiag;         values[1]  = diag;
      indexes[0] = globalDim-2;     indexes[1] = globalDim-1; 
    }
    else {                                    // Middle rows
      numEntries = 3;
      values[0]  = offDiag;         values[1]  = diag;          values[2]  = offDiag;
      indexes[0] = rowIndex-1;      indexes[1] = rowIndex;      indexes[2] = rowIndex+1; 
    }
    TEST_FOR_EXCEPT( 0!=A_epetra->InsertGlobalValues(rowIndex,numEntries,values,indexes) );
  }

  // (B.4) Finish the construction of the Epetra_CrsMatrix
  TEST_FOR_EXCEPT( 0!=A_epetra->FillComplete() );

  //
  // (C) Wrap the above created (Epetra_CrsMatrix) Epetra_Operator object created above
  // into a Thyra::EpetraLinearOp object to turn it into a Thyra::LinearOpBase object
  //

  RefCountPtr<Thyra::LinearOpBase<double> > A = rcp(new Thyra::EpetraLinearOp(A_epetra));

  //
  // (D) Finally return the Thyra-wrapped Epetra matrix object
  //
  
  return A;

  // Note that when this function returns the returned
  // RefCountPtr-wrapped Thyra::LinearOpBase object will own all of the
  // Epetra objects that went into its construction and these objects
  // will stay around until all of the RefCountPtr objects to the
  // allocated Thyra::LinearOpBase object are removed and destruction
  // occurs!

} // end createTridiagLinearOp()
