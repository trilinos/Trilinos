// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

Teuchos::RCP<Epetra_Operator>
createTridiagEpetraLinearOp(
  const int      globalDim
#ifdef HAVE_MPI
  ,MPI_Comm      mpiComm
#endif
  ,const double  diagScale
  ,const bool    verbose
  ,std::ostream  &out
  )
{
  using Teuchos::RCP;
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
  RCP<Epetra_CrsMatrix> A_epetra = rcp(new Epetra_CrsMatrix(::Copy,epetra_map,3));
  // Note that Epetra_CrsMatrix is *not* a handle object so have to use RCP above.
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
    TEUCHOS_TEST_FOR_EXCEPT( 0!=A_epetra->InsertGlobalValues(rowIndex,numEntries,values,indexes) );
  }

  // (B.4) Finish the construction of the Epetra_CrsMatrix
  TEUCHOS_TEST_FOR_EXCEPT( 0!=A_epetra->FillComplete() );

  // Return the Epetra_Operator object
  return A_epetra;

  // Note that when this function returns the returned RCP-wrapped
  // Epetra_Operator object will own all of the Epetra objects that went into
  // its construction and these objects will stay around until all of the
  // RCP objects to the allocated Epetra_Operator object are removed
  // and destruction occurs!

} // end createTridiagLinearOp()
