// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Galeri_VbrMatrices.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Teuchos_ParameterList.hpp"

using namespace Galeri;

// =========== //
// main driver //
// =========== //

int main(int argc, char* argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // pointer to the map to be created
  Epetra_Map*            Map;
  // pointer to the CrsMatrix to be created
  Epetra_CrsMatrix*      CrsMatrix;
  // pointer to the VbrMatrix to be created
  Epetra_VbrMatrix*      VbrMatrix;
  // container for parameters
  Teuchos::ParameterList GaleriList;
  // here we specify the global dimension of the problem
  GaleriList.set("n", 2 * Comm.NumProc());

  try
  {
    // Creates a simple linear map; for more details on the map creation
    // refer to the documentation
    Map = CreateMap("Linear", Comm, GaleriList);

    // Creates a diagonal matrix with 1's on the diagonal
    CrsMatrix = CreateCrsMatrix("Diag", Map, GaleriList);

    // at this point we can create the VbrMatrix. Each block of the Vbr will
    // be 2 x 2. (Any positive value is accepted, like 1.)
    VbrMatrix = CreateVbrMatrix(CrsMatrix, 2);

    // print out the matrix
    cout << *VbrMatrix;

    // To created objects must be free'd using delete
    delete Map;
    delete CrsMatrix;
    delete VbrMatrix;
  }
  catch (Exception& rhs)
  {
    if (Comm.MyPID() == 0)
    {
      cerr << "Caught exception: ";
      rhs.Print();
    }
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
}
