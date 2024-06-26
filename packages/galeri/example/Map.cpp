// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Galeri_Maps.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Teuchos_ParameterList.hpp"

using namespace Galeri;

int main(int argc, char* argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // Creates an Epetra_Map corresponding to a 2D Cartesian grid
  // on the unit square. For parallel runs, the nodes are divided into
  // strips, so that the total number of subdomains is Comm.NumProc() x 1.
  
  // Pointer to the object to be created
  Epetra_Map* Map = 0; 
  // Type of the object
  std::string MapType = "Cartesian2D";
  // Container for parameters
  Teuchos::ParameterList GaleriList;
  GaleriList.set("nx", 2 * Comm.NumProc()); 
  GaleriList.set("ny", 2);
  GaleriList.set("mx", Comm.NumProc());
  GaleriList.set("my", 1);

  try
  {
    // Creation of the map
#ifndef GALERI_TEST_USE_LONGLONG_GO
    Map = CreateMap("Cartesian2D", Comm, GaleriList);
#else
    Map = CreateMap64("Cartesian2D", Comm, GaleriList);
#endif
    // print out the map
    cout << *Map;

    // To created object must be free'd using delete
    delete Map;
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
