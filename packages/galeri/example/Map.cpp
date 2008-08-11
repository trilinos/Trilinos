// @HEADER
// ************************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
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
//
// Questions about Galeri? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
//
// ************************************************************************
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

int main(int argv, char* argc[])
{
#ifdef HAVE_MPI
  MPI_Init(&argv, &argc);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  bool verbose = (Comm.MyPID() == 0);

  // Creates an Epetra_Map corresponding to a 2D Cartesian grid
  // on the unit square. For parallel runs, the nodes are divided into
  // strips, so that the total number of subdomains is Comm.NumProc() x 1.
  
  // Pointer to the object to be created
  Epetra_Map* Map = 0; 
  // Type of the object
  string MapType = "Cartesian2D";
  // Container for parameters
  Teuchos::ParameterList GaleriList;
  GaleriList.set("nx", 2 * Comm.NumProc()); 
  GaleriList.set("ny", 2);
  GaleriList.set("mx", Comm.NumProc());
  GaleriList.set("my", 1);

  try
  {
    // Creation of the map
    Map = CreateMap("Cartesian2D", Comm, GaleriList);

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
