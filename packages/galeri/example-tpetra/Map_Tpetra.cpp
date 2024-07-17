// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Galeri_XpetraMaps.hpp"

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_ParameterList.hpp"

#ifndef GALERI_TEST_USE_LONGLONG_GO
#define GO int
#else
#define GO long long
#endif

using namespace Galeri;

int main(int argc, char *argv[])
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  typedef Tpetra::Map<int, GO> Tpetra_Map;

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  // Create comm
  RCP<const Teuchos::Comm<int>> comm = Teuchos::DefaultComm<int>::getComm();

  // Creates an Tpetra_Map corresponding to a 2D Cartesian grid
  // on the unit square. For parallel runs, the nodes are divided into
  // strips, so that the total number of subdomains is comm->getSize() x 1.

  // Type of the object
  std::string mapType = "Cartesian2D";
  // Container for parameters
  Teuchos::ParameterList galeriList;
  galeriList.set("nx", 2 * comm->getSize());
  galeriList.set("ny", 2);
  galeriList.set("mx", comm->getSize());
  galeriList.set("my", 1);

  try
  {
    // Creation of the map
    auto map = Galeri::Xpetra::CreateMap<int, GO, Tpetra_Map>(mapType, comm, galeriList);

    // print out the map
    auto out = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
    map->describe(*out, Teuchos::EVerbosityLevel::VERB_EXTREME);
  }
  catch (Exception &rhs)
  {
    if (comm->getRank() == 0)
    {
      cerr << "Caught exception: ";
      rhs.Print();

#ifdef HAVE_MPI
      MPI_Finalize();
#endif
      return (EXIT_FAILURE);
    }
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return (EXIT_SUCCESS);
}
