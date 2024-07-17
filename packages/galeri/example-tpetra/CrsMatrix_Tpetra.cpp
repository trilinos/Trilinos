// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Galeri_XpetraMaps.hpp"
#include "Galeri_MatrixTraits.hpp"
#include "Galeri_XpetraMatrixTypes.hpp"

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_ParameterList.hpp"

#ifndef GALERI_TEST_USE_LONGLONG_GO
#define GO int
#else
#define GO long long
#endif
#define Scalar int
#define LO int
#define Node Tpetra::KokkosClassic::DefaultNode::DefaultNodeType

using namespace Galeri;

int main(int argc, char *argv[])
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  typedef Tpetra::Map<LO, GO, Node> Tpetra_Map;
  typedef Tpetra::CrsMatrix<Scalar, LO, GO, Node> Tpetra_CrsMatrix;

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  // Create comm
  RCP<const Teuchos::Comm<int>> comm = Teuchos::DefaultComm<int>::getComm();

  std::string mapType = "Cartesian2D";
  auto mapParameters = Teuchos::ParameterList("Tpetra::Map");
  // here we specify the global dimension of the problem
  int nx = 5 * comm->getSize();
  int ny = 5 * comm->getSize();
  mapParameters.set("nx", nx);
  mapParameters.set("ny", ny);

  try
  {
    // Creation of the map
    auto map = RCP{Galeri::Xpetra::CreateMap<Scalar, GO, Tpetra_Map>(mapType, comm, mapParameters)};

    // Creates a diagonal matrix with 1's on the diagonal
    auto matrix = Galeri::Xpetra::BigStar2D<Scalar, LO, GO, Tpetra_Map, Tpetra_CrsMatrix>(map, nx, ny, 20, -8, -8, -8, -8, 2, 2, 2, 2, 1, 1, 1, 1);

    // print out the matrix
    auto out = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
    matrix->describe(*out, Teuchos::EVerbosityLevel::VERB_EXTREME);
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
