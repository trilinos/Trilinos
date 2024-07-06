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
#include "Galeri_XpetraProblemFactory.hpp"

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

int main(int argc, char* argv[])
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  typedef Tpetra::Map<LO, GO, Node> Tpetra_Map;
  typedef Tpetra::CrsMatrix<Scalar, LO, GO, Node> Tpetra_CrsMatrix;
  typedef Tpetra::MultiVector<Scalar, LO, GO, Node> Tpetra_MultiVector;
  typedef Teuchos::ScalarTraits<Scalar> ScalarTraits;

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  // Create comm
  RCP<const Teuchos::Comm<int>> comm = Teuchos::DefaultComm<int>::getComm();

  // Here we create the linear problem
  //
  //   Matrix * LHS = RHS
  //
  // with Matrix arising from a 5-point formula discretization.

  std::string mapType = "Cartesian2D";
  auto mapParameters = Teuchos::ParameterList("Tpetra::Map");
  // dimension of the problem is nx x ny
  mapParameters.set("nx", 10 * comm->getSize());
  mapParameters.set("ny", 10);
  // total number of processors is mx x my
  mapParameters.set("mx", comm->getSize());
  mapParameters.set("my", 1);

  auto out = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));

  try
  {
    // Creation of the map
    auto map = RCP{Galeri::Xpetra::CreateMap<Scalar, GO, Tpetra_Map>(mapType, comm, mapParameters)};

    // Creation of linear problem
    auto problem = Galeri::Xpetra::BuildProblem<Scalar, LO, GO, Tpetra_Map, Tpetra_CrsMatrix, Tpetra_MultiVector>("Laplace2D", map, mapParameters);

    // Build Matrix and MultiVectors
    auto matrix = problem->BuildMatrix();
    auto LHS = rcp(new Tpetra_MultiVector(matrix->getDomainMap(), 1));
    auto RHS = rcp(new Tpetra_MultiVector(matrix->getRangeMap(), 1));
    auto ExactSolution = rcp(new Tpetra_MultiVector(matrix->getDomainMap(), 1));

    ExactSolution->randomize(0, 100);
    LHS->putScalar(ScalarTraits::zero());

    matrix->apply(*ExactSolution, *RHS);

    matrix->describe(*out, Teuchos::EVerbosityLevel::VERB_EXTREME);
    LHS->describe(*out, Teuchos::EVerbosityLevel::VERB_EXTREME);
    RHS->describe(*out, Teuchos::EVerbosityLevel::VERB_EXTREME);
    ExactSolution->describe(*out, Teuchos::EVerbosityLevel::VERB_EXTREME);

    // at this point any LinearSolver can be used which understands the Tpetra objects. For example: Amesos2 or Ifpack2
  }
  catch (Galeri::Exception &rhs)
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
