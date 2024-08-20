// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
//
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/**
  \file   SimpleSolve.cpp
  \author Eric Bavier <etbavie@sandia.gov>
  \date   Sat Jul 17 10:35:39 2010

  \brief  Simple example of Amesos2 usage.

  This example solves a simple sparse system of linear equations using the
  Amesos2 interface to the Superlu solver.
  */

#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_VerboseObject.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include "Amesos2.hpp"
#include "Amesos2_Version.hpp"

int main(int argc, char *argv[]) {
  Tpetra::ScopeGuard mpiSession(&argc, &argv);

  typedef Tpetra::MultiVector<> MV;
  typedef MV::scalar_type Scalar;
  typedef MV::local_ordinal_type LO;
  typedef MV::global_ordinal_type GO;

  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef Tpetra::CrsMatrix<Scalar, LO, GO> MAT;

  using std::endl;
  using Teuchos::tuple;
  using Tpetra::global_size_t;

  std::ostream &out              = std::cout;
  RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(rcpFromRef(out));

  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  size_t myRank                       = comm->getRank();

  out << "Amesos2 stand-alone test" << endl
      << endl;

  const size_t numVectors = 1;

  int numGlobalElements                         = 1000;
  RCP<const Tpetra::Map<> > map                 = Tpetra::createUniformContigMap<LO, GO>(numGlobalElements, comm);
  const size_t numMyElements                    = map->getLocalNumElements();
  Teuchos::ArrayView<const GO> myGlobalElements = map->getLocalElementList();

  RCP<MAT> A = Tpetra::createCrsMatrix<Scalar>(map, 3);

  // 1D Laplace
  for (size_t i = 0; i < numMyElements; i++) {
    if (myGlobalElements[i] == 0) {
      A->insertGlobalValues(myGlobalElements[i],
                            tuple<GO>(myGlobalElements[i], myGlobalElements[i] + 1),
                            tuple<Scalar>(2.0, -1.0));
    } else if (myGlobalElements[i] == numGlobalElements - 1) {
      A->insertGlobalValues(myGlobalElements[i],
                            tuple<GO>(myGlobalElements[i] - 1, myGlobalElements[i]),
                            tuple<Scalar>(-1.0, 2.0));
    } else {
      A->insertGlobalValues(myGlobalElements[i],
                            tuple<GO>(myGlobalElements[i] - 1, myGlobalElements[i], myGlobalElements[i] + 1),
                            tuple<Scalar>(-1.0, 2.0, -1.0));
    }
  }
  A->fillComplete();

  /* Create X */
  RCP<MV> X = rcp(new MV(map, numVectors));
  Teuchos::ScalarTraits<Scalar>::seedrandom(846930886);
  X->randomize();

  /* Print X norm */
  {
    Teuchos::Array<ST::magnitudeType> norms(1);
    X->norm2(norms);
    if (myRank == 0)
      *fos << "||X_true|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << endl;
  }

  /* Create B  */
  RCP<MV> B = rcp(new MV(map, numVectors));
  A->apply(*X, *B, Teuchos::NO_TRANS, (Scalar)1.0, (Scalar)0.0);

  /* Reset X */
  X->putScalar((Scalar)0.0);

  // Create solver interface to Superlu through Amesos::Factory
  RCP<Amesos2::Solver<MAT, MV> > solver = Amesos2::create<MAT, MV>("Superlu", A, X, B);

  // Solve
  solver->symbolicFactorization().numericFactorization().solve();

  /* Print X norm */
  {
    Teuchos::Array<ST::magnitudeType> norms(1);
    X->norm2(norms);
    if (myRank == 0)
      *fos << "||X_directSolve|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << endl;
  }

  //   /* Print the solution */
  //   RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(rcpFromRef(out));

  //   *fos << "Solution :" << endl;
  //   X->describe(*fos,Teuchos::VERB_EXTREME);
  //   *fos << endl;

  // We are done.
  return 0;
}
