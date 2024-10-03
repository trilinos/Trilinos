// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Assert.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_FancyOStream.hpp>

#define XPETRA_ENABLED  // TODO

#include "MueLu_MatrixFactory.hpp"
#include "MueLu_MatrixTypes.hpp"

#include "Xpetra_Map.hpp"
#include "Xpetra_CrsMatrix.hpp"

#ifdef HAVE_XPETRA_TPETRA
#include "Xpetra_TpetraMap.hpp"
#include "Xpetra_TpetraCrsMatrix.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsMatrix.hpp"
#endif

#ifdef HAVE_XPETRA_EPETRA
#include "Xpetra_EpetraMap.hpp"
#include "Xpetra_EpetraCrsMatrix.hpp"
#endif

// Define data types
typedef double SC;
typedef int LO;
typedef int GO;

using Teuchos::RCP;

// This is a test to see if the gallery is working with allowed matrices types.

int main(int argc, char* argv[]) {
  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, &blackhole);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  GO nx = 4;
  GO ny = 4;
  GO nz = 4;
  std::string matrixType("Laplace1D");

  Teuchos::ParameterList pl;
  GO numGlobalElements = nx * ny;
  LO indexBase         = 0;

  Teuchos::ParameterList matrixList;
  matrixList.set("nx", nx);
  matrixList.set("ny", ny);
  matrixList.set("nz", nz);

  RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

#ifdef HAVE_XPETRA_TPETRA
  // Tpetra::CrsMatrix
  {
    RCP<const Tpetra::Map<LO, GO> > map   = rcp(new Tpetra::Map<LO, GO>(numGlobalElements, indexBase, comm));
    RCP<Tpetra::CrsMatrix<SC, LO, GO> > A = MueLu::Gallery::CreateCrsMatrix<SC, LO, GO, Tpetra::Map<LO, GO>, Tpetra::CrsMatrix<SC, LO, GO> >(matrixType, map, matrixList);
    A->describe(*out, Teuchos::VERB_EXTREME);
  }
#endif

#ifdef HAVE_XPETRA_TPETRA
  // Xpetra::TpetraCrsMatrix
  {
    RCP<const Xpetra::TpetraMap<LO, GO> > map   = rcp(new Xpetra::TpetraMap<LO, GO>(numGlobalElements, indexBase, comm));
    RCP<Xpetra::TpetraCrsMatrix<SC, LO, GO> > A = MueLu::Gallery::CreateCrsMatrix<SC, LO, GO, Xpetra::TpetraMap<LO, GO>, Xpetra::TpetraCrsMatrix<SC, LO, GO> >(matrixType, map, matrixList);
    A->describe(*out, Teuchos::VERB_EXTREME);
  }
#endif

#ifdef HAVE_XPETRA_EPETRA
  // Xpetra::EpetraCrsMatrix
  {
    RCP<const Xpetra::EpetraMap> map = rcp(new Xpetra::EpetraMap(numGlobalElements, indexBase, comm));
    RCP<Xpetra::EpetraCrsMatrix> A   = MueLu::Gallery::CreateCrsMatrix<SC, LO, GO, Xpetra::EpetraMap, Xpetra::EpetraCrsMatrix>(matrixType, map, matrixList);
    A->describe(*out, Teuchos::VERB_EXTREME);
  }
#endif

#ifdef HAVE_XPETRA_TPETRA
  // Xpetra::CrsMatrix (Tpetra)
  {
    RCP<const Xpetra::Map<LO, GO> > map   = rcp(new Xpetra::TpetraMap<LO, GO>(numGlobalElements, indexBase, comm));
    RCP<Xpetra::CrsMatrix<SC, LO, GO> > A = MueLu::Gallery::CreateCrsMatrix<SC, LO, GO, Xpetra::Map<LO, GO>, Xpetra::CrsMatrix<SC, LO, GO> >(matrixType, map, matrixList);
    A->describe(*out, Teuchos::VERB_EXTREME);
  }
#endif

#ifdef HAVE_XPETRA_EPETRA
  // Xpetra::CrsMatrix (Epetra)
  {
    RCP<const Xpetra::Map<LO, GO> > map   = rcp(new Xpetra::EpetraMap(numGlobalElements, indexBase, comm));
    RCP<Xpetra::CrsMatrix<SC, LO, GO> > A = MueLu::Gallery::CreateCrsMatrix<SC, LO, GO, Xpetra::Map<LO, GO>, Xpetra::CrsMatrix<SC, LO, GO> >(matrixType, map, matrixList);
    A->describe(*out, Teuchos::VERB_EXTREME);
  }
#endif

#ifdef HAVE_XPETRA_TPETRA
  // Xpetra::Matrix (Tpetra)
  {
    RCP<const Xpetra::Map<LO, GO> > map = rcp(new Xpetra::TpetraMap<LO, GO>(numGlobalElements, indexBase, comm));
    RCP<Xpetra::Matrix<SC, LO, GO> > A  = MueLu::Gallery::CreateCrsMatrix<SC, LO, GO, Xpetra::Map<LO, GO>, Xpetra::Matrix<SC, LO, GO> >(matrixType, map, matrixList);
    A->describe(*out, Teuchos::VERB_EXTREME);
  }
#endif

#ifdef HAVE_XPETRA_EPETRA
  // Xpetra::Matrix (Epetra)
  {
    RCP<const Xpetra::Map<LO, GO> > map = rcp(new Xpetra::EpetraMap(numGlobalElements, indexBase, comm));
    RCP<Xpetra::Matrix<SC, LO, GO> > A  = MueLu::Gallery::CreateCrsMatrix<SC, LO, GO, Xpetra::Map<LO, GO>, Xpetra::Matrix<SC, LO, GO> >(matrixType, map, matrixList);
    A->describe(*out, Teuchos::VERB_EXTREME);
  }
#endif

  return EXIT_SUCCESS;
}
