// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*
 * Structure2D_epetra.cpp
 *
 *  Created on: Oct 24, 2011
 *      Author: wiesner
 */

#include <unistd.h>
#include <iostream>
#include <fstream>
#include <sstream>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

// Epetra
#include <EpetraExt_CrsMatrixIn.h>
#include <EpetraExt_VectorIn.h>
#include <EpetraExt_MultiVectorIn.h>

#include <Epetra_LinearProblem.h>

// AztecOO
#include <AztecOO.h>

// Xpetra
#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>

// MueLu
#include <MueLu.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_ParameterListInterpreter.hpp>

#if defined(HAVE_MUELU_EPETRA)
#include <MueLu_EpetraOperator.hpp>

// prescribe types
// run plain Epetra
typedef double Scalar;
typedef int LocalOrdinal;
typedef int GlobalOrdinal;
typedef Xpetra::EpetraNode Node;
#endif

/*!
 *  2d structural mechanics example for Epetra
 *
 *  (Nearly) Symmetric problem (except of Dirichlet boundaries) solved with AMG solver using a
 *  3 level multigrid with smoothed aggregation transfer operators.
 *
 */

int main(int argc, char* argv[]) {
#if defined(HAVE_MUELU_EPETRA)
#include "MueLu_UseShortNames.hpp"
  using Teuchos::RCP;
  using Teuchos::rcp;

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, &blackhole);

  bool success = false;
  try {
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
    RCP<Teuchos::FancyOStream> out      = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    out->setOutputToRootOnly(0);

#ifndef HAVE_XPETRA_INT_LONG_LONG
    *out << "Warning: scaling test was not compiled with long long int support" << std::endl;
#endif

    // =========================================================================
    // Parameters initialization
    // =========================================================================
    Teuchos::CommandLineProcessor clp(false);

    std::string xmlFileName = "xml/muelu_ParameterList.xml";
    clp.setOption("xml", &xmlFileName, "read parameters from a file [default = 'xml/muelu_ParameterList.xml']");

    int globalNumDofs = 0;  // 7020;
    clp.setOption("globalNumDofs", &globalNumDofs, "global number of degrees of freedom [has to be set by user, default = 0 -> error]");
    int nDofsPerNode = 1;
    clp.setOption("nDofsPerNode", &nDofsPerNode, "number of degrees of freedom per node [has to be set by user, default = 1]");
    int nProcs             = comm->getSize();
    std::string dsolveType = "cg";
    clp.setOption("solver", &dsolveType, "solve type: (none | cg | gmres | standalone) [default = cg]");
    double dtol = 1e-12;
    clp.setOption("tol", &dtol, "solver convergence tolerance [default = 1e-12]");
    std::string problemFile = "stru2d";
    clp.setOption("problem", &problemFile, "string for problem file (e.g. 'stru2d' expects 'stru2d_A.txt', 'stru2d_b.txt' and 'stru2d_ns.txt')");
    std::string coordsFile = "";
    clp.setOption("coordinates", &coordsFile, "file name containing coordinates in matrix market format");

    switch (clp.parse(argc, argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
    }

    if (globalNumDofs == 0) {
      std::cout << "Please specify '--globalNumDofs'! Simulation cannot run without that parameter correctly set" << std::endl;
      return EXIT_FAILURE;
    }

    int nLocalDofs     = (int)globalNumDofs / nProcs;
    nLocalDofs         = nLocalDofs - (nLocalDofs % nDofsPerNode);
    int nCumulatedDofs = 0;
    MueLu_sumAll(comm, nLocalDofs, nCumulatedDofs);

    if (comm->getRank() == nProcs - 1) {
      nLocalDofs += globalNumDofs - nCumulatedDofs;
    }

    // read in problem
    Epetra_Map emap(globalNumDofs, nLocalDofs, 0, *Xpetra::toEpetra(comm));
    Epetra_CrsMatrix* ptrA    = 0;
    Epetra_Vector* ptrf       = 0;
    Epetra_MultiVector* ptrNS = 0;

    std::cout << "Reading matrix market file" << std::endl;

    std::stringstream ssA, ssB, ssNS;
    ssA << problemFile << "_A.txt";
    ssB << problemFile << "_b.txt";
    ssNS << problemFile << "_ns.txt";
    std::string fileA  = ssA.str();
    std::string fileB  = ssB.str();
    std::string fileNS = ssNS.str();
    EpetraExt::MatrixMarketFileToCrsMatrix(fileA.c_str(), emap, emap, emap, ptrA);
    EpetraExt::MatrixMarketFileToVector(fileB.c_str(), emap, ptrf);
    EpetraExt::MatrixMarketFileToMultiVector(fileNS.c_str(), emap, ptrNS);
    RCP<Epetra_CrsMatrix> epA    = Teuchos::rcp(ptrA);
    RCP<Epetra_Vector> epB       = Teuchos::rcp(ptrf);
    RCP<Epetra_MultiVector> epNS = Teuchos::rcp(ptrNS);

    // read in coordinates
    RCP<MultiVector> xCoords = Teuchos::null;
    if (coordsFile != "") {
      Epetra_MultiVector* ptrcoords = 0;
      Epetra_Map coords_emap(globalNumDofs / nDofsPerNode, nLocalDofs / nDofsPerNode, 0, *Xpetra::toEpetra(comm));
      EpetraExt::MatrixMarketFileToMultiVector(coordsFile.c_str(), coords_emap, ptrcoords);
      RCP<Epetra_MultiVector> epCoords = Teuchos::rcp(ptrcoords);
      xCoords                          = Teuchos::rcp(new Xpetra::EpetraMultiVectorT<int, Node>(epCoords));
    }

    // Epetra_CrsMatrix -> Xpetra::Matrix
    RCP<CrsMatrix> exA       = Teuchos::rcp(new Xpetra::EpetraCrsMatrixT<int, Node>(epA));
    RCP<CrsMatrixWrap> crsOp = Teuchos::rcp(new CrsMatrixWrap(exA));
    RCP<Matrix> Op           = Teuchos::rcp_dynamic_cast<Matrix>(crsOp);
    Op->SetFixedBlockSize(nDofsPerNode);

    RCP<MultiVector> xNS = Teuchos::rcp(new Xpetra::EpetraMultiVectorT<int, Node>(epNS));

    // Epetra_Map -> Xpetra::Map
    const RCP<const Map> map = Xpetra::toXpetra<GO, Node>(emap);

    ParameterListInterpreter mueLuFactory(xmlFileName, *comm);
    RCP<Hierarchy> H         = mueLuFactory.CreateHierarchy();
    RCP<MueLu::Level> Finest = H->GetLevel(0);
    Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);
    Finest->Set("A", Op);
    Finest->Set("Nullspace", xNS);
    if (xCoords != Teuchos::null) Finest->Set("Coordinates", xCoords);

    mueLuFactory.SetupHierarchy(*H);

#ifdef HAVE_MUELU_AZTECOO

    H->IsPreconditioner(true);
    MueLu::EpetraOperator mueluPrec(H);  // Wrap MueLu preconditioner into an Epetra Operator

    // create a solution vector
    RCP<Epetra_Vector> epX = rcp(new Epetra_Vector(epA->RowMap()));
    epX->PutScalar((Scalar)0.0);

    Epetra_LinearProblem eProblem(epA.get(), epX.get(), epB.get());

    // AMG as preconditioner within AztecOO
    AztecOO solver(eProblem);
    solver.SetPrecOperator(&mueluPrec);
    if (dsolveType == "cg")
      solver.SetAztecOption(AZ_solver, AZ_cg);
    else if (dsolveType == "gmres")
      solver.SetAztecOption(AZ_solver, AZ_gmres);
    else {  // use fix point method instead
      solver.SetAztecOption(AZ_solver, AZ_fixed_pt);
    }
    solver.SetAztecOption(AZ_output, 1);

    solver.Iterate(500, dtol);

    {  // TODO: simplify this
      RCP<Vector> mueluX = rcp(new Xpetra::EpetraVectorT<int, Node>(epX));
      RCP<Vector> mueluB = rcp(new Xpetra::EpetraVectorT<int, Node>(epB));
      // Print relative residual norm
      Teuchos::ScalarTraits<SC>::magnitudeType residualNorms = Utilities::ResidualNorm(*Op, *mueluX, *mueluB)[0];
      if (comm->getRank() == 0)
        std::cout << "||Residual|| = " << residualNorms << std::endl;
    }
#endif  // HAVE_MUELU_AZTECOO

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
#else
  return EXIT_SUCCESS;
#endif  // #ifdef defined(HAVE_MUELU_EPETRA) and defined(HAVE_MUELU_SERIAL)
}
