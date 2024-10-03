// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*
 * Navier2D_epetra.cpp
 *
 *  Created on: Mar 26, 2011
 *      Author: wiesner
 */

#include <unistd.h>
#include <iostream>
#include <fstream>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

// Epetra
#include <EpetraExt_CrsMatrixIn.h>
#include <EpetraExt_VectorIn.h>
#include <EpetraExt_MultiVectorIn.h>

// Xpetra
#include <Xpetra_Map.hpp>
#include <Xpetra_EpetraMap.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Parameters.hpp>
#include <Xpetra_MapExtractorFactory.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_StridedMapFactory.hpp>
#include <Xpetra_StridedMap.hpp>

// MueLu
#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Memory.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_UncoupledAggregationFactory.hpp"
#include "MueLu_PgPFactory.hpp"
#include "MueLu_GenericRFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_BlockedRAPFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_PreDropFunctionConstVal.hpp"
#include "MueLu_NullspaceFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_DirectSolver.hpp"
#include "MueLu_EpetraOperator.hpp"
#include "MueLu_SubBlockAFactory.hpp"
#include "MueLu_BlockedPFactory.hpp"
#include "MueLu_BlockedGaussSeidelSmoother.hpp"
#include "MueLu_CoarseMapFactory.hpp"
#include "MueLu_AmalgamationFactory.hpp"
#include <MueLu_ParameterListInterpreter.hpp>

#include <Epetra_LinearProblem.h>
#include <AztecOO.h>

#include "Navier2D_Helpers.h"

/*!
 *  2d Navier Stokes example (for Epetra)
 *
 *  using block matrices
 */

int main(int argc, char* argv[]) {
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
  typedef double Scalar;
  typedef int LocalOrdinal;
  typedef int GlobalOrdinal;
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef Xpetra::EpetraNode Node;
#include "MueLu_UseShortNames.hpp"

  using Teuchos::RCP;
  using Teuchos::rcp;
  using namespace MueLuTests;
  using namespace Teuchos;

  oblackholestream blackhole;
  GlobalMPISession mpiSession(&argc, &argv, &blackhole);

  bool success = false;
  bool verbose = true;
  try {
    // default parameters
    std::string xmlFile = "myXML.xml";

    // Note: use --help to list available options.
    CommandLineProcessor clp(false);
    clp.setOption("xml", &xmlFile, "xml file with solver parameters for a 2x2 blocked NS example");

    switch (clp.parse(argc, argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
    }

    RCP<const Comm<int> > comm = DefaultComm<int>::getComm();
    RCP<FancyOStream> out      = fancyOStream(rcpFromRef(std::cout));
    out->setOutputToRootOnly(0);
    *out << MueLu::MemUtils::PrintMemoryUsage() << std::endl;

    // Timing
    Time myTime("global");
    TimeMonitor MM(myTime);

    GO maxCoarseSize = 1;  // FIXME clp doesn't like long long int

    int globalNumDofs = 1500;  // used for the maps
    int nDofsPerNode  = 3;     // used for generating the fine level null-space

    // build strided maps
    // striding information: 2 velocity dofs and 1 pressure dof = 3 dofs per node
    std::vector<size_t> stridingInfo;
    stridingInfo.push_back(2);
    stridingInfo.push_back(1);

    /////////////////////////////////////// build strided maps
    // build strided maps:
    // xstridedfullmap: full map (velocity and pressure dof gids), continous
    // xstridedvelmap: only velocity dof gid maps (i.e. 0,1,3,4,6,7...)
    // xstridedpremap: only pressure dof gid maps (i.e. 2,5,8,...)
    Xpetra::UnderlyingLib lib             = Xpetra::UseEpetra;
    RCP<const StridedMap> xstridedfullmap = StridedMapFactory::Build(lib, globalNumDofs, 0, stridingInfo, comm, -1);
    RCP<const StridedMap> xstridedvelmap  = StridedMapFactory::Build(xstridedfullmap, 0);
    RCP<const StridedMap> xstridedpremap  = StridedMapFactory::Build(xstridedfullmap, 1);

    /////////////////////////////////////// transform Xpetra::Map objects to Epetra
    // this is needed for AztecOO
    const RCP<const Epetra_Map> fullmap = rcpFromRef(Xpetra::toEpetra(*xstridedfullmap));
    RCP<const Epetra_Map> velmap        = rcpFromRef(Xpetra::toEpetra(*xstridedvelmap));
    RCP<const Epetra_Map> premap        = rcpFromRef(Xpetra::toEpetra(*xstridedpremap));

    /////////////////////////////////////// import problem matrix and RHS from files (-> Epetra)

    // read in problem
    Epetra_CrsMatrix* ptrA    = 0;
    Epetra_Vector* ptrf       = 0;
    Epetra_MultiVector* ptrNS = 0;

    *out << "Reading matrix market file" << std::endl;
    EpetraExt::MatrixMarketFileToCrsMatrix("A_re1000_5932.txt", *fullmap, *fullmap, *fullmap, ptrA);
    EpetraExt::MatrixMarketFileToVector("b_re1000_5932.txt", *fullmap, ptrf);
    RCP<Epetra_CrsMatrix> epA    = rcp(ptrA);
    RCP<Epetra_Vector> epv       = rcp(ptrf);
    RCP<Epetra_MultiVector> epNS = rcp(ptrNS);

    /////////////////////////////////////// split system into 2x2 block system

    *out << "Split matrix into 2x2 block matrix" << std::endl;

    // split fullA into A11,..., A22
    RCP<Epetra_CrsMatrix> A11;
    RCP<Epetra_CrsMatrix> A12;
    RCP<Epetra_CrsMatrix> A21;
    RCP<Epetra_CrsMatrix> A22;

    if (SplitMatrix2x2(epA, *velmap, *premap, A11, A12, A21, A22) == false)
      *out << "Problem with splitting matrix" << std::endl;

    /////////////////////////////////////// transform Epetra objects to Xpetra (needed for MueLu)

    // build Xpetra objects from Epetra_CrsMatrix objects
    RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > xA11 = rcp(new Xpetra::EpetraCrsMatrixT<GO, Node>(A11));
    RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > xA12 = rcp(new Xpetra::EpetraCrsMatrixT<GO, Node>(A12));
    RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > xA21 = rcp(new Xpetra::EpetraCrsMatrixT<GO, Node>(A21));
    RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > xA22 = rcp(new Xpetra::EpetraCrsMatrixT<GO, Node>(A22));

    /////////////////////////////////////// generate MapExtractor object

    std::vector<RCP<const Xpetra::Map<LO, GO, Node> > > xmaps;
    xmaps.push_back(xstridedvelmap);
    xmaps.push_back(xstridedpremap);

    RCP<const Xpetra::MapExtractor<Scalar, LO, GO, Node> > map_extractor = Xpetra::MapExtractorFactory<Scalar, LO, GO, Node>::Build(xstridedfullmap->getMap(), xmaps);

    /////////////////////////////////////// build blocked transfer operator
    // using the map extractor
    RCP<Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node> > bOp = rcp(new Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node>(map_extractor, map_extractor, 10));
    bOp->setMatrix(0, 0, Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(xA11)));
    bOp->setMatrix(0, 1, Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(xA12)));
    bOp->setMatrix(1, 0, Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(xA21)));
    bOp->setMatrix(1, 1, Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(xA22)));

    bOp->fillComplete();

    //////////////////////////////////////// prepare setup
    ParameterListInterpreter mueLuFactory(xmlFile, *comm);

    RCP<Hierarchy> H = mueLuFactory.CreateHierarchy();
    H->setDefaultVerbLevel(VERB_HIGH);
    H->SetMaxCoarseSize(maxCoarseSize);

    RCP<MueLu::Level> Finest = H->GetLevel(0);
    Finest->setDefaultVerbLevel(VERB_HIGH);
    Finest->Set("A", rcp_dynamic_cast<Matrix>(bOp));

    ////////////////////////////////////////// prepare null space for A11
    RCP<MultiVector> nullspace11 = MultiVectorFactory::Build(xstridedvelmap, 2);  // this is a 2D standard null space

    for (int i = 0; i < nDofsPerNode - 1; ++i) {
      ArrayRCP<Scalar> nsValues = nullspace11->getDataNonConst(i);
      int numBlocks             = nsValues.size() / (nDofsPerNode - 1);
      for (int j = 0; j < numBlocks; ++j) {
        nsValues[j * (nDofsPerNode - 1) + i] = 1.0;
      }
    }

    Finest->Set("Nullspace1", nullspace11);

    ////////////////////////////////////////// prepare null space for A22
    RCP<MultiVector> nullspace22 = MultiVectorFactory::Build(xstridedpremap, 1);  // this is a 2D standard null space
    ArrayRCP<Scalar> nsValues22  = nullspace22->getDataNonConst(0);
    for (int j = 0; j < nsValues22.size(); ++j) {
      nsValues22[j] = 1.0;
    }

    Finest->Set("Nullspace2", nullspace22);

    /////////////////////////////////// BEGIN setup

    mueLuFactory.SetupHierarchy(*H);

    ///////////////////////////////////// END setup

    *out << std::endl;

    RCP<MultiVector> xLsg = MultiVectorFactory::Build(xstridedfullmap, 1);

    // Use AMG directly as an iterative method
    {
      xLsg->putScalar((SC)0.0);

      // Epetra_Vector -> Xpetra::Vector
      RCP<Vector> xRhs = rcp(new Xpetra::EpetraVectorT<int, Node>(epv));

      // calculate initial (absolute) residual
      Array<ScalarTraits<SC>::magnitudeType> norms(1);
      xRhs->norm2(norms);
      *out << "||x_0|| = " << norms[0] << std::endl;

      // apply ten multigrid iterations
      H->Iterate(*xRhs, *xLsg, 100);

      // calculate and print residual
      RCP<MultiVector> xTmp = MultiVectorFactory::Build(xstridedfullmap, 1);
      bOp->apply(*xLsg, *xTmp, NO_TRANS, (SC)1.0, (SC)0.0);
      xRhs->update((SC)-1.0, *xTmp, (SC)1.0);
      xRhs->norm2(norms);
      *out << "||r|| = " << norms[0] << std::endl;
    }

    // TODO: don't forget to add Aztec as prerequisite in CMakeLists.txt!
    //
    // Solve Ax = b using AMG as a preconditioner in AztecOO
    //
    {
      RCP<Epetra_Vector> X = rcp(new Epetra_Vector(epv->Map()));
      X->PutScalar(0.0);
      Epetra_LinearProblem epetraProblem(epA.get(), X.get(), epv.get());

      AztecOO aztecSolver(epetraProblem);
      aztecSolver.SetAztecOption(AZ_solver, AZ_gmres);

      MueLu::EpetraOperator aztecPrec(H);
      aztecSolver.SetPrecOperator(&aztecPrec);

      int maxIts = 50;
      double tol = 1e-8;

      aztecSolver.Iterate(maxIts, tol);
    }

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
#else
  std::cout << "Epetra (and/or EpetraExt) are not available. Skip test." << std::endl;
  return EXIT_SUCCESS;
#endif
}
