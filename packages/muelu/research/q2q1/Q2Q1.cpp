// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <cstdio>
#include <unistd.h>
#include <iostream>

#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_ArrayRCP.hpp>

#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_MapExtractorFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_StridedMapFactory.hpp>
#include <Xpetra_VectorFactory.hpp>

#include <MueLu.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_BaseClass.hpp>

#include <MueLu_BlockedDirectSolver.hpp>
#include <MueLu_BlockedPFactory.hpp>
#include <MueLu_BlockedRAPFactory.hpp>
#include <MueLu_CoalesceDropFactory.hpp>
#include <MueLu_ConstraintFactory.hpp>
#include <MueLu_EminPFactory.hpp>
#include <MueLu_FactoryManager.hpp>
#include <MueLu_FilteredAFactory.hpp>
#include <MueLu_GenericRFactory.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_PatternFactory.hpp>
#include <MueLu_Q2Q1PFactory.hpp>
#include <MueLu_Q2Q1uPFactory.hpp>
#include <MueLu_SmootherFactory.hpp>
#include <MueLu_SubBlockAFactory.hpp>

#include <MueLu_Utilities.hpp>

#include <MueLu_UseDefaultTypes.hpp>
#include "MueLu_SmootherFactory.hpp"
#include <MueLu_Ifpack2Smoother.hpp>
#include "MueLu_TrilinosSmoother.hpp"

#ifdef HAVE_MUELU_BELOS
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosXpetraAdapter.hpp>     // => This header defines Belos::XpetraOp
#include <BelosMueLuAdapter.hpp>      // => This header defines Belos::MueLuOp
#endif

#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "MatrixMarket_Tpetra.hpp"
#ifdef HAVE_TEKO
#endif
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_TpetraVectorSpace.hpp"
#include "Teko_Utilities.hpp"
#include "MueLu_CreateTpetraPreconditioner.hpp"
#ifdef out
#include "Teko_InverseFactory.hpp"
#include "Teko_InverseLibrary.hpp"
#include "Teko_StridedEpetraOperator.hpp"
#include "Teko_EpetraBlockPreconditioner.hpp"
#include "Teko_LSCPreconditionerFactory.hpp"
#include "Teko_InvLSCStrategy.hpp"
#include "Teko_SIMPLEPreconditionerFactory.hpp"
#endif

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_LinearOpWithSolveFactoryBase.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_PreconditionerFactoryHelpers.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_DefaultPreconditioner.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_VectorBase.hpp"
#include <Thyra_LinearOpWithSolveBase.hpp>
#include <Thyra_VectorBase.hpp>
#include <Thyra_SolveSupportTypes.hpp>
#include <Thyra_MueLuPreconditionerFactory.hpp>
#include <Thyra_MueLuTpetraQ2Q1PreconditionerFactory.hpp>
#include "Stratimikos_MueluTpetraHelpers.hpp"
#include <Thyra_Ifpack2PreconditionerFactory.hpp>



int main(int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  using Teuchos::Array;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::null;
  using Teuchos::as;
  using Tpetra::MatrixMarket::Reader;
  using Thyra::tpetraVectorSpace;


  // =========================================================================
  // MPI initialization using Teuchos
  // =========================================================================
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);

  bool success = true;
  bool verbose = true;
  try {
    RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    TEUCHOS_TEST_FOR_EXCEPTION(comm->getSize() > 1, MueLu::Exceptions::RuntimeError,
        "For now, Q2Q1 only works in serial. "
        "We are working on the parallel implementation.");

    // =========================================================================
    // Parameter initialization
    // =========================================================================
    Teuchos::CommandLineProcessor clp(false);

    Xpetra::Parameters xpetraParameters(clp);

    std::string xmlFileName  = "driver.xml";   clp.setOption("xml",      &xmlFileName,   "read parameters from a file [default = 'driver.xml']");
    double      tol          = 1e-12;          clp.setOption("tol",      &tol,           "solver convergence tolerance");
    int         n            = 9;              clp.setOption("n",        &n,             "problem size (1D)");
    int         maxLevels    = 2;              clp.setOption("nlevels",  &maxLevels,     "max num levels");
    std::string type         = "structured";   clp.setOption("type",     &type,          "structured/unstructured");

    switch (clp.parse(argc, argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:          break;
    }

    // Read data from files
    typedef Tpetra::CrsMatrix<SC,LO,GO>           TP_Crs;
    typedef Tpetra::Operator<SC,LO,GO>            TP_Op;
    typedef Tpetra::MultiVector<SC,LO,GO,NO>      TP_Mv;
    typedef Tpetra::Map<LO,GO,NO>                 TP_Map;
    typedef Thyra::TpetraVectorSpace<SC,LO,GO,NO> THTP_Vs;
    RCP<NO>    node =Tpetra::DefaultPlatform::getDefaultPlatform ().getNode ();

    RCP<TP_Op> A11=Reader<TP_Crs>::readSparseFile("Q2Q1_9x9_A.mm", comm,node);
    RCP<TP_Op> A21=Reader<TP_Crs>::readSparseFile("Q2Q1_9x9_B.mm", comm,node);
    RCP<TP_Op> A12=Reader<TP_Crs>::readSparseFile("Q2Q1_9x9_Bt.mm",comm,node);

    RCP<const TP_Map > cmap = A11->getRangeMap();

    RCP<TP_Mv > Vcoords=Reader<TP_Crs>::readDenseFile("VelCoords9x9.mm",comm,node,cmap);
    RCP<TP_Mv > Pcoords=Reader<TP_Crs>::readDenseFile("PresCoords9x9.mm",comm,node,cmap);


    Teuchos::ArrayRCP<const SC> slop =Utils2::ReadMultiVector("p2vMap9x9.mm",
                             Xpetra::toXpetra(A21->getRangeMap()))->getData(0);
    Teuchos::ArrayRCP<LO> p2vMap(n*n);
    for (int i = 0; i < n*n; i++)  p2vMap[i] = (LO) slop[i];

    // Convert matrices to Teko/Thyra operators
   
    RCP<const THTP_Vs > domain11=tpetraVectorSpace<SC>(A11->getDomainMap());
    RCP<const THTP_Vs > domain12=tpetraVectorSpace<SC>(A12->getDomainMap());
    RCP<const THTP_Vs > domain21=tpetraVectorSpace<SC>(A21->getDomainMap());

    RCP<const THTP_Vs > range11 =tpetraVectorSpace<SC>(A11->getRangeMap());
    RCP<const THTP_Vs > range12 =tpetraVectorSpace<SC>(A12->getRangeMap());
    RCP<const THTP_Vs > range21 =tpetraVectorSpace<SC>(A21->getRangeMap());

    Teko::LinearOp thA11 = Thyra::tpetraLinearOp<double>(range11,domain11,A11);
    Teko::LinearOp thA12 = Thyra::tpetraLinearOp<double>(range12,domain12,A12);
    Teko::LinearOp thA21 = Thyra::tpetraLinearOp<double>(range21,domain21,A21);

    // Bang together the parameter list. Right now, all the MueLu details is
    // hardwired in the MueLu-TpetraQ2Q1 classes. We probably want to switch
    // things so that several of these hardwired features could be modified
    // via parameter lists.

    RCP<Teuchos::ParameterList> StratList = Teuchos::rcp(new Teuchos::ParameterList);
    Teuchos::ParameterList & Q2Q1List     = StratList->sublist("Preconditioner Types").sublist("MueLu-TpetraQ2Q1");
    Teuchos::ParameterList & LSolveTypes  = StratList->sublist("Linear Solver Types");
    Teuchos::ParameterList & BelosList    = LSolveTypes.sublist("Belos");
    Teuchos::ParameterList & BelosSolTypes= BelosList.sublist("Solver Types");
    Teuchos::ParameterList & GmresDetails = BelosSolTypes.sublist("Block GMRES");

    StratList->set("Linear Solver Type","Belos");
    StratList->set("Preconditioner Type","MueLu-TpetraQ2Q1");
    BelosList.set("Solver Type", "Block GMRES");
    GmresDetails.set("Maximum Iterations",    20);
    GmresDetails.set("Convergence Tolerance", 1e-12);
    GmresDetails.set("Verbosity",Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
    GmresDetails.set("Output Frequency",      1);
    GmresDetails.set("Output Style",          1);
 
    Q2Q1List.set("Velcoords",   Vcoords);
    Q2Q1List.set("Prescoords",  Pcoords);
    Q2Q1List.set("p2vMap"   ,   p2vMap);
    Q2Q1List.set("A11"   ,   thA11);
    Q2Q1List.set("A12"   ,   thA12);
    Q2Q1List.set("A21"   ,   thA21);

    // Stratimikos vodou
    
    typedef Thyra::PreconditionerFactoryBase<SC>         Base;
    typedef Thyra::Ifpack2PreconditionerFactory<TP_Crs > Impl;
    typedef Thyra::LinearOpWithSolveFactoryBase<SC>      LOWSFB;
    typedef Thyra::LinearOpWithSolveBase<SC>             LOWSB;
    typedef Thyra::MultiVectorBase<SC>                   TH_Mvb;

    Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;  

    Thyra::addMueLuToStratimikosBuilder(linearSolverBuilder);    
    Stratimikos::enableMueLuTpetraQ2Q1<LO,GO,NO>(linearSolverBuilder,"MueLu-TpetraQ2Q1");

    linearSolverBuilder.setParameterList( StratList );
    RCP<const LOWSFB > lowsFactory = Thyra::createLinearSolveStrategy(linearSolverBuilder);
    RCP<LOWSB        > nsA = lowsFactory->createOp();

    // I've hack together a big matrix that does not use strided maps by simply
    // reading the data again. Normally, this would be supplied by Eric Cyr
    // and would be Teko operators.
   
    int NumEle = A12->getRangeMap()->getNodeNumElements()+A21->getRangeMap()->getNodeNumElements();
    RCP<const TP_Map > FullMap= Utils::Map2TpetraMap(*(MapFactory::createUniformContigMap(Xpetra::UseTpetra, NumEle, comm)));

    RCP<TP_Op> BigMatrix =Tpetra::MatrixMarket::Reader<TP_Crs>::readSparseFile(
                "BigA9x9.mm", FullMap,FullMap,FullMap,FullMap, true,true,false);

    const RCP<Thyra::LinearOpBase<SC> > ThBigA = Thyra::createLinearOp(BigMatrix);
    Thyra::initializeOp<SC>(*lowsFactory, ThBigA, nsA.ptr());

    RCP<TP_Mv> TpetX = Tpetra::createVector<SC>(FullMap);
    RCP<TP_Mv> TpetB = Tpetra::createVector<SC>(FullMap);
    TpetX->randomize();
    BigMatrix->apply(*TpetX,*TpetB);
    TpetX->putScalar(0.0);

    RCP<TH_Mvb > Stratx = Thyra::createMultiVector( TpetX);
    RCP<TH_Mvb > Stratb = Thyra::createMultiVector( TpetB);

    Thyra::SolveStatus<SC> solveStatus = Thyra::solve(*nsA, Thyra::NOTRANS, *Stratb, Stratx.ptr());

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
