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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
/*
 * PgPtest.cpp
 *
 *  Created on: 24.09.2011
 *      Author: tobias
 */

#include <iostream>
#include <unistd.h>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>

#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_CrsOperator.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_Parameters.hpp>

#include "MueLu_Hierarchy.hpp"
#include "MueLu_UCAggregationFactory.hpp"
#include "MueLu_PgPFactory.hpp"
#include "MueLu_GenericRFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_IfpackSmoother.hpp"
#include "MueLu_Ifpack2Smoother.hpp"
#include "MueLu_SmootherFactory.hpp"

#include "MueLu_AmesosSmoother.hpp"
#include "MueLu_Amesos2Smoother.hpp"
#include "MueLu_Utilities.hpp"

#include "MueLu_NoFactory.hpp"

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraMatrixFactory.hpp>

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"
/**********************************************************************************/

// Belos
//#ifdef HAVE_MUELU_BELOS
//#include "BelosConfigDefs.hpp"
//#include "BelosLinearProblem.hpp"
//#include "BelosBlockCGSolMgr.hpp"
//#include "BelosMueLuAdapter.hpp" // this header defines Belos::MueLuOp()
//#endif

int main(int argc, char *argv[]) {
  using Teuchos::RCP;
  using Teuchos::rcp;

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  /**********************************************************************************/
  /* SET TEST PARAMETERS                                                            */
  /**********************************************************************************/
  // Note: use --help to list available options.
  Teuchos::CommandLineProcessor clp(false);

  // Default is Laplace1D with nx = 8748.
  // It's a nice size for 1D and perfect aggregation. (6561=3^8)
  //Nice size for 1D and perfect aggregation on small numbers of processors. (8748=4*3^7)
  Galeri::Xpetra::Parameters<GO> matrixParameters(clp, 8748); // manage parameters of the test case
  Xpetra::Parameters xpetraParameters(clp);             // manage parameters of xpetra

  // custom parameters
  LO maxLevels = 5;
  LO its=10;
  std::string coarseSolver="ifpack2";
  int pauseForDebugger=0;
  clp.setOption("maxLevels",&maxLevels,"maximum number of levels allowed");
  clp.setOption("its",&its,"number of multigrid cycles");
  clp.setOption("coarseSolver",&coarseSolver,"amesos2 or ifpack2 (Tpetra specific. Ignored for Epetra)");
  clp.setOption("debug",&pauseForDebugger,"pause to attach debugger");

  switch (clp.parse(argc,argv)) {
  case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
  case Teuchos::CommandLineProcessor::PARSE_ERROR:
  case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
  case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
  }

  matrixParameters.check();
  xpetraParameters.check();
  // TODO: check custom parameters

  if (comm->getRank() == 0) {
    std::cout << xpetraParameters << matrixParameters;
    // TODO: print custom parameters
  }

  if (pauseForDebugger) {
    Utils::PauseForDebugger();
  }

  /**********************************************************************************/
  /* CREATE INITIAL MATRIX                                                          */
  /**********************************************************************************/
  const RCP<const Map> map = MapFactory::Build(xpetraParameters.GetLib(), matrixParameters.GetNumGlobalElements(), 0, comm);
  RCP<Operator> Op = Galeri::Xpetra::CreateCrsMatrix<SC, LO, GO, Map, CrsOperator>(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList()); //TODO: Operator vs. CrsOperator
  /**********************************************************************************/
  /*                                                                                */
  /**********************************************************************************/

  // dump matrix to file
  //std::string fileName = "Amat.mm";
  //Utils::Write(fileName,Op);

  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map,1);
  nullSpace->putScalar( (SC) 1.0);
  Teuchos::Array<ST::magnitudeType> norms(1);
  nullSpace->norm1(norms);
  if (comm->getRank() == 0)
    std::cout << "||NS|| = " << norms[0] << std::endl;

  RCP<MueLu::Hierarchy<SC,LO,GO,NO,LMO> > H = rcp( new Hierarchy() );
  H->setDefaultVerbLevel(Teuchos::VERB_HIGH);
  H->SetMaxCoarseSize(1);
  RCP<MueLu::Level> Finest = H->GetLevel();
  Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);
  Finest->Request("A",MueLu::NoFactory::get());
  Finest->Set("A",Op,MueLu::NoFactory::get());

  RCP<UCAggregationFactory> UCAggFact = rcp(new UCAggregationFactory());
  UCAggFact->SetMinNodesPerAggregate(3);
  UCAggFact->SetMaxNeighAlreadySelected(0);
  UCAggFact->SetOrdering(MueLu::AggOptions::NATURAL);
  UCAggFact->SetPhase3AggCreation(0.5);

  RCP<TentativePFactory> TentPFact = rcp(new TentativePFactory(UCAggFact));

  RCP<PgPFactory>        Pfact = Teuchos::rcp( new PgPFactory(TentPFact) );
  //Pfact->SetMinimizationMode(MueLu::PgPFactory::L2NORM);
  //Pfact->SetDampingFactor(0.);
  RCP<RFactory>          Rfact = rcp( new GenericRFactory(Pfact) );
  //RCP<GenericPRFactory>  PRfact = rcp( new GenericPRFactory(Pfact,Rfact));
  RCP<RAPFactory>        Acfact = rcp( new RAPFactory(Pfact,Rfact) );

  RCP<SmootherPrototype> smooProto;
  Teuchos::ParameterList ifpackList;

  ifpackList.set("relaxation: sweeps", (LO) 1);
  ifpackList.set("relaxation: damping factor", (SC) 1.0);
  /*
  ifpackList.set("type", "Chebyshev");
  ifpackList.set("chebyshev: degree", (int) 1);
  ifpackList.set("chebyshev: max eigenvalue", (double) 2.0);
  ifpackList.set("chebyshev: min eigenvalue", (double) 1.0);
  ifpackList.set("chebyshev: zero starting solution", false);
  */
  if (xpetraParameters.GetLib() == Xpetra::UseEpetra) {
#ifdef HAVE_MUELU_IFPACK
    ifpackList.set("relaxation: type", "symmetric Gauss-Seidel");
    smooProto = rcp( new IfpackSmoother("point relaxation stand-alone",ifpackList) );
#endif
  } else if (xpetraParameters.GetLib() == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_IFPACK2
      ifpackList.set("relaxation: type", "Symmetric Gauss-Seidel");
    smooProto = rcp( new Ifpack2Smoother("RELAXATION",ifpackList) );
#endif
  }
  if (smooProto == Teuchos::null) {
    throw(MueLu::Exceptions::RuntimeError("main: smoother error"));
  }

  RCP<SmootherFactory> SmooFact = rcp( new SmootherFactory(smooProto) );
  Acfact->setVerbLevel(Teuchos::VERB_HIGH);

  Teuchos::ParameterList status;
  status = H->FullPopulate(*Pfact,*Rfact,*Acfact,*SmooFact,0,maxLevels);
  //RCP<MueLu::Level> coarseLevel = H.GetLevel(1);
  //RCP<Operator> P = coarseLevel->template Get< RCP<Operator> >("P");
  //fileName = "Pfinal.mm";
  //Utils::Write(fileName,P);
  if (comm->getRank() == 0) {
    std::cout  << "======================\n Multigrid statistics \n======================" << std::endl;
    status.print(std::cout,Teuchos::ParameterList::PrintOptions().indent(2));
  }

  //FIXME we should be able to just call smoother->SetNIts(50) ... but right now an exception gets thrown

  RCP<SmootherPrototype> coarseProto;
  if (xpetraParameters.GetLib() == Xpetra::UseEpetra) {
#ifdef HAVE_MUELU_AMESOS
    if (comm->getRank() == 0) std::cout << "CoarseGrid: AMESOS" << std::endl;
    Teuchos::ParameterList amesosList;
    amesosList.set("PrintTiming",true);
    coarseProto = rcp( new AmesosSmoother("Amesos_Klu",amesosList) );
    //#elif HAVE_MUELU_IFPACK...
#endif
  } else if (xpetraParameters.GetLib() == Xpetra::UseTpetra) {
    if (coarseSolver=="amesos2") {
#ifdef HAVE_MUELU_AMESOS2
      if (comm->getRank() == 0) std::cout << "CoarseGrid: AMESOS2" << std::endl;
      Teuchos::ParameterList paramList; //unused
      coarseProto = rcp( new Amesos2Smoother("Superlu", paramList) );
#else
      std::cout  << "AMESOS2 not available (try --coarseSolver=ifpack2)" << std::endl;
      return EXIT_FAILURE;
#endif // HAVE_MUELU_AMESOS2
    } else if(coarseSolver=="ifpack2") {
#if defined(HAVE_MUELU_IFPACK2)
        if (comm->getRank() == 0) std::cout << "CoarseGrid: IFPACK2" << std::endl;
        Teuchos::ParameterList ifpack2List;
        ifpack2List.set("fact: ilut level-of-fill",99); // TODO ??
        ifpack2List.set("fact: drop tolerance", 0);
        ifpack2List.set("fact: absolute threshold", 0);
        ifpack2List.set("fact: relative threshold", 0);
        coarseProto = rcp( new Ifpack2Smoother("ILUT",ifpack2List) );
#else
        std::cout  << "IFPACK2 not available (try --coarseSolver=amesos2)" << std::endl;
        return EXIT_FAILURE;
#endif
    } else {
      std::cout  << "Unknow coarse grid solver (try  --coarseSolver=ifpack2 or --coarseSolver=amesos2)" << std::endl;
      return EXIT_FAILURE;
    }

  }
  if (coarseProto == Teuchos::null) {
    throw(MueLu::Exceptions::RuntimeError("main: coarse smoother error"));
  }

  SmootherFactory coarseSolveFact(coarseProto);
  H->SetCoarsestSolver(coarseSolveFact,MueLu::PRE);

  // Define RHS
  RCP<MultiVector> X = MultiVectorFactory::Build(map,1);
  RCP<MultiVector> RHS = MultiVectorFactory::Build(map,1);

  //X->setSeed(846930886);
  //X->randomize();
  X->putScalar(1.0);
  X->norm2(norms);
  if (comm->getRank() == 0)
    std::cout << "||X_true|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

  Op->apply(*X,*RHS,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);

  // Use AMG directly as an iterative method
  {
    X->putScalar( (SC) 0.0);

    H->Iterate(*RHS,its,*X);

    X->norm2(norms);
    if (comm->getRank() == 0)
      std::cout << "||X_" << std::setprecision(2) << its << "|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;
  }

  RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  out->setOutputToRootOnly(0);

  for (int i=0; i<H->GetNumLevels(); i++) {
    RCP<Level> l = H->GetLevel(i);
    *out << std::endl << "Level " << i << std::endl;
    l->print(*out);
  }


  return EXIT_SUCCESS;

}


