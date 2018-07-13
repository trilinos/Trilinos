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
#include <iostream>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>

#include "MueLu_Hierarchy.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_IfpackSmoother.hpp"
#include "MueLu_Ifpack2Smoother.hpp"
#include "MueLu_GenericPRFactory.hpp"

#include "MueLu_MergedSmoother.hpp"
#include "MueLu_AmesosSmoother.hpp"
#include "MueLu_Amesos2Smoother.hpp"
#include "MueLu_Utilities.hpp"

#include <Xpetra_Map.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Parameters.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>

#include "MueLu_UseDefaultTypes.hpp"
#include <unistd.h>
/**********************************************************************************/

// Belos
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosMueLuAdapter.hpp" // this header defines Belos::MueLuOp()

namespace MueLuTests {
#include "MueLu_UseShortNames.hpp"

  RCP<SmootherPrototype> gimmeGaussSeidelProto(Xpetra::UnderlyingLib lib) {

    RCP<SmootherPrototype> smooProto;
    Teuchos::ParameterList ifpackList;
    ifpackList.set("relaxation: sweeps", (LO) 1);
    ifpackList.set("relaxation: damping factor", (SC) 1.0);

    if (lib == Xpetra::UseEpetra) {
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_IFPACK)
      ifpackList.set("relaxation: type", "symmetric Gauss-Seidel");
      smooProto = rcp( new IfpackSmoother("point relaxation stand-alone",ifpackList) );
#endif
    } else if (lib == Xpetra::UseTpetra) {
#if defined(HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_IFPACK2)
      ifpackList.set("relaxation: type", "Symmetric Gauss-Seidel");
      smooProto = rcp( new Ifpack2Smoother("RELAXATION", ifpackList) );
#endif
    }
    if (smooProto == Teuchos::null) {
      throw(MueLu::Exceptions::RuntimeError("gimmeGaussSeidelSmoother: smoother error"));
    }

    return smooProto;
  }

  RCP<SmootherPrototype> gimmeCoarseProto(Xpetra::UnderlyingLib lib, const std::string& coarseSolver, int rank) {
    RCP<SmootherPrototype> coarseProto;
    if (lib == Xpetra::UseEpetra) {
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_AMESOS)
      if (rank == 0) std::cout << "CoarseGrid: AMESOS" << std::endl;
      Teuchos::ParameterList amesosList;
      amesosList.set("PrintTiming",true);
      coarseProto = rcp( new AmesosSmoother("Amesos_Klu",amesosList) );
      //#elif...
#endif
    } else if (lib == Xpetra::UseTpetra) {
      if (coarseSolver=="amesos2") {
#if defined(HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_AMESOS2)
        if (rank == 0) std::cout << "CoarseGrid: AMESOS2" << std::endl;
        Teuchos::ParameterList paramList; //unused
        coarseProto = rcp( new Amesos2Smoother("Superlu", paramList) );
#else
        std::cout  << "AMESOS2 not available (try --coarseSolver=ifpack2)" << std::endl;
        return Teuchos::null; // TODO test for exception //EXIT_FAILURE;
#endif // HAVE_MUELU_TPETRA && HAVE_MUELU_AMESOS2
      } else if(coarseSolver=="ifpack2") {
#if defined(HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_IFPACK2)
        if (rank == 0) std::cout << "CoarseGrid: IFPACK2" << std::endl;
        Teuchos::ParameterList ifpack2List;
        ifpack2List.set("fact: ilut level-of-fill",99); // TODO ??
        ifpack2List.set("fact: drop tolerance", 0);
        ifpack2List.set("fact: absolute threshold", 0);
        ifpack2List.set("fact: relative threshold", 0);
        coarseProto = rcp( new Ifpack2Smoother("ILUT",ifpack2List) );
#else
        std::cout  << "IFPACK2 not available (try --coarseSolver=amesos2)" << std::endl;
        //TODO        TEUCHOS_TEST_FOR_EXCEPTION
        return Teuchos::null;
#endif
      } else {
        std::cout  << "Unknow coarse grid solver (try  --coarseSolver=ifpack2 or --coarseSolver=amesos2)" << std::endl;
        return Teuchos::null;
      }

    }
    if (coarseProto == Teuchos::null) {
      throw(MueLu::Exceptions::RuntimeError("main: coarse smoother error"));
    }

    return coarseProto;
  }


  RCP<SmootherPrototype> gimmeMergedSmoother(int nSmoothers, Xpetra::UnderlyingLib lib, const std::string& coarseSolver, int rank) {
    ArrayRCP<RCP<SmootherPrototype> > smootherList(nSmoothers);

    for (int i=0; i<nSmoothers; i++)
      smootherList[i] = gimmeGaussSeidelProto(lib);

    return rcp (new MergedSmoother(smootherList));
    //verbose mode: return rcp (new MergedSmoother(smootherList, true));
  }

}

int main(int argc, char *argv[]) {
#include "MueLu_UseShortNames.hpp"

  using namespace MueLuTests;

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
  int nSmoothers=2;
  LO maxLevels = 3;
  LO its=10;
  std::string coarseSolver="ifpack2";
  int pauseForDebugger=0;
  clp.setOption("nSmoothers",&nSmoothers,"number of Gauss-Seidel smoothers in the MergedSmootehrs");
  clp.setOption("maxLevels",&maxLevels,"maximum number of levels allowed. If 1, then a MergedSmoother is used on the coarse grid");
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
    matrixParameters.print();
    xpetraParameters.print();
    // TODO: print custom parameters
  }

  if (pauseForDebugger) {
    Utils::PauseForDebugger();
  }

  /**********************************************************************************/
  /* CREATE INITIAL MATRIX                                                          */
  /**********************************************************************************/
  const RCP<const Map> map = MapFactory::Build(xpetraParameters.GetLib(), matrixParameters.GetNumGlobalElements(), 0, comm);
  RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
    Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList()); //TODO: Matrix vs. CrsMatrixWrap
  RCP<Matrix> Op = Pr->BuildMatrix();
  /**********************************************************************************/
  /*                                                                                */
  /**********************************************************************************/

  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map,1);
  nullSpace->putScalar( (SC) 1.0);
  Teuchos::Array<Teuchos::ScalarTraits<SC>::magnitudeType> norms(1);
  nullSpace->norm1(norms);
  if (comm->getRank() == 0)
    std::cout << "||NS|| = " << norms[0] << std::endl;

  RCP<MueLu::Hierarchy<SC,LO,GO,NO> > H = rcp( new Hierarchy() );
  H->setDefaultVerbLevel(Teuchos::VERB_HIGH);
  RCP<MueLu::Level> Finest = rcp( new MueLu::Level() );
  Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);

  Finest->Set("A",Op);
  Finest->Set("Nullspace",nullSpace);
  Finest->Request("Nullspace"); //FIXME putting this in to avoid error until Merge needs business
  //FIXME is implemented

  Finest->Set("NullSpace",nullSpace);
  H->SetLevel(Finest);

  RCP<CoupledAggregationFactory> CoupledAggFact = rcp(new CoupledAggregationFactory());
  RCP<TentativePFactory>    TentPFact = rcp(new TentativePFactory(CoupledAggFact));

  RCP<SaPFactory>       Pfact = rcp( new SaPFactory(TentPFact) );
  RCP<Factory>         Rfact = rcp( new TransPFactory() );
  RCP<GenericPRFactory> PRfact = rcp( new GenericPRFactory(Pfact,Rfact));
  RCP<RAPFactory>       Acfact = rcp( new RAPFactory() );

  RCP<SmootherPrototype> smooProto = gimmeMergedSmoother(nSmoothers, xpetraParameters.GetLib(), coarseSolver, comm->getRank());
  RCP<SmootherFactory> SmooFact = rcp( new SmootherFactory(smooProto) );
  Acfact->setVerbLevel(Teuchos::VERB_HIGH);

  Teuchos::ParameterList status;
  status = H->FullPopulate(PRfact,Acfact,SmooFact,0,maxLevels);
  if (comm->getRank() == 0) {
    std::cout  << "======================\n Multigrid statistics \n======================" << std::endl;
    status.print(std::cout,Teuchos::ParameterList::PrintOptions().indent(2));
  }

  RCP<SmootherPrototype> coarseProto;
  if (maxLevels != 1) {
    coarseProto = gimmeCoarseProto(xpetraParameters.GetLib(), coarseSolver, comm->getRank());
  } else {
    coarseProto = gimmeMergedSmoother(nSmoothers, xpetraParameters.GetLib(), coarseSolver, comm->getRank());
  }
  if (coarseProto == Teuchos::null) return EXIT_FAILURE;
  SmootherFactory coarseSolveFact(coarseProto);
  H->SetCoarsestSolver(coarseSolveFact,MueLu::PRE);

  // Define RHS
  RCP<MultiVector> X = MultiVectorFactory::Build(map,1);
  RCP<MultiVector> RHS = MultiVectorFactory::Build(map,1);

  X->setSeed(846930886);
  X->randomize();
  X->norm2(norms);
  if (comm->getRank() == 0)
    std::cout << "||X_true|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

  Op->apply(*X,*RHS,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);

  // Use AMG directly as an iterative method
  {
    X->putScalar( (SC) 0.0);

    H->PrintResidualHistory(true);
    H->Iterate(*RHS,*X,its);

    X->norm2(norms);
    if (comm->getRank() == 0)
      std::cout << "||X_" << std::setprecision(2) << its << "|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;
  }

  return EXIT_SUCCESS;

}
