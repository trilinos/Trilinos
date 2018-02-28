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
#include <stdexcept>

// Teuchos
#include <Teuchos_StandardCatchMacros.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>
#include <Galeri_XpetraMaps.hpp>

#include <Kokkos_DefaultNode.hpp> // For Epetra only runs this points to FakeKokkos in Xpetra

#include "Xpetra_ConfigDefs.hpp"
#include <Xpetra_Map.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_MatrixMatrix.hpp>

#include <MueLu.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_MLParameterListInterpreter.hpp>
#include <MueLu_CreateXpetraPreconditioner.hpp>

// Belos
#ifdef HAVE_MUELU_BELOS
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosSolverFactory.hpp"
#include "BelosXpetraAdapter.hpp"
#include "BelosMueLuAdapter.hpp"
#include "BelosPseudoBlockCGSolMgr.hpp"
#endif

using Teuchos::RCP;

//----------------------------------------------------------------------------------------------------------
//
// This example demonstrates how to use MueLu in a fashion that looks like ML's LevelWrap
//
// In this example, we suppose that the user provides a fine matrix A and P & R operators.  These are
// given to MueLu which then forms the next finest matrix and makes a hierarchy from there.
//
//----------------------------------------------------------------------------------------------------------

const std::string thickSeparator = "==========================================================================================================================";
const std::string thinSeparator  = "--------------------------------------------------------------------------------------------------------------------------";

const std::string prefSeparator = "=====================================";

namespace MueLuExamples {

#ifdef HAVE_MUELU_BELOS
  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void solve_system_hierarchy(Xpetra::UnderlyingLib& lib,
                              Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>>&   A,
                              Teuchos::RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>>&   X,
                              Teuchos::RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>>&   B,
                              Teuchos::RCP<MueLu::Hierarchy<Scalar,LocalOrdinal,GlobalOrdinal,Node>>& H,
                              Teuchos::RCP<Teuchos::ParameterList>& SList) {
#include "MueLu_UseShortNames.hpp"
    using Teuchos::rcp;

    typedef Xpetra::MultiVector <Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;
    typedef Belos::OperatorT<MV>                                         OP;

    // Construct a Belos LinearProblem object
    RCP<OP> belosOp   = rcp(new Belos::XpetraOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>(A));
    RCP<OP> belosPrec = rcp(new Belos::MueLuOp <Scalar,LocalOrdinal,GlobalOrdinal,Node>(H));

    RCP<Belos::LinearProblem<Scalar, MV, OP> > belosProblem =
        rcp(new Belos::LinearProblem<Scalar, MV, OP>(belosOp, X, B));
    belosProblem->setRightPrec(belosPrec);
    belosProblem->setProblem(X,B);
    Belos::SolverFactory<Scalar, MV, OP> BelosFactory;

    // We need to register our manager - or do we want these auto registered?
    // This change happens due to implementing the DII system for Belos with
    // auto registration of all the solver.
    Belos::Impl::registerSolverSubclassForTypes<
      Belos::PseudoBlockCGSolMgr<SC, MV, OP>, SC, MV, OP> ("PSEUDOBLOCK CG");

    RCP<Belos::SolverManager<Scalar, MV, OP> > BelosSolver =
        BelosFactory.create(std::string("CG"), SList);
    BelosSolver->setProblem(belosProblem);
    Belos::ReturnType result = BelosSolver->solve();
    TEUCHOS_TEST_FOR_EXCEPTION(result == Belos::Unconverged, std::runtime_error, "Belos failed to converge");
  }

  // --------------------------------------------------------------------------------------
  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void solve_system_list(Xpetra::UnderlyingLib& lib,
                         Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>>& A,
                         Teuchos::RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>>& X,
                         Teuchos::RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>>& B,
                         Teuchos::ParameterList& MueLuList,
                         Teuchos::RCP<Teuchos::ParameterList>& SList) {
#include "MueLu_UseShortNames.hpp"
    using Teuchos::rcp;

    Teuchos::RCP<MueLu::Hierarchy<Scalar,LocalOrdinal,GlobalOrdinal,Node> > H =
        MueLu::CreateXpetraPreconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node>(
            A, MueLuList, Teuchos::null, Teuchos::null);

    typedef Xpetra::MultiVector <Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;
    typedef Belos::OperatorT<MV>                                         OP;

    // Construct a Belos LinearProblem object
    RCP<OP> belosOp   = rcp(new Belos::XpetraOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>(A));
    RCP<OP> belosPrec = rcp(new Belos::MueLuOp <Scalar,LocalOrdinal,GlobalOrdinal,Node>(H));

    RCP<Belos::LinearProblem<Scalar, MV, OP> > belosProblem =
        rcp(new Belos::LinearProblem<Scalar, MV, OP>(belosOp, X, B));
    belosProblem->setRightPrec(belosPrec);
    belosProblem->setProblem(X,B);

    Belos::SolverFactory<SC, MV, OP> BelosFactory;

    // We need to register our manager - or do we want these auto registered?
    // This change happens due to implementing the DII system for Belos with
    // auto registration of all the solver.
    Belos::Impl::registerSolverSubclassForTypes<
      Belos::PseudoBlockCGSolMgr<SC, MV, OP>, SC, MV, OP> ("PSEUDOBLOCK CG");


    Teuchos::RCP<Belos::SolverManager<SC, MV, OP> > BelosSolver = BelosFactory.create(std::string("CG"), SList);
    BelosSolver->setProblem(belosProblem);
    Belos::ReturnType result = BelosSolver->solve();
    TEUCHOS_TEST_FOR_EXCEPTION(result == Belos::Unconverged, std::runtime_error, "Belos failed to converge");
  }
#endif


  // --------------------------------------------------------------------------------------
  // This routine generate's the user's original A matrix and nullspace
  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void generate_user_matrix_and_nullspace(std::string& matrixType,
                                          Xpetra::UnderlyingLib& lib,
                                          Teuchos::ParameterList& galeriList,
                                          Teuchos::RCP<const Teuchos::Comm<int>>& comm,
                                          Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>>&      A,
                                          Teuchos::RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>>& nullspace) {
#include "MueLu_UseShortNames.hpp"

    RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    Teuchos::FancyOStream& out = *fancy;

    RCP<const Map>   map;
    RCP<MultiVector> coordinates;
    if (matrixType == "Laplace1D") {
      map = Galeri::Xpetra::CreateMap<LO, GO, Node>(lib, "Cartesian1D", comm, galeriList);
      coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,MultiVector>("1D", map, galeriList);

    } else if (matrixType == "Laplace2D" || matrixType == "Star2D" || matrixType == "BigStar2D" || matrixType == "Elasticity2D") {
      map = Galeri::Xpetra::CreateMap<LO, GO, Node>(lib, "Cartesian2D", comm, galeriList);
      coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,MultiVector>("2D", map, galeriList);

    } else if (matrixType == "Laplace3D" || matrixType == "Brick3D" || matrixType == "Elasticity3D") {
      map = Galeri::Xpetra::CreateMap<LO, GO, Node>(lib, "Cartesian3D", comm, galeriList);
      coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,MultiVector>("3D", map, galeriList);
    }

    // Expand map to do multiple DOF per node for block problems
    if (matrixType == "Elasticity2D" || matrixType == "Elasticity3D")
      map = Xpetra::MapFactory<LO,GO,Node>::Build(map, (matrixType == "Elasticity2D" ? 2 : 3));

    out << "Processor subdomains in x direction: " << galeriList.get<GO>("mx") << std::endl
        << "Processor subdomains in y direction: " << galeriList.get<GO>("my") << std::endl
        << "Processor subdomains in z direction: " << galeriList.get<GO>("mz") << std::endl
        << "========================================================" << std::endl;

    RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
        Galeri::Xpetra::BuildProblem<SC,LO,GO,Map,CrsMatrixWrap,MultiVector>(matrixType, map, galeriList);

    A = Pr->BuildMatrix();

    if (matrixType == "Elasticity2D" || matrixType == "Elasticity3D") {
      nullspace = Pr->BuildNullspace();
      A->SetFixedBlockSize((matrixType == "Elasticity2D") ? 2 : 3);
    }
  }
}


// --------------------------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib lib, int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>
  using Teuchos::TimeMonitor;

  bool success = false;
  bool verbose = true;
  try {
#if defined(HAVE_MUELU_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT)
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    Teuchos::FancyOStream& out = *fancy;

    typedef Teuchos::ScalarTraits<SC> STS;

    // =========================================================================
    // Parameters initialization
    // =========================================================================
    GO nx = 100, ny = 100, nz = 100;
    Galeri::Xpetra::Parameters<GO> galeriParameters(clp, nx, ny, nz, "Laplace2D"); // manage parameters of the test case
    Xpetra::Parameters             xpetraParameters(clp);                          // manage parameters of Xpetra

    switch (clp.parse(argc, argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
    }

    Teuchos::ParameterList galeriList = galeriParameters.GetParameterList();
    out << thickSeparator << std::endl << xpetraParameters << galeriParameters;

    // =========================================================================
    // Problem construction
    // =========================================================================
    RCP<const Map>   map;
    RCP<Matrix> A,P,R, Ac;
    RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal,Node> > nullspace;
    std::string matrixType = galeriParameters.GetMatrixType();
    MueLuExamples::generate_user_matrix_and_nullspace(matrixType,lib,galeriList,comm,A,nullspace);
    map=A->getRowMap();

    // =========================================================================
    // Setups and solves
    // =========================================================================
    RCP<Vector> X = VectorFactory::Build(map);
    RCP<Vector> B = VectorFactory::Build(map);
    B->setSeed(846930886);
    B->randomize();
    RCP<TimeMonitor> tm;

#ifdef HAVE_MUELU_BELOS
    // Belos Options
    RCP<Teuchos::ParameterList> SList = rcp(new Teuchos::ParameterList );
    SList->set("Verbosity",Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
    SList->set("Output Frequency",10);
    SList->set("Output Style",Belos::Brief);
    SList->set("Maximum Iterations",200);
    SList->set("Convergence Tolerance",1e-10);
#endif


    // =========================================================================
    // Solve #1 (standard MueLu)
    // =========================================================================
    out << thickSeparator << std::endl;
    out << prefSeparator << " Solve 1: Standard "<< prefSeparator <<std::endl;
    {
      // Use an ML-style parameter list for variety
      Teuchos::ParameterList MLList;
      MLList.set("ML output", 10);
      MLList.set("coarse: type","Amesos-Superlu");
#ifdef HAVE_AMESOS2_KLU2
      MLList.set("coarse: type","Amesos-KLU");
#endif
      MLParameterListInterpreter mueLuFactory(MLList);
      RCP<Hierarchy> H = mueLuFactory.CreateHierarchy();
      Teuchos::RCP<FactoryManagerBase> LevelFactory = mueLuFactory.GetFactoryManager(1);
      H->setlib(lib);
      H->AddNewLevel();
      H->GetLevel(1)->Keep("Nullspace",LevelFactory->GetFactory("Nullspace").get());
      H->GetLevel(0)->Set("A", A);
      mueLuFactory.SetupHierarchy(*H);

#ifdef HAVE_MUELU_BELOS
      // Solve
      MueLuExamples::solve_system_hierarchy<Scalar,LocalOrdinal,GlobalOrdinal,Node>(lib,A,X,B,H,SList);
#endif

      // Extract R,P & Ac for LevelWrap Usage
      H->GetLevel(1)->Get("R",R);
      H->GetLevel(1)->Get("P",P);
      H->GetLevel(1)->Get("A",Ac);

      // extract coarse level null space from level 1 that we have to inject for the next runs...
      nullspace = H->GetLevel(1)->template Get<RCP<MultiVector> >("Nullspace",LevelFactory->GetFactory("Nullspace").get());
    }
    out << thickSeparator << std::endl;

    // =========================================================================
    // Solve #2 (level wrap, the long way, using pre-done Ac)
    // =========================================================================
    out << thickSeparator << std::endl;
    out << prefSeparator << " Solve 2: LevelWrap, Long Way, P, R, Ac "<< prefSeparator <<std::endl;
    {
      // Start w/ an ML-style parameter list
      Teuchos::ParameterList MLList;
      MLList.set("ML output", 10);
      MLList.set("coarse: type","Amesos-Superlu");
#ifdef HAVE_AMESOS2_KLU2
      MLList.set("coarse: type","Amesos-KLU");
#endif
      FactoryManager M1;
      M1.SetFactory("A",        MueLu::NoFactory::getRCP());
      M1.SetFactory("P",        MueLu::NoFactory::getRCP());
      M1.SetFactory("R",        MueLu::NoFactory::getRCP());

      MLParameterListInterpreter mueLuFactory(MLList);
      mueLuFactory.AddFactoryManager(1, 1, Teuchos::rcpFromRef(M1));
      RCP<Hierarchy> H = mueLuFactory.CreateHierarchy();
      H->setlib(lib);
      H->GetLevel(0)->Set("A", A);
      H->AddNewLevel();
      H->GetLevel(1)->Set("R", R);
      H->GetLevel(1)->Set("P", P);
      H->GetLevel(1)->Set("A", Ac);
      H->GetLevel(1)->Set("Nullspace", nullspace);
      mueLuFactory.SetupHierarchy(*H);

#ifdef HAVE_MUELU_BELOS
      MueLuExamples::solve_system_hierarchy(lib,A,X,B,H,SList);
#endif

    }
    out << thickSeparator << std::endl;


    // =========================================================================
    // Solve #3 (level wrap, the long way, using P, R and nullspace)
    // =========================================================================
    out << thickSeparator << std::endl;
    out << prefSeparator << " Solve 3: LevelWrap, Long Way, P, R "<< prefSeparator <<std::endl;
    {

      // Start w/ an ML-style parameter list
      Teuchos::ParameterList MLList;
      MLList.set("ML output", 10);
      MLList.set("coarse: type","Amesos-Superlu");
#ifdef HAVE_AMESOS2_KLU2
      MLList.set("coarse: type","Amesos-KLU");
#endif
      FactoryManager M1;
      M1.SetFactory("P",        MueLu::NoFactory::getRCP());
      M1.SetFactory("R",        MueLu::NoFactory::getRCP());

      MLParameterListInterpreter mueLuFactory(MLList);
      mueLuFactory.AddFactoryManager(1, 1, Teuchos::rcpFromRef(M1));
      RCP<Hierarchy> H = mueLuFactory.CreateHierarchy();
      H->setlib(lib);
      H->GetLevel(0)->Set("A", A);
      H->AddNewLevel();
      H->GetLevel(1)->Set("R", R);
      H->GetLevel(1)->Set("P", P);
      H->GetLevel(1)->Set("Nullspace", nullspace);
      mueLuFactory.SetupHierarchy(*H);
#ifdef HAVE_MUELU_BELOS
      MueLuExamples::solve_system_hierarchy(lib,A,X,B,H,SList);
#endif

    }
    out << thickSeparator << std::endl;

    // =========================================================================
    // Solve #4 (level wrap, the fast way, everything)
    // =========================================================================
    out << thickSeparator << std::endl;
    out << prefSeparator << " Solve 4: LevelWrap, Fast Way, P, R, Ac "<< prefSeparator <<std::endl;
    {
      Teuchos::ParameterList MueLuList, level1;
      level1.set("A",Ac);
      level1.set("R",R);
      level1.set("P",P);
      level1.set("Nullspace",nullspace);
      MueLuList.set("level 1",level1);
      MueLuList.set("verbosity","high");
      MueLuList.set("coarse: max size",100);
      MueLuList.set("max levels",4);
#ifdef HAVE_MUELU_BELOS
      MueLuExamples::solve_system_list(lib,A,X,B,MueLuList,SList);
#endif
    }

    // =========================================================================
    // Solve #5 (level wrap, the fast way, P, R + Nullspace)
    // =========================================================================
    out << thickSeparator << std::endl;
    out << prefSeparator << " Solve 5: LevelWrap, Fast Way, P, R "<< prefSeparator <<std::endl;
    {
      Teuchos::ParameterList MueLuList, level1;
      level1.set("R",R);
      level1.set("P",P);
      level1.set("Nullspace",nullspace);
      MueLuList.set("level 1",level1);
      MueLuList.set("verbosity","high");
      MueLuList.set("coarse: max size",100);
#ifdef HAVE_MUELU_BELOS
      MueLuExamples::solve_system_list(lib,A,X,B,MueLuList,SList);
#endif
    }


    // =========================================================================
    // Solve #6 (level wrap, the fast way, P only, explicit transpose)
    // =========================================================================
    out << thickSeparator << std::endl;
    out << prefSeparator << " Solve 6: LevelWrap, Fast Way, P only, explicit transpose "<< prefSeparator <<std::endl;
    {
      Teuchos::ParameterList MueLuList, level1;
      level1.set("P",P);
      level1.set("Nullspace",nullspace);
      MueLuList.set("level 1",level1);
      MueLuList.set("verbosity","high");
      MueLuList.set("coarse: max size",100);
      MueLuList.set("transpose: use implicit",false);
      MueLuList.set("max levels",4);
#ifdef HAVE_MUELU_BELOS
      MueLuExamples::solve_system_list(lib,A,X,B,MueLuList,SList);
#endif
    }


    // =========================================================================
    // Solve #7 (level wrap, the fast way, P only, implicit transpose)
    // =========================================================================
    out << thickSeparator << std::endl;
    out << prefSeparator << " Solve 7: LevelWrap, Fast Way, P only, implicit transpose "<< prefSeparator <<std::endl;
    {
      Teuchos::ParameterList MueLuList, level1;
      level1.set("P",P);
      level1.set("Nullspace",nullspace);
      MueLuList.set("level 1",level1);
      MueLuList.set("verbosity","high");
      MueLuList.set("coarse: max size",100);
      MueLuList.set("transpose: use implicit",true);
      MueLuList.set("max levels",2);
#ifdef HAVE_MUELU_BELOS
      MueLuExamples::solve_system_list(lib,A,X,B,MueLuList,SList);
#endif
    }
#endif
    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}



//- -- --------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char *argv[]) {
  return Automatic_Test_ETI(argc,argv);
}


