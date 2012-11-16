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
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include <Xpetra_MultiVectorFactory.hpp>

#include "MueLu_SaPFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_UCAggregationFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_Utilities.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace MueLuTests {

  TEUCHOS_UNIT_TEST(SaPFactory, Test0)
  {
    out << "version: " << MueLu::Version() << std::endl;

    RCP<SaPFactory> sapFactory = rcp(new SaPFactory);
    TEST_EQUALITY(sapFactory != Teuchos::null, true);

    out << *sapFactory << std::endl;

  }

  TEUCHOS_UNIT_TEST(SaPFactory, GetSetMethods)
  {
    out << "version: " << MueLu::Version() << std::endl;

    RCP<SaPFactory> sapFactory = rcp(new SaPFactory);
    sapFactory->SetDampingFactor( (Scalar)4/3 );
    TEST_EQUALITY(((Scalar)4/3) == sapFactory->GetDampingFactor(), true);
    sapFactory->SetDiagonalView("roomWithAView");
    TEST_EQUALITY( sapFactory->GetDiagonalView(), "roomWithAView");

  } //GetSetMethods

#if defined(HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_EPETRAEXT) && defined(HAVE_MUELU_IFPACK2)
  TEUCHOS_UNIT_TEST(SaPFactory, EpetraVsTpetra)
  {
    out << "version: " << MueLu::Version() << std::endl;
    out << "Compare results of Epetra and Tpetra" << std::endl;
    out << "for 3 level AMG solver using smoothed aggregation with" << std::endl;
    out << "one SGS sweep on each multigrid level as pre- and postsmoother" << std::endl;

    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    Teuchos::Array<ST::magnitudeType> results(2);

    // run test only on 1 proc
    if(comm->getSize() == 1)
    {
        Xpetra::UnderlyingLib lib = Xpetra::UseEpetra;

        // run Epetra and Tpetra test
        for (int run = 0; run < 2; run++) //TODO: create a subfunction instead or Tuple of UnderlyingLib
        {
            if (run == 0) lib = Xpetra::UseEpetra;
            else lib = Xpetra::UseTpetra;

            // generate problem
            LO maxLevels = 3;
            LO its=10;
            LO nEle = 63;
            const RCP<const Map> map = MapFactory::Build(lib, nEle, 0, comm);
            Teuchos::ParameterList matrixParameters;
            matrixParameters.set("nx",nEle);

            RCP<Matrix> Op = Galeri::Xpetra::CreateCrsMatrix<SC, LO, GO, Map, CrsMatrixWrap>("Laplace1D", map, matrixParameters);

            // build nullspace
            RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map,1);
            nullSpace->putScalar( (SC) 1.0);
            Teuchos::Array<ST::magnitudeType> norms(1);
            nullSpace->norm1(norms);
            if (comm->getRank() == 0)
              out << "||NS|| = " << norms[0] << std::endl;

            // fill hierarchy
            RCP<Hierarchy> H = rcp( new Hierarchy() );
            H->setDefaultVerbLevel(Teuchos::VERB_HIGH);

            RCP<Level> Finest = H->GetLevel();
            Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);
            Finest->Set("A",Op);                      // set fine level matrix
            Finest->Set("Nullspace",nullSpace);       // set null space information for finest level

            // define transfer operators
            RCP<UCAggregationFactory> UCAggFact = rcp(new UCAggregationFactory());
            UCAggFact->SetMinNodesPerAggregate(3);
            UCAggFact->SetMaxNeighAlreadySelected(0);
            UCAggFact->SetOrdering(MueLu::AggOptions::NATURAL);
            UCAggFact->SetPhase3AggCreation(0.5);

            RCP<TentativePFactory> Ptentfact = rcp(new TentativePFactory());
            RCP<SaPFactory>        Pfact = rcp( new SaPFactory());
            RCP<FactoryBase2>      Rfact = rcp( new TransPFactory() );
            RCP<RAPFactory>        Acfact = rcp( new RAPFactory() );
            H->SetMaxCoarseSize(1);

            // setup smoothers
            Teuchos::ParameterList smootherParamList;
            smootherParamList.set("relaxation: type", "Symmetric Gauss-Seidel");
            smootherParamList.set("relaxation: sweeps", (LO) 1);
            smootherParamList.set("relaxation: damping factor", (SC) 1.0);
            RCP<SmootherPrototype> smooProto = rcp( new TrilinosSmoother("RELAXATION", smootherParamList) );
            RCP<SmootherFactory> SmooFact = rcp( new SmootherFactory(smooProto) );
            Acfact->setVerbLevel(Teuchos::VERB_HIGH);

            RCP<SmootherFactory> coarseSolveFact = rcp(new SmootherFactory(smooProto, Teuchos::null));

            FactoryManager M;
            M.SetFactory("P", Pfact);
            M.SetFactory("R", Rfact);
            M.SetFactory("A", Acfact);
            M.SetFactory("Ptent", Ptentfact);
            M.SetFactory("Aggregates", UCAggFact);
            M.SetFactory("Smoother", SmooFact);
            M.SetFactory("CoarseSolver", coarseSolveFact);

            H->Setup(M, 0, maxLevels);

            // test some basic multigrid data
            RCP<Level> coarseLevel = H->GetLevel(1);
            TEST_EQUALITY(coarseLevel->IsRequested("A",MueLu::NoFactory::get()), false);
            TEST_EQUALITY(coarseLevel->IsRequested("P",MueLu::NoFactory::get()), false);
            TEST_EQUALITY(coarseLevel->IsRequested("PreSmoother",MueLu::NoFactory::get()), false);
            TEST_EQUALITY(coarseLevel->IsRequested("PostSmoother",MueLu::NoFactory::get()), false);
            TEST_EQUALITY(coarseLevel->IsRequested("R",MueLu::NoFactory::get()), false);
            TEST_EQUALITY(coarseLevel->IsAvailable("A",MueLu::NoFactory::get()), true);
            TEST_EQUALITY(coarseLevel->IsAvailable("P",MueLu::NoFactory::get()), true);
            TEST_EQUALITY(coarseLevel->IsAvailable("PreSmoother",MueLu::NoFactory::get()), true);
            TEST_EQUALITY(coarseLevel->IsAvailable("PostSmoother",MueLu::NoFactory::get()), true);
            TEST_EQUALITY(coarseLevel->IsAvailable("R",MueLu::NoFactory::get()), true);
            TEST_EQUALITY(coarseLevel->GetKeepFlag("A",MueLu::NoFactory::get()), MueLu::Final);
            TEST_EQUALITY(coarseLevel->GetKeepFlag("P",MueLu::NoFactory::get()), MueLu::Final);
            TEST_EQUALITY(coarseLevel->GetKeepFlag("PreSmoother",MueLu::NoFactory::get()), MueLu::Final);
            TEST_EQUALITY(coarseLevel->GetKeepFlag("PostSmoother",MueLu::NoFactory::get()), MueLu::Final);
            TEST_EQUALITY(coarseLevel->GetKeepFlag("R",MueLu::NoFactory::get()), MueLu::Final);
            TEST_EQUALITY(coarseLevel->IsRequested("P",Pfact.get()), false);
            TEST_EQUALITY(coarseLevel->IsRequested("P",Ptentfact.get()), false);
            TEST_EQUALITY(coarseLevel->IsRequested("PreSmoother",SmooFact.get()), false);
            TEST_EQUALITY(coarseLevel->IsRequested("PostSmoother",SmooFact.get()), false);
            TEST_EQUALITY(coarseLevel->IsRequested("R",Rfact.get()), false);
            TEST_EQUALITY(coarseLevel->IsRequested("A",Acfact.get()), false);
            TEST_EQUALITY(coarseLevel->IsAvailable("P",Pfact.get()), false);
            TEST_EQUALITY(coarseLevel->IsAvailable("P",Ptentfact.get()), false);
            TEST_EQUALITY(coarseLevel->IsAvailable("PreSmoother",SmooFact.get()), false);
            TEST_EQUALITY(coarseLevel->IsAvailable("PostSmoother",SmooFact.get()), false);
            TEST_EQUALITY(coarseLevel->IsAvailable("R",Rfact.get()), false);
            TEST_EQUALITY(coarseLevel->IsAvailable("A",Acfact.get()), false);
            TEST_EQUALITY(coarseLevel->GetKeepFlag("P",Pfact.get()), 0);
            TEST_EQUALITY(coarseLevel->GetKeepFlag("P",Ptentfact.get()), 0);
            TEST_EQUALITY(coarseLevel->GetKeepFlag("PreSmoother",SmooFact.get()), 0);
            TEST_EQUALITY(coarseLevel->GetKeepFlag("PostSmoother",SmooFact.get()), 0);
            TEST_EQUALITY(coarseLevel->GetKeepFlag("R",Rfact.get()), 0);
            TEST_EQUALITY(coarseLevel->GetKeepFlag("A",Acfact.get()), 0);
            RCP<Matrix> P1 = coarseLevel->Get< RCP<Matrix> >("P");
            RCP<Matrix> R1 = coarseLevel->Get< RCP<Matrix> >("R");
            TEST_EQUALITY(P1->getGlobalNumRows(), 63);
            TEST_EQUALITY(P1->getGlobalNumCols(), 21);
            TEST_EQUALITY(R1->getGlobalNumRows(), 21);
            TEST_EQUALITY(R1->getGlobalNumCols(), 63);
            RCP<Level> coarseLevel2 = H->GetLevel(2);
            TEST_EQUALITY(coarseLevel2->IsRequested("A",MueLu::NoFactory::get()), false);
            TEST_EQUALITY(coarseLevel2->IsRequested("P",MueLu::NoFactory::get()), false);
            TEST_EQUALITY(coarseLevel2->IsRequested("R",MueLu::NoFactory::get()), false);
            TEST_EQUALITY(coarseLevel2->IsRequested("PreSmoother",MueLu::NoFactory::get()), false);
            TEST_EQUALITY(coarseLevel2->IsRequested("PostSmoother",MueLu::NoFactory::get()), false);
            TEST_EQUALITY(coarseLevel2->IsAvailable("A",MueLu::NoFactory::get()), true);
            TEST_EQUALITY(coarseLevel2->IsAvailable("P",MueLu::NoFactory::get()), true);
            TEST_EQUALITY(coarseLevel2->IsAvailable("PreSmoother",MueLu::NoFactory::get()), true);
            TEST_EQUALITY(coarseLevel2->IsAvailable("PostSmoother",MueLu::NoFactory::get()), false);
            TEST_EQUALITY(coarseLevel2->IsAvailable("R",MueLu::NoFactory::get()), true);
            TEST_EQUALITY(coarseLevel2->GetKeepFlag("A",MueLu::NoFactory::get()), MueLu::Final);
            TEST_EQUALITY(coarseLevel2->GetKeepFlag("P",MueLu::NoFactory::get()), MueLu::Final);
            TEST_EQUALITY(coarseLevel2->GetKeepFlag("PreSmoother",MueLu::NoFactory::get()), MueLu::Final);
            TEST_EQUALITY(coarseLevel2->GetKeepFlag("PostSmoother",MueLu::NoFactory::get()), 0);
            TEST_EQUALITY(coarseLevel2->GetKeepFlag("R",MueLu::NoFactory::get()), MueLu::Final);
            TEST_EQUALITY(coarseLevel2->IsRequested("P",Pfact.get()), false);
            TEST_EQUALITY(coarseLevel2->IsRequested("P",Ptentfact.get()), false);
            TEST_EQUALITY(coarseLevel2->IsRequested("R",Rfact.get()), false);
            TEST_EQUALITY(coarseLevel2->IsAvailable("P",Pfact.get()), false);
            TEST_EQUALITY(coarseLevel2->IsAvailable("P",Ptentfact.get()), false);
            TEST_EQUALITY(coarseLevel2->IsAvailable("PreSmoother",SmooFact.get()), false);
            TEST_EQUALITY(coarseLevel2->IsAvailable("PostSmoother",SmooFact.get()), false);
            TEST_EQUALITY(coarseLevel2->IsAvailable("R",Rfact.get()), false);
            TEST_EQUALITY(coarseLevel2->GetKeepFlag("P",Pfact.get()), 0);
            TEST_EQUALITY(coarseLevel2->GetKeepFlag("P",Ptentfact.get()), 0);
            TEST_EQUALITY(coarseLevel2->GetKeepFlag("PreSmoother",SmooFact.get()), 0);
            TEST_EQUALITY(coarseLevel2->GetKeepFlag("PostSmoother",SmooFact.get()), 0);
            TEST_EQUALITY(coarseLevel2->GetKeepFlag("R",Rfact.get()), 0);
            RCP<Matrix> P2 = coarseLevel2->Get< RCP<Matrix> >("P");
            RCP<Matrix> R2 = coarseLevel2->Get< RCP<Matrix> >("R");
            TEST_EQUALITY(P2->getGlobalNumRows(), 21);
            TEST_EQUALITY(P2->getGlobalNumCols(), 7);
            TEST_EQUALITY(R2->getGlobalNumRows(), 7);
            TEST_EQUALITY(R2->getGlobalNumCols(), 21);

            Teuchos::RCP<Xpetra::Matrix<Scalar,LO,GO> > PtentTPtent = MueLu::Utils<Scalar,LO,GO>::TwoMatrixMultiply(P1,true,P1,false);
            TEST_EQUALITY(PtentTPtent->getGlobalMaxNumRowEntries()-3<1e-12, true);
            TEST_EQUALITY(P1->getGlobalMaxNumRowEntries()-2<1e-12, true);
            TEST_EQUALITY(P2->getGlobalMaxNumRowEntries()-2<1e-12, true);

            // Define RHS
            RCP<MultiVector> X = MultiVectorFactory::Build(map,1);
            RCP<MultiVector> RHS = MultiVectorFactory::Build(map,1);

            X->putScalar(1.0);
            X->norm2(norms);
            if (comm->getRank() == 0)
              out << "||X_true|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

            Op->apply(*X,*RHS,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);

            // Use AMG directly as an iterative method
            {
              X->putScalar( (SC) 0.0);

              H->Iterate(*RHS,its,*X);

              X->norm2(norms);
              if (comm->getRank() == 0)
                out << "||X_" << std::setprecision(2) << its << "|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;
              results[run] = norms[0];
            }
        }

        TEST_EQUALITY(results[0] - results[1] < 1e-10, true); // check results of EPETRA vs TPETRA
    } // comm->getSize == 1

  } //SaPFactory_EpetraVsTpetra
#endif

}//namespace MueLuTests

