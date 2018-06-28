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
/*
 * PgPFactory.cpp
 *
 *  Created on: 13.10.2011
 *      Author: tobias
 */

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>

#include <MueLu_TestHelpers.hpp>
#include <MueLu_Version.hpp>

#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_MatrixMatrix.hpp>

#include <MueLu_TentativePFactory.hpp>
#include <MueLu_PgPFactory.hpp>
#include <MueLu_GenericRFactory.hpp>
#include <MueLu_TransPFactory.hpp>
#include <MueLu_TrilinosSmoother.hpp>
#include <MueLu_CoupledAggregationFactory.hpp>
#include <MueLu_RAPFactory.hpp>
#include <MueLu_SmootherFactory.hpp>
#include <MueLu_Utilities.hpp>

namespace MueLuTests {

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PgPFactory, Test0, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    out << "version: " << MueLu::Version() << std::endl;

    RCP<PgPFactory> pgpFactory = rcp(new PgPFactory);
    TEST_EQUALITY(pgpFactory != Teuchos::null, true);

    out << *pgpFactory << std::endl;
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PgPFactory, nonsymExample, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   if defined(IFPACK) && defined(IFPACK2)
#   include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    out << "Skipping test because not all required packages are enabled (Ifpack, Ifpack2)." << std::endl;
    return;
    out << "version: " << MueLu::Version() << std::endl;
    out << "Test PgPFactory within" << std::endl;
    out << "level AMG solver using Petrov Galerkin smoothed aggregation with" << std::endl;
    out << "one SGS sweep on each multigrid level as pre- and postsmoother" << std::endl;

    typedef typename Teuchos::ScalarTraits<SC>::magnitudeType magnitude_type;

    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    Teuchos::Array<magnitude_type> results(2);

    // used Xpetra lib (for maps and smoothers)
    Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();

    // generate problem
    LocalOrdinal maxLevels = 3;
    LocalOrdinal its=10;
    LocalOrdinal nEle = 63;
    const RCP<const Map> map = MapFactory::Build(lib, nEle, 0, comm);
    Teuchos::ParameterList matrixParameters;
    matrixParameters.set("nx",nEle);

    // create nonsymmetric tridiagonal matrix
    Scalar epsilon = 1e-3;
    RCP<Matrix> Op = Galeri::Xpetra::TriDiag<SC,LocalOrdinal,GlobalOrdinal,Map,CrsMatrixWrap>(map, nEle, 1.0, 1.0-epsilon, epsilon);

    // build nullspace
    RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map,1);
    nullSpace->putScalar( (SC) 1.0);
    Teuchos::Array<magnitude_type> norms(1);
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
    RCP<CoupledAggregationFactory> CoupledAggFact = rcp(new CoupledAggregationFactory());
    CoupledAggFact->SetMinNodesPerAggregate(3);
    CoupledAggFact->SetMaxNeighAlreadySelected(0);
    CoupledAggFact->SetOrdering("natural");
    CoupledAggFact->SetPhase3AggCreation(0.5);

    RCP<TentativePFactory> Ptentfact = rcp(new TentativePFactory());
    RCP<PgPFactory>        Pfact = rcp( new PgPFactory());
    RCP<Factory>      Rfact = rcp( new GenericRFactory() );
    RCP<RAPFactory>        Acfact = rcp( new RAPFactory() );
    H->SetMaxCoarseSize(1);

    // setup smoothers
    Teuchos::ParameterList smootherParamList;
    smootherParamList.set("relaxation: type", "Symmetric Gauss-Seidel");
    smootherParamList.set("relaxation: sweeps", (LocalOrdinal) 1);
    smootherParamList.set("relaxation: damping factor", (SC) 1.0);
    RCP<SmootherPrototype> smooProto = rcp( new TrilinosSmoother("RELAXATION", smootherParamList) );
    RCP<SmootherFactory> SmooFact = rcp( new SmootherFactory(smooProto) );
    Acfact->setVerbLevel(Teuchos::VERB_HIGH);

    RCP<SmootherFactory> coarseSolveFact = rcp(new SmootherFactory(smooProto, Teuchos::null));

    FactoryManager M;
    M.SetKokkosRefactor(false);
    M.SetFactory("P", Pfact);
    M.SetFactory("R", Rfact);
    M.SetFactory("A", Acfact);
    M.SetFactory("Ptent", Ptentfact);
    M.SetFactory("Aggregates", CoupledAggFact);
    M.SetFactory("Smoother", SmooFact);
    M.SetFactory("CoarseSolver", coarseSolveFact);

    H->Setup(M, 0, maxLevels);

    // test some basic multigrid data
    RCP<Level> coarseLevel = H->GetLevel(1);
    coarseLevel->print(out);
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
    TEST_EQUALITY(coarseLevel->IsAvailable("P",Pfact.get()), false);
    TEST_EQUALITY(coarseLevel->IsAvailable("P",Ptentfact.get()), false);
    TEST_EQUALITY(coarseLevel->IsAvailable("PreSmoother",SmooFact.get()), false);
    TEST_EQUALITY(coarseLevel->IsAvailable("PostSmoother",SmooFact.get()), false);
    TEST_EQUALITY(coarseLevel->IsAvailable("R",Rfact.get()), false);
    TEST_EQUALITY(coarseLevel->GetKeepFlag("P",Pfact.get()), 0);
    TEST_EQUALITY(coarseLevel->GetKeepFlag("P",Ptentfact.get()), 0);
    TEST_EQUALITY(coarseLevel->GetKeepFlag("PreSmoother",SmooFact.get()), 0);
    TEST_EQUALITY(coarseLevel->GetKeepFlag("PostSmoother",SmooFact.get()), 0);
    TEST_EQUALITY(coarseLevel->GetKeepFlag("R",Rfact.get()), 0);

    RCP<Matrix> P1 = coarseLevel->Get< RCP<Matrix> >("P");
    RCP<Matrix> R1 = coarseLevel->Get< RCP<Matrix> >("R");
    RCP<Level> coarseLevel2 = H->GetLevel(2);
    coarseLevel2->print(out);
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
    TEST_EQUALITY(coarseLevel->GetKeepFlag("A",MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(coarseLevel->GetKeepFlag("P",MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(coarseLevel->GetKeepFlag("PreSmoother",MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(coarseLevel->GetKeepFlag("PostSmoother",MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(coarseLevel->GetKeepFlag("R",MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(coarseLevel2->IsRequested("P",Pfact.get()), false);
    TEST_EQUALITY(coarseLevel2->IsRequested("P",Ptentfact.get()), false);
    TEST_EQUALITY(coarseLevel2->IsRequested("PreSmoother",SmooFact.get()), false);
    TEST_EQUALITY(coarseLevel2->IsRequested("PostSmoother",SmooFact.get()), false);
    TEST_EQUALITY(coarseLevel2->IsRequested("R",Rfact.get()), false);
    TEST_EQUALITY(coarseLevel2->IsAvailable("P",Pfact.get()), false);
    TEST_EQUALITY(coarseLevel2->IsAvailable("P",Ptentfact.get()), false);
    TEST_EQUALITY(coarseLevel2->IsAvailable("PreSmoother",SmooFact.get()), false);
    TEST_EQUALITY(coarseLevel2->IsAvailable("PostSmoother",SmooFact.get()), false);
    TEST_EQUALITY(coarseLevel2->IsAvailable("R",Rfact.get()), false);
    TEST_EQUALITY(coarseLevel->GetKeepFlag("P",Pfact.get()), 0);
    TEST_EQUALITY(coarseLevel->GetKeepFlag("P",Ptentfact.get()), 0);
    TEST_EQUALITY(coarseLevel->GetKeepFlag("PreSmoother",SmooFact.get()), 0);
    TEST_EQUALITY(coarseLevel->GetKeepFlag("PostSmoother",SmooFact.get()), 0);
    TEST_EQUALITY(coarseLevel->GetKeepFlag("R",Rfact.get()), 0);
    RCP<Matrix> P2 = coarseLevel2->Get< RCP<Matrix> >("P");
    RCP<Matrix> R2 = coarseLevel2->Get< RCP<Matrix> >("R");

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

      H->Iterate(*RHS,*X,its);

      X->norm2(norms);
      if (comm->getRank() == 0)
        out << "||X_" << std::setprecision(2) << its << "|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;
    }
#   else
    out << "Skipping test because some required packages are not enabled (Ifpack, Ifpack2)." << std::endl;
#   endif
  } //nonsymExample


#if 0
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PgPFactory, NonStandardMaps, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
    Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();

    // generate problem
    LocalOrdinal maxLevels = 3;
    //LocalOrdinal its=10;
    GlobalOrdinal nEle = 63;
    GlobalOrdinal nIndexBase = 10;
    const RCP<const Map> map = MapFactory::Build(lib, nEle, nIndexBase, comm);

    RCP<CrsMatrixWrap> mtx = Galeri::Xpetra::MatrixTraits<Map,CrsMatrixWrap>::Build(map, 3);

    LocalOrdinal NumMyElements = map->getNodeNumElements();
    Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = map->getNodeElementList();
    GlobalOrdinal NumGlobalElements = map->getGlobalNumElements();
    assert(NumGlobalElements == nEle);

    GlobalOrdinal NumEntries;
    LocalOrdinal nnz=2;
    std::vector<Scalar> Values(nnz);
    std::vector<GlobalOrdinal> Indices(nnz);

    Scalar a = 2.0;
    Scalar b = -1.9;
    Scalar c = -0.1;

    for (LocalOrdinal i = 0; i < NumMyElements; ++i)
    {
      if (MyGlobalElements[i] == nIndexBase)
      {
        // off-diagonal for first row
        Indices[0] = nIndexBase;
        NumEntries = 1;
        Values[0] = c;
      }
      else if (MyGlobalElements[i] == nIndexBase + NumGlobalElements - 1)
      {
        // off-diagonal for last row
        Indices[0] = nIndexBase + NumGlobalElements - 2;
        NumEntries = 1;
        Values[0] = b;
      }
      else
      {
        // off-diagonal for internal row
        Indices[0] = MyGlobalElements[i] - 1;
        Values[1] = b;
        Indices[1] = MyGlobalElements[i] + 1;
        Values[0] = c;
        NumEntries = 2;
      }

      // put the off-diagonal entries
      // Xpetra wants ArrayViews (sigh)
      Teuchos::ArrayView<Scalar> av(&Values[0],NumEntries);
      Teuchos::ArrayView<GlobalOrdinal> iv(&Indices[0],NumEntries);
      mtx->insertGlobalValues(MyGlobalElements[i], iv, av);

      // Put in the diagonal entry
      mtx->insertGlobalValues(MyGlobalElements[i],
          Teuchos::tuple<GlobalOrdinal>(MyGlobalElements[i]),
          Teuchos::tuple<Scalar>(a) );

    } //for (LocalOrdinal i = 0; i < NumMyElements; ++i)


    mtx->fillComplete(map,map);

    std::cout << map->getIndexBase() << std::endl;

    RCP<Matrix> Op = Teuchos::rcp_dynamic_cast<Matrix>(mtx);

    // build nullspace
    RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map,1);
    nullSpace->putScalar( (SC) 1.0);

    // fill hierarchy
    RCP<Hierarchy> H = rcp( new Hierarchy() );
    H->setDefaultVerbLevel(Teuchos::VERB_HIGH);
    RCP<Level> Finest = H->GetLevel(); // first associate level with hierarchy (for defaultFactoryHandler!)

    Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);
    Finest->Set("A",Op);                      // set fine level matrix
    Finest->Set("Nullspace",nullSpace);       // set null space information for finest level

    // define transfer operators
    RCP<CoupledAggregationFactory> CoupledAggFact = rcp(new CoupledAggregationFactory());
    CoupledAggFact->SetMinNodesPerAggregate(3);
    CoupledAggFact->SetMaxNeighAlreadySelected(0);
    CoupledAggFact->SetOrdering("natural");
    CoupledAggFact->SetPhase3AggCreation(0.5);

    RCP<TentativePFactory> Ptentfact = rcp(new TentativePFactory());
    RCP<PgPFactory> Pfact = rcp(new PgPFactory());
    RCP<Factory>          Rfact = rcp( new GenericRFactory() );
    RCP<RAPFactory>        Acfact = rcp( new RAPFactory() );
    H->SetMaxCoarseSize(1);

    // setup smoothers
    Teuchos::ParameterList smootherParamList;
    smootherParamList.set("relaxation: type", "Symmetric Gauss-Seidel");
    smootherParamList.set("relaxation: sweeps", (LocalOrdinal) 1);
    smootherParamList.set("relaxation: damping factor", (SC) 1.0);
    RCP<SmootherPrototype> smooProto = rcp( new TrilinosSmoother("RELAXATION", smootherParamList) );
    RCP<SmootherFactory> SmooFact = rcp( new SmootherFactory(smooProto) );
    Acfact->setVerbLevel(Teuchos::VERB_HIGH);

    FactoryManager M;
    M.SetKokkosRefactor(false);
    M.SetFactory("P", Pfact);
    M.SetFactory("R", Rfact);
    M.SetFactory("A", Acfact);
    M.SetFactory("Ptent", Ptentfact);
    M.SetFactory("Aggregates", CoupledAggFact);
    M.SetFactory("Smoother", SmooFact);

    H->Setup(M, 0, maxLevels);

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
    TEST_EQUALITY(coarseLevel->IsRequested("PreSmoother",SmooFact.get()), false);
    TEST_EQUALITY(coarseLevel->IsRequested("PostSmoother",SmooFact.get()), false);
    TEST_EQUALITY(coarseLevel->IsRequested("R",Rfact.get()), false);
    TEST_EQUALITY(coarseLevel->IsAvailable("P",Pfact.get()), false);
    TEST_EQUALITY(coarseLevel->IsAvailable("PreSmoother",SmooFact.get()), false);
    TEST_EQUALITY(coarseLevel->IsAvailable("PostSmoother",SmooFact.get()), false);
    TEST_EQUALITY(coarseLevel->IsAvailable("R",Rfact.get()), false);
    TEST_EQUALITY(coarseLevel->GetKeepFlag("P",Pfact.get()), 0);
    TEST_EQUALITY(coarseLevel->GetKeepFlag("PreSmoother",SmooFact.get()), 0);
    TEST_EQUALITY(coarseLevel->GetKeepFlag("PostSmoother",SmooFact.get()), 0);
    TEST_EQUALITY(coarseLevel->GetKeepFlag("R",Rfact.get()), 0);
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
    TEST_EQUALITY(coarseLevel2->IsRequested("R",Rfact.get()), false);
    TEST_EQUALITY(coarseLevel2->IsAvailable("P",Pfact.get()), false);
    TEST_EQUALITY(coarseLevel2->IsAvailable("PreSmoother",SmooFact.get()), false);
    TEST_EQUALITY(coarseLevel2->IsAvailable("PostSmoother",SmooFact.get()), false);
    TEST_EQUALITY(coarseLevel2->IsAvailable("R",Rfact.get()), false);
    TEST_EQUALITY(coarseLevel2->GetKeepFlag("P",Pfact.get()), 0);
    TEST_EQUALITY(coarseLevel2->GetKeepFlag("PreSmoother",SmooFact.get()), 0);
    TEST_EQUALITY(coarseLevel2->GetKeepFlag("PostSmoother",SmooFact.get()), 0);
    TEST_EQUALITY(coarseLevel2->GetKeepFlag("R",Rfact.get()), 0);

  }
#endif

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PgPFactory, MinimizationModes, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    out << "version: " << MueLu::Version() << std::endl;
    out << "Test PgPFactory (minimization modes)" << std::endl;

    RCP<PgPFactory>        Pfact = rcp( new PgPFactory());
    TEST_EQUALITY(Pfact != Teuchos::null, true);
    Pfact->SetMinimizationMode(MueLu::ANORM);
    TEST_EQUALITY(Pfact->GetMinimizationMode(), MueLu::ANORM);
    Pfact->SetMinimizationMode(MueLu::L2NORM);
    TEST_EQUALITY(Pfact->GetMinimizationMode(), MueLu::L2NORM);
    Pfact->SetMinimizationMode(MueLu::DINVANORM);
    TEST_EQUALITY(Pfact->GetMinimizationMode(), MueLu::DINVANORM);

  }

#if 0 // TODO check me
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PgPFactory, ColumnBasedOmegas, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    out << "version: " << MueLu::Version() << std::endl;
    out << "Test PgPFactory (column based omegas)" << std::endl;


    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    // used Xpetra lib (for maps and smoothers)
    Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();

    // generate problem
    LocalOrdinal maxLevels = 3;
    LocalOrdinal its=10;
    LocalOrdinal nEle = 63;
    const RCP<const Map> map = MapFactory::Build(lib, nEle, 0, comm);
    Teuchos::ParameterList matrixParameters;
    matrixParameters.set("nx",nEle);

    // create nonsymmetric tridiagonal matrix
    RCP<Matrix> Op = Galeri::Xpetra::TriDiag<SC,LocalOrdinal,GlobalOrdinal,Map,CrsMatrixWrap>(map, nEle, 2.0, -1.0, -1.0);

    // build nullspace
    RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map,1);
    nullSpace->putScalar( (SC) 1.0);
    Teuchos::Array<Teuchos::ScalarTraits<SC>::magnitudeType> norms(1);
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
    RCP<CoupledAggregationFactory> CoupledAggFact = rcp(new CoupledAggregationFactory());
    CoupledAggFact->SetMinNodesPerAggregate(3);
    CoupledAggFact->SetMaxNeighAlreadySelected(0);
    CoupledAggFact->SetOrdering("natural");
    CoupledAggFact->SetPhase3AggCreation(0.5);

    RCP<TentativePFactory> Ptentfact = rcp(new TentativePFactory(CoupledAggFact));
    RCP<PgPFactory>        Pfact = rcp( new PgPFactory(Ptentfact));
    RCP<Factory>          Rfact = rcp( new GenericRFactory(Pfact) );
    RCP<RAPFactory>        Acfact = rcp( new RAPFactory() );
    H->SetMaxCoarseSize(1);

    // keep column based omegas
    Finest->Keep("ColBasedOmega",Pfact.get());

    // setup smoothers
    Teuchos::ParameterList smootherParamList;
    smootherParamList.set("relaxation: type", "Symmetric Gauss-Seidel");
    smootherParamList.set("relaxation: sweeps", (LocalOrdinal) 1);
    smootherParamList.set("relaxation: damping factor", (SC) 1.0);
    RCP<SmootherPrototype> smooProto = rcp( new TrilinosSmoother("RELAXATION", smootherParamList) );
    RCP<SmootherFactory> SmooFact = rcp( new SmootherFactory(smooProto) );
    Acfact->setVerbLevel(Teuchos::VERB_HIGH);

    Teuchos::ParameterList status;
    status = H->FullPopulate(*Pfact,*Rfact, *Acfact,*SmooFact,0,maxLevels);

    SmootherFactory coarseSolveFact(smooProto);
    H->SetCoarsestSolver(coarseSolveFact,MueLu::PRE);

    // test some basic multigrid data
    RCP<Level> coarseLevel = H->GetLevel(1);
    coarseLevel->print(out);
    TEST_EQUALITY(coarseLevel->IsRequested("A",MueLu::NoFactory::get()), false);

    RCP<Matrix> P1 = coarseLevel->Get< RCP<Matrix> >("P");
    RCP<Matrix> R1 = coarseLevel->Get< RCP<Matrix> >("R");
    RCP<Level> coarseLevel2 = H->GetLevel(2);
    coarseLevel2->print(out);
    TEST_EQUALITY(coarseLevel2->IsRequested("A",MueLu::NoFactory::get()), false);
    RCP<Matrix> P2 = coarseLevel2->Get< RCP<Matrix> >("P");
    RCP<Matrix> R2 = coarseLevel2->Get< RCP<Matrix> >("R");

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

      H->Iterate(*RHS,*X,its);

      X->norm2(norms);
      if (comm->getRank() == 0)
        out << "||X_" << std::setprecision(2) << its << "|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;
    }

    //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
    RCP<Level> l0 = H->GetLevel(0);
    TEST_EQUALITY(l0->IsRequested("A",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l0->IsRequested("P",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l0->IsRequested("PreSmoother",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l0->IsRequested("PostSmoother",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l0->IsRequested("R",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l0->IsAvailable("A",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l0->IsAvailable("P",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l0->IsAvailable("PreSmoother",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l0->IsAvailable("PostSmoother",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l0->IsAvailable("R",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l0->IsAvailable("ColBasedOmega",Pfact.get()), false);
    TEST_EQUALITY(l0->IsRequested("P",Pfact.get()), false);
    TEST_EQUALITY(l0->GetKeepFlag("A",MueLu::NoFactory::get()), MueLu::UserData);
    TEST_EQUALITY(l0->GetKeepFlag("ColBasedOmega",Pfact.get()), MueLu::Keep);
    RCP<Level> l1 = H->GetLevel(1);
    TEST_EQUALITY(l1->IsRequested("A",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l1->IsRequested("P",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l1->IsRequested("PreSmoother",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l1->IsRequested("PostSmoother",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l1->IsRequested("R",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l1->IsAvailable("A",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l1->IsAvailable("P",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l1->IsAvailable("PreSmoother",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l1->IsAvailable("PostSmoother",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l1->IsAvailable("R",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l1->IsAvailable("ColBasedOmega",Pfact.get()), true);
    TEST_EQUALITY(l1->IsRequested("P",Pfact.get()), false);
    TEST_EQUALITY(l1->GetKeepFlag("A",MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(l1->GetKeepFlag("P",MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(l1->GetKeepFlag("R",MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(l1->GetKeepFlag("PreSmoother",MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(l1->GetKeepFlag("PostSmoother",MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(l1->GetKeepFlag("ColBasedOmega",Pfact.get()), MueLu::Keep);
    RCP<Level> l2 = H->GetLevel(2);
    TEST_EQUALITY(l2->IsRequested("A",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l2->IsRequested("P",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l2->IsRequested("PreSmoother",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l2->IsRequested("PostSmoother",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l2->IsRequested("R",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l2->IsAvailable("A",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l2->IsAvailable("P",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l2->IsAvailable("PreSmoother",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l2->IsAvailable("PostSmoother",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l2->IsAvailable("R",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l2->IsAvailable("ColBasedOmega",Pfact.get()), true);
    TEST_EQUALITY(l2->IsRequested("P",Pfact.get()), false);
    TEST_EQUALITY(l2->GetKeepFlag("A",MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(l2->GetKeepFlag("P",MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(l2->GetKeepFlag("R",MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(l2->GetKeepFlag("PreSmoother",MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(l2->GetKeepFlag("ColBasedOmega",Pfact.get()), MueLu::Keep);
  }
#endif

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PgPFactory, ReUseOmegas, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
#   if !defined(HAVE_MUELU_IFPACK2)
    MUELU_TESTING_DO_NOT_TEST(Xpetra::UseTpetra,"Ifpack2");
#   endif
    out << "version: " << MueLu::Version() << std::endl;
    out << "Test PgPFactory (reuse row based omegas for restriction operator)" << std::endl;

    typedef typename Teuchos::ScalarTraits<SC>::magnitudeType magnitude_type;

    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    // used Xpetra lib (for maps and smoothers)
    Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();

    // generate problem
    LocalOrdinal maxLevels = 3;
    LocalOrdinal its       = 10;
    LocalOrdinal nEle      = 63;
    const RCP<const Map> map = MapFactory::Build(lib, nEle, 0, comm);
    Teuchos::ParameterList matrixParameters;
    matrixParameters.set("nx", nEle);

    // create nonsymmetric tridiagonal matrix
    RCP<Matrix> Op = Galeri::Xpetra::TriDiag<SC,LocalOrdinal,GlobalOrdinal,Map,CrsMatrixWrap>(map, nEle, 2.0, -1.0, -1.0);

    GO nx = nEle;
    Teuchos::ParameterList galeriList;
    galeriList.set("nx", nx);
    RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<Scalar,LocalOrdinal,GlobalOrdinal,Map,RealValuedMultiVector>("1D", map, galeriList);

    // build nullspace
    RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map,1);
    nullSpace->putScalar( (SC) 1.0);
    Teuchos::Array<magnitude_type> norms(1);
    nullSpace->norm1(norms);
    if (comm->getRank() == 0)
      out << "||NS|| = " << norms[0] << std::endl;

    // fill hierarchy
    RCP<Hierarchy> H = rcp( new Hierarchy() );
    H->setDefaultVerbLevel(Teuchos::VERB_HIGH);

    RCP<Level> Finest = H->GetLevel();
    Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);
    Finest->Set("A",         Op);              // set fine level matrix
    Finest->Set("Nullspace", nullSpace);       // set null space information for finest level
    Finest->Set("Coordinates", coordinates);   // set coordinates for finest level

    // define transfer operators
    RCP<CoupledAggregationFactory> CoupledAggFact = rcp(new CoupledAggregationFactory());
    CoupledAggFact->SetMinNodesPerAggregate(3);
    CoupledAggFact->SetMaxNeighAlreadySelected(0);
    CoupledAggFact->SetOrdering("natural");
    CoupledAggFact->SetPhase3AggCreation(0.5);

    RCP<TentativePFactory> Ptentfact = rcp(new TentativePFactory());
    RCP<PgPFactory>        Pfact     = rcp(new PgPFactory());
    RCP<Factory>           Rfact     = rcp(new GenericRFactory());
    RCP<RAPFactory>        Acfact    = rcp(new RAPFactory());
    H->SetMaxCoarseSize(1);

    Pfact->ReUseDampingParameters(true);

    // setup smoothers
    Teuchos::ParameterList smootherParamList;
    smootherParamList.set("relaxation: type",           "Symmetric Gauss-Seidel");
    smootherParamList.set("relaxation: sweeps",         (LocalOrdinal) 1);
    smootherParamList.set("relaxation: damping factor", (SC) 1.0);
    RCP<SmootherPrototype> smooProto = rcp(new TrilinosSmoother("RELAXATION", smootherParamList));
    RCP<SmootherFactory>   SmooFact  = rcp(new SmootherFactory(smooProto));
    Acfact->setVerbLevel(Teuchos::VERB_HIGH);

    RCP<SmootherFactory> coarseSolveFact = rcp(new SmootherFactory(smooProto, Teuchos::null));

    FactoryManager M;
    M.SetKokkosRefactor(false);
    M.SetFactory("P",             Pfact);
    M.SetFactory("R",             Rfact);
    M.SetFactory("A",             Acfact);
    M.SetFactory("Ptent",         Ptentfact);
    M.SetFactory("Aggregates",    CoupledAggFact);
    M.SetFactory("Smoother",      SmooFact);
    M.SetFactory("CoarseSolver",  coarseSolveFact);

    H->Setup(M, 0, maxLevels);

    // test some basic multigrid data
    RCP<Level> coarseLevel = H->GetLevel(1);
    coarseLevel->print(out);
    TEST_EQUALITY(coarseLevel->IsRequested("A", MueLu::NoFactory::get()), false);

    RCP<Matrix> P1 = coarseLevel->Get< RCP<Matrix> >("P");
    RCP<Matrix> R1 = coarseLevel->Get< RCP<Matrix> >("R");
    RCP<Level> coarseLevel2 = H->GetLevel(2);
    coarseLevel2->print(out);
    TEST_EQUALITY(coarseLevel2->IsRequested("A",MueLu::NoFactory::get()), false);
    RCP<Matrix> P2 = coarseLevel2->Get< RCP<Matrix> >("P");
    RCP<Matrix> R2 = coarseLevel2->Get< RCP<Matrix> >("R");

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

      H->Iterate(*RHS,*X,its);

      X->norm2(norms);
      if (comm->getRank() == 0)
        out << "||X_" << std::setprecision(2) << its << "|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;
    }

    //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
    RCP<Level> l0 = H->GetLevel(0);
    TEST_EQUALITY(l0->IsRequested("A",            MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l0->IsRequested("P",            MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l0->IsRequested("PreSmoother",  MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l0->IsRequested("PostSmoother", MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l0->IsRequested("R",            MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l0->IsRequested("RowBasedOmega",Pfact.get()),             false);
    TEST_EQUALITY(l0->IsAvailable("A",            MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l0->IsAvailable("P",            MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l0->IsAvailable("PreSmoother",  MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l0->IsAvailable("PostSmoother", MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l0->IsAvailable("RowBasedOmega",Pfact.get()),             false);
    TEST_EQUALITY(l0->IsAvailable("R",            MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l0->IsRequested("P",            Pfact.get()),             false);
    TEST_EQUALITY(l0->GetKeepFlag("A",            MueLu::NoFactory::get()), MueLu::UserData);
    RCP<Level> l1 = H->GetLevel(1);
    TEST_EQUALITY(l1->IsRequested("A",            MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l1->IsRequested("P",            MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l1->IsRequested("PreSmoother",  MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l1->IsRequested("PostSmoother", MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l1->IsRequested("R",            MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l1->IsRequested("RowBasedOmega",Pfact.get()),             false);
    TEST_EQUALITY(l1->IsAvailable("A",            MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l1->IsAvailable("P",            MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l1->IsAvailable("PreSmoother",  MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l1->IsAvailable("PostSmoother", MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l1->IsAvailable("R",            MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l1->IsAvailable("RowBasedOmega",Pfact.get()),             false);
    TEST_EQUALITY(l1->IsRequested("P",            Pfact.get()),             false);
    TEST_EQUALITY(l1->GetKeepFlag("A",            MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(l1->GetKeepFlag("P",            MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(l1->GetKeepFlag("R",            MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(l1->GetKeepFlag("PreSmoother",  MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(l1->GetKeepFlag("PostSmoother", MueLu::NoFactory::get()), MueLu::Final);
    RCP<Level> l2 = H->GetLevel(2);
    TEST_EQUALITY(l2->IsRequested("A",            MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l2->IsRequested("P",            MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l2->IsRequested("PreSmoother",  MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l2->IsRequested("PostSmoother", MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l2->IsRequested("R",            MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l2->IsRequested("RowBasedOmega",Pfact.get()),             false);
    TEST_EQUALITY(l2->IsAvailable("A",            MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l2->IsAvailable("P",            MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l2->IsAvailable("PreSmoother",  MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l2->IsAvailable("PostSmoother", MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l2->IsAvailable("R",            MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l2->IsAvailable("RowBasedOmega",Pfact.get()),             false);
    TEST_EQUALITY(l2->IsRequested("P",            Pfact.get()),             false);
    TEST_EQUALITY(l2->GetKeepFlag("A",            MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(l2->GetKeepFlag("P",            MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(l2->GetKeepFlag("R",            MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(l2->GetKeepFlag("PreSmoother",  MueLu::NoFactory::get()), MueLu::Final);
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PgPFactory, ReUseOmegasTransP, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
    // reuse row based omegas in PgPFactory activated but not used, since TransPFactory is set as restriction factory
    // check if RowBasedOmega is not stored in Level!
#   include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
#   if !defined(HAVE_MUELU_AMESOS2) || !defined(HAVE_MUELU_IFPACK2)
    MUELU_TESTING_DO_NOT_TEST(Xpetra::UseTpetra,"Amesos2, Ifpack2");
#   endif
    out << "version: " << MueLu::Version() << std::endl;
    out << "Test PgPFactory (reuse row based omegas for restriction operator)" << std::endl;

    typedef typename Teuchos::ScalarTraits<SC>::magnitudeType magnitude_type;

    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    // used Xpetra lib (for maps and smoothers)
    Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();

    // generate problem
    LocalOrdinal maxLevels = 3;
    LocalOrdinal its=10;
    LocalOrdinal nEle = 63;
    const RCP<const Map> map = MapFactory::Build(lib, nEle, 0, comm);
    Teuchos::ParameterList matrixParameters;
    matrixParameters.set("nx",nEle);

    // create nonsymmetric tridiagonal matrix
    RCP<Matrix> Op = Galeri::Xpetra::TriDiag<SC,LocalOrdinal,GlobalOrdinal,Map,CrsMatrixWrap>(map, nEle, 2.0, -1.0, -1.0);

    GO nx = nEle;
    Teuchos::ParameterList galeriList;
    galeriList.set("nx", nx);
    RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<Scalar,LocalOrdinal,GlobalOrdinal,Map,RealValuedMultiVector>("1D", map, galeriList);

    // build nullspace
    RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map,1);
    nullSpace->putScalar( (SC) 1.0);
    Teuchos::Array<magnitude_type> norms(1);
    nullSpace->norm1(norms);
    if (comm->getRank() == 0)
      out << "||NS|| = " << norms[0] << std::endl;

    // fill hierarchy
    RCP<Hierarchy> H = rcp( new Hierarchy() );
    H->setDefaultVerbLevel(Teuchos::VERB_HIGH);

    RCP<Level> Finest = H->GetLevel();
    Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);
    Finest->Set("A", Op);                      // set fine level matrix
    Finest->Set("Nullspace", nullSpace);       // set null space information for finest level
    Finest->Set("Coordinates", coordinates);   // set coordinates for finest level

    // define transfer operators
    RCP<CoupledAggregationFactory> CoupledAggFact = rcp(new CoupledAggregationFactory());
    CoupledAggFact->SetMinNodesPerAggregate(3);
    CoupledAggFact->SetMaxNeighAlreadySelected(0);
    CoupledAggFact->SetOrdering("natural");
    CoupledAggFact->SetPhase3AggCreation(0.5);

    RCP<TentativePFactory> Ptentfact = rcp(new TentativePFactory());
    RCP<PgPFactory>        Pfact = rcp( new PgPFactory());
    RCP<Factory>          Rfact = rcp( new TransPFactory() );
    RCP<RAPFactory>        Acfact = rcp( new RAPFactory() );
    H->SetMaxCoarseSize(1);

    Pfact->ReUseDampingParameters(true);

    // setup smoothers
    Teuchos::ParameterList smootherParamList;
    smootherParamList.set("relaxation: type", "Symmetric Gauss-Seidel");
    smootherParamList.set("relaxation: sweeps", (LocalOrdinal) 1);
    smootherParamList.set("relaxation: damping factor", (SC) 1.0);
    RCP<SmootherPrototype> smooProto = rcp( new TrilinosSmoother("RELAXATION", smootherParamList) );
    RCP<SmootherFactory> SmooFact = rcp( new SmootherFactory(smooProto) );
    Acfact->setVerbLevel(Teuchos::VERB_HIGH);

    RCP<SmootherFactory> coarseSolveFact = rcp(new SmootherFactory(smooProto, Teuchos::null));

    FactoryManager M;
    M.SetKokkosRefactor(false);
    M.SetFactory("P", Pfact);
    M.SetFactory("R", Rfact);
    M.SetFactory("A", Acfact);
    M.SetFactory("Ptent", Ptentfact);
    M.SetFactory("Aggregates", CoupledAggFact);
    M.SetFactory("Smoother", SmooFact);
    M.SetFactory("CoarseSolver", coarseSolveFact);

    H->Setup(M, 0, maxLevels);

    // test some basic multigrid data
    RCP<Level> coarseLevel = H->GetLevel(1);
    coarseLevel->print(out);
    TEST_EQUALITY(coarseLevel->IsRequested("A",MueLu::NoFactory::get()), false);

    RCP<Matrix> P1 = coarseLevel->Get< RCP<Matrix> >("P");
    RCP<Matrix> R1 = coarseLevel->Get< RCP<Matrix> >("R");
    RCP<Level> coarseLevel2 = H->GetLevel(2);
    coarseLevel2->print(out);
    TEST_EQUALITY(coarseLevel2->IsRequested("A",MueLu::NoFactory::get()), false);
    RCP<Matrix> P2 = coarseLevel2->Get< RCP<Matrix> >("P");
    RCP<Matrix> R2 = coarseLevel2->Get< RCP<Matrix> >("R");

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

      H->Iterate(*RHS,*X,its);

      X->norm2(norms);
      if (comm->getRank() == 0)
        out << "||X_" << std::setprecision(2) << its << "|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;
    }

    //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
    RCP<Level> l0 = H->GetLevel(0);
    TEST_EQUALITY(l0->IsRequested("A",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l0->IsRequested("P",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l0->IsRequested("PreSmoother",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l0->IsRequested("PostSmoother",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l0->IsRequested("R",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l0->IsRequested("RowBasedOmega",Pfact.get()), false);
    TEST_EQUALITY(l0->IsAvailable("A",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l0->IsAvailable("P",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l0->IsAvailable("PreSmoother",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l0->IsAvailable("PostSmoother",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l0->IsAvailable("RowBasedOmega",Pfact.get()), false);
    TEST_EQUALITY(l0->IsAvailable("R",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l0->IsRequested("P",Pfact.get()), false);
    TEST_EQUALITY(l0->GetKeepFlag("A",MueLu::NoFactory::get()), MueLu::UserData);
    RCP<Level> l1 = H->GetLevel(1);
    TEST_EQUALITY(l1->IsRequested("A",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l1->IsRequested("P",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l1->IsRequested("PreSmoother",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l1->IsRequested("PostSmoother",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l1->IsRequested("R",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l1->IsRequested("RowBasedOmega",Pfact.get()), false);
    TEST_EQUALITY(l1->IsAvailable("A",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l1->IsAvailable("P",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l1->IsAvailable("PreSmoother",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l1->IsAvailable("PostSmoother",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l1->IsAvailable("R",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l1->IsAvailable("RowBasedOmega",Pfact.get()), false);
    TEST_EQUALITY(l1->IsRequested("P",Pfact.get()), false);
    TEST_EQUALITY(l1->GetKeepFlag("A",MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(l1->GetKeepFlag("P",MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(l1->GetKeepFlag("R",MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(l1->GetKeepFlag("PreSmoother",MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(l1->GetKeepFlag("PostSmoother",MueLu::NoFactory::get()), MueLu::Final);
    RCP<Level> l2 = H->GetLevel(2);
    TEST_EQUALITY(l2->IsRequested("A",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l2->IsRequested("P",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l2->IsRequested("PreSmoother",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l2->IsRequested("PostSmoother",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l2->IsRequested("R",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l2->IsRequested("RowBasedOmega",Pfact.get()), false);
    TEST_EQUALITY(l2->IsAvailable("A",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l2->IsAvailable("P",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l2->IsAvailable("PreSmoother",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l2->IsAvailable("PostSmoother",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(l2->IsAvailable("R",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(l2->IsAvailable("RowBasedOmega",Pfact.get()), false);
    TEST_EQUALITY(l2->IsRequested("P",Pfact.get()), false);
    TEST_EQUALITY(l2->GetKeepFlag("A",MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(l2->GetKeepFlag("P",MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(l2->GetKeepFlag("R",MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(l2->GetKeepFlag("PreSmoother",MueLu::NoFactory::get()), MueLu::Final);
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PgPFactory, EpetraVsTpetra, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   if defined(HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT) && defined(HAVE_MUELU_IFPACK) && defined(HAVE_MUELU_IFPACK2)
#   include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
    out << "version: " << MueLu::Version() << std::endl;
    out << "Compare results of Epetra and Tpetra" << std::endl;
    out << "for 3 level AMG solver using Petrov Galerkin smoothed aggregation with" << std::endl;
    out << "one SGS sweep on each multigrid level as pre- and postsmoother" << std::endl;

    MUELU_TESTING_LIMIT_EPETRA_SCOPE_TPETRA_IS_DEFAULT(Scalar,GlobalOrdinal,Node);

    typedef typename Teuchos::ScalarTraits<SC>::magnitudeType magnitude_type;

    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    Teuchos::Array<magnitude_type> results(2);

    // run test only on 1 procs
    // then we can check shape of transfer operators
    // furthermore slightly different results in parallel for Epetra and Tpetra (due to smoothers?)
    if(comm->getSize() == 1)
    {
      Xpetra::UnderlyingLib lib = Xpetra::UseEpetra;

      // run Epetra and Tpetra test
      for (int run = 0; run < 2; run++) //TODO: create a subfunction instead or Tuple of UnderlyingLib
      {
        if (run == 0) lib = Xpetra::UseEpetra;
        else lib = Xpetra::UseTpetra;

        // generate problem
        LocalOrdinal maxLevels = 3;
        LocalOrdinal its=10;
        GlobalOrdinal nEle = 63;
        const RCP<const Map> map = MapFactory::Build(lib, nEle, 0, comm);
        Teuchos::ParameterList matrixParameters;
        matrixParameters.set("nx",nEle);

        RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
          Galeri::Xpetra::BuildProblem<SC, LocalOrdinal, GlobalOrdinal, Map, CrsMatrixWrap, MultiVector>("Laplace1D", map, matrixParameters);
        RCP<Matrix> Op = Pr->BuildMatrix();
        RCP<RealValuedMultiVector> coordinates = Pr->BuildCoords();

        // build nullspace
        RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map,1);
        nullSpace->putScalar( (SC) 1.0);
        Teuchos::Array<magnitude_type> norms(1);
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
        // Finest->Set("Coordinates", coordinates);  // set coordinates for finest level

        // define transfer operators
        RCP<CoupledAggregationFactory> CoupledAggFact = rcp(new CoupledAggregationFactory());
        CoupledAggFact->SetMinNodesPerAggregate(3);
        CoupledAggFact->SetMaxNeighAlreadySelected(0);
        CoupledAggFact->SetOrdering("natural");
        CoupledAggFact->SetPhase3AggCreation(0.5);

        RCP<TentativePFactory> Ptentfact = rcp(new TentativePFactory());
        RCP<PgPFactory>         Pfact = rcp( new PgPFactory());
        RCP<Factory>           Rfact = rcp( new GenericRFactory() );
        RCP<RAPFactory>        Acfact = rcp( new RAPFactory() );
        H->SetMaxCoarseSize(1);

        // setup smoothers
        Teuchos::ParameterList smootherParamList;
        smootherParamList.set("relaxation: type", "Symmetric Gauss-Seidel");
        smootherParamList.set("relaxation: sweeps", (LocalOrdinal) 1);
        smootherParamList.set("relaxation: damping factor", (SC) 1.0);
        RCP<SmootherPrototype> smooProto = rcp( new TrilinosSmoother("RELAXATION", smootherParamList) );
        RCP<SmootherFactory> SmooFact = rcp( new SmootherFactory(smooProto) );
        Acfact->setVerbLevel(Teuchos::VERB_HIGH);

        RCP<SmootherFactory> coarseSolveFact = rcp(new SmootherFactory(smooProto, Teuchos::null));

        FactoryManager M;
        M.SetKokkosRefactor(false);
        M.SetFactory("P", Pfact);
        M.SetFactory("R", Rfact);
        M.SetFactory("A", Acfact);
        M.SetFactory("Ptent", Ptentfact);
        M.SetFactory("Aggregates", CoupledAggFact);
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
        TEST_EQUALITY(coarseLevel->IsRequested("P",Pfact.get()), false);
        TEST_EQUALITY(coarseLevel->GetKeepFlag("A",MueLu::NoFactory::get()), MueLu::Final);
        TEST_EQUALITY(coarseLevel->GetKeepFlag("P",MueLu::NoFactory::get()), MueLu::Final);
        TEST_EQUALITY(coarseLevel->GetKeepFlag("PreSmoother",MueLu::NoFactory::get()), MueLu::Final);
        TEST_EQUALITY(coarseLevel->GetKeepFlag("PostSmoother",MueLu::NoFactory::get()), MueLu::Final);
        TEST_EQUALITY(coarseLevel->GetKeepFlag("R",MueLu::NoFactory::get()), MueLu::Final);
        TEST_EQUALITY(coarseLevel->IsRequested("P",Ptentfact.get()), false);
        TEST_EQUALITY(coarseLevel->IsRequested("PreSmoother",SmooFact.get()), false);
        TEST_EQUALITY(coarseLevel->IsRequested("PostSmoother",SmooFact.get()), false);
        TEST_EQUALITY(coarseLevel->IsRequested("R",Rfact.get()), false);
        TEST_EQUALITY(coarseLevel->IsAvailable("P",Pfact.get()), false);
        TEST_EQUALITY(coarseLevel->IsAvailable("P",Ptentfact.get()), false);
        TEST_EQUALITY(coarseLevel->IsAvailable("PreSmoother",SmooFact.get()), false);
        TEST_EQUALITY(coarseLevel->IsAvailable("PostSmoother",SmooFact.get()), false);
        TEST_EQUALITY(coarseLevel->IsAvailable("R",Rfact.get()), false);
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

        Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > PtentTPtent = Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*P1,true,*P1,false,out);
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

          H->Iterate(*RHS,*X,its);

          X->norm2(norms);
          if (comm->getRank() == 0)
            out << "||X_" << std::setprecision(2) << its << "|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;
          results[run] = norms[0];
        }
      }

      std::cout << results[0] << " VS " << results[1] << std::endl;
      TEST_EQUALITY(results[0] - results[1] < 1e-10, true); // check results of EPETRA vs TPETRA
    } // comm->getSize == 1
#   else
    out << "Skipping test because some required packages are not enabled (Tpetra, Epetra, EpetraExt, Ifpack, Ifpack2)." << std::endl;
#   endif

  } //EpetraVsTpetra


#   define MUELU_ETI_GROUP(Scalar, LO, GO, Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PgPFactory, Test0, Scalar, LO, GO, Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PgPFactory, nonsymExample, Scalar, LO, GO, Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PgPFactory, MinimizationModes, Scalar, LO, GO, Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PgPFactory, ReUseOmegas, Scalar, LO, GO, Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PgPFactory, ReUseOmegasTransP, Scalar, LO, GO, Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PgPFactory, EpetraVsTpetra, Scalar, LO, GO, Node)

# include <MueLu_ETI_4arg.hpp>

}//namespace MueLuTests




