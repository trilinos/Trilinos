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
 * GenericRFactory.cpp
 *
 *  Created on: 20.09.2011
 *      Author: tobias
 */

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>

#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Version.hpp"
#include "MueLu_Needs.hpp"

#include "MueLu_PgPFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_GenericRFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_CoupledAggregationFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_RAPFactory.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

//TODO: remove this
#ifdef HAVE_MUELU_EPETRAEXT
#include "EpetraExt_RowMatrixOut.h"
#include "Xpetra_EpetraCrsMatrix.hpp"
#endif

namespace MueLuTests {


  //this macro declares the unit-test-class:
  TEUCHOS_UNIT_TEST(GenericRFactory, Constructor)
  {
    //we are now in a class method declared by the above macro, and
    //that method has these input arguments:
    //Teuchos::FancyOStream& out, bool& success
    out << "version: " << MueLu::Version() << std::endl;

    //TEST_EQUALITY(needs != Teuchos::null, true);
  }

  TEUCHOS_UNIT_TEST(GenericRFactory, SymmetricProblem)
  {
    out << "version: " << MueLu::Version() << std::endl;
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    // generate problem
    LO maxLevels = 3;
    LO nEle = 63;
    const RCP<const Map> map = MapFactory::Build(TestHelpers::Parameters::getLib(), nEle, 0, comm);
    Teuchos::ParameterList matrixParameters;
    matrixParameters.set("nx",nEle);

    RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
        Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>("Laplace1D", map, matrixParameters);
    RCP<Matrix> Op = Pr->BuildMatrix();

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
    RCP<CoupledAggregationFactory> CoupledAggFact = rcp(new CoupledAggregationFactory());
    CoupledAggFact->SetMinNodesPerAggregate(3);
    CoupledAggFact->SetMaxNeighAlreadySelected(0);
    CoupledAggFact->SetOrdering(MueLu::AggOptions::NATURAL);
    CoupledAggFact->SetPhase3AggCreation(0.5);

    RCP<SaPFactory>         Pfact = rcp( new SaPFactory());
    RCP<Factory>           Rfact = rcp( new GenericRFactory() );
    H->SetMaxCoarseSize(1);

    // setup smoothers
    Teuchos::ParameterList smootherParamList;
    smootherParamList.set("relaxation: type", "Symmetric Gauss-Seidel");
    smootherParamList.set("relaxation: sweeps", (LO) 1);
    smootherParamList.set("relaxation: damping factor", (SC) 1.0);
    RCP<SmootherPrototype> smooProto = rcp( new TrilinosSmoother("RELAXATION", smootherParamList) );
    RCP<SmootherFactory> SmooFact = rcp( new SmootherFactory(smooProto) );
    //Acfact->setVerbLevel(Teuchos::VERB_HIGH);

    RCP<SmootherFactory> coarseSolveFact = rcp(new SmootherFactory(smooProto, Teuchos::null));

    FactoryManager M;
    M.SetFactory("P", Pfact);
    M.SetFactory("R", Rfact);
    M.SetFactory("Aggregates", CoupledAggFact);
    M.SetFactory("Smoother", SmooFact);
    M.SetFactory("CoarseSolver", coarseSolveFact);

    H->Setup(M, 0, maxLevels);

    RCP<Level> coarseLevel = H->GetLevel(1);
    RCP<Matrix> P1 = coarseLevel->Get< RCP<Matrix> >("P");
    RCP<Matrix> R1 = coarseLevel->Get< RCP<Matrix> >("R");
    RCP<Level> coarseLevel2 = H->GetLevel(2);
    RCP<Matrix> P2 = coarseLevel2->Get< RCP<Matrix> >("P");
    RCP<Matrix> R2 = coarseLevel2->Get< RCP<Matrix> >("R");

    TEST_EQUALITY(Finest->IsAvailable("PreSmoother"), true);
    TEST_EQUALITY(Finest->IsAvailable("PostSmoother"), true);
    TEST_EQUALITY(coarseLevel->IsAvailable("PreSmoother"), true);
    TEST_EQUALITY(coarseLevel->IsAvailable("PostSmoother"), true);
    TEST_EQUALITY(coarseLevel2->IsAvailable("PreSmoother"), true);
    TEST_EQUALITY(coarseLevel2->IsAvailable("PostSmoother"), false);

    // test some basic multgrid data
    TEST_EQUALITY(P1->getGlobalNumEntries(), R1->getGlobalNumEntries());
    TEST_EQUALITY(P1->getGlobalNumRows(), R1->getGlobalNumCols());
    TEST_EQUALITY(P1->getGlobalNumCols(), R1->getGlobalNumRows());
    TEST_EQUALITY(P2->getGlobalNumEntries(), R2->getGlobalNumEntries());
    TEST_EQUALITY(P2->getGlobalNumRows(), R2->getGlobalNumCols());
    TEST_EQUALITY(P2->getGlobalNumCols(), R2->getGlobalNumRows());


    //RCP<Teuchos::FancyOStream> fos = getFancyOStream(Teuchos::rcpFromRef(cout));

    // since A is chosen symmetric, it is P^T = R
    // check P^T * P = R * P
    // note: the Epetra matrix-matrix multiplication using implicit transpose is buggy in parallel case
    //       (for multiplication of a square matrix with a rectangular matrix)
    //       however it seems to work for two rectangular matrices
    Teuchos::RCP<Xpetra::Matrix<Scalar,LO,GO> > RP = MueLu::Utils<Scalar,LO,GO>::Multiply(*R1,false,*P1,false);
    Teuchos::RCP<Xpetra::Matrix<Scalar,LO,GO> > PtP = MueLu::Utils<Scalar,LO,GO>::Multiply(*P1,true,*P1,false);

    RCP<Vector> x = VectorFactory::Build(RP->getDomainMap());
    RCP<Vector> bRP  = VectorFactory::Build(RP->getRangeMap());
    RCP<Vector> bPtP = VectorFactory::Build(PtP->getRangeMap());

    x->randomize();
    RP->apply(*x,*bRP);
    PtP->apply(*x,*bPtP);

    TEST_EQUALITY(bRP->norm1() - bPtP->norm1() < 1e-12, true);

    Teuchos::RCP<Xpetra::Matrix<Scalar,LO,GO> > RP2 = MueLu::Utils<Scalar,LO,GO>::Multiply(*R2,false,*P2,false);
    Teuchos::RCP<Xpetra::Matrix<Scalar,LO,GO> > PtP2 = MueLu::Utils<Scalar,LO,GO>::Multiply(*P2,true,*P2,false);

    x = VectorFactory::Build(RP2->getDomainMap());
    bRP  = VectorFactory::Build(RP2->getRangeMap());
    bPtP = VectorFactory::Build(PtP2->getRangeMap());

    x->randomize();
    RP2->apply(*x,*bRP);
    PtP2->apply(*x,*bPtP);

    TEST_EQUALITY(bRP->norm1() - bPtP->norm1() < 1e-12, true);


    //R1->describe(*fos,Teuchos::VERB_EXTREME);

    /*RCP<CrsMatrixWrap> crsP1 = rcp_dynamic_cast<CrsMatrixWrap>(P1);
    RCP<CrsMatrix> crsMat = crsP1->getCrsMatrix();
    RCP<Xpetra::EpetraCrsMatrix> epcrsMat = rcp_dynamic_cast<Xpetra::EpetraCrsMatrix>(crsMat);
    RCP<Epetra_CrsMatrix> epMat = epcrsMat->getEpetra_CrsMatrixNonConst();
    EpetraExt::RowMatrixToMatrixMarketFile( "Test.mat", *epMat);*/

    //P1->describe(*fos,Teuchos::VERB_EXTREME);

    //R1->getRangeMap()->describe(*fos,Teuchos::VERB_EXTREME);
    //P1->describe(*fos,Teuchos::VERB_EXTREME);
    //R1->describe(*fos,Teuchos::VERB_EXTREME);



  }

  // check Hierarchy::Setup routine with GenericRFactory as restriction factory
  TEUCHOS_UNIT_TEST(GenericRFactory, GenericRSetup)
  {
    out << "version: " << MueLu::Version() << std::endl;

    for (int i=1; i<5; i++) {
      // generate problem
      RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
      RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO, LMO>::Build1DPoisson(50*i*comm->getSize());

      // Multigrid Hierarchy
      Hierarchy H(A);
      H.setVerbLevel(Teuchos::VERB_HIGH);

      // build nullspace
      RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(),1);
      nullSpace->putScalar( (SC) 1.0);
      Teuchos::Array<ST::magnitudeType> norms(1);
      nullSpace->norm1(norms);
      if (comm->getRank() == 0)
        out << "||NS|| = " << norms[0] << std::endl;

      RCP<PgPFactory>         Pfact = rcp( new PgPFactory());
      RCP<Factory>           Rfact = rcp( new GenericRFactory() );
      RCP<RAPFactory>        Acfact = rcp( new RAPFactory() );

      // setup smoothers
      Teuchos::ParameterList smootherParamList;
      smootherParamList.set("relaxation: type", "Symmetric Gauss-Seidel");
      smootherParamList.set("relaxation: sweeps", (LO) 1);
      smootherParamList.set("relaxation: damping factor", (SC) 1.0);
      RCP<SmootherPrototype> smooProto = rcp( new TrilinosSmoother("RELAXATION", smootherParamList) );
      RCP<SmootherFactory> SmooFact = rcp( new SmootherFactory(smooProto) );
      Acfact->setVerbLevel(Teuchos::VERB_HIGH);

      // Multigrid setup phase (using default parameters)
      FactoryManager M0; // how to build aggregates and smoother of the first level
      M0.SetFactory("A", Acfact);
      M0.SetFactory("P", Pfact);
      M0.SetFactory("R", Rfact);
      M0.SetFactory("Smoother", SmooFact);
      M0.SetFactory("CoarseSolver", SmooFact);

      FactoryManager M1; // first coarse level (Plain aggregation)
      M1.SetFactory("A", Acfact);
      M1.SetFactory("P", Pfact);
      M1.SetFactory("R", Rfact);
      M1.SetFactory("Smoother", SmooFact);
      M1.SetFactory("CoarseSolver", SmooFact);

      FactoryManager M2; // last level (SA)
      M2.SetFactory("A", Acfact);
      M2.SetFactory("P", Pfact);
      M2.SetFactory("R", Rfact);
      M2.SetFactory("Smoother", SmooFact);
      M2.SetFactory("CoarseSolver", SmooFact);

      bool bIsLastLevel = false;
      if(!bIsLastLevel) bIsLastLevel = H.Setup(0, Teuchos::null,  ptrInArg(M0), ptrInArg(M1));
      if(!bIsLastLevel) bIsLastLevel = H.Setup(1, ptrInArg(M0), ptrInArg(M1), ptrInArg(M2));
      if(!bIsLastLevel) bIsLastLevel = H.Setup(2, ptrInArg(M1), ptrInArg(M2), Teuchos::null );

      RCP<Level> l0 = H.GetLevel(0);
      RCP<Level> l1;
      RCP<Level> l2;

      if(H.GetNumLevels() > 1) l1 = H.GetLevel(1);
      if(H.GetNumLevels() > 2) l2 = H.GetLevel(2);

      /*RCP<Teuchos::FancyOStream> stdout = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      l0->print(*stdout,Teuchos::VERB_EXTREME);
      if(l1 != Teuchos::null) l1->print(*stdout,Teuchos::VERB_EXTREME);
      if(l2 != Teuchos::null) l2->print(*stdout,Teuchos::VERB_EXTREME);*/

      TEST_EQUALITY(l0->IsAvailable("PreSmoother",  MueLu::NoFactory::get()), true);
      TEST_EQUALITY(l0->IsAvailable("PostSmoother", MueLu::NoFactory::get()), true);
      TEST_EQUALITY(l0->IsAvailable("A",            MueLu::NoFactory::get()), true);
      TEST_EQUALITY(l0->GetKeepFlag("PreSmoother",  MueLu::NoFactory::get()), MueLu::Final);
      TEST_EQUALITY(l0->GetKeepFlag("A",            MueLu::NoFactory::get()), MueLu::UserData);
      TEST_EQUALITY(l0->IsAvailable("P",            MueLu::NoFactory::get()), false);
      TEST_EQUALITY(l0->IsAvailable("R",            MueLu::NoFactory::get()), false);
      TEST_EQUALITY(l0->IsRequested("P",            MueLu::NoFactory::get()), false);
      TEST_EQUALITY(l0->IsRequested("R",            MueLu::NoFactory::get()), false);

      if(l1 != Teuchos::null) {
        TEST_EQUALITY(l1->IsAvailable("PreSmoother",  MueLu::NoFactory::get()), true);
        TEST_EQUALITY(l1->IsAvailable("PostSmoother", MueLu::NoFactory::get()), true);
        TEST_EQUALITY(l1->IsAvailable("P",            MueLu::NoFactory::get()), true);
        TEST_EQUALITY(l1->IsAvailable("R",            MueLu::NoFactory::get()), true);
        TEST_EQUALITY(l1->IsAvailable("A",            MueLu::NoFactory::get()), true);
        TEST_EQUALITY(l1->GetKeepFlag("A",            MueLu::NoFactory::get()), MueLu::Final);
        TEST_EQUALITY(l1->GetKeepFlag("P",            MueLu::NoFactory::get()), MueLu::Final);
        TEST_EQUALITY(l1->GetKeepFlag("R",            MueLu::NoFactory::get()), MueLu::Final);
        TEST_EQUALITY(l1->GetKeepFlag("PostSmoother", MueLu::NoFactory::get()), MueLu::Final);
        TEST_EQUALITY(l1->IsRequested("Graph",        MueLu::NoFactory::get()), false);
        TEST_EQUALITY(l1->IsRequested("Aggregates",   MueLu::NoFactory::get()), false);
        TEST_EQUALITY(l1->IsRequested("Nullspace",    MueLu::NoFactory::get()), false);
      }
      if(l2 != Teuchos::null) {
        TEST_EQUALITY(l2->IsAvailable("PreSmoother",  MueLu::NoFactory::get()), true);
        TEST_EQUALITY(l2->IsAvailable("PostSmoother", MueLu::NoFactory::get()), true);
        TEST_EQUALITY(l2->IsAvailable("P",            MueLu::NoFactory::get()), true);
        TEST_EQUALITY(l2->IsAvailable("R",            MueLu::NoFactory::get()), true);
        TEST_EQUALITY(l2->IsAvailable("A",            MueLu::NoFactory::get()), true);
        TEST_EQUALITY(l2->IsRequested("Graph",        MueLu::NoFactory::get()), false);
        TEST_EQUALITY(l2->IsRequested("Aggregates",   MueLu::NoFactory::get()), false);
        TEST_EQUALITY(l2->IsRequested("Nullspace",    MueLu::NoFactory::get()), false);
        TEST_EQUALITY(l2->GetKeepFlag("A",            MueLu::NoFactory::get()), MueLu::Final);
        TEST_EQUALITY(l2->GetKeepFlag("P",            MueLu::NoFactory::get()), MueLu::Final);
        TEST_EQUALITY(l2->GetKeepFlag("R",            MueLu::NoFactory::get()), MueLu::Final);
        TEST_EQUALITY(l2->GetKeepFlag("PreSmoother",  MueLu::NoFactory::get()), MueLu::Final);
        TEST_EQUALITY(l2->GetKeepFlag("PostSmoother", MueLu::NoFactory::get()), MueLu::Final);
      }
    } // end for i=1..5

  }

}//namespace MueLuTests


