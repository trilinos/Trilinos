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
#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_NotayAggregationFactory.hpp"
#include "MueLu_Types.hpp"

namespace MueLuTests {

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(NotayAggregation, InitialAggregation1D, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);

    //    using TST                   = Teuchos::ScalarTraits<SC>;
    //    using magnitude_type        = typename TST::magnitudeType;
    //    using TMT                   = Teuchos::ScalarTraits<magnitude_type>;
    //    using real_type             = typename TST::coordinateType;
    //    using RealValuedMultiVector = Xpetra::MultiVector<real_type,LO,GO,NO>;
    using test_factory          = TestHelpers::TestFactory<SC, LO, GO, NO>;

    out << "version: " << MueLu::Version() << std::endl;
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    int rank = comm->getRank();
    int numproc = comm->getSize();
    RCP<const Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(16*comm->getSize());
    const int numRows = static_cast<int>(A->getNodeNumRows());
    RCP<Aggregates> aggregates = rcp(new Aggregates(A->getMap()));
    RCP<NotayAggregationFactory> NAF = rcp(new NotayAggregationFactory());
    std::vector<unsigned> aggStat(numRows, MueLu::READY);
    LO numUnaggregatedNodes = numRows, numDirichletNodes = 0;

    Array<LO> orderingVector(numRows);
    for (LO i = 0; i < numRows; i++) {
      orderingVector[i] = i;
    }

    Teuchos::ParameterList params;
    params.set("aggregation: pairwise: tie threshold",1e-6);
    NAF->BuildInitialAggregates(params, A, orderingVector(),
                                Teuchos::ScalarTraits<SC>::magnitude(10.0),
                                *aggregates, aggStat, numUnaggregatedNodes, numDirichletNodes);

    auto v2a = aggregates->GetVertex2AggId()->getData(0);
    Teuchos::ArrayRCP<LO> sizes = aggregates->ComputeAggregateSizes();
    std::cout << "p=" << rank << " | aggregate sizes="
              << sizes.view(0, sizes.size()) << std::endl;

#if 0
    printf("[%d] Aggregates: ",rank);
    for(int i=0; i<(int)v2a.size(); i++)
      printf("%d(%d) ",i,v2a[i]);
    printf("\n");
#endif

    //    TEST_EQUALITY(numUnaggregatedNodes, 0);
    if(numproc == 1) {
      TEST_EQUALITY(aggregates->GetNumAggregates(),7);
    }
    else {
      TEST_EQUALITY(aggregates->GetNumAggregates(),8);
    }

    // On proc 0 gets picked as a Dirichlet (so it does not get aggregated) and the last aggregate is a singleton
    // All the other ranks wind up with strict pairs
    int expected;
    for(int i=0; i<(int)sizes.size(); i++) {
      if(numproc == 1) expected = 2;
      else expected = ((rank == 0 || rank == 3) && i == (int)sizes.size() -1) ? 1 : 2;
      TEST_EQUALITY(sizes[i], expected);
    }

  } // InitialAggregation1D



  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(NotayAggregation, InitialAggregation2D, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);

    //    using TST                   = Teuchos::ScalarTraits<SC>;
    //    using magnitude_type        = typename TST::magnitudeType;
    //    using TMT                   = Teuchos::ScalarTraits<magnitude_type>;
    //    using real_type             = typename TST::coordinateType;
    //    using RealValuedMultiVector = Xpetra::MultiVector<real_type,LO,GO,NO>;
    using test_factory          = TestHelpers::TestFactory<SC, LO, GO, NO>;

    out << "version: " << MueLu::Version() << std::endl;
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    int rank = comm->getRank();
    int numproc = comm->getSize();
    Xpetra::UnderlyingLib lib = TestHelpers::Parameters::getLib();

    // Do a 2D Star2D Matrix with up/down as a "weak connection"
    Teuchos::ParameterList mp;
    mp.set("matrixType","Star2D");
    const int nx = 3;
    mp.set("nx",(GO)nx);
    mp.set("ny",(GO)nx*comm->getSize());
    mp.set("a",2.2);
    mp.set("b",-0.1);  mp.set("c",-0.1);
    mp.set("d",-1.0); mp.set("e",-1.0);
    mp.set("z1",0.0);  mp.set("z2",0.0);  mp.set("z3",0.0);  mp.set("z4",0.0);
    RCP<const Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::BuildMatrix(mp,lib);
    const int numRows = static_cast<int>(A->getNodeNumRows());

    RCP<Aggregates> aggregates = rcp(new Aggregates(A->getMap()));
    RCP<NotayAggregationFactory> NAF = rcp(new NotayAggregationFactory());
    std::vector<unsigned> aggStat(numRows, MueLu::READY);
    LO numUnaggregatedNodes = numRows, numDirichletNodes = 0;

    Array<LO> orderingVector(numRows);
    for (LO i = 0; i < numRows; i++) {
      orderingVector[i] = i;
    }

    Teuchos::ParameterList params;
    params.set("aggregation: pairwise: tie threshold",1e-6);
    NAF->BuildInitialAggregates(params, A, orderingVector(),
                                Teuchos::ScalarTraits<SC>::magnitude(4.1),
                                *aggregates, aggStat, numUnaggregatedNodes, numDirichletNodes);

    auto v2a = aggregates->GetVertex2AggId()->getData(0);
    Teuchos::ArrayRCP<LO> sizes = aggregates->ComputeAggregateSizes();

    TEST_EQUALITY(numUnaggregatedNodes, 0);

    // For this problem, the four corners will be detected as "Dirichlet" and ignored, this will
    // generate some number of singletons
    Teuchos::Array<int> expected;
    if(numproc == 1) {
      expected = Teuchos::Array<int>({2, 1, 1, 1});
    } else if(rank == 0 || rank == numproc - 1) {
      expected = Teuchos::Array<int>({2, 2, 2, 1});
    } else {
      expected = Teuchos::Array<int>({2, 2, 2, 1, 1, 1});
    }

    TEST_EQUALITY(sizes.size(), expected.size());


    for(int i = 0; i < static_cast<int>(sizes.size()); i++) {
      TEST_EQUALITY(sizes[i], expected[i]);
    }

#if 0
    printf("[%d] Aggregates: ",rank);
    for(int i=0; i<(int)v2a.size(); i++)
      printf("%d(%d) ",i,v2a[i]);
    printf("\n");
#endif

  } // InitialAggregation2D


  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(NotayAggregation, IntermediateProlongator2D, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);

    //    using TST                   = Teuchos::ScalarTraits<SC>;
    //    using magnitude_type        = typename TST::magnitudeType;
    //    using TMT                   = Teuchos::ScalarTraits<magnitude_type>;
    //    using real_type             = typename TST::coordinateType;
    //    using RealValuedMultiVector = Xpetra::MultiVector<real_type,LO,GO,NO>;
    using test_factory          = TestHelpers::TestFactory<SC, LO, GO, NO>;

    out << "version: " << MueLu::Version() << std::endl;
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib = TestHelpers::Parameters::getLib();

    // Do a 2D Star2D Matrix with up/down as a "weak connection"
    Teuchos::ParameterList mp;
    mp.set("matrixType","Star2D");
    const int nx = 3;
    mp.set("nx",(GO)nx);
    mp.set("ny",(GO)nx*comm->getSize());
    mp.set("a",2.2);
    mp.set("b",-0.1);  mp.set("c",-0.1);
    mp.set("d",-1.0); mp.set("e",-1.0);
    mp.set("z1",0.0);  mp.set("z2",0.0);  mp.set("z3",0.0);  mp.set("z4",0.0);
    RCP<const Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::BuildMatrix(mp,lib);
    const int numRows = static_cast<int>(A->getNodeNumRows());

    RCP<Aggregates> aggregates = rcp(new Aggregates(A->getMap()));
    RCP<NotayAggregationFactory> NAF = rcp(new NotayAggregationFactory());
    std::vector<unsigned> aggStat(numRows, MueLu::READY);
    LO numUnaggregatedNodes = numRows, numDirichletNodes = 0;
    typename Matrix::local_matrix_type intermediateP;

    Array<LO> orderingVector(numRows);
    for (LO i = 0; i < numRows; i++) {
      orderingVector[i] = i;
    }

    Teuchos::ParameterList params;
    params.set("aggregation: pairwise: tie threshold",1e-6);
    NAF->BuildInitialAggregates(params, A, orderingVector(),
                                Teuchos::ScalarTraits<SC>::magnitude(4.1),
                                *aggregates, aggStat, numUnaggregatedNodes, numDirichletNodes);

    out << "numDirichletNodes=" << numDirichletNodes << std::endl;

    auto v2a = aggregates->GetVertex2AggId()->getData(0);
    Teuchos::ArrayRCP<LO> sizes = aggregates->ComputeAggregateSizes();

    TEST_EQUALITY(numUnaggregatedNodes, 0);

    NAF->BuildIntermediateProlongator(A->getNodeNumRows(), numDirichletNodes,
                                      aggregates->GetNumAggregates(),
                                      v2a.view(0, A->getNodeNumRows()),
                                      intermediateP);

  } // IntermediateProlongator2D


  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(NotayAggregation, CoarseLocalMatrix2D, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);

    using test_factory          = TestHelpers::TestFactory<SC, LO, GO, NO>;

    out << "version: " << MueLu::Version() << std::endl;
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib = TestHelpers::Parameters::getLib();

    // Do a 2D Star2D Matrix with up/down as a "weak connection"
    Teuchos::ParameterList mp;
    mp.set("matrixType","Star2D");
    const GO nx = 3;
    mp.set("nx", nx);
    mp.set("ny", nx*comm->getSize());
    mp.set("a",  2.2);
    mp.set("b", -0.1); mp.set("c", -0.1);
    mp.set("d", -1.0); mp.set("e", -1.0);
    mp.set("z1", 0.0); mp.set("z2", 0.0);  mp.set("z3", 0.0);  mp.set("z4", 0.0);
    RCP<const Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::BuildMatrix(mp,lib);
    const int numRows = static_cast<int>(A->getNodeNumRows());

    RCP<Aggregates> aggregates = rcp(new Aggregates(A->getMap()));
    RCP<NotayAggregationFactory> NAF = rcp(new NotayAggregationFactory());
    std::vector<unsigned> aggStat(numRows, MueLu::READY);
    LO numUnaggregatedNodes = numRows, numDirichletNodes = 0;
    typename Matrix::local_matrix_type intermediateP;
    typename Matrix::local_matrix_type coarseA = A->getLocalMatrix();

    Array<LO> orderingVector(numRows);
    for (LO i = 0; i < numRows; i++) {
      orderingVector[i] = i;
    }

    Teuchos::ParameterList params;
    params.set("aggregation: pairwise: tie threshold",1e-6);
    NAF->BuildInitialAggregates(params, A, orderingVector(),
                                Teuchos::ScalarTraits<SC>::magnitude(4.1),
                                *aggregates, aggStat, numUnaggregatedNodes, numDirichletNodes);

    auto v2a = aggregates->GetVertex2AggId()->getData(0);
    Teuchos::ArrayRCP<LO> sizes = aggregates->ComputeAggregateSizes();

    TEST_EQUALITY(numUnaggregatedNodes, 0);

    NAF->BuildOnRankLocalMatrix(A->getLocalMatrix(), coarseA);
    NAF->BuildIntermediateProlongator(A->getNodeNumRows(), numDirichletNodes,
                                      aggregates->GetNumAggregates(),
                                      v2a.view(0, A->getNodeNumRows()),
                                      intermediateP);
    NAF->BuildCoarseLocalMatrix(intermediateP, coarseA);

  } // CoarseLocalMatrix2D


  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(NotayAggregation, BuildNotayAggregates2D, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);

    using test_factory          = TestHelpers::TestFactory<SC, LO, GO, NO>;

    out << "version: " << MueLu::Version() << std::endl;
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib = TestHelpers::Parameters::getLib();

    // Do a 2D Star2D Matrix with up/down as a "weak connection"
    Teuchos::ParameterList mp;
    mp.set("matrixType","Star2D");
    const int nx = 9;
    mp.set("nx",(GO)nx);
    mp.set("ny",(GO)nx*comm->getSize());
    mp.set("a",2.2);
    mp.set("b",-0.1);  mp.set("c",-0.1);
    mp.set("d",-1.0); mp.set("e",-1.0);
    mp.set("z1",0.0);  mp.set("z2",0.0);  mp.set("z3",0.0);  mp.set("z4",0.0);
    RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::BuildMatrix(mp,lib);

    // Use default ordering
    {
      MueLu::Level currentLevel;
      test_factory::createSingleLevelHierarchy(currentLevel);
      currentLevel.SetFactoryManager(Teuchos::null);
      currentLevel.Request("A");
      currentLevel.Set("A", A);
      
      RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
      RCP<CoalesceDropFactory> dropFact  = rcp(new CoalesceDropFactory());
      dropFact->SetFactory("UnAmalgamationInfo", amalgFact);
      
      RCP<NotayAggregationFactory> NAF   = rcp(new NotayAggregationFactory());
      NAF->SetFactory("Graph", dropFact);
      NAF->SetFactory("DofsPerNode", dropFact);
      NAF->SetParameter("aggregation: pairwise: size", Teuchos::ParameterEntry(4));
      NAF->SetParameter("aggregation: Dirichlet threshold", Teuchos::ParameterEntry(4.1));
      
      currentLevel.Request("Aggregates", NAF.get());
      currentLevel.Request(*NAF);
      NAF->Build(currentLevel);
      out << "Pairwise aggregation factory is done" << std::endl;
      RCP<Aggregates> aggregates = currentLevel.Get<RCP<Aggregates> >("Aggregates", NAF.get());
      out << "Recovered aggregates from muelu level" << std::endl;
      currentLevel.Release("Aggregates", NAF.get());
      std::cout << "Testing pairwise aggregates" << std::endl;
      
      auto v2a = aggregates->GetVertex2AggId()->getData(0);
      Teuchos::ArrayRCP<LO> sizes = aggregates->ComputeAggregateSizes();
      
      TEST_EQUALITY(aggregates->AggregatesCrossProcessors(),false);
    }

    // Use random ordering
    {
      MueLu::Level currentLevel;
      test_factory::createSingleLevelHierarchy(currentLevel);
      currentLevel.SetFactoryManager(Teuchos::null);
      currentLevel.Request("A");
      currentLevel.Set("A", A);
      
      RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
      RCP<CoalesceDropFactory> dropFact  = rcp(new CoalesceDropFactory());
      dropFact->SetFactory("UnAmalgamationInfo", amalgFact);
      
      RCP<NotayAggregationFactory> NAF   = rcp(new NotayAggregationFactory());
      NAF->SetFactory("Graph", dropFact);
      NAF->SetFactory("DofsPerNode", dropFact);
      NAF->SetParameter("aggregation: pairwise: size", Teuchos::ParameterEntry(4));
      NAF->SetParameter("aggregation: Dirichlet threshold", Teuchos::ParameterEntry(4.1));
      NAF->SetParameter("aggregation: ordering", Teuchos::ParameterEntry(std::string("random")));
      
      currentLevel.Request("Aggregates", NAF.get());
      currentLevel.Request(*NAF);
      NAF->Build(currentLevel);
      out << "Pairwise aggregation factory is done" << std::endl;
      RCP<Aggregates> aggregates = currentLevel.Get<RCP<Aggregates> >("Aggregates", NAF.get());
      out << "Recovered aggregates from muelu level" << std::endl;
      currentLevel.Release("Aggregates", NAF.get());
      std::cout << "Testing pairwise aggregates" << std::endl;
      
      auto v2a = aggregates->GetVertex2AggId()->getData(0);
      Teuchos::ArrayRCP<LO> sizes = aggregates->ComputeAggregateSizes();
      
      TEST_EQUALITY(aggregates->AggregatesCrossProcessors(),false);
    }

#if defined(HAVE_MUELU_KOKKOS_REFACTOR)
    // Use Cuthill-McKee if we have it (requires Kokkos refactor)
   {
      MueLu::Level currentLevel;
      test_factory::createSingleLevelHierarchy(currentLevel);
      currentLevel.SetFactoryManager(Teuchos::null);
      currentLevel.Request("A");
      currentLevel.Set("A", A);
      
      RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
      RCP<CoalesceDropFactory> dropFact  = rcp(new CoalesceDropFactory());
      dropFact->SetFactory("UnAmalgamationInfo", amalgFact);
      
      RCP<NotayAggregationFactory> NAF   = rcp(new NotayAggregationFactory());
      NAF->SetFactory("Graph", dropFact);
      NAF->SetFactory("DofsPerNode", dropFact);
      NAF->SetParameter("aggregation: pairwise: size", Teuchos::ParameterEntry(4));
      NAF->SetParameter("aggregation: Dirichlet threshold", Teuchos::ParameterEntry(4.1));
      NAF->SetParameter("aggregation: ordering", Teuchos::ParameterEntry(std::string("cuthill-mckee")));
      
      currentLevel.Request("Aggregates", NAF.get());
      currentLevel.Request(*NAF);
      NAF->Build(currentLevel);
      out << "Pairwise aggregation factory is done" << std::endl;
      RCP<Aggregates> aggregates = currentLevel.Get<RCP<Aggregates> >("Aggregates", NAF.get());
      out << "Recovered aggregates from muelu level" << std::endl;
      currentLevel.Release("Aggregates", NAF.get());
      std::cout << "Testing pairwise aggregates" << std::endl;
      
      auto v2a = aggregates->GetVertex2AggId()->getData(0);
      Teuchos::ArrayRCP<LO> sizes = aggregates->ComputeAggregateSizes();
      
      TEST_EQUALITY(aggregates->AggregatesCrossProcessors(),false);
    }

#endif 

  } // BuildNotayAggregates2D

#  define MUELU_ETI_GROUP(Scalar, LO, GO, Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(NotayAggregation,InitialAggregation1D,Scalar,LO,GO,Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(NotayAggregation,InitialAggregation2D,Scalar,LO,GO,Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(NotayAggregation,IntermediateProlongator2D,Scalar,LO,GO,Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(NotayAggregation,CoarseLocalMatrix2D,Scalar,LO,GO,Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(NotayAggregation,BuildNotayAggregates2D,Scalar,LO,GO,Node)


#include <MueLu_ETI_4arg.hpp>


} // namespace MueLuTests
