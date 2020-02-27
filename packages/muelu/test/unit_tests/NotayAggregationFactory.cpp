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

#include "MueLu_NotayAggregationFactory.hpp"
#include "MueLu_Aggregates.hpp"
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
    RCP<Aggregates> aggregates = rcp(new Aggregates(A->getMap()));    
    RCP<NotayAggregationFactory> NAF = rcp(new NotayAggregationFactory());
    std::vector<unsigned> aggStat(A->getMap()->getNodeNumElements(),MueLu::READY);
    LO numUnaggregatedNodes;

    Teuchos::ParameterList params;
    params.set("aggregation: Dirichlet threshold",10.0);
    NAF->Build_InitialAggregation(params,A,*aggregates,aggStat,numUnaggregatedNodes);

    auto v2a = aggregates->GetVertex2AggId()->getData(0);
    Teuchos::ArrayRCP<LO> sizes = aggregates->ComputeAggregateSizes();

#if 0
    printf("[%d] Aggregates: ",rank);
    for(int i=0; i<(int)v2a.size(); i++)
      printf("%d(%d) ",i,v2a[i]);
    printf("\n");
#endif 

    //    TEST_EQUALITY(numUnaggregatedNodes, 0);
    TEST_EQUALITY(aggregates->AggregatesCrossProcessors(),false);
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
      else expected = (rank == 0 && i == (int)sizes.size() -1) ? 1 : 2;
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
#if 0
    int rank = comm->getRank();
#endif
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
                                                                                
    RCP<Aggregates> aggregates = rcp(new Aggregates(A->getMap()));    
    RCP<NotayAggregationFactory> NAF = rcp(new NotayAggregationFactory());
    std::vector<unsigned> aggStat(A->getMap()->getNodeNumElements(),MueLu::READY);
    LO numUnaggregatedNodes;

    Teuchos::ParameterList params;
    params.set("aggregation: Dirichlet threshold",4.1);
    NAF->Build_InitialAggregation(params,A,*aggregates,aggStat,numUnaggregatedNodes);

    auto v2a = aggregates->GetVertex2AggId()->getData(0);
    Teuchos::ArrayRCP<LO> sizes = aggregates->ComputeAggregateSizes();
    
    TEST_EQUALITY(numUnaggregatedNodes, 0);
    TEST_EQUALITY(aggregates->AggregatesCrossProcessors(),false);

    // For this problem, the four corners will be detected as "Dirichlet" and ignored, this will
    // generate some number of singletons
    int expected;
    for(int i=0; i<(int)sizes.size(); i++) {
      if(numproc == 1) expected = (i==0) ? 2 : 1;
      else expected = (i==sizes.size()-1) ? 1 : 2;
      TEST_EQUALITY(sizes[i], expected);
    }

#if 0
    printf("[%d] Aggregates: ",rank);
    for(int i=0; i<(int)v2a.size(); i++)
      printf("%d(%d) ",i,v2a[i]);
    printf("\n");
#endif 

  } // InitialAggregation2D

#  define MUELU_ETI_GROUP(Scalar, LO, GO, Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(NotayAggregation,InitialAggregation1D,Scalar,LO,GO,Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(NotayAggregation,InitialAggregation2D,Scalar,LO,GO,Node)


#include <MueLu_ETI_4arg.hpp>


} // namespace MueLuTests
