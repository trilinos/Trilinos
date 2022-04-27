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

#include "MueLu_CoupledAggregationFactory.hpp"

namespace MueLuTests {

#include "MueLu_UseShortNamesOrdinal.hpp"

  using std::string; //?? TODO

  //TODO: should go in the Aggregates class
  template <class LocalOrdinal,
           class GlobalOrdinal,
           class Node>
             void printAggregates(MueLu::Aggregates<LocalOrdinal, GlobalOrdinal, Node>& aggregates, Teuchos::FancyOStream& out) {
               RCP<LOVector> Final_ = LOVectorFactory::Build( aggregates.GetVertex2AggId()->getMap() );

               ArrayRCP<LO> Final = Final_->getDataNonConst(0);
               ArrayRCP<const LO> vertex2AggId = aggregates.GetVertex2AggId()->getData(0);
               ArrayRCP<const LO> procWinner   = aggregates.GetProcWinner()->getData(0);

               for (size_t i=0; i<aggregates.GetVertex2AggId()->getMap()->getLocalNumElements(); i++)
                 Final[i] = vertex2AggId[i] + procWinner[i]*1000;

               out << *Final_ << std::endl;
             }

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoupledAggregationFactory, Constructor, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
    MUELU_TEST_EPETRA_ONLY_FOR_DOUBLE_AND_INT(Scalar, LocalOrdinal, GlobalOrdinal) {

      out << "version: " << MueLu::Version() << std::endl;
      RCP<CoupledAggregationFactory> aggFact = rcp(new CoupledAggregationFactory());
      TEST_EQUALITY(aggFact != Teuchos::null, true);
    }
  } // Constructor

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoupledAggregationFactory, Build, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
    MUELU_TEST_EPETRA_ONLY_FOR_DOUBLE_AND_INT(Scalar, LocalOrdinal, GlobalOrdinal) {
      //    typedef double Scalar;
#include "MueLu_UseShortNames.hpp"

      out << "version: " << MueLu::Version() << std::endl;

      RCP<Matrix> Op = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(16);
      RCP<Graph> graph = rcp(new Graph(Op->getCrsGraph(), "someGraphLabel"));

      {

        CoupledAggregationFactory aggFact;
        aggFact.SetPrintFlag(6);
        aggFact.SetMinNodesPerAggregate(2);
        aggFact.SetMaxNeighAlreadySelected(5);
        aggFact.SetOrdering("natural");
        aggFact.SetPhase3AggCreation(0.5);

        RCP<Aggregates> aggregates;

        aggregates = aggFact.Build(*graph);
        printAggregates(*aggregates, out);
      }

      {

        CoupledAggregationFactory aggFact;
        aggFact.SetPrintFlag(6);
        aggFact.SetMinNodesPerAggregate(2);
        aggFact.SetMaxNeighAlreadySelected(5);
        aggFact.SetOrdering(MueLu::AggOptions::RANDOM);
        aggFact.SetPhase3AggCreation(0.5);

        RCP<Aggregates> aggregates;

        aggregates = aggFact.Build(*graph);
        printAggregates(*aggregates, out);
      }

      {

        CoupledAggregationFactory aggFact;
        aggFact.SetPrintFlag(6);
        aggFact.SetMinNodesPerAggregate(2);
        aggFact.SetMaxNeighAlreadySelected(5);
        aggFact.SetOrdering(MueLu::AggOptions::GRAPH);
        aggFact.SetPhase3AggCreation(0.5);

        RCP<Aggregates> aggregates;

        aggregates = aggFact.Build(*graph);
        printAggregates(*aggregates, out);
      }
    }
  } // Build

  //
  // INSTANTIATIONS
  //

  typedef double Scalar;                             // Scalar is not relevant for this test
  typedef KokkosClassic::DefaultNode::DefaultNodeType Node; // Kokkos Node is not relevant for this test

  typedef long int LongInt;                          // macros dislike parameters with space...
#ifdef HAVE_XPETRA_INT_LONG_LONG
  typedef long long int LongLongInt;
#endif

#define UNIT_TEST_GROUP_4(SC, LO, GO, NO)                          \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoupledAggregationFactory, Constructor, SC, LO, GO, NO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoupledAggregationFactory, Build,       SC, LO, GO, NO)

#define UNIT_TEST_GROUP_2(LO, GO)                                       \
  UNIT_TEST_GROUP_4(Scalar, LO, GO, Node)

  UNIT_TEST_GROUP_2(int, int)
    UNIT_TEST_GROUP_2(int, LongInt)
    UNIT_TEST_GROUP_2(LongInt, LongInt)

#ifdef HAVE_XPETRA_INT_LONG_LONG
    UNIT_TEST_GROUP_2(LongInt, LongLongInt)
#endif

} // namespace <anonymous>
