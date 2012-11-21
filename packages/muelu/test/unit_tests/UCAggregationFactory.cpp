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
#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_UCAggregationFactory.hpp"

namespace {

  using std::string; //?? TODO

  //TODO: should go in the Aggregates class
  template <class LocalOrdinal,
            class GlobalOrdinal,
            class Node,
            class LocalMatOps>
  void printAggregates(MueLu::Aggregates<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>& aggregates, Teuchos::FancyOStream& out) {
#include "MueLu_UseShortNamesOrdinal.hpp"
    RCP<LOVector> Final_ = LOVectorFactory::Build( aggregates.GetVertex2AggId()->getMap() );

    ArrayRCP<LO> Final = Final_->getDataNonConst(0);
    ArrayRCP<const LO> vertex2AggId = aggregates.GetVertex2AggId()->getData(0);
    ArrayRCP<const LO> procWinner   = aggregates.GetProcWinner()->getData(0);

    for (size_t i=0; i<aggregates.GetVertex2AggId()->getMap()->getNodeNumElements(); i++)
      Final[i] = vertex2AggId[i] + procWinner[i]*1000;

    out << *Final_ << std::endl;
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_5_DECL(UCAggregationFactory, Constructor, Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps)
  {
    MUELU_TEST_EPETRA_ONLY_FOR_DOUBLE_AND_INT(Scalar, LocalOrdinal, GlobalOrdinal) {
#include "MueLu_UseShortNames.hpp"

      out << "version: " << MueLu::Version() << std::endl;
      RCP<UCAggregationFactory> aggFact = rcp(new UCAggregationFactory());
      TEST_EQUALITY(aggFact != Teuchos::null, true);
    }
  } // Constructor

  TEUCHOS_UNIT_TEST_TEMPLATE_5_DECL(UCAggregationFactory, Build, Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps)
  {
    MUELU_TEST_EPETRA_ONLY_FOR_DOUBLE_AND_INT(Scalar, LocalOrdinal, GlobalOrdinal) {
      //    typedef double Scalar;
#include "MueLu_UseShortNames.hpp"

      out << "version: " << MueLu::Version() << std::endl;

      RCP<Matrix> Op = TestHelpers::TestFactory<SC, LO, GO, NO, LMO>::Build1DPoisson(16);
      RCP<Graph> graph = rcp(new Graph(Op->getCrsGraph(), "someGraphLabel"));

      {

        UCAggregationFactory aggFact;
        aggFact.SetPrintFlag(6);
        aggFact.SetMinNodesPerAggregate(2);
        aggFact.SetMaxNeighAlreadySelected(5);
        aggFact.SetOrdering(MueLu::AggOptions::NATURAL);
        aggFact.SetPhase3AggCreation(0.5);

        RCP<Aggregates> aggregates;

        aggregates = aggFact.Build(*graph);
        printAggregates(*aggregates, out);
      }

      {

        UCAggregationFactory aggFact;
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

        UCAggregationFactory aggFact;
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
  typedef Kokkos::DefaultNode::DefaultNodeType Node; // Kokkos Node is not relevant for this test

  typedef long int LongInt;                          // macros dislike parameters with space...
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
  typedef long long int LongLongInt;
#endif

#define UNIT_TEST_GROUP_5(SC, LO, GO, NO, LMO)                          \
  TEUCHOS_UNIT_TEST_TEMPLATE_5_INSTANT(UCAggregationFactory, Constructor, SC, LO, GO, NO, LMO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_5_INSTANT(UCAggregationFactory, Build,       SC, LO, GO, NO, LMO)

#define UNIT_TEST_GROUP_2(LO, GO)                                       \
  typedef Kokkos::DefaultKernels<Scalar,LO,Node>::SparseOps LMO ## LO;  \
  UNIT_TEST_GROUP_5(Scalar, LO, GO, Node, LMO ## LO)

  UNIT_TEST_GROUP_2(int, int)
  UNIT_TEST_GROUP_2(int, LongInt)
  UNIT_TEST_GROUP_2(LongInt, LongInt)

#ifdef HAVE_TEUCHOS_LONG_LONG_INT
  UNIT_TEST_GROUP_2(LongInt, LongLongInt)
#endif

} // namespace <anonymous>
