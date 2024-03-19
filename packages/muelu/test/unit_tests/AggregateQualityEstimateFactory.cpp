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
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include "MueLu_config.hpp"

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"
#include "MueLu_UncoupledAggregationFactory.hpp"
#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_CoarseMapFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include <MueLu_CoalesceDropFactory.hpp>

#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_IO.hpp>

#include "MueLu_Utilities.hpp"
#include "MueLu_AggregateQualityEstimateFactory.hpp"

namespace MueLuTests {

template <typename LocalOrdinal>
std::vector<LocalOrdinal> flip_agg_horizontal(const std::vector<LocalOrdinal>& indices, LocalOrdinal nx) {
  std::vector<LocalOrdinal> flipped;

  for (LocalOrdinal idx : indices) {
    LocalOrdinal x = idx % nx;
    LocalOrdinal y = idx / nx;

    flipped.push_back(nx * y - x);
  }

  std::sort(flipped.begin(), flipped.end());
  return flipped;
}

template <typename LocalOrdinal>
std::vector<LocalOrdinal> flip_agg_vertical(const std::vector<LocalOrdinal>& indices, LocalOrdinal nx) {
  std::vector<LocalOrdinal> flipped;

  for (LocalOrdinal idx : indices) {
    LocalOrdinal x = idx % nx;
    LocalOrdinal y = idx / nx;

    flipped.push_back(-nx * y + x);
  }

  std::sort(flipped.begin(), flipped.end());
  return flipped;
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(AggregateQualityEstimateFactory, Constructor, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;
  typedef Xpetra::MultiVector<MT, LO, GO, Node> MultiVectorDouble;

  out << "version: " << MueLu::Version() << std::endl;

  typedef MueLu::AggregateQualityEstimateFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> AggregateQualityEstimateFactory;

  RCP<AggregateQualityEstimateFactory> aggQualityEstimateFactory = rcp(new AggregateQualityEstimateFactory);
  TEST_EQUALITY(aggQualityEstimateFactory != Teuchos::null, true);

}  // Constructor test

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(AggregateQualityEstimateFactory, Poisson2D, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  out << "version: " << MueLu::Version() << std::endl;

  typedef MueLu::AggregateQualityEstimateFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> AggregateQualityEstimateFactory;
  typedef typename Teuchos::ScalarTraits<Scalar> TST;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;
  typedef Xpetra::MultiVector<MT, LO, GO, NO> MultiVectorDouble;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();

  Level level;
  TestHelpers::TestFactory<Scalar, LO, GO, NO>::createSingleLevelHierarchy(level);

  GO nx          = 20 * comm->getSize();
  GO ny          = nx;
  RCP<Matrix> Op = TestHelpers::TestFactory<Scalar, LO, GO, NO>::Build2DPoisson(nx, ny);
  level.Set("A", Op);

  AggregateQualityEstimateFactory aggQualityEstimateFactory;
  std::cout << *(aggQualityEstimateFactory.GetValidParameterList()) << std::endl;
  aggQualityEstimateFactory.SetParameter("aggregate qualities: check symmetry", Teuchos::ParameterEntry(false));
  aggQualityEstimateFactory.SetParameter("aggregate qualities: good aggregate threshold", Teuchos::ParameterEntry(100.0));
  aggQualityEstimateFactory.SetParameter("aggregate qualities: file output", Teuchos::ParameterEntry(false));

  level.Request("AggregateQualities", &aggQualityEstimateFactory);
  level.Request(aggQualityEstimateFactory);

  out << "Getting aggregate qualities...\n\n";

  RCP<MultiVectorDouble> aggQualities = level.Get<RCP<MultiVectorDouble>>("AggregateQualities", &aggQualityEstimateFactory);

  out << "Testing aggregate qualities to make sure all aggregates are of good quality...\n\n";

  ArrayRCP<const MT> aggQualitiesLocalData = aggQualities->getData(0);

  for (size_t i = 0; i < aggQualities->getLocalLength(); ++i) {
    out << "Aggregate " << i << ": " << aggQualitiesLocalData[i] << "\n";
    TEST_COMPARE(aggQualitiesLocalData[i], >, 0.0);
    TEST_COMPARE(aggQualitiesLocalData[i], <=, 20.0);
  }

}  // Poisson 2D test

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(AggregateQualityEstimateFactory, AnisotropicDiffusion2D, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;
  typedef Xpetra::MultiVector<MT, LO, GO, Node> MultiVectorDouble;
  typedef MueLu::AggregateQualityEstimateFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> AggregateQualityEstimateFactory;

  out << "version: " << MueLu::Version() << std::endl;

  typedef typename Teuchos::ScalarTraits<Scalar> TST;

  // Don't test for complex or ordinal - matrix reader won't work
  if (TST::isComplex || TST::isOrdinal) {
    success = true;
    return;
  }

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();

  LO nx = (LO)20;

  LO num_nodes = (LO)nx * nx;
  RCP<Matrix> A;
  RCP<Map> map_for_read = MapFactory::Build(TestHelpers::Parameters::getLib(), num_nodes, 0, comm);
  try {
    A = Xpetra::IO<SC, LO, GO, NO>::Read("TestMatrices/aniso2dx.mat", map_for_read);
  } catch (...) {
    // Sometimes the matrix reader just fails.
    return;
  };

  std::vector<std::string> test_matrices = {"TestMatrices/aniso2dx.mat", "TestMatrices/aniso2dy.mat", "TestMatrices/iso2d.mat"};
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;

  const std::vector<LO> AGG0_NODES  = {0, 1};
  const std::vector<LO> AGG1_NODES  = {0, nx};
  const std::vector<LO> AGG2_NODES  = {0, nx, nx + 1};
  const std::vector<LO> AGG3_NODES  = {0, 1, nx, nx + 1};
  const std::vector<LO> AGG4_NODES  = {0, 1, 2};
  const std::vector<LO> AGG5_NODES  = {0, nx, 2 * nx};
  const std::vector<LO> AGG6_NODES  = {0, nx, nx + 1, 2 * nx};
  const std::vector<LO> AGG7_NODES  = {1, nx, nx + 1, nx + 2};
  const std::vector<LO> AGG8_NODES  = {1, nx, nx + 1, nx + 2, 2 * nx + 1};
  const std::vector<LO> AGG9_NODES  = {0, 1, 2, nx, nx + 1, nx + 2};
  const std::vector<LO> AGG10_NODES = {0, 1, nx, nx + 1, 2 * nx, 2 * nx + 1};
  const std::vector<LO> AGG11_NODES = {0, 1, 2, nx, nx + 1, nx + 2, 2 * nx, 2 * nx + 1, 2 * nx + 2};
  const std::vector<LO> AGG12_NODES = {0, 1, nx, nx + 1, nx + 2};
  const std::vector<LO> AGG13_NODES = {0, 1, nx, nx + 1, 2 * nx};
  const std::vector<LO> AGG14_NODES = {0, 1, 2, nx, nx + 1, nx + 2, 2 * nx};
  const std::vector<LO> AGG15_NODES = {0, 1, 2, nx, nx + 1, nx + 2, 2 * nx + 1};
  const std::vector<LO> AGG16_NODES = {0, 1, 2, nx, nx + 1, 2 * nx, 2 * nx + 1};
  const std::vector<LO> AGG17_NODES = {0, 1, nx, nx + 1, nx + 2, 2 * nx, 2 * nx + 1};
  const std::vector<LO> AGG18_NODES = {1, nx, nx + 1, nx + 2, 2 * nx + 1, 2 * nx + 2};

  const std::vector<std::vector<LO>> AGG_TYPES = {AGG0_NODES, AGG1_NODES, AGG2_NODES,
                                                  AGG3_NODES, AGG4_NODES, AGG5_NODES,
                                                  AGG6_NODES, AGG7_NODES, AGG8_NODES,
                                                  AGG9_NODES, AGG10_NODES, AGG11_NODES,
                                                  AGG12_NODES, AGG13_NODES, AGG14_NODES,
                                                  AGG15_NODES, AGG16_NODES, AGG17_NODES,
                                                  AGG18_NODES};

  const std::vector<MT> AGG_QUALITIES = {2002., 2.002, 2670.000811,
                                         2001.9998, 4003.999998, 4.004,
                                         3003.666728, 4003.999198, 4003.998397,
                                         4003.999198, 2001.999599, 4003.999198,
                                         4004.666568, 2403.067729, 4886.458381,
                                         4003.998397, 4005.666904, 4004.665757,
                                         4004.665841};

  for (int i = 0; i < 3; i++) {
    out << "Aggregating Matrix " << i << "..." << std::endl;

    Level level;
    TestHelpers::TestFactory<Scalar, LO, GO, NO>::createSingleLevelHierarchy(level);

    RCP<Matrix> A_agg = Xpetra::IO<SC, LO, GO, NO>::Read(test_matrices.at(i), map_for_read);

    level.Set("A", A_agg);

    RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory);
    RCP<CoalesceDropFactory> dropFact  = rcp(new CoalesceDropFactory);
    dropFact->SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.1));
    dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

    RCP<UncoupledAggregationFactory> aggFactory = rcp(new UncoupledAggregationFactory);
    aggFactory->SetFactory("Graph", dropFact);

    out << "Coalesce Drop: " << dropFact.get() << std::endl;
    out << "Aggregation: " << aggFactory.get() << std::endl;

    level.Request("Aggregates", aggFactory.get());

    out << "Building aggregates..." << std::endl;

    aggFactory->Build(level);

    out << "Built aggregates! Getting them...\n"
        << std::endl;

    RCP<Aggregates> aggs = level.Get<RCP<Aggregates>>("Aggregates", aggFactory.get());

    out << "Got aggregates! Computing quality estimate...\n"
        << std::endl;

    level.print(out, Teuchos::VERB_EXTREME);

    level.Set("A", A);

    RCP<CoarseMapFactory> coarseMapFact = rcp(new CoarseMapFactory);
    coarseMapFact->SetFactory("Aggregates", aggFactory);

    AggregateQualityEstimateFactory aggQualityEstimateFactory;
    aggQualityEstimateFactory.SetFactory("Aggregates", aggFactory);
    aggQualityEstimateFactory.SetFactory("CoarseMap", coarseMapFact);

    aggQualityEstimateFactory.SetParameter("aggregate qualities: check symmetry", Teuchos::ParameterEntry(false));
    aggQualityEstimateFactory.SetParameter("aggregate qualities: good aggregate threshold", Teuchos::ParameterEntry(100.0));
    aggQualityEstimateFactory.SetParameter("aggregate qualities: file output", Teuchos::ParameterEntry(false));

    level.Request("AggregateQualities", &aggQualityEstimateFactory);
    // level.Request(aggQualityEstimateFactory);

    aggQualityEstimateFactory.Build(level);

    RCP<MultiVectorDouble> aggQualities = level.Get<RCP<MultiVectorDouble>>("AggregateQualities", &aggQualityEstimateFactory);

    ArrayRCP<const MT> aggQualitiesLocalData = aggQualities->getData(0);

    ArrayRCP<LO> aggSortedVertices, aggsToIndices, aggSizes;
    AggregateQualityEstimateFactory::ConvertAggregatesData(aggs, aggSortedVertices, aggsToIndices, aggSizes);

    for (size_t j = 0; j < aggQualities->getLocalLength(); ++j) {
      std::vector<LO> nodes;
      for (LO k = 0; k < aggSizes[j]; ++k) {
        nodes.push_back(aggSortedVertices[aggsToIndices[j] + k]);
      }

      std::sort(nodes.begin(), nodes.end());

      bool onBoundary = false;

      for (size_t k = 0; k < nodes.size(); ++k) {
        if (nodes[k] % nx == 0 || nodes[k] % nx == nx - 1 || nodes[k] / nx == 0 || nodes[k] / nx == nx - 1) {
          onBoundary = true;
          break;
        }
      }

      if (onBoundary) continue;

      bool assert_performed = false;

      for (size_t agg_id = 0; agg_id < AGG_TYPES.size(); ++agg_id) {
        if (AGG_TYPES[agg_id].size() != nodes.size()) continue;

        const std::vector<LO>& unflipped_agg = AGG_TYPES[agg_id];

        for (int flip_id = 0; flip_id < 4; ++flip_id) {
          std::vector<LO> flipped_agg;

          switch (flip_id) {
            case 0:
              flipped_agg = unflipped_agg;
              break;
            case 1:
              flipped_agg = flip_agg_horizontal(unflipped_agg, nx);
              break;
            case 2:
              flipped_agg = flip_agg_vertical(unflipped_agg, nx);
              break;
            case 3:
              flipped_agg = flip_agg_horizontal(flip_agg_vertical(unflipped_agg, nx), nx);
              break;
          }

          LO difference = nodes[0] - flipped_agg[0];

          bool aggFound = true;

          for (size_t k = 1; k < nodes.size(); ++k) {
            if (difference != nodes[k] - flipped_agg[k]) {
              aggFound = false;
              break;
            }

            if (nodes[k] % nx == 0 || nodes[k] % nx == nx - 1 || nodes[k] / nx == 0 || nodes[k] / nx == nx - 1) {
              aggFound = false;
              break;
            }
          }

          if (!aggFound) continue;

          TEST_FLOATING_EQUALITY(aggQualitiesLocalData[j], AGG_QUALITIES[agg_id], 1e-3);

          assert_performed = true;
          break;
        }

        if (assert_performed) break;
      }
    }
  }

}  // Anisotropic Diffusion 2D test

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(AggregateQualityEstimateFactory, ConvectionDiffusion2D, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;
  typedef Xpetra::MultiVector<MT, LO, GO, Node> MultiVectorDouble;
  typedef MueLu::AggregateQualityEstimateFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> AggregateQualityEstimateFactory;

  out << "version: " << MueLu::Version() << std::endl;

  typedef typename Teuchos::ScalarTraits<Scalar> TST;

  // Don't test for complex - matrix reader won't work
  if (TST::isComplex) {
    success = true;
    return;
  }

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();

  LO nx = (LO)20;

  LO num_nodes          = (LO)nx * nx;
  RCP<Map> map_for_read = MapFactory::Build(TestHelpers::Parameters::getLib(), num_nodes, 0, comm);
  RCP<Matrix> A         = Xpetra::IO<SC, LO, GO, NO>::Read("TestMatrices/cd2dx.mat", map_for_read);

  std::vector<std::string> test_matrices = {"TestMatrices/cd2dx.mat", "TestMatrices/cd2dy.mat", "TestMatrices/aniso2dx.mat", "TestMatrices/aniso2dy.mat", "TestMatrices/iso2d.mat"};

  const std::vector<LO> AGG0_NODES  = {0, 1};
  const std::vector<LO> AGG1_NODES  = {0, nx};
  const std::vector<LO> AGG2_NODES  = {0, nx, nx + 1};
  const std::vector<LO> AGG3_NODES  = {0, 1, nx, nx + 1};
  const std::vector<LO> AGG4_NODES  = {0, 1, 2};
  const std::vector<LO> AGG5_NODES  = {0, nx, 2 * nx};
  const std::vector<LO> AGG6_NODES  = {0, nx, nx + 1, 2 * nx};
  const std::vector<LO> AGG7_NODES  = {1, nx, nx + 1, nx + 2};
  const std::vector<LO> AGG8_NODES  = {1, nx, nx + 1, nx + 2, 2 * nx + 1};
  const std::vector<LO> AGG9_NODES  = {0, 1, 2, nx, nx + 1, nx + 2};
  const std::vector<LO> AGG10_NODES = {0, 1, nx, nx + 1, 2 * nx, 2 * nx + 1};
  const std::vector<LO> AGG11_NODES = {0, 1, 2, nx, nx + 1, nx + 2, 2 * nx, 2 * nx + 1, 2 * nx + 2};
  const std::vector<LO> AGG12_NODES = {0, 1, nx, nx + 1, nx + 2};
  const std::vector<LO> AGG13_NODES = {0, 1, nx, nx + 1, 2 * nx};
  const std::vector<LO> AGG14_NODES = {0, 1, 2, nx, nx + 1, nx + 2, 2 * nx};
  const std::vector<LO> AGG15_NODES = {0, 1, 2, nx, nx + 1, nx + 2, 2 * nx + 1};
  const std::vector<LO> AGG16_NODES = {0, 1, 2, nx, nx + 1, 2 * nx, 2 * nx + 1};
  const std::vector<LO> AGG17_NODES = {0, 1, nx, nx + 1, nx + 2, 2 * nx, 2 * nx + 1};
  const std::vector<LO> AGG18_NODES = {1, nx, nx + 1, nx + 2, 2 * nx + 1, 2 * nx + 2};

  const std::vector<std::vector<LO>> AGG_TYPES = {AGG0_NODES, AGG1_NODES, AGG2_NODES,
                                                  AGG3_NODES, AGG4_NODES, AGG5_NODES,
                                                  AGG6_NODES, AGG7_NODES, AGG8_NODES,
                                                  AGG9_NODES, AGG10_NODES, AGG11_NODES,
                                                  AGG12_NODES, AGG13_NODES, AGG14_NODES,
                                                  AGG15_NODES, AGG16_NODES, AGG17_NODES,
                                                  AGG18_NODES};

  const std::vector<MT> AGG_QUALITIES = {2.008365, 480.190476, 640.925507,
                                         480.190465, 4.016730, 960.380952,
                                         960.380906, 720.956363, 960.380860,
                                         480.190453, 960.380906, 960.380860,
                                         576.902089, 961.050518, 962.055240,
                                         961.050432, 1172.649129, 960.380860,
                                         961.050785};

  for (size_t i = 0; i < test_matrices.size(); i++) {
    out << "Aggregating Matrix " << i << "..." << std::endl;

    Level level;
    TestHelpers::TestFactory<Scalar, LO, GO, NO>::createSingleLevelHierarchy(level);

    RCP<Matrix> A_agg = Xpetra::IO<SC, LO, GO, NO>::Read(test_matrices.at(i), map_for_read);

    level.Set("A", A_agg);

    RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory);
    RCP<CoalesceDropFactory> dropFact  = rcp(new CoalesceDropFactory);
    dropFact->SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.1));
    dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

    RCP<UncoupledAggregationFactory> aggFactory = rcp(new UncoupledAggregationFactory);
    aggFactory->SetFactory("Graph", dropFact);

    out << "Coalesce Drop: " << dropFact.get() << std::endl;
    out << "Aggregation: " << aggFactory.get() << std::endl;

    level.Request("Aggregates", aggFactory.get());

    out << "Building aggregates..." << std::endl;

    aggFactory->Build(level);

    out << "Built aggregates! Getting them...\n"
        << std::endl;

    RCP<Aggregates> aggs = level.Get<RCP<Aggregates>>("Aggregates", aggFactory.get());

    out << "Got aggregates! Computing quality estimate...\n"
        << std::endl;

    level.print(out, Teuchos::VERB_EXTREME);

    level.Set("A", A);

    RCP<CoarseMapFactory> coarseMapFact = rcp(new CoarseMapFactory);
    coarseMapFact->SetFactory("Aggregates", aggFactory);

    AggregateQualityEstimateFactory aggQualityEstimateFactory;
    aggQualityEstimateFactory.SetFactory("Aggregates", aggFactory);
    aggQualityEstimateFactory.SetFactory("CoarseMap", coarseMapFact);

    aggQualityEstimateFactory.SetParameter("aggregate qualities: check symmetry", Teuchos::ParameterEntry(true));
    aggQualityEstimateFactory.SetParameter("aggregate qualities: good aggregate threshold", Teuchos::ParameterEntry(100.0));
    aggQualityEstimateFactory.SetParameter("aggregate qualities: file output", Teuchos::ParameterEntry(false));

    level.Request("AggregateQualities", &aggQualityEstimateFactory);
    // level.Request(aggQualityEstimateFactory);

    aggQualityEstimateFactory.Build(level);

    RCP<MultiVectorDouble> aggQualities = level.Get<RCP<MultiVectorDouble>>("AggregateQualities", &aggQualityEstimateFactory);

    ArrayRCP<const MT> aggQualitiesLocalData = aggQualities->getData(0);

    ArrayRCP<LO> aggSortedVertices, aggsToIndices, aggSizes;
    AggregateQualityEstimateFactory::ConvertAggregatesData(aggs, aggSortedVertices, aggsToIndices, aggSizes);

    for (size_t j = 0; j < aggQualities->getLocalLength(); ++j) {
      std::vector<LO> nodes;
      for (LO k = 0; k < aggSizes[j]; ++k) {
        nodes.push_back(aggSortedVertices[aggsToIndices[j] + k]);
      }

      std::sort(nodes.begin(), nodes.end());

      bool onBoundary = false;

      for (size_t k = 0; k < nodes.size(); ++k) {
        if (nodes[k] % nx == 0 || nodes[k] % nx == nx - 1 || nodes[k] / nx == 0 || nodes[k] / nx == nx - 1) {
          onBoundary = true;
          break;
        }
      }

      if (onBoundary) continue;

      bool assert_performed = false;

      for (size_t agg_id = 0; agg_id < AGG_TYPES.size(); ++agg_id) {
        if (AGG_TYPES[agg_id].size() != nodes.size()) continue;

        const std::vector<LO>& unflipped_agg = AGG_TYPES[agg_id];

        for (int flip_id = 0; flip_id < 4; ++flip_id) {
          std::vector<LO> flipped_agg;

          switch (flip_id) {
            case 0:
              flipped_agg = unflipped_agg;
              break;
            case 1:
              flipped_agg = flip_agg_horizontal(unflipped_agg, nx);
              break;
            case 2:
              flipped_agg = flip_agg_vertical(unflipped_agg, nx);
              break;
            case 3:
              flipped_agg = flip_agg_horizontal(flip_agg_vertical(unflipped_agg, nx), nx);
              break;
          }

          LO difference = nodes[0] - flipped_agg[0];

          bool aggFound = true;

          for (size_t k = 1; k < nodes.size(); ++k) {
            if (difference != nodes[k] - flipped_agg[k]) {
              aggFound = false;
              break;
            }

            if (nodes[k] % nx == 0 || nodes[k] % nx == nx - 1 || nodes[k] / nx == 0 || nodes[k] / nx == nx - 1) {
              aggFound = false;
              break;
            }
          }

          if (!aggFound) continue;

          TEST_FLOATING_EQUALITY(aggQualitiesLocalData[j], AGG_QUALITIES[agg_id], 1e-3);

          assert_performed = true;
          break;
        }

        if (assert_performed) break;
      }
    }
  }

}  // Convection Diffusion 2D test

#define MUELU_ETI_GROUP(Scalar, LO, GO, Node)                                                              \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(AggregateQualityEstimateFactory, Constructor, Scalar, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(AggregateQualityEstimateFactory, Poisson2D, Scalar, LO, GO, Node)
//  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(AggregateQualityEstimateFactory,AnisotropicDiffusion2D,Scalar,LO,GO,Node)

//  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(AggregateQualityEstimateFactory,ConvectionDiffusion2D,Scalar,LO,GO,Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
