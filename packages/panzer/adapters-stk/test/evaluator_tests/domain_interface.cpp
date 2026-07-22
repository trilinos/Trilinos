// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Panzer
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_WorksetFactory.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_WorksetDescriptor.hpp"
#include "Panzer_WorksetNeeds.hpp"
#include "Panzer_IntegrationDescriptor.hpp"
#include "Panzer_BasisDescriptor.hpp"
#include "Panzer_Evaluator_DomainInterface.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_Evaluator_WithBaseImpl.hpp"

// Teuchos
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_UnitTestHarness.hpp"

class MockEvaluator : public panzer::EvaluatorWithBaseImpl<panzer::Traits> {
  int expectedStartIndex_;
  int expectedEndIndex_;

public:

  void evaluateFields(panzer::Traits::EvalData workset) override
  {
    // Temporary vars needed by test macros
    auto& out = std::cout;
    int success = 0;

    TEST_EQUALITY(this->cellStartIndex(workset), expectedStartIndex_);
    TEST_EQUALITY(this->cellEndIndex(workset), expectedEndIndex_);
    TEST_EQUALITY(success, 0);
  }

  void setExpectedIndices(const int& startIndex,
                          const int& endIndex)
  {
    expectedStartIndex_ = startIndex;
    expectedEndIndex_ = endIndex;
  }
};


TEUCHOS_UNIT_TEST(domain_interface, base)
{
  using namespace Teuchos;
  using namespace panzer;
  using namespace panzer_stk;

  RCP<MpiComm<int>> comm = rcp(new MpiComm<int>(MPI_COMM_WORLD));

  const int myRank = comm->getRank();
  const int numProcs = comm->getSize();
  TEUCHOS_ASSERT(numProcs == 2);

  Teuchos::RCP<panzer_stk::STK_Interface> mesh;
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Blocks",1);
    pl->set("Y Blocks",1);
    pl->set("X Elements",3);
    pl->set("Y Elements",4);
    panzer_stk::SquareQuadMeshFactory factory;
    factory.setParameterList(pl);
    mesh = factory.buildUncommitedMesh(MPI_COMM_WORLD);
    factory.completeMeshConstruction(*mesh,MPI_COMM_WORLD);
  }

  IntegrationDescriptor integrationDescriptor(2,IntegrationDescriptor::VOLUME);
  BasisDescriptor basisDescriptor(1,"HGrad");
  WorksetNeeds worksetNeeds;
  worksetNeeds.addBasis(basisDescriptor);
  worksetNeeds.addIntegrator(integrationDescriptor);

  WorksetDescriptor worksetDescriptor("eblock-0_0", WorksetSizeType::ALL_ELEMENTS,true,false);

  WorksetFactory worksetFactory(mesh);

  auto worksets = worksetFactory.getWorksets(worksetDescriptor,
                                             worksetNeeds);

  TEST_EQUALITY(worksets->size(),1);

  const auto& workset = (*worksets)[0];

  MockEvaluator e;

  // These values are known by construction of the
  // QuadMeshFactory. Changing the parallel decomposition in the mesh
  // facotry may break these values.

  if (myRank == 0) {

    e.setDomain(DomainEvaluator::OWNED);
    e.setExpectedIndices(0,8);
    e.evaluateFields(workset);

    e.setDomain(DomainEvaluator::GHOST);
    e.setExpectedIndices(8,12);
    e.evaluateFields(workset);

    e.setDomain(DomainEvaluator::REAL);
    e.setExpectedIndices(0,12);
    e.evaluateFields(workset);

    e.setDomain(DomainEvaluator::VIRTUAL);
    e.setExpectedIndices(12,20);
    e.evaluateFields(workset);

    e.setDomain(DomainEvaluator::EXTERNAL);
    e.setExpectedIndices(8,20);
    e.evaluateFields(workset);

    e.setDomain(DomainEvaluator::ALL);
    e.setExpectedIndices(0,20);
    e.evaluateFields(workset);

  } else if (myRank == 1) {

    e.setDomain(DomainEvaluator::OWNED);
    e.setExpectedIndices(0,4);
    e.evaluateFields(workset);

    e.setDomain(DomainEvaluator::GHOST);
    e.setExpectedIndices(4,8);
    e.evaluateFields(workset);

    e.setDomain(DomainEvaluator::REAL);
    e.setExpectedIndices(0,8);
    e.evaluateFields(workset);

    e.setDomain(DomainEvaluator::VIRTUAL);
    e.setExpectedIndices(8,14);
    e.evaluateFields(workset);

    e.setDomain(DomainEvaluator::EXTERNAL);
    e.setExpectedIndices(4,14);
    e.evaluateFields(workset);

    e.setDomain(DomainEvaluator::ALL);
    e.setExpectedIndices(0,14);
    e.evaluateFields(workset);

  }

}
