// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Panzer_STK_Version.hpp"
#include "PanzerAdaptersSTK_config.hpp"

#include "Panzer_WorksetContainer.hpp"
#include "Panzer_IntrepidBasisFactory.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_OrientationsInterface.hpp"

#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_STK_WorksetFactory.hpp"

namespace panzer {
namespace {

void
testApplyOrientations(const bool by_container,
                      const bool by_factory,
                      const bool use_needs,
                      const bool should_fail,
                      Teuchos::FancyOStream & out,
                      bool & success)
{
  // In this test we want to create a simple 2 cell 2D mesh
  // of quads which share a center face.
  // If orientations are not applied, then the HDiv 1 basis will have an orientation (vector)
  // on the interface that points outward from the cell center.
  // If orientations are applied, then both vectors will point in the same direction.

  auto comm = Teuchos::DefaultComm<int>::getComm();
  const int rank = comm->getRank();

  // Only works if we have 2 processors
  TEST_EQUALITY(comm->getSize(), 2);

  // So that I don't spell the name wrong later
  const std::string element_block = "eblock-0_0";

  // Basis has a single vector per face
  const panzer::BasisDescriptor basis(1, "HDiv");

  // Integrator has a single point per face
  const panzer::IntegrationDescriptor integrator(0, panzer::IntegrationDescriptor::SURFACE);

  // Create mesh
  Teuchos::RCP<panzer_stk::STK_Interface> mesh;
  {
    // Each cell is 2x2 (same as reference cell)
    auto pl = Teuchos::rcp(new Teuchos::ParameterList);
    pl->set("X Elements",2);
    pl->set("Y Elements",1);
    pl->set("X0",-2.0);
    pl->set("Xf", 2.0);
    pl->set("Y0",-1.0);
    pl->set("Yf", 1.0);

    panzer_stk::SquareQuadMeshFactory factory;
    factory.setParameterList(pl);
    mesh = factory.buildMesh(MPI_COMM_WORLD);
  }

  // Need a dof manager to act as a global indexer
  Teuchos::RCP<panzer::DOFManager> dof_manager;
  {

    // Build a connectivity manager for the given mesh
    const auto conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    // Initialize the dof manager with the conn manager
    dof_manager = Teuchos::rcp(new panzer::DOFManager(conn_manager,MPI_COMM_WORLD));

    // Build an intrepid basis and a related field pattern for seeding the DOFManager
    {
      const auto & cell_topology = *mesh->getCellTopology(element_block);
      auto intrepid_basis = panzer::createIntrepid2Basis<PHX::Device::execution_space,double,double>(basis.getType(),basis.getOrder(), cell_topology);
      auto field_pattern = Teuchos::rcp(new panzer::Intrepid2FieldPattern(intrepid_basis));

      // Add arbitrary field to dof manager to lock in the field pattern
      dof_manager->addField(element_block, "B", field_pattern);
    }

    // Finalize the dof manager
    dof_manager->buildGlobalUnknowns();
  }

  // Create a needs map
  std::map<std::string, WorksetNeeds> needs_map;
  needs_map[element_block].addBasis(basis);
  needs_map[element_block].addIntegrator(integrator);

  // Build a workset container
  Teuchos::RCP<WorksetContainer> workset_container;
  {
    // Build workset factory
    auto factory = Teuchos::rcp(new panzer_stk::WorksetFactory(mesh));

    if(by_factory)
      factory->setOrientationsInterface(Teuchos::rcp(new OrientationsInterface(dof_manager)));

    if(use_needs)
      workset_container = Teuchos::rcp(new WorksetContainer(factory, needs_map));
    else
      workset_container = Teuchos::rcp(new WorksetContainer(factory, {}));

    if(by_container)
      workset_container->setGlobalIndexer(dof_manager);
  }

  // Grab the one and only workset
  auto worksets = workset_container->getWorksets(WorksetDescriptor(element_block, WorksetSizeType::ALL_ELEMENTS, true, by_container));
  TEST_EQUALITY(worksets->size(), 1);

  const auto & workset = (*worksets)[0];
  const auto & basis_values = workset.getBasisValues(basis,integrator);

  // Check to see if the basis values actually has orientations applied
  if(by_factory or by_container){
    TEST_ASSERT(basis_values.orientationsApplied());
  } else {
    TEST_ASSERT(not basis_values.orientationsApplied());
  }

  // Grab the basis values vector for the given integration points
  const auto basis_vector_device = basis_values.getVectorBasisValues(false).get_static_view();
  auto basis_vector = Kokkos::create_mirror_view(basis_vector_device);
  Kokkos::deep_copy(basis_vector,basis_vector_device);

  // Note that this is for a partitioned workset, which will have 1 owned cells, 1 ghost cell, and 3 virtual cells for a total of 5 cells
  TEST_EQUALITY(basis_vector.extent_int(0),5); // cells
  TEST_EQUALITY(basis_vector.extent_int(1),4); // basis
  TEST_EQUALITY(basis_vector.extent_int(2),4); // points
  TEST_EQUALITY(basis_vector.extent_int(3),2); // vector dim

  // Print out the vector values to see what they are
  out << "Basis Values:\n";
  for(int c=0; c<2; ++c)
    for(int b=0; b<4; ++b)
      for(int p=0; p<4; ++p)
        out << c << ", " << b << ", " << p << ": ["<<basis_vector(c,b,p,0)<<", "<<basis_vector(c,b,p,1)<<"]\n";

  // We are going to make the interface between cell 0 (owned) and 1 (ghost)
  // The basis function associated with that interface is (cell,basis) (0,1) and (1,3)
  // Points are aligned with basis, so (cell,point) are also (0,1) and (1,3)
  // Note that on rank 1, the cells are reversed

  double normal_0[2] = {0}, normal_1[2] = {0};
  if(rank==0){
    normal_0[0] = basis_vector(0,1,1,0);
    normal_0[1] = basis_vector(0,1,1,1);
    normal_1[0] = basis_vector(1,3,3,0);
    normal_1[1] = basis_vector(1,3,3,1);
  } else {
    normal_0[0] = basis_vector(0,3,1,0);
    normal_0[1] = basis_vector(0,3,1,1);
    normal_1[0] = basis_vector(1,1,3,0);
    normal_1[1] = basis_vector(1,1,3,1);
  }

  const double vec[2] = {normal_0[0]+normal_1[0], normal_0[1]+normal_1[1]};

  if(should_fail){
    // The vectors oppose each other (1,0)+(-1,0), so the magnitude should be 0
    TEUCHOS_ASSERT(vec[0]*vec[0]+vec[1]*vec[1] < 1.e-8);
  } else {
    // The vectors are the same (1,0)+(1,0), so the magnitude should be 2
    TEUCHOS_ASSERT(vec[0]*vec[0]+vec[1]*vec[1] - 4 < 1.e-8);
  }

}

}
}

TEUCHOS_UNIT_TEST(apply_orientations, none)
{
  panzer::testApplyOrientations(false, false, true, true, out, success);
}

TEUCHOS_UNIT_TEST(apply_orientations, by_container)
{
  panzer::testApplyOrientations(true, false, true, false, out, success);
}

TEUCHOS_UNIT_TEST(apply_orientations, by_factory)
{
  panzer::testApplyOrientations(false, true, false, false, out, success);
}

TEUCHOS_UNIT_TEST(apply_orientations, by_both)
{
  panzer::testApplyOrientations(true, true, false, false, out, success);
}
