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

using Teuchos::RCP;
using Teuchos::rcp;

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Panzer_STK_Version.hpp"
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_CubeHexMeshFactory.hpp"
#include "Panzer_Workset_Builder.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_GlobalData.hpp"
#include "Panzer_BC.hpp"

#include "user_app_EquationSetFactory.hpp"

namespace panzer {

void testInitialization(const Teuchos::RCP<Teuchos::ParameterList>& ipb,
                        std::vector<panzer::BC>& bcs, const int integration_order);

void testIpMatch(const panzer::WorksetDetails& d0, const panzer::WorksetDetails& d1,
                 const index_t num_cells, Teuchos::FancyOStream& out, bool& success);

/*
namespace {
std::string prws (const panzer::Workset& w) {
  std::stringstream ss;
  ss << w;
  for (size_t i = 0; i < w.numDetails(); ++i) {
    ss << "details " << i << ":\n";
    const panzer::WorksetDetails& d = w.details(i);
    for (size_t j = 0; j < d.int_rules.size(); ++j) {
      const panzer::IntegrationValues2<double>& ir = *d.int_rules[j];
      const int num_ip = ir.ip_coordinates.extent_int(1);
      const size_t num_dim = ir.ip_coordinates.extent(2);
      ss << "int_rule " << j << ":\n";
      ss << "cub_points:\n";
      for (int ip = 0; ip < num_ip; ++ip) {
        for (size_t dim = 0; dim < num_dim; ++dim)
          ss << " " << ir.cub_points(ip, dim);
        ss << "\n";
      }
      ss << "side_cub_points:\n";
      for (int ip = 0; ip < num_ip; ++ip) {
        for (size_t dim = 0; dim < ir.side_cub_points.extent(1); ++dim)
          ss << " " << ir.side_cub_points(ip, dim);
        ss << "\n";
      }
      ss << "cub_weights:\n";
      for (int ip = 0; ip < num_ip; ++ip)
        ss << " " << ir.cub_weights(ip);
      ss << "\n";
      for (size_t cell = 0; cell < w.num_cells; ++cell) {
        ss << "cell " << cell << ":\n";
        ss << "ip_coordinates:\n";
        for (int ip = 0; ip < num_ip; ++ip) {
          for (size_t dim = 0; dim < num_dim; ++dim)
            ss << " " << ir.ip_coordinates(cell, ip, dim);
          ss << "\n";
        }
        ss << "jac_det:\n";
        for (int ip = 0; ip < num_ip; ++ip)
          ss << " " << ir.jac_det(cell, ip);
        ss << "\n";
        ss << "weighted_measure:\n";
        for (int ip = 0; ip < num_ip; ++ip)
          ss << " " << ir.weighted_measure(cell, ip);
        ss << "\n";
        ss << "jac:\n";
        for (int ip = 0; ip < num_ip; ++ip) {
          for (size_t d0 = 0; d0 < num_dim; ++d0) {
            for (size_t d1 = 0; d1 < num_dim; ++d1)
              ss << " " << ir.jac(cell, ip, d0, d1);
            ss << " |";
          }
          ss << "\n";
        }
        ss << "jac_inv:\n";
        for (int ip = 0; ip < num_ip; ++ip) {
          for (size_t d0 = 0; d0 < num_dim; ++d0) {
            for (size_t d1 = 0; d1 < num_dim; ++d1)
              ss << " " << ir.jac_inv(cell, ip, d0, d1);
            ss << " |";
          }
          ss << "\n";
        }
      }
    }
  }
  return ss.str();
}
}
*/

TEUCHOS_UNIT_TEST(workset_builder, stk_edge)
{

  RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
  pl->set("X Blocks",2);
  pl->set("Y Blocks",1);
  pl->set("Z Blocks",1);
  pl->set("X Elements",2);
  pl->set("Y Elements",2);
  pl->set("Z Elements",1);
  pl->set("Build Interface Sidesets", true);

  int myRank=0;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  panzer_stk::CubeHexMeshFactory factory;
  factory.setParameterList(pl);
  RCP<panzer_stk::STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
  mesh->writeToExodus("test.exo");

  std::vector<std::string> element_blocks;
  mesh->getElementBlockNames(element_blocks);
  const std::size_t workset_size = 20;

  // Use a high integration order to test ordering of the fields in
  // IntegrationValues2.
  const int default_integration_order = 4;
  Teuchos::RCP<Teuchos::ParameterList> ipb = Teuchos::parameterList("Physics Blocks");
  std::vector<panzer::BC> bcs;
  testInitialization(ipb, bcs, default_integration_order);

  // build physics blocks
  //////////////////////////////////////////////////////////////
  std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physicsBlocks;
  {
    Teuchos::RCP<user_app::MyFactory> eqset_factory = Teuchos::rcp(new user_app::MyFactory);

    std::map<std::string,std::string> block_ids_to_physics_ids;
    block_ids_to_physics_ids["eblock-0_0_0"] = "test physics";
    block_ids_to_physics_ids["eblock-1_0_0"] = "test physics";

    std::map<std::string,Teuchos::RCP<const shards::CellTopology> > block_ids_to_cell_topo;
    block_ids_to_cell_topo["eblock-0_0_0"] = mesh->getCellTopology("eblock-0_0_0");
    block_ids_to_cell_topo["eblock-1_0_0"] = mesh->getCellTopology("eblock-1_0_0");

    Teuchos::RCP<panzer::GlobalData> gd = panzer::createGlobalData();

    panzer::buildPhysicsBlocks(block_ids_to_physics_ids,
                               block_ids_to_cell_topo,
                               ipb,
                               default_integration_order,
                               workset_size,
                               eqset_factory,
                               gd,
                               false,
                               physicsBlocks);
  }

  const int eb_idxs[2][2] = {{0,1}, {1,0}};
  for (int ebi = 0; ebi < 2; ++ebi) {
    std::string sideset = "vertical_0";
    Teuchos::RCP<const panzer::PhysicsBlock> pb_a = panzer::findPhysicsBlock(element_blocks[eb_idxs[ebi][0]], physicsBlocks);
    Teuchos::RCP<const panzer::PhysicsBlock> pb_b = panzer::findPhysicsBlock(element_blocks[eb_idxs[ebi][1]], physicsBlocks);
    Teuchos::RCP<std::map<unsigned,panzer::Workset> > worksets = panzer_stk::buildBCWorksets(
      *mesh,
      pb_a->getWorksetNeeds(),pb_a->elementBlockID(),
      pb_b->getWorksetNeeds(),pb_b->elementBlockID(),
      sideset);

    TEST_EQUALITY(worksets->size(), 1);
    panzer::Workset& workset = worksets->begin()->second;
    //std::cout << "prws(workset) for element block index " << ebi << ":\n" << prws(workset) << "\n";
    TEST_EQUALITY(workset.numDetails(), 2);
    testIpMatch(workset.details(0), workset.details(1), workset.num_cells, out, success);
  }
}

void testInitialization(const Teuchos::RCP<Teuchos::ParameterList>& ipb,
                        std::vector<panzer::BC>& bcs, const int integration_order)
{
  // Physics block
  Teuchos::ParameterList& physics_block = ipb->sublist("test physics");
  {
    Teuchos::ParameterList& p = physics_block.sublist("a");
    p.set("Type","Energy");
    p.set("Prefix","");
    p.set("Model ID","solid");
    p.set("Basis Type","HGrad");
    p.set("Basis Order",2);
    p.set("Integration Order",integration_order);
  }
  {
    Teuchos::ParameterList& p = physics_block.sublist("b");
    p.set("Type","Energy");
    p.set("Prefix","ION_");
    p.set("Model ID","ion solid");
    p.set("Basis Type","HGrad");
    p.set("Basis Order",1);
    p.set("Integration Order",integration_order);
  }

  {
    std::size_t bc_id = 0;
    panzer::BCType neumann = BCT_Dirichlet;
    std::string sideset_id = "left";
    std::string element_block_id = "eblock-0_0_0";
    std::string dof_name = "UX";
    std::string strategy = "Constant";
    double value = 5.0;
    Teuchos::ParameterList p;
    p.set("Value",value);
    panzer::BC bc(bc_id, neumann, sideset_id, element_block_id, dof_name,
                  strategy, p);
    bcs.push_back(bc);
  }
  {
    std::size_t bc_id = 0;
    panzer::BCType neumann = BCT_Dirichlet;
    std::string sideset_id = "right";
    std::string element_block_id = "eblock-1_0_0";
    std::string dof_name = "UX";
    std::string strategy = "Constant";
    double value = 5.0;
    Teuchos::ParameterList p;
    p.set("Value",value);
    panzer::BC bc(bc_id, neumann, sideset_id, element_block_id, dof_name,
                  strategy, p);
    bcs.push_back(bc);
  }
  {
    std::size_t bc_id = 0;
    panzer::BCType neumann = BCT_Dirichlet;
    std::string sideset_id = "top";
    std::string element_block_id = "eblock-1_0_0";
    std::string dof_name = "UX";
    std::string strategy = "Constant";
    double value = 5.0;
    Teuchos::ParameterList p;
    p.set("Value",value);
    panzer::BC bc(bc_id, neumann, sideset_id, element_block_id, dof_name,
                  strategy, p);
    bcs.push_back(bc);
  }
}

void testIpMatch(const panzer::WorksetDetails& d0, const panzer::WorksetDetails& d1,
                 const index_t num_cells, Teuchos::FancyOStream& out, bool& success)
{
  TEST_EQUALITY(d0.int_rules.size(), d1.int_rules.size());
  for (std::size_t iri = 0; iri < d0.int_rules.size(); ++iri) {
    const std::size_t num_ip = d0.int_rules[iri]->cub_points.extent(0),
      num_dim = d0.int_rules[iri]->cub_points.extent(1);

    auto d0_rules_view = d0.int_rules[iri]->ip_coordinates.get_view();
    auto d0_rules_h = Kokkos::create_mirror_view(d0_rules_view);
    Kokkos::deep_copy(d0_rules_h, d0_rules_view);

    auto d1_rules_view = d1.int_rules[iri]->ip_coordinates.get_view();
    auto d1_rules_h = Kokkos::create_mirror_view(d1_rules_view);
    Kokkos::deep_copy(d1_rules_h, d1_rules_view);

    for (index_t cell = 0; cell < num_cells; ++cell)
      for (std::size_t ip = 0; ip < num_ip; ++ip)
        for (std::size_t dim = 0; dim < num_dim; ++dim)
          TEST_EQUALITY(d0_rules_h(cell, ip, dim), d1_rules_h(cell, ip, dim))
  }
}

}
