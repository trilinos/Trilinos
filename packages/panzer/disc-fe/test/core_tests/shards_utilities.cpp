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

#include "Panzer_Shards_Utilities.hpp"

using std::cout;
using std::endl;
using panzer::getLocalSideIndexFromGlobalNodeList;

namespace panzer {
  TEUCHOS_UNIT_TEST(shards_utilities, quad_test)
  {

    const shards::CellTopology
      quad(shards::getCellTopologyData<shards::Quadrilateral<4> >());
    const CellTopologyData* quad_data = quad.getCellTopologyData();

    //cout << quad << endl;

    unsigned local_face = 0;
    std::vector<unsigned> element_gids(quad.getNodeCount());
    const unsigned num_nodes_per_side =
      quad_data->side[0].topology->node_count;
    std::vector<unsigned> edge_gids(num_nodes_per_side);

    /* GIDs

          2 *----* 8
            |    |
            |    |
         10 *----* 30

    */
    element_gids[0] = 10;
    element_gids[1] = 30;
    element_gids[2] = 8;
    element_gids[3] = 2;

    edge_gids[0] = 8;
    edge_gids[1] = 2;
    local_face = getLocalSideIndexFromGlobalNodeList(element_gids, edge_gids,
						     quad);
    TEST_EQUALITY(local_face, 2);

    edge_gids[0] = 10;
    edge_gids[1] = 30;
    local_face = getLocalSideIndexFromGlobalNodeList(element_gids, edge_gids,
						     quad);
    TEST_EQUALITY(local_face, 0);

    edge_gids[0] = 2;
    edge_gids[1] = 10;
    local_face = getLocalSideIndexFromGlobalNodeList(element_gids, edge_gids,
						     quad);
    TEST_EQUALITY(local_face, 3);

    edge_gids[0] = 30;
    edge_gids[1] = 8;
    local_face = getLocalSideIndexFromGlobalNodeList(element_gids, edge_gids,
						     quad);
    TEST_EQUALITY(local_face, 1);

    // Edge with order reversed
    edge_gids[0] = 8;
    edge_gids[1] = 30;
    local_face = getLocalSideIndexFromGlobalNodeList(element_gids, edge_gids,
						     quad);
    TEST_EQUALITY(local_face, 1);

    // Give it an incorrect edge and make sure it throws
    bool caught_throw = false;
    try {

      edge_gids[0] = 8;
      edge_gids[1] = 10;
      local_face = getLocalSideIndexFromGlobalNodeList(element_gids, edge_gids,
						       quad);

    }
    catch (const std::runtime_error& ) {
      caught_throw = true;
      //cout << e.what() << std::endl;
    }

    TEST_ASSERT(caught_throw);

  }

  TEUCHOS_UNIT_TEST(shards_utilities, hex_test)
  {

    const shards::CellTopology
      hex(shards::getCellTopologyData<shards::Hexahedron<8> >());
    const CellTopologyData* hex_data = hex.getCellTopologyData();

    //cout << hex << endl;

    unsigned local_face = 0;
    std::vector<unsigned> element_gids(hex.getNodeCount());
    const unsigned num_nodes_per_side =
      hex_data->side[0].topology->node_count;
    std::vector<unsigned> edge_gids(num_nodes_per_side);

    // local to global mapping for element
    // 0 --> 10
    // 1 --> 18
    // 2 --> 5
    // 3 --> 7
    // 4 --> 1
    // 5 --> 30
    // 6 --> 20
    // 7 --> 17

    element_gids[0] = 10;
    element_gids[1] = 18;
    element_gids[2] = 5;
    element_gids[3] = 7;
    element_gids[4] = 1;
    element_gids[5] = 30;
    element_gids[6] = 20;
    element_gids[7] = 17;

    edge_gids[0] = 10;
    edge_gids[1] = 18;
    edge_gids[2] = 30;
    edge_gids[3] = 1;
    local_face = getLocalSideIndexFromGlobalNodeList(element_gids, edge_gids,
						     hex);
    TEST_EQUALITY(local_face, 0);

    edge_gids[0] = 18;
    edge_gids[1] = 5;
    edge_gids[2] = 20;
    edge_gids[3] = 30;
    local_face = getLocalSideIndexFromGlobalNodeList(element_gids, edge_gids,
						     hex);
    TEST_EQUALITY(local_face, 1);

    edge_gids[0] = 5;
    edge_gids[1] = 7;
    edge_gids[2] = 17;
    edge_gids[3] = 20;
    local_face = getLocalSideIndexFromGlobalNodeList(element_gids, edge_gids,
						     hex);
    TEST_EQUALITY(local_face, 2);

    edge_gids[0] = 10;
    edge_gids[1] = 1;
    edge_gids[2] = 17;
    edge_gids[3] = 7;
    local_face = getLocalSideIndexFromGlobalNodeList(element_gids, edge_gids,
						     hex);
    TEST_EQUALITY(local_face, 3);

    edge_gids[0] = 10;
    edge_gids[1] = 7;
    edge_gids[2] = 5;
    edge_gids[3] = 18;
    local_face = getLocalSideIndexFromGlobalNodeList(element_gids, edge_gids,
						     hex);
    TEST_EQUALITY(local_face, 4);

    edge_gids[0] = 1;
    edge_gids[1] = 30;
    edge_gids[2] = 20;
    edge_gids[3] = 17;
    local_face = getLocalSideIndexFromGlobalNodeList(element_gids, edge_gids,
						     hex);
    TEST_EQUALITY(local_face, 5);

    // Valid edge with bad ordering
    edge_gids[0] = 30;
    edge_gids[1] = 20;
    edge_gids[2] = 1;
    edge_gids[3] = 17;
    local_face = getLocalSideIndexFromGlobalNodeList(element_gids, edge_gids,
						     hex);
    TEST_EQUALITY(local_face, 5);

    // Give it an incorrect edge and make sure it throws
    bool caught_throw = false;
    try {

      edge_gids[0] = 30;
      edge_gids[1] = 7;
      edge_gids[2] = 18;
      edge_gids[3] = 17;
      local_face = getLocalSideIndexFromGlobalNodeList(element_gids, edge_gids,
						       hex);

    }
    catch (const std::runtime_error& ) {
      caught_throw = true;
      //cout << e.what() << std::endl;
    }

    TEST_ASSERT(caught_throw);

  }

}
