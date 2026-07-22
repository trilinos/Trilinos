// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_STK_CheckSidesetOverlap.hpp"
#include <vector>
#include <algorithm>

namespace panzer_stk {

  /// Returns true if the sidesets overlap.
  bool checkSidesetOverlap(const std::string& side_a_name,
                           const std::string& side_b_name,
                           const panzer_stk::STK_Interface& mesh) {

    const bool print_debug = false;

    // Get the locally owned nodes of sideset a
    std::vector<stk::mesh::EntityId> gids_a;
    {
      std::vector<stk::mesh::Entity> nodes_a;
      stk::mesh::Part* part_a = mesh.getSideset(side_a_name);
      TEUCHOS_TEST_FOR_EXCEPTION(part_a==nullptr,std::runtime_error,
                                 "panzer::checkSidesetOverlap: Unknown side set name \"" << side_a_name << "\"");
      stk::mesh::Selector selector_a = *part_a & mesh.getMetaData()->locally_owned_part();
      const bool sort_by_gid = true;
      stk::mesh::get_selected_entities(selector_a,mesh.getBulkData()->buckets(mesh.getNodeRank()),nodes_a,sort_by_gid);
      // convert the entities to global ids
      gids_a.resize(nodes_a.size());
      size_t i = 0;
      for (auto&& node : nodes_a) {
        gids_a[i] = mesh.getBulkData()->identifier(node);
        ++i;
      }
    }

    // Get all nodes of sideset b (including nodes from all mpi processes)
    std::vector<stk::mesh::EntityId> gids_b;
    {
      std::vector<stk::mesh::Entity> nodes_b;
      stk::mesh::Part* part_b = mesh.getSideset(side_b_name);
      TEUCHOS_TEST_FOR_EXCEPTION(part_b==nullptr,std::runtime_error,
                                 "panzer::checkSidesetOverlap: Unknown side set name \"" << side_b_name << "\"");
      stk::mesh::Selector selector_b = *part_b;
      const bool sort_by_gid = true;
      stk::mesh::get_selected_entities(selector_b,mesh.getBulkData()->buckets(mesh.getNodeRank()),nodes_b,sort_by_gid);
      // convert the entities to global ids
      gids_b.resize(nodes_b.size());
      size_t i = 0;
      for (auto&& node : nodes_b) {
        gids_b[i] = mesh.getBulkData()->identifier(node);
        ++i;
      }
    }

    // Sort the element gids so we can use binary search
    std::sort(gids_b.begin(),gids_b.end());

    if (print_debug) {
      Teuchos::FancyOStream os(Teuchos::rcpFromRef(std::cout));
      os.setShowProcRank(true);
      os << std::endl;
      os << "gids_a.size()=" << gids_a.size() << std::endl;
      for (auto&& gid : gids_a)
        os << "gid_a=" << gid << std::endl;
      os << "gids_b.size()=" << gids_b.size() << std::endl;
      for (auto&& gid : gids_b)
        os << "gid_b=" << gid << std::endl;
    }

    // Search for each node in a in b
    // 0 = no overlap, 1 = overlap
    // We use int for MPI communication
    int has_local_overlap = 0;
    for (auto&& a : gids_a) {
      if (std::binary_search(gids_b.begin(),gids_b.end(),a)) {
        has_local_overlap = 1;
        break;
      }
    }
    int has_overlap = 0;
    Teuchos::reduceAll(*mesh.getComm(),Teuchos::REDUCE_SUM,1,&has_local_overlap,&has_overlap);
    if (has_overlap == 0)
      return false;

    return true;
  }
}
