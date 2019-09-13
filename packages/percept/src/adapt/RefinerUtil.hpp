// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_RefinerUtil_hpp
#define adapt_RefinerUtil_hpp

#include <iostream>
#include <stdexcept>
#include <string>
#include <sstream>
#include <cmath>
#include <utility>
#include <math.h>
#include <map>
#include <set>
#include <vector>

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyData.h>

#include <percept/stk_mesh.hpp>
#include <percept/PerceptMesh.hpp>
#include <percept/ProgressMeter.hpp>
#include <adapt/UniformRefinerPattern.hpp>
#include <adapt/Refiner.hpp>
#include <adapt/NodeRegistry.hpp>

#include <adapt/SubDimCell.hpp>

#include <adapt/RefinementInfoByType.hpp>


namespace percept {


class RefinerUtil
{
public:

  static BlockNamesType
  getBlockNames(const std::string& block_name, unsigned proc_rank, percept::PerceptMesh& eMesh, const std::string& geomFile="");

  static BlockNamesType
  correctBlockNamesForPartPartConsistency(percept::PerceptMesh& eMesh, BlockNamesType& blocks, const std::string& geomFile);

  static void 
  remove_existing_nodes(PerceptMesh& eMesh);
  static void 
  add_new_nodes(PerceptMesh& eMesh);
  static void 
  collect_locally_owned_entities(PerceptMesh& eMesh,  stk::mesh::EntityRank rank,  std::vector<stk::mesh::Entity>& elements);
  static void 
  get_parent_entity_and_id(PerceptMesh& eMesh, stk::mesh::EntityRank rank, stk::mesh::Entity& element, stk::mesh::Entity& parent_elem, bool debug);
  static void 
  find_or_set_parent_child_relation(PerceptMesh& eMesh, stk::mesh::EntityRank rank, stk::mesh::Entity& element, stk::mesh::Entity& parent_elem, size_t& i_ft, std::vector<stk::mesh::Entity>& ft_new_elements);

  static void
  rebuild_family_tree(PerceptMesh& eMesh, bool debug=false);

  /// NodeRegistry
  static void rebuild_node_registry(PerceptMesh& eMesh, NodeRegistry& nodeRegistry, bool initNR = true, PerceptMesh *eMeshNR = 0, NodeRegistry *compareNR=0, bool skipEmpty = true);
  static void save_node_registry(PerceptMesh& eMesh, NodeRegistry& nodeRegistry, const std::string& msg, bool doComm=true);
};

}

#endif
