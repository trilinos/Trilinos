// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include "stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp"
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_io/WriteMesh.hpp>
#include <stk_tools/mesh_tools/DisconnectBlocks.hpp>
#include <stk_tools/mesh_tools/DisconnectBlocksImpl.hpp>
#include <stk_tools/mesh_tools/DetectHingesImpl.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_util/environment/WallTime.hpp>
#include "stk_unit_test_utils/getOption.h"
#include "stk_util/parallel/ParallelReduce.hpp"
#include "stk_util/util/SortAndUnique.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/base/FEMHelpers.hpp"
#include "stk_mesh/base/EntityLess.hpp"
#include "stk_mesh/base/Comm.hpp"

#include <string>
#include <algorithm>

namespace {

#define EPS 0.15

stk::mesh::Part & create_part(stk::mesh::MetaData& meta, const stk::topology topology, const std::string & blockName, int64_t blockId)
{
  stk::mesh::Part& part = meta.declare_part_with_topology(blockName, topology);
  stk::io::put_io_part_attribute(part);
  meta.set_part_id(part, blockId);
  return part;
}

//void print_node_count(stk::mesh::BulkData& bulk, const std::string str)
//{
//  stk::mesh::EntityVector nodes;
//  bulk.get_entities(stk::topology::NODE_RANK, bulk.mesh_meta_data().universal_part(), nodes);
//
//  std::cout << str << std::endl;
//  std::cout << "p:" << bulk.parallel_rank() << " node vec size: " << nodes.size() << std::endl;
//}

bool is_new_owner(const stk::mesh::BulkData& bulk, const stk::mesh::Entity& elem, const stk::mesh::Entity& node)
{
  stk::mesh::EntityVector entityVec(bulk.begin_elements(node), bulk.begin_elements(node)+bulk.num_elements(node));
  std::sort(entityVec.begin(), entityVec.end(), stk::mesh::EntityLess(bulk));

  return (entityVec[0] == elem) ? true : false;
}

void distribute_mesh(stk::mesh::BulkData& bulk, const stk::mesh::EntityIdProcVec& idProcVec)
{
  stk::mesh::EntityProcVec procVec;
  if(bulk.parallel_rank() == 0) {
    for(const stk::mesh::EntityIdProc& idProc : idProcVec) {

      stk::mesh::Entity elem = bulk.get_entity(stk::topology::ELEM_RANK, idProc.first);
      ThrowRequire(bulk.is_valid(elem));
      if(idProc.second != bulk.parallel_rank()) {
        procVec.push_back(std::make_pair(elem, idProc.second));
//        std::cout << "p:" << bulk.parallel_rank() << " element: " << bulk.entity_key(elem) << std::endl;
      }

      stk::mesh::Entity const * elemNodes = bulk.begin_nodes(elem);
      for(unsigned i = 0; i < bulk.num_nodes(elem); i++) {
        stk::mesh::Entity node = elemNodes[i];
        ThrowRequire(bulk.is_valid(node));
        if(idProc.second != bulk.parallel_rank() && is_new_owner(bulk, elem, node)) {
          procVec.push_back(std::make_pair(node, idProc.second));
//          std::cout << "p:" << bulk.parallel_rank() << " node: " << bulk.entity_key(node) << std::endl;
        }
      }
    }
  }
  bulk.change_entity_owner(procVec);

//  print_node_count(bulk, "distribute_mesh");
}

stk::mesh::PartVector setup_mesh_1block_1quad(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  std::string meshDesc = "0,1,QUAD_4_2D,1,2,4,3,block_1";
  std::vector<double> coordinates = { 0,0, 1,0, 0,1, 1,1 };
  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1};
}

stk::mesh::PartVector setup_mesh_1block_2quad(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  std::string meshDesc = "0,1,QUAD_4_2D,1,2,4,3,block_1\n"
                         "0,2,QUAD_4_2D,2,5,6,4,block_1";
  std::vector<double> coordinates = { 0,0, 1,0, 0,1, 1,1, 2,0, 2,1 };
  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1};
}

stk::mesh::PartVector setup_mesh_2block_2quad(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);
  std::string meshDesc = "0,1,QUAD_4_2D,1,2,4,3,block_1\n"
                         "0,2,QUAD_4_2D,2,5,6,4,block_2";
  std::vector<double> coordinates = { 0,0, 1,0, 0,1, 1,1, 2,0, 2,1 };
  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1,&block2};
}

stk::mesh::PartVector setup_mesh_1block_2quad_1node_hinge(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  std::string meshDesc = "0,1,QUAD_4_2D,1,2,4,3,block_1\n"
                         "0,2,QUAD_4_2D,7,5,6,4,block_1";
  std::vector<double> coordinates = { 0,0, (1-EPS),0, 0,1, 1,1, 2,0, 2,1, (1+EPS),0 };
  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1};
}

stk::mesh::PartVector setup_mesh_2block_2quad_1node_hinge(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);
  std::string meshDesc = "0,1,QUAD_4_2D,1,2,4,3,block_1\n"
                         "0,2,QUAD_4_2D,7,5,6,4,block_2";
  std::vector<double> coordinates = { 0,0, (1-EPS),0, 0,1, 1,1, 2,0, 2,1, (1+EPS),0 };
  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1,&block2};
}

stk::mesh::PartVector setup_mesh_1block_2quad_2hinge(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  std::string meshDesc = "0,1,QUAD_4_2D,3,4,1,2,block_1\n"
                         "0,2,QUAD_4_2D,2,6,4,5,block_1";
  std::vector<double> coordinates = { 0,2, 2,1, 1,2, 2,3, 3,2, 4,2 };
  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1};
}

stk::mesh::PartVector setup_mesh_2block_2quad_2hinge(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);
  std::string meshDesc = "0,1,QUAD_4_2D,3,4,1,2,block_1\n"
                         "0,2,QUAD_4_2D,2,6,4,5,block_2";
  std::vector<double> coordinates = { 0,2, 2,1, 1,2, 2,3, 3,2, 4,2 };
  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1,&block2};
}

stk::mesh::PartVector setup_mesh_1block_3quad_1hinge(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  std::string meshDesc = "0,1,QUAD_4_2D,1,2,4,3,block_1\n"
                         "0,2,QUAD_4_2D,2,5,6,4,block_1\n"
                         "0,3,QUAD_4_2D,4,7,8,9,block_1";
  std::vector<double> coordinates = { 0,0, 1,0, 0,1, 1,1, 2,0, 2,(1-EPS), 2,(1+EPS), 2,2, 1,2 };
  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1};
}

stk::mesh::PartVector setup_mesh_1block_3quad_1hinge_linear_stack(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  std::string meshDesc = "0,1,QUAD_4_2D,1,2,4,3,block_1\n"
                         "0,2,QUAD_4_2D,2,5,6,4,block_1\n"
                         "0,3,QUAD_4_2D,7,8,9,6,block_1";
  std::vector<double> coordinates = { 0,0, 1,0, 0,1, 1,1, (2-EPS),0, 2,1, (2+EPS),0, 3,0, 3,1 };
  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1};
}

stk::mesh::PartVector setup_mesh_3block_3quad_1hinge_linear_stack(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);
  stk::mesh::Part & block3 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_3", 3);
  std::string meshDesc = "0,1,QUAD_4_2D,1,2,4,3,block_1\n"
                         "0,2,QUAD_4_2D,2,5,6,4,block_2\n"
                         "0,3,QUAD_4_2D,7,8,9,6,block_3";
  std::vector<double> coordinates = { 0,0, 1,0, 0,1, 1,1, (2-EPS),0, 2,1, (2+EPS),0, 3,0, 3,1 };
  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1,&block2,&block3};
}

stk::mesh::PartVector setup_mesh_1block_4quad_bowtie_1hinge(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  std::string meshDesc = "0,1,QUAD_4_2D,1,2,4,3,block_1\n"
                         "0,2,QUAD_4_2D,5,6,7,4,block_1\n"
                         "0,3,QUAD_4_2D,4,8,9,10,block_1\n"
                         "0,4,QUAD_4_2D,4,11,12,13,block_1";
  std::vector<double> coordinates = { 0,0, (1-EPS),0, 0,(1-EPS), 1,1, (1+EPS),0,
                                      2,0, 2,(1-EPS), 2,(1+EPS), 2,2, (1+EPS),2, (1-EPS),2, 0,2, 0,(1+EPS) };
  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1};
}

stk::mesh::PartVector setup_mesh_1block_4quad_2hinge(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  std::string meshDesc = "0,1,QUAD_4_2D,1,2,4,3,block_1\n"
                         "0,2,QUAD_4_2D,2,5,6,4,block_1\n"
                         "0,3,QUAD_4_2D,6,7,8,9,block_1\n"
                         "0,4,QUAD_4_2D,3,9,8,10,block_1";
  std::vector<double> coordinates = { 0,0, 1,0, 0,1, 1,(1-EPS), 2,0, 2,1, 2,2, 1,2, 1,(1+EPS), 0,2 };
  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1};
}

stk::mesh::PartVector setup_mesh_4block_4quad_2hinge(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);
  stk::mesh::Part & block3 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_3", 3);
  stk::mesh::Part & block4 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_4", 4);
  std::string meshDesc = "0,1,QUAD_4_2D,1,2,4,3,block_1\n"
                         "0,2,QUAD_4_2D,2,5,6,4,block_2\n"
                         "0,3,QUAD_4_2D,6,7,8,9,block_3\n"
                         "0,4,QUAD_4_2D,3,9,8,10,block_4";
  std::vector<double> coordinates = { 0,0, 1,0, 0,1, 1,(1-EPS), 2,0, 2,1, 2,2, 1,2, 1,(1+EPS), 0,2 };
  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1,&block2,&block3,&block4};
}

stk::mesh::PartVector setup_mesh_1block_4quad_4hinge(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  std::string meshDesc = "0,1,QUAD_4_2D,1,2,4,3,block_1\n"
                         "0,2,QUAD_4_2D,2,5,6,7,block_1\n"
                         "0,3,QUAD_4_2D,6,9,10,8,block_1\n"
                         "0,4,QUAD_4_2D,3,11,10,12,block_1";
  std::vector<double> coordinates = { 0,0, 1,0, 0,1, 1,1, 2,0, 2,1, 1,1, 1,1, 2,2, 1,2, 1,1, 0,2 };
  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1};
}

stk::mesh::PartVector setup_mesh_4block_4quad_4hinge(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);
  stk::mesh::Part & block3 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_3", 3);
  stk::mesh::Part & block4 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_4", 4);
  std::string meshDesc = "0,1,QUAD_4_2D,1,2,4,3,block_1\n"
                         "0,2,QUAD_4_2D,2,5,6,7,block_2\n"
                         "0,3,QUAD_4_2D,6,9,10,8,block_3\n"
                         "0,4,QUAD_4_2D,3,11,10,12,block_4";
  std::vector<double> coordinates = { 0,0, 1,0, 0,1, 1,1, 2,0, 2,1, 1,1, 1,1, 2,2, 1,2, 1,1, 0,2 };
  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1,&block2,&block3,&block4};
}

stk::mesh::PartVector setup_mesh_1block_4quad_pacman(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  std::string meshDesc = "0,1,QUAD_4_2D,1,2,4,3,block_1\n"
                         "0,2,QUAD_4_2D,2,5,6,4,block_1\n"
                         "0,3,QUAD_4_2D,4,7,8,9,block_1\n"
                         "0,4,QUAD_4_2D,3,4,9,10,block_1";
  std::vector<double> coordinates = { 0,0, 1,0, 0,1, 1,1, 2,0, 2,(1-EPS), 2,(1+EPS), 2,2, 1,2, 0,2 };
  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1};
}

stk::mesh::PartVector setup_mesh_4block_4quad_pacman(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 1);
  stk::mesh::Part & block3 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_3", 1);
  stk::mesh::Part & block4 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_4", 1);
  std::string meshDesc = "0,1,QUAD_4_2D,1,2,4,3,block_1\n"
                         "0,2,QUAD_4_2D,2,5,6,4,block_2\n"
                         "0,3,QUAD_4_2D,4,7,8,9,block_3\n"
                         "0,4,QUAD_4_2D,3,4,9,10,block_4";
  std::vector<double> coordinates = { 0,0, 1,0, 0,1, 1,1, 2,0, 2,(1-EPS), 2,(1+EPS), 2,2, 1,2, 0,2 };
  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1,&block2,&block3,&block4};
}

stk::mesh::PartVector setup_mesh_1block_4quad_1hinge(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  std::string meshDesc = "0,1,QUAD_4_2D,1,2,4,3,block_1\n"
                         "0,2,QUAD_4_2D,5,6,7,4,block_1\n"
                         "0,3,QUAD_4_2D,7,8,9,4,block_1\n"
                         "0,4,QUAD_4_2D,4,9,10,11,block_1";
  std::vector<double> coordinates = { 0,0, (1-EPS),0, 0,(1-EPS), 1,1, (1+EPS),0, 2,0, 2,1, 2,2, 1,2, 0,2, 0,(1+EPS) };
  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1};
}

stk::mesh::PartVector setup_mesh_4block_4quad_1hinge(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);
  stk::mesh::Part & block3 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_3", 3);
  stk::mesh::Part & block4 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_4", 4);
  std::string meshDesc = "0,1,QUAD_4_2D,1,2,4,3,block_1\n"
                         "0,2,QUAD_4_2D,5,6,7,4,block_2\n"
                         "0,3,QUAD_4_2D,7,8,9,4,block_3\n"
                         "0,4,QUAD_4_2D,4,9,10,11,block_4";
  std::vector<double> coordinates = { 0,0, (1-EPS),0, 0,(1-EPS), 1,1, (1+EPS),0, 2,0, 2,1, 2,2, 1,2, 0,2, 0,(1+EPS) };
  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1,&block2,&block3,&block4};
}

stk::mesh::PartVector setup_mesh_2block_3quad_2tri_1hinge(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::TRI_3_2D, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);
  std::string meshDesc = "0,1,TRI_3_2D,1,2,4,block_1\n"
                         "0,2,TRI_3_2D,1,4,3,block_1\n"
                         "0,3,QUAD_4_2D,5,6,7,4,block_2\n"
                         "0,4,QUAD_4_2D,7,8,9,4,block_2\n"
                         "0,5,QUAD_4_2D,4,9,10,11,block_2";
  std::vector<double> coordinates = { 0,0, (1-EPS),0, 0,(1-EPS), 1,1, (1+EPS),0, 2,0, 2,1, 2,2, 1,2, 0,2, 0,(1+EPS) };
  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1,&block2};
}

stk::mesh::PartVector setup_mesh_5block_3quad_2tri_1hinge(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::TRI_3_2D, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::TRI_3_2D, "block_2", 2);
  stk::mesh::Part & block3 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_3", 3);
  stk::mesh::Part & block4 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_4", 4);
  stk::mesh::Part & block5 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_5", 5);
  std::string meshDesc = "0,1,TRI_3_2D,1,2,4,block_1\n"
                         "0,2,TRI_3_2D,1,4,3,block_2\n"
                         "0,3,QUAD_4_2D,5,6,7,4,block_3\n"
                         "0,4,QUAD_4_2D,7,8,9,4,block_4\n"
                         "0,5,QUAD_4_2D,4,9,10,11,block_5";
  std::vector<double> coordinates = { 0,0, (1-EPS),0, 0,(1-EPS), 1,1, (1+EPS),0, 2,0, 2,1, 2,2, 1,2, 0,2, 0,(1+EPS) };
  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1,&block2,&block3,&block4,&block5};
}

stk::mesh::PartVector setup_mesh_1block_1hex(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_1", 1);
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1";
  std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1};

  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1};
}

stk::mesh::PartVector setup_mesh_1block_2hex(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_1", 1);
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                         "0,2,HEX_8,5,6,7,8,9,10,11,12,block_1";
  std::vector<double> coordinates = {
    0,0,0, 1,0,0, 1,1,0, 0,1,0,
    0,0,1, 1,0,1, 1,1,1, 0,1,1,
    0,0,2, 1,0,2, 1,1,2, 0,1,2
  };

  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1};
}

stk::mesh::PartVector setup_mesh_1block_2hex_1node_hinge(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_1", 1);
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                         "0,2,HEX_8,5,9,10,11,12,13,14,15,block_1";
  std::vector<double> coordinates = {
    0,0,0, 1,0,0, 1,1,0, 0,1,0,
    0,0,1, 1,0,1, 1,1,1, 0,1,1,
    1,0,(1+EPS), 1,1,(1+EPS), 0,1,(1+EPS), 0,0,2,
    1,0,2, 1,1,2, 0,1,2
  };

  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1};
}

stk::mesh::PartVector setup_mesh_1block_2hex_2node_hinge(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_1", 1);
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                         "0,2,HEX_8,2,9,4,10,11,12,13,14,block_1";
  std::vector<double> coordinates = {
    0,0,0, 1,0,0, 0.5,0.5,0, 0,1,0,
    0,0,1, 1,0,1, 0.5,0.5,1, 0,1,1,
    1,1,0, (0.5+EPS),(0.5+EPS),0, 1,EPS,1, 1,1,1, 0,(1+EPS),1, (0.5+EPS),(0.5+EPS),1
  };

  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1};
}

stk::mesh::PartVector setup_mesh_1block_3hex_1node_hinge(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_1", 1);
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                         "0,2,HEX_8,5,6,7,8,9,10,11,12,block_1\n"
                         "0,3,HEX_8,13,14,15,16,6,17,18,19,block_1";
  std::vector<double> coordinates = {
    0,0,0, 1,0,0, 1,1,0, 0,1,0,
    0,0,1, 1,0,1, 1,1,1, 0,1,1,
    0,0,2, 1,0,2, 1,1,2, 0,1,2,
    (1+EPS),0,0, 2,0,0, 2,1,0, (1+EPS),1,0, 2,0,1, 2,1,1, (1+EPS),1,1
  };

  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1};
}

stk::mesh::PartVector setup_mesh_1block_2hex_face_test(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_1", 1);
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                         "0,2,HEX_8,2,9,10,3,6,11,12,7,block_1";
  std::vector<double> coordinates = {
    0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1,
    2,0,0, 2,1,0, 2,0,1, 2,1,1
  };

  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1};
}


stk::mesh::PartVector setup_mesh_1block_8hex_flower_1node_hinge(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_1", 1);
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,7,block_1\n"
                         "0,3,HEX_8,16,17,18,19,20,7,21,22,block_1\n"
                         "0,4,HEX_8,23,24,25,26,7,27,28,29,block_1\n"
                         "0,5,HEX_8,30,31,7,32,33,34,35,36,block_1\n"
                         "0,6,HEX_8,37,38,39,7,40,41,42,43,block_1\n"
                         "0,7,HEX_8,7,44,45,46,47,48,49,50,block_1\n"
                         "0,8,HEX_8,51,7,52,53,54,55,56,57,block_1";
  std::vector<double> coordinates = {
    0,0,0, (1-EPS),0,0, (1-EPS),(1-EPS),0, 0,(1-EPS),0, 0,0,(1-EPS), (1-EPS),0,(1-EPS), 1,1,1, 0,(1-EPS),(1-EPS),
    (1+EPS),0,0, 2,0,0, 2,(1-EPS),0, (1+EPS),(1-EPS),0, (1+EPS),0,(1-EPS), 2,0,(1-EPS), 2,(1-EPS),(1-EPS),
    0,(1+EPS),0, (1-EPS),(1+EPS),0, (1-EPS),2,0, 0,2,0, 0,(1+EPS),(1-EPS), (1-EPS),2,(1-EPS), 0,2,(1-EPS),
    (1+EPS),(1+EPS),0, 2,(1+EPS),0, 2,2,0, (1+EPS),2,0, 2,(1+EPS),(1-EPS), 2,2,(1-EPS), (1+EPS),2,(1-EPS),
    0,0,(1+EPS), (1-EPS),0,(1+EPS), 0,(1-EPS),(1+EPS), 0,0,2, (1-EPS),0,2, (1-EPS),(1-EPS),2, 0,(1-EPS),2,
    (1+EPS),0,(1+EPS), 2,0,(1+EPS), 2,(1-EPS),(1+EPS), (1+EPS),0,2, 2,0,2, 2,(1-EPS),2, (1+EPS),(1-EPS),2,
    2,(1+EPS),(1+EPS), 2,2,(1+EPS), (1+EPS),2,(1+EPS), (1+EPS),(1+EPS),2, 2,(1+EPS),2, 2,2,2, (1+EPS),2,2,
    0,(1+EPS),(1+EPS), (1-EPS),2,(1+EPS), 0,2,(1+EPS), 0,(1+EPS),2, (1-EPS),(1+EPS),2, (1-EPS),2,2, 0,2,2
  };

  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1};
}


stk::mesh::PartVector setup_mesh_1block_2tet_1node_hinge(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::TET_4, "block_1", 1);
  std::string meshDesc = "0,1,TET_4,1,2,3,4,block_1\n"
                         "0,2,TET_4,3,5,6,7,block_1";
  std::vector<double> coordinates = {
    0,0,0, 1,0,0, 0,1,0, 0,0,1,
    1,1,0, 0,2,0, 0,1,1
  };

  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1};
}

stk::mesh::PartVector setup_mesh_1block_2hex_1edge_hinge(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_1", 1);
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                         "0,2,HEX_8,4,9,10,11,8,12,13,14,block_1";
  std::vector<double> coordinates = {
    0,0,0, 0,0,1, 1,0,1, 1,0,0, 0,1,0, 0,1,1, 1,1,1, 1,1,0,
    2,0,0, 2,0,-1, 1,0,-1, 2,1,0, 2,1,-1, 1,1,-1
  };

  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1};
}

stk::mesh::PartVector setup_mesh_2block_2hex_1edge_hinge(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_2", 2);
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                         "0,2,HEX_8,4,9,10,11,8,12,13,14,block_2";
  std::vector<double> coordinates = {
    0,0,0, 0,0,1, 1,0,1, 1,0,0, 0,1,0, 0,1,1, 1,1,1, 1,1,0,
    2,0,0, 2,0,-1, 1,0,-1, 2,1,0, 2,1,-1, 1,1,-1
  };

  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1,&block2};
}

stk::mesh::PartVector setup_mesh_1block_3hex_1edge_hinge(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_1", 1);
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                         "0,2,HEX_8,2,9,10,3,6,11,12,7,block_1\n"
                         "0,3,HEX_8,13,14,15,16,8,7,17,18,block_1";
  std::vector<double> coordinates = {
    0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1,
    2,0,0, 2,1,0, 2,0,1, 2,1,1,
    0,(1+EPS),0, 1,(1+EPS),0, 1,2,0, 0,2,0, 1,2,1, 0,2,1,
  };

  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1};
}

stk::mesh::PartVector setup_mesh_3block_3hex_1edge_hinge(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_2", 2);
  stk::mesh::Part & block3 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_3", 3);
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                         "0,2,HEX_8,2,9,10,3,6,11,12,7,block_2\n"
                         "0,3,HEX_8,13,14,15,16,8,7,17,18,block_3";
  std::vector<double> coordinates = {
    0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1,
    2,0,0, 2,1,0, 2,0,1, 2,1,1,
    0,(1+EPS),0, 1,(1+EPS),0, 1,2,0, 0,2,0, 1,2,1, 0,2,1,
  };

  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1,&block2,&block3};
}

stk::mesh::PartVector setup_mesh_1block_3hex_1node_hinge_1edge_hinge(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_1", 1);
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                         "0,2,HEX_8,4,9,10,11,8,12,13,14,block_1\n"
                         "0,3,HEX_8,8,15,16,17,18,19,20,21,block_1";
  std::vector<double> coordinates = {
    0,0,0, 0,0,1, 1,0,1, 1,0,0, 0,1,0, 0,1,1, 1,1,1, 1,1,0,
    2,0,0, 2,0,-1, 1,0,-1, 2,1,0, 2,1,-1, 1,1,-1,
    2,(1+EPS),0, 2,(1+EPS),-1, 1,(1+EPS),-1, 1,2,0, 2,2,0, 2,2,-1, 1,2,-1
  };

  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1};
}

stk::mesh::PartVector setup_mesh_3block_3hex_1node_hinge_1edge_hinge(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_2", 21);
  stk::mesh::Part & block3 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_3", 3);
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                         "0,2,HEX_8,4,9,10,11,8,12,13,14,block_2\n"
                         "0,3,HEX_8,8,15,16,17,18,19,20,21,block_3";
  std::vector<double> coordinates = {
    0,0,0, 0,0,1, 1,0,1, 1,0,0, 0,1,0, 0,1,1, 1,1,1, 1,1,0,
    2,0,0, 2,0,-1, 1,0,-1, 2,1,0, 2,1,-1, 1,1,-1,
    2,(1+EPS),0, 2,(1+EPS),-1, 1,(1+EPS),-1, 1,2,0, 2,2,0, 2,2,-1, 1,2,-1
  };

  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1,&block2,&block3};
}

stk::mesh::PartVector setup_mesh_1block_3hex_1node_hinge_1edge_hinge2(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_1", 1);
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                         "0,2,HEX_8,4,9,10,11,8,12,13,14,block_1\n"
                         "0,3,HEX_8,15,12,16,17,18,19,20,21,block_1";
  std::vector<double> coordinates = {
    0,0,0, 0,0,1, 1,0,1, 1,0,0, 0,1,0, 0,1,1, 1,1,1, 1,1,0,
    2,0,0, 2,0,-1, 1,0,-1, 2,1,0, 2,1,-1, 1,1,-1,
    1,(1+EPS),0, 2,(1+EPS),-1, 1,(1+EPS),-1, 1,2,0, 2,2,0, 2,2,-1, 1,2,-1
  };

  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1};
}

stk::mesh::PartVector setup_mesh_3block_3hex_1node_hinge_1edge_hinge2(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_2", 2);
  stk::mesh::Part & block3 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_3", 3);
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                         "0,2,HEX_8,4,9,10,11,8,12,13,14,block_2\n"
                         "0,3,HEX_8,15,12,16,17,18,19,20,21,block_3";
  std::vector<double> coordinates = {
    0,0,0, 0,0,1, 1,0,1, 1,0,0, 0,1,0, 0,1,1, 1,1,1, 1,1,0,
    2,0,0, 2,0,-1, 1,0,-1, 2,1,0, 2,1,-1, 1,1,-1,
    1,(1+EPS),0, 2,(1+EPS),-1, 1,(1+EPS),-1, 1,2,0, 2,2,0, 2,2,-1, 1,2,-1
  };

  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1,&block2,&block3};
}

stk::mesh::PartVector setup_mesh_1block_3hex_1node_hinge_1edge_hinge3(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_1", 1);
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                         "0,2,HEX_8,4,9,10,11,8,12,13,14,block_1\n"
                         "0,3,HEX_8,15,16,13,17,18,19,20,21,block_1";
  std::vector<double> coordinates = {
    0,0,0, 0,0,1, 1,0,1, 1,0,0, 0,1,0, 0,1,1, 1,1,1, 1,1,0,
    2,0,0, 2,0,-1, 1,0,-1, 2,1,0, 2,1,-1, 1,1,-1,
    1,(1+EPS),0, 2,(1+EPS),0, 1,(1+EPS),-1, 1,2,0, 2,2,0, 2,2,-1, 1,2,-1
  };

  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1};
}

stk::mesh::PartVector setup_mesh_3block_3hex_1node_hinge_1edge_hinge3(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_2", 2);
  stk::mesh::Part & block3 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_3", 3);
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                         "0,2,HEX_8,4,9,10,11,8,12,13,14,block_2\n"
                         "0,3,HEX_8,15,16,13,17,18,19,20,21,block_3";
  std::vector<double> coordinates = {
    0,0,0, 0,0,1, 1,0,1, 1,0,0, 0,1,0, 0,1,1, 1,1,1, 1,1,0,
    2,0,0, 2,0,-1, 1,0,-1, 2,1,0, 2,1,-1, 1,1,-1,
    1,(1+EPS),0, 2,(1+EPS),0, 1,(1+EPS),-1, 1,2,0, 2,2,0, 2,2,-1, 1,2,-1
  };

  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1,&block2,&block3};
}

stk::mesh::PartVector setup_mesh_1block_4hex_bowtie_1edge_hinge(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_1", 1);
  std::string meshDesc = "0,1,HEX_8,1,2,4,3,14,15,17,16,block_1\n"
                         "0,2,HEX_8,5,6,7,4,18,19,20,17,block_1\n"
                         "0,3,HEX_8,4,8,9,10,17,21,22,23,block_1\n"
                         "0,4,HEX_8,4,11,12,13,17,24,25,26,block_1";
  std::vector<double> coordinates = {
    0,0,0, (1-EPS),0,0, 0,(1-EPS),0, 1,1,0, (1+EPS),0,0, 2,0,0, 2,(1-EPS),0, 2,(1+EPS),0,
    2,2,0, (1+EPS),2,0, (1-EPS),2,0, 0,2,0, 0,(1+EPS),0,
    0,0,1, (1-EPS),0,1, 0,(1-EPS),1, 1,1,1, (1+EPS),0,1, 2,0,1, 2,(1-EPS),1, 2,(1+EPS),1,
    2,2,1, (1+EPS),2,1, (1-EPS),2,1, 0,2,1, 0,(1+EPS),1
  };

  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1};
}

stk::mesh::PartVector setup_mesh_4block_4hex_bowtie_1edge_hinge(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_2", 2);
  stk::mesh::Part & block3 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_3", 3);
  stk::mesh::Part & block4 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_4", 4);
  std::string meshDesc = "0,1,HEX_8,1,2,4,3,14,15,17,16,block_1\n"
                         "0,2,HEX_8,5,6,7,4,18,19,20,17,block_2\n"
                         "0,3,HEX_8,4,8,9,10,17,21,22,23,block_3\n"
                         "0,4,HEX_8,4,11,12,13,17,24,25,26,block_4";
  std::vector<double> coordinates = {
    0,0,0, (1-EPS),0,0, 0,(1-EPS),0, 1,1,0, (1+EPS),0,0, 2,0,0, 2,(1-EPS),0, 2,(1+EPS),0,
    2,2,0, (1+EPS),2,0, (1-EPS),2,0, 0,2,0, 0,(1+EPS),0,
    0,0,1, (1-EPS),0,1, 0,(1-EPS),1, 1,1,1, (1+EPS),0,1, 2,0,1, 2,(1-EPS),1, 2,(1+EPS),1,
    2,2,1, (1+EPS),2,1, (1-EPS),2,1, 0,2,1, 0,(1+EPS),1
  };

  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1,&block2,&block3,&block4};
}

stk::mesh::PartVector setup_mesh_1block_two_by_two_hex_2edge_hinge(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_1", 1);
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                         "0,2,HEX_8,4,9,10,11,8,12,13,14,block_1\n"
                         "0,3,HEX_8,5,6,7,8,15,16,17,18,block_1\n"
                         "0,4,HEX_8,8,12,13,14,18,19,20,21,block_1";
  std::vector<double> coordinates = {
    0,0,0, 0,0,1, 1,0,1, 1,0,0, 0,1,0, 0,1,1, 1,1,1, 1,1,0,
    2,0,0, 2,0,-1, 1,0,-1, 2,1,0, 2,1,-1, 1,1,-1,
    0,2,0, 0,2,1, 1,2,1, 1,2,0,
    2,2,0, 2,2,-1, 1,2,-1,
  };

  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1};
}

stk::mesh::PartVector setup_mesh_4block_two_by_two_hex_2edge_hinge(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_2", 2);
  stk::mesh::Part & block3 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_3", 3);
  stk::mesh::Part & block4 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_4", 4);
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                         "0,2,HEX_8,4,9,10,11,8,12,13,14,block_2\n"
                         "0,3,HEX_8,5,6,7,8,15,16,17,18,block_3\n"
                         "0,4,HEX_8,8,12,13,14,18,19,20,21,block_4";
  std::vector<double> coordinates = {
    0,0,0, 0,0,1, 1,0,1, 1,0,0, 0,1,0, 0,1,1, 1,1,1, 1,1,0,
    2,0,0, 2,0,-1, 1,0,-1, 2,1,0, 2,1,-1, 1,1,-1,
    0,2,0, 0,2,1, 1,2,1, 1,2,0,
    2,2,0, 2,2,-1, 1,2,-1,
  };

  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1,&block2,&block3,&block4};
}

stk::mesh::PartVector setup_mesh_1block_four_hex_one_edge_one_node_hinge(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_1", 1);
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,10,11,12,13,block_1\n"
                         "0,2,HEX_8,1,4,5,6,10,13,14,15,block_1\n"
                         "0,3,HEX_8,1,7,8,9,10,16,17,18,block_1\n"
                         "0,4,HEX_8,13,19,20,21,22,23,24,25,block_1";
  std::vector<double> coordinates = {
    0,0,0, 1,EPS,0, 1,1,0, 0,1,0, -1,1,0, -1,0,0, 0,-1,0, 1,-1,0, 1,-EPS,0,
    0,0,1, 1,EPS,1, 1,1,1, 0,1,1, -1,1,1, -1,0,1, 0,-1,1, 1,-1,1, 1,-EPS,1,
    0,2,1, -1,2,1, -1,(1+EPS),1, 0,1,2, 0,2,2, -1,2,2, -1,1,2
  };

  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1};
}

stk::mesh::PartVector setup_mesh_4block_four_hex_one_edge_one_node_hinge(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_2", 2);
  stk::mesh::Part & block3 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_3", 3);
  stk::mesh::Part & block4 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_4", 4);
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,10,11,12,13,block_1\n"
                         "0,2,HEX_8,1,4,5,6,10,13,14,15,block_2\n"
                         "0,3,HEX_8,1,7,8,9,10,16,17,18,block_3\n"
                         "0,4,HEX_8,13,19,20,21,22,23,24,25,block_4";
  std::vector<double> coordinates = {
    0,0,0, 1,EPS,0, 1,1,0, 0,1,0, -1,1,0, -1,0,0, 0,-1,0, 1,-1,0, 1,-EPS,0,
    0,0,1, 1,EPS,1, 1,1,1, 0,1,1, -1,1,1, -1,0,1, 0,-1,1, 1,-1,1, 1,-EPS,1,
    0,2,1, -1,2,1, -1,(1+EPS),1, 0,1,2, 0,2,2, -1,2,2, -1,1,2
  };

  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1,&block2,&block3,&block4};
}

stk::mesh::PartVector setup_mesh_1block_four_hex_2node_hinge(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_1", 1);
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,7,8,9,10,block_1\n"
                         "0,2,HEX_8,1,4,5,6,7,10,11,12,block_1\n"
                         "0,3,HEX_8,7,13,14,15,16,17,18,19,block_1\n"
                         "0,4,HEX_8,10,20,21,22,23,24,25,26,block_1";
  std::vector<double> coordinates = {
    1,1,0, 2,1,0, 2,2,0, 1,2,0, 0,2,0, 0,1,0,
    1,1,1, 2,1,1, 2,2,1, 1,2,1, 0,2,1, 0,1,1,
    0,(1-EPS),1, 0,0,1, 1,0,1, 1,1,2, 0,1,2, 0,0,2, 1,0,2,
    1,3,1, 0,3,1, 0,(2+EPS),1, 1,2,2, 1,3,2, 0,3,2, 0,2,2
  };

  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1};
}

stk::mesh::PartVector setup_mesh_4block_four_hex_2node_hinge(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_2", 2);
  stk::mesh::Part & block3 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_3", 3);
  stk::mesh::Part & block4 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_4", 4);
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,7,8,9,10,block_1\n"
                         "0,2,HEX_8,1,4,5,6,7,10,11,12,block_2\n"
                         "0,3,HEX_8,7,13,14,15,16,17,18,19,block_3\n"
                         "0,4,HEX_8,10,20,21,22,23,24,25,26,block_4";
  std::vector<double> coordinates = {
    1,1,0, 2,1,0, 2,2,0, 1,2,0, 0,2,0, 0,1,0,
    1,1,1, 2,1,1, 2,2,1, 1,2,1, 0,2,1, 0,1,1,
    0,(1-EPS),1, 0,0,1, 1,0,1, 1,1,2, 0,1,2, 0,0,2, 1,0,2,
    1,3,1, 0,3,1, 0,(2+EPS),1, 1,2,2, 1,3,2, 0,3,2, 0,2,2
  };

  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1,&block2,&block3,&block4};
}

stk::mesh::PartVector setup_mesh_1block_four_hex_2node_one_edge_hinge_manual(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_1", 1);
  std::string meshDesc = "2,1,HEX_8,1,2,3,4,7,8,9,10,block_1\n"
                         "2,2,HEX_8,27,28,5,6,7,10,11,12,block_1\n"
                         "0,3,HEX_8,7,13,14,15,16,17,18,19,block_1\n"
                         "1,4,HEX_8,10,20,21,22,23,24,25,26,block_1";
  std::vector<double> coordinates = {
    (1+EPS),1,0, 2,1,0, 2,2,0, (1+EPS),2,0, 0,2,0, 0,1,0,
    1,1,1, 2,1,1, 2,2,1, 1,2,1, 0,2,1, 0,1,1,
    0,(1-EPS),1, 0,0,1, 1,0,1, 1,1,2, 0,1,2, 0,0,2, 1,0,2,
    1,3,1, 0,3,1, 0,(2+EPS),1, 1,2,2, 1,3,2, 0,3,2, 0,2,2,
    (1-EPS),1,0, (1-EPS),2,0
  };

  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1};
}

stk::mesh::PartVector setup_mesh_1block_four_hex_2node_one_edge_hinge(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_1", 1);
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,7,8,9,10,block_1\n"
                         "0,2,HEX_8,27,28,5,6,7,10,11,12,block_1\n"
                         "0,3,HEX_8,7,13,14,15,16,17,18,19,block_1\n"
                         "0,4,HEX_8,10,20,21,22,23,24,25,26,block_1";
  std::vector<double> coordinates = {
    (1+EPS),1,0, 2,1,0, 2,2,0, (1+EPS),2,0, 0,2,0, 0,1,0,
    1,1,1, 2,1,1, 2,2,1, 1,2,1, 0,2,1, 0,1,1,
    0,(1-EPS),1, 0,0,1, 1,0,1, 1,1,2, 0,1,2, 0,0,2, 1,0,2,
    1,3,1, 0,3,1, 0,(2+EPS),1, 1,2,2, 1,3,2, 0,3,2, 0,2,2,
    (1-EPS),1,0, (1-EPS),2,0
  };

  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1};
}

stk::mesh::PartVector setup_mesh_1block_eight_tri_1node_hinge(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::TRI_3_2D, "block_1", 1);
  std::string meshDesc = "0,1,TRI_3_2D,1,2,7,block_1\n"
                         "0,2,TRI_3_2D,2,3,7,block_1\n"
                         "0,3,TRI_3_2D,4,7,6,block_1\n"
                         "0,4,TRI_3_2D,5,8,7,block_1\n"
                         "0,5,TRI_3_2D,6,7,9,block_1\n"
                         "0,6,TRI_3_2D,7,8,13,block_1\n"
                         "0,7,TRI_3_2D,7,11,10,block_1\n"
                         "0,8,TRI_3_2D,7,12,11,block_1";
  std::vector<double> coordinates = {
    0,-EPS, 1,0, 2,-EPS,
    0,EPS, 2,EPS, 0,1,
    1,1, 2,1, 0,(2-EPS),
    0,(2+EPS), 1,2, 2,(2+EPS),
    2,(2-EPS)
  };

  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1};
}

stk::mesh::PartVector setup_mesh_4block_4quad_bowtie_1hinge(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);
  stk::mesh::Part & block3 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_3", 3);
  stk::mesh::Part & block4 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_4", 4);
  std::string meshDesc = "0,1,QUAD_4_2D,1,2,4,3,block_1\n"
                         "0,2,QUAD_4_2D,5,6,7,4,block_2\n"
                         "0,3,QUAD_4_2D,4,8,9,10,block_3\n"
                         "0,4,QUAD_4_2D,4,11,12,13,block_4";
  std::vector<double> coordinates = {
    0,0, (1-EPS),0, 0,(1-EPS), 1,1, (1+EPS),0,
    2,0, 2,(1-EPS), 2,(1+EPS), 2,2, (1+EPS),2, (1-EPS),2, 0,2, 0,(1+EPS)
  };
  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1,&block2,&block3,&block4};
}

stk::mesh::PartVector setup_mesh_3block_3quad_1hinge(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);
  stk::mesh::Part & block3 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_3", 3);
  std::string meshDesc = "0,1,QUAD_4_2D,1,2,4,3,block_1\n"
                         "0,2,QUAD_4_2D,2,5,6,4,block_2\n"
                         "0,3,QUAD_4_2D,4,7,8,9,block_3";
  std::vector<double> coordinates = { 0,0, 1,0, 0,1, 1,1, 2,0, 2,(1-EPS), 2,(1+EPS), 2,2, 1,2 };
  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  return {&block1,&block2,&block3};
}


void print_hinge_info(const stk::mesh::BulkData& bulk,
                      const stk::tools::impl::HingeNodeVector& hingeNodes,
                      const stk::tools::impl::HingeEdgeVector& hingeEdges)
{
  std::ostringstream os;
  if(hingeNodes.size() > 0) {
    os << "PRINTING HINGE NODES on Proc " << bulk.parallel_rank() << " : " << std::endl;
    for(const stk::tools::impl::HingeNode& node : hingeNodes) {
      if(hinge_node_is_locally_owned(bulk, node)) {
        os << "\tHinge node id: " << bulk.identifier(node.get_node()) << std::endl;
      }
    }
  }
  if(hingeEdges.size() > 0) {
    os << "PRINTING HINGE EDGES on Proc " << bulk.parallel_rank() << " : " << std::endl;
    for(const stk::tools::impl::HingeEdge& edge : hingeEdges) {
      if(hinge_edge_is_locally_owned(bulk, edge)) {
        os << "\tHinge edge ids: " << bulk.identifier(edge.first.get_node())
           << ", " << bulk.identifier(edge.second.get_node()) << std::endl;
      }
    }
  }

  for(int i = 0; i < bulk.parallel_size(); i++) {
    if(i == bulk.parallel_rank()) {
      std::cout << os.str() << std::endl;
    }
    MPI_Barrier(bulk.parallel());
  }
}

void output_mesh(stk::mesh::BulkData & bulk)
{
  std::string writeOutput = stk::unit_test_util::get_option("--output", "off");
  const std::string fileName = std::string(::testing::UnitTest::GetInstance()->current_test_info()->name()) + ".g";
  if (writeOutput == "on") {
    stk::io::write_mesh(fileName, bulk);
  }
}

bool is_debug()
{
  return stk::unit_test_util::has_option("--debug");
}

// Common Decompositions
void two_elements_decomposition(stk::mesh::BulkData& bulk)
{
  if(bulk.parallel_size() == 2) {
    stk::mesh::EntityIdProcVec idProcVec{ {2u,1} };
    distribute_mesh(bulk, idProcVec);
  }
}
void three_elements_decomposition(stk::mesh::BulkData& bulk)
{
  if(bulk.parallel_size() == 2) {
    stk::mesh::EntityIdProcVec idProcVec{ {3u,1} };
    distribute_mesh(bulk, idProcVec);
  }
  else if(bulk.parallel_size() == 3) {
    stk::mesh::EntityIdProcVec idProcVec{ {2u,1}, {3u,2} };
    distribute_mesh(bulk, idProcVec);
  }
}

void four_elements_decomposition(stk::mesh::BulkData& bulk)
{
  if(bulk.parallel_size() == 2) {
    stk::mesh::EntityIdProcVec idProcVec{ {2u,1}, {3u,1} };
    distribute_mesh(bulk, idProcVec);
  }
  else if(bulk.parallel_size() == 3) {
    stk::mesh::EntityIdProcVec idProcVec{ {3u,2} };
    distribute_mesh(bulk, idProcVec);
  }
  else if(bulk.parallel_size() == 4) {
    stk::mesh::EntityIdProcVec idProcVec{ {2u,1}, {3u,2}, {4u,3} };
    distribute_mesh(bulk, idProcVec);
  }
}

void four_elements_decomposition2(stk::mesh::BulkData& bulk)
{
  if(bulk.parallel_size() == 2) {
    stk::mesh::EntityIdProcVec idProcVec{ {1u,1}, {3u,1} };
    distribute_mesh(bulk, idProcVec);
  }
  else if(bulk.parallel_size() == 3) {
    stk::mesh::EntityIdProcVec idProcVec{ {2u,1}, {3u,2} };
    distribute_mesh(bulk, idProcVec);
  }
  else if(bulk.parallel_size() == 4) {
    stk::mesh::EntityIdProcVec idProcVec{ {2u,1}, {3u,2}, {4u,3} };
    distribute_mesh(bulk, idProcVec);
  }
}

void five_elements_decomposition(stk::mesh::BulkData& bulk)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 5)
    return;

  if(bulk.parallel_size() == 2) {
    stk::mesh::EntityIdProcVec idProcVec{ {1u,1}, {3u,1} };
    distribute_mesh(bulk, idProcVec);
  }
  else if(bulk.parallel_size() == 3) {
    stk::mesh::EntityIdProcVec idProcVec{ {2u,1}, {3u,2} };
    distribute_mesh(bulk, idProcVec);
  }
  else if(bulk.parallel_size() == 4) {
    stk::mesh::EntityIdProcVec idProcVec{ {2u,1}, {3u,2}, {4u,3} };
    distribute_mesh(bulk, idProcVec);
  }
  else if(bulk.parallel_size() == 5) {
    stk::mesh::EntityIdProcVec idProcVec{ {2u,1}, {3u,2}, {4u,3}, {5u,4} };
    distribute_mesh(bulk, idProcVec);
  }
}

void verify_test_run(stk::ParallelMachine pm, int locallyRanTest)
{
  int globallyRanTest;
  stk::all_reduce_max(pm, &locallyRanTest, &globallyRanTest, 1);
  EXPECT_EQ(1, globallyRanTest);
}

template<typename T>
void test_two_element_one_hinge_grouping(const stk::mesh::BulkData& bulk, const T& hinge)
{
  stk::mesh::Entity elem1 = bulk.get_entity(stk::topology::ELEMENT_RANK, 1u);
  stk::mesh::Entity elem2 = bulk.get_entity(stk::topology::ELEMENT_RANK, 2u);
  stk::tools::impl::HingeGroupVector groupings = stk::tools::impl::get_convex_groupings(bulk, hinge);
  int ranTest = 0;

  if(!groupings.empty()) {
    ranTest = 1;
    EXPECT_EQ(groupings.size(), 2u);
    EXPECT_EQ(groupings[0].size(), 1u);
    EXPECT_EQ(groupings[1].size(), 1u);
    if(groupings[0][0] == elem1) {
      EXPECT_EQ(groupings[1][0], elem2);
    }
    if(groupings[0][0] == elem2) {
      EXPECT_EQ(groupings[1][0], elem1);
    }
  }

  verify_test_run(bulk.parallel(), ranTest);
}

stk::mesh::EntityVector get_entities_from_id_range(const stk::mesh::BulkData& bulk, stk::topology::rank_t rank, unsigned count)
{
  stk::mesh::EntityVector returnVec;
  for(unsigned i = 1; i <= count; i++)
  {
    stk::mesh::Entity elem = bulk.get_entity(rank, i);
    returnVec.push_back(elem);
  }
  return returnVec;
}

stk::mesh::EntityVector get_nodes_from_id_range(const stk::mesh::BulkData& bulk, unsigned count)
{
  return get_entities_from_id_range(bulk, stk::topology::NODE_RANK, count);
}

stk::mesh::EntityVector get_elements_from_id_range(const stk::mesh::BulkData& bulk, unsigned count)
{
  return get_entities_from_id_range(bulk, stk::topology::ELEMENT_RANK, count);
}

//void snip_and_get_remaining_hinge_count(stk::mesh::BulkData& bulk,
//                                        stk::mesh::EntityVector& elemVec,
//                                        stk::mesh::EntityVector& nodeVec,
//                                        std::pair<unsigned,unsigned>& hingeCounts)
//{
//  stk::tools::impl::snip_all_hinges_between_blocks(bulk);
//
//  stk::mesh::get_entities(bulk, stk::topology::ELEMENT_RANK, elemVec);
//  stk::mesh::get_entities(bulk, stk::topology::NODE_RANK, nodeVec);
//
//  hingeCounts = stk::tools::impl::get_hinge_count(bulk);
//}

//std::pair<unsigned,unsigned> get_locally_owned_elem_node_pair(const stk::mesh::BulkData& bulk)
//{
//  const stk::mesh::BucketVector& elemBuckets = bulk.get_buckets(stk::topology::ELEMENT_RANK, bulk.mesh_meta_data().locally_owned_part());
//  const stk::mesh::BucketVector& nodeBuckets = bulk.get_buckets(stk::topology::NODE_RANK, bulk.mesh_meta_data().locally_owned_part());
//
//  unsigned numElems = 0;
//  unsigned numNodes = 0;
//
//  for(const stk::mesh::Bucket* elemBucket : elemBuckets) {
//    numElems += elemBucket->size();
//  }
//
//  for(const stk::mesh::Bucket* nodeBucket : nodeBuckets) {
//    numNodes += nodeBucket->size();
//  }
//
//  return std::make_pair(numElems, numNodes);
//}

std::pair<unsigned,unsigned> get_reduced_entity_counts(const stk::mesh::BulkData& bulk)
{
  // node edge face elem
  std::vector<size_t> reducedEntityCounts;

  stk::mesh::comm_mesh_counts(bulk, reducedEntityCounts);

  return std::make_pair(reducedEntityCounts[stk::topology::ELEMENT_RANK], reducedEntityCounts[stk::topology::NODE_RANK]);
}

} //namespace

//Start of Tests

TEST(DetectHinge2D, EmptyMesh)
{
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  std::pair<unsigned, unsigned> hingeCount = stk::tools::impl::get_hinge_count(bulk);
  EXPECT_EQ(0u, hingeCount.first);
  EXPECT_EQ(0u, hingeCount.second);
}

TEST(DetectHinge2D, SingleBlockNoHinge)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1)
    return;
  stk::mesh::MetaData meta(2);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_1quad(bulk);
  std::pair<unsigned, unsigned> hingeCount = stk::tools::impl::get_hinge_count(bulk);
  EXPECT_EQ(0u, hingeCount.first);
  EXPECT_EQ(0u, hingeCount.second);
  output_mesh(bulk);
}

TEST(DetectHinge2D, SingleBlockTwoElementsNoHinge)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 2)
    return;
  stk::mesh::MetaData meta(2);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_2quad(bulk);
  two_elements_decomposition(bulk);

  std::pair<unsigned, unsigned> hingeCount = stk::tools::impl::get_hinge_count(bulk);
  EXPECT_EQ(0u, hingeCount.first);
  EXPECT_EQ(0u, hingeCount.second);
  output_mesh(bulk);
}

TEST(DetectHinge2D, SingleBlockTwoElementsOneNodeHinge)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 2)
    return;
  stk::mesh::MetaData meta(2);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_2quad_1node_hinge(bulk);
  two_elements_decomposition(bulk);

  std::pair<unsigned, unsigned> hingeCount = stk::tools::impl::get_hinge_count(bulk);
  EXPECT_EQ(1u, hingeCount.first);
  EXPECT_EQ(0u, hingeCount.second);
  output_mesh(bulk);
}

TEST(DetectHinge2D, SingleBlockTwoElementsTwoNodeHinges)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 2)
    return;
  stk::mesh::MetaData meta(2);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_2quad_2hinge(bulk);
  two_elements_decomposition(bulk);

  std::pair<unsigned, unsigned> hingeCount = stk::tools::impl::get_hinge_count(bulk);
  EXPECT_EQ(2u, hingeCount.first);
  EXPECT_EQ(0u, hingeCount.second);
  output_mesh(bulk);
}

TEST(DetectHinge2D, AreNodesPartOfASide)
{
  stk::mesh::MetaData meta(2);
  stk::mesh::BulkData bulk(meta, MPI_COMM_SELF);
  setup_mesh_1block_1quad(bulk);
  stk::mesh::EntityVector elems = get_elements_from_id_range(bulk, 1);
  stk::mesh::EntityVector nodes = get_nodes_from_id_range(bulk, 4);
  EXPECT_TRUE(stk::tools::impl::common_nodes_are_part_of_a_side(bulk, elems, stk::mesh::EntityVector{nodes[0],nodes[1]}));
  EXPECT_TRUE(stk::tools::impl::common_nodes_are_part_of_a_side(bulk, elems, stk::mesh::EntityVector{nodes[1],nodes[3]}));
  EXPECT_TRUE(stk::tools::impl::common_nodes_are_part_of_a_side(bulk, elems, stk::mesh::EntityVector{nodes[3],nodes[2]}));
  EXPECT_TRUE(stk::tools::impl::common_nodes_are_part_of_a_side(bulk, elems, stk::mesh::EntityVector{nodes[2],nodes[0]}));
  EXPECT_FALSE(stk::tools::impl::common_nodes_are_part_of_a_side(bulk, elems, stk::mesh::EntityVector{nodes[0],nodes[3]}));
  EXPECT_FALSE(stk::tools::impl::common_nodes_are_part_of_a_side(bulk, elems, stk::mesh::EntityVector{nodes[1],nodes[2]}));
}

TEST(DetectHinge2D, SingleBlockThreeElementsOneHinge_Decomp1)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 3)
    return;
  stk::mesh::MetaData meta(2);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_3quad_1hinge(bulk);

  three_elements_decomposition(bulk);

  std::pair<unsigned, unsigned> hingeCount = stk::tools::impl::get_hinge_count(bulk);
  EXPECT_EQ(1u, hingeCount.first);
  EXPECT_EQ(0u, hingeCount.second);
  output_mesh(bulk);
}

TEST(DetectHinge2D, SingleBlockThreeElementsOneHinge_Decomp2)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 3)
    return;
  stk::mesh::MetaData meta(2);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_3quad_1hinge(bulk);

  if(bulk.parallel_size() == 2) {
    stk::mesh::EntityIdProcVec idProcVec{ {2u,1}, {3u,1} };
    distribute_mesh(bulk, idProcVec);
  }
  else if(bulk.parallel_size() == 3) {
    stk::mesh::EntityIdProcVec idProcVec{ {2u,1}, {3u,2} };
    distribute_mesh(bulk, idProcVec);
  }

  std::pair<unsigned, unsigned> hingeCount = stk::tools::impl::get_hinge_count(bulk);
  EXPECT_EQ(1u, hingeCount.first);
  EXPECT_EQ(0u, hingeCount.second);
  output_mesh(bulk);
}

TEST(DetectHinge2D, SingleBlockThreeElementsOneHinge_LinearStack)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 3)
    return;
  stk::mesh::MetaData meta(2);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_3quad_1hinge_linear_stack(bulk);

  three_elements_decomposition(bulk);

  std::pair<unsigned, unsigned> hingeCount = stk::tools::impl::get_hinge_count(bulk);
  EXPECT_EQ(1u, hingeCount.first);
  EXPECT_EQ(0u, hingeCount.second);
  output_mesh(bulk);
}

TEST(DetectHinge2D, SingleBlockFourElementsBowtie_Decomp1)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 4)
    return;
  stk::mesh::MetaData meta(2);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_4quad_bowtie_1hinge(bulk);

  if(bulk.parallel_size() == 2) {
    stk::mesh::EntityIdProcVec idProcVec{ {3u,1} };
    distribute_mesh(bulk, idProcVec);
  }
  else if(bulk.parallel_size() == 3) {
    stk::mesh::EntityIdProcVec idProcVec{ {2u,1}, {3u,2} };
    distribute_mesh(bulk, idProcVec);
  }
  else if(bulk.parallel_size() == 4) {
    stk::mesh::EntityIdProcVec idProcVec{ {2u,1}, {3u,2}, {4u,3} };
    distribute_mesh(bulk, idProcVec);
  }

  std::pair<unsigned, unsigned> hingeCount = stk::tools::impl::get_hinge_count(bulk);
  EXPECT_EQ(1u, hingeCount.first);
  EXPECT_EQ(0u, hingeCount.second);
  output_mesh(bulk);
}

TEST(DetectHinge2D, SingleBlockFourElementsBowtie_Decomp2)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 4)
    return;
  stk::mesh::MetaData meta(2);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_4quad_bowtie_1hinge(bulk);

  if(bulk.parallel_size() == 2) {
    stk::mesh::EntityIdProcVec idProcVec{ {3u,1}, {4u,1} };
    distribute_mesh(bulk, idProcVec);
  }
  else if(bulk.parallel_size() == 3) {
    stk::mesh::EntityIdProcVec idProcVec{ {3u,2} };
    distribute_mesh(bulk, idProcVec);
  }
  else if(bulk.parallel_size() == 4) {
    stk::mesh::EntityIdProcVec idProcVec{ {2u,1}, {3u,2}, {4u,3} };
    distribute_mesh(bulk, idProcVec);
  }

  std::pair<unsigned, unsigned> hingeCount = stk::tools::impl::get_hinge_count(bulk);
  EXPECT_EQ(1u, hingeCount.first);
  EXPECT_EQ(0u, hingeCount.second);
  output_mesh(bulk);
}

TEST(DetectHinge2D, SingleBlockFourElementsTwoHinge_Decomp1)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 4)
    return;
  stk::mesh::MetaData meta(2);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_4quad_2hinge(bulk);

  if(bulk.parallel_size() == 2) {
    stk::mesh::EntityIdProcVec idProcVec{ {3u,1}, {4u,1} };
    distribute_mesh(bulk, idProcVec);
  }
  else if(bulk.parallel_size() == 3) {
    stk::mesh::EntityIdProcVec idProcVec{ {2u,1}, {3u,2} };
    distribute_mesh(bulk, idProcVec);
  }
  else if(bulk.parallel_size() == 4) {
    stk::mesh::EntityIdProcVec idProcVec{ {2u,1}, {3u,2}, {4u,3} };
    distribute_mesh(bulk, idProcVec);
  }

  std::pair<unsigned, unsigned> hingeCount = stk::tools::impl::get_hinge_count(bulk);
  EXPECT_EQ(2u, hingeCount.first);
  EXPECT_EQ(0u, hingeCount.second);
  output_mesh(bulk);
}

TEST(DetectHinge2D, SingleBlockFourElementsTwoHinge_Decomp2)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 4)
    return;
  stk::mesh::MetaData meta(2);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_4quad_2hinge(bulk);

  four_elements_decomposition(bulk);

  std::pair<unsigned, unsigned> hingeCount = stk::tools::impl::get_hinge_count(bulk);
  EXPECT_EQ(2u, hingeCount.first);
  EXPECT_EQ(0u, hingeCount.second);
  output_mesh(bulk);
}

TEST(DetectHinge2D, SingleBlockFourElementsFourHinge_Decomp1)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 4)
    return;
  stk::mesh::MetaData meta(2);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_4quad_4hinge(bulk);

  four_elements_decomposition(bulk);

  std::pair<unsigned, unsigned> hingeCount = stk::tools::impl::get_hinge_count(bulk);
  EXPECT_EQ(4u, hingeCount.first);
  EXPECT_EQ(0u, hingeCount.second);
  output_mesh(bulk);
}

TEST(DetectHinge2D, SingleBlockFourElementsFourHinge_Decomp2)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 4)
    return;
  stk::mesh::MetaData meta(2);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_4quad_4hinge(bulk);

  if(bulk.parallel_size() == 2) {
    stk::mesh::EntityIdProcVec idProcVec{ {1u,1}, {3u,1} };
    distribute_mesh(bulk, idProcVec);
  }
  else if(bulk.parallel_size() == 3) {
    stk::mesh::EntityIdProcVec idProcVec{ {3u,2} };
    distribute_mesh(bulk, idProcVec);
  }
  else if(bulk.parallel_size() == 4) {
    stk::mesh::EntityIdProcVec idProcVec{ {2u,1}, {3u,2}, {4u,3} };
    distribute_mesh(bulk, idProcVec);
  }

  std::pair<unsigned, unsigned> hingeCount = stk::tools::impl::get_hinge_count(bulk);
  EXPECT_EQ(4u, hingeCount.first);
  EXPECT_EQ(0u, hingeCount.second);
  output_mesh(bulk);
}

TEST(DetectHinge2D, SingleBlockFourElementsNoHingePacman)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 4)
    return;
  stk::mesh::MetaData meta(2);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_4quad_pacman(bulk);

  four_elements_decomposition2(bulk);

  std::pair<unsigned, unsigned> hingeCount = stk::tools::impl::get_hinge_count(bulk);
  EXPECT_EQ(0u, hingeCount.first);
  EXPECT_EQ(0u, hingeCount.second);
  output_mesh(bulk);
}

TEST(DetectHinge2D, SingleBlockFourElementsOneHinge)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 4)
    return;
  stk::mesh::MetaData meta(2);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_4quad_1hinge(bulk);

  four_elements_decomposition2(bulk);

  std::pair<unsigned, unsigned> hingeCount = stk::tools::impl::get_hinge_count(bulk);
  EXPECT_EQ(1u, hingeCount.first);
  EXPECT_EQ(0u, hingeCount.second);
  output_mesh(bulk);
}

TEST(DetectHinge2D, TwoBlockFiveElementsOneHinge)
{
  stk::mesh::MetaData meta(2);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_2block_3quad_2tri_1hinge(bulk);
  five_elements_decomposition(bulk);

  std::pair<unsigned, unsigned> hingeCount = stk::tools::impl::get_hinge_count(bulk);
  EXPECT_EQ(1u, hingeCount.first);
  EXPECT_EQ(0u, hingeCount.second);
  output_mesh(bulk);
}

TEST(DetectHinge3D, SingleBlockNoHinge)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1)
    return;
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_1hex(bulk);
  std::pair<unsigned, unsigned> hingeCount = stk::tools::impl::get_hinge_count(bulk);
  EXPECT_EQ(0u, hingeCount.first);
  EXPECT_EQ(0u, hingeCount.second);
  output_mesh(bulk);
}

TEST(DetectHinge3D, SingleBlockTwoElementsNoHinge)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 2)
    return;
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_2hex(bulk);

  if(bulk.parallel_size() == 2) {
    stk::mesh::EntityIdProcVec idProcVec{ {2u,1} };
    distribute_mesh(bulk, idProcVec);
  }

  std::pair<unsigned, unsigned> hingeCount = stk::tools::impl::get_hinge_count(bulk);
  EXPECT_EQ(0u, hingeCount.first);
  EXPECT_EQ(0u, hingeCount.second);
  output_mesh(bulk);
}

TEST(DetectHinge3D, SingleBlockTwoElementsNoHingeFaceTest)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1)
    return;
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_2hex_face_test(bulk);

  stk::mesh::Entity node = bulk.get_entity(stk::topology::NODE_RANK, 7u);
  stk::tools::impl::PairwiseSideInfoVector infoVec = stk::tools::impl::get_hinge_info_vec(bulk, node);

  EXPECT_EQ(1u, get_side_count(infoVec));
  output_mesh(bulk);
}

TEST(DetectHinge3D, SingleBlockTwoElementsOneNodeHinge)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 2)
    return;
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_2hex_1node_hinge(bulk);

  if(bulk.parallel_size() == 2) {
    stk::mesh::EntityIdProcVec idProcVec{ {2u,1} };
    distribute_mesh(bulk, idProcVec);
  }

  std::pair<unsigned, unsigned> hingeCount = stk::tools::impl::get_hinge_count(bulk);
  EXPECT_EQ(1u, hingeCount.first);
  EXPECT_EQ(0u, hingeCount.second);
  output_mesh(bulk);
}

TEST(DetectHinge3D, SingleBlockTwoElementsTwoNodeHinge)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 2)
    return;
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_2hex_2node_hinge(bulk);

  if(bulk.parallel_size() == 2) {
    stk::mesh::EntityIdProcVec idProcVec{ {2u,1} };
    distribute_mesh(bulk, idProcVec);
  }

  std::pair<unsigned, unsigned> hingeCount = stk::tools::impl::get_hinge_count(bulk);
  EXPECT_EQ(2u, hingeCount.first);
  EXPECT_EQ(0u, hingeCount.second);
  output_mesh(bulk);
}

TEST(DetectHinge3D, SingleBlockThreeElementsOneNodeHinge)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 3)
    return;
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_3hex_1node_hinge(bulk);

  three_elements_decomposition(bulk);

  std::pair<unsigned, unsigned> hingeCount = stk::tools::impl::get_hinge_count(bulk);
  EXPECT_EQ(1u, hingeCount.first);
  EXPECT_EQ(0u, hingeCount.second);
  output_mesh(bulk);
}


TEST(DetectHinge3D, SingleBlockEightElementsOneNodeHingeFlower)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 8)
    return;
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_8hex_flower_1node_hinge(bulk);

  if(bulk.parallel_size() == 2) {
    stk::mesh::EntityIdProcVec idProcVec{ {5u,1}, {6u,1}, {7u,1}, {8u,1} };
    distribute_mesh(bulk, idProcVec);
  }
  else if(bulk.parallel_size() == 3) {
    stk::mesh::EntityIdProcVec idProcVec{ {3u,1}, {4u,1}, {5u,1}, {6u,2}, {7u,2}, {8u,2} };
    distribute_mesh(bulk, idProcVec);
  }
  else if(bulk.parallel_size() == 4) {
    stk::mesh::EntityIdProcVec idProcVec{ {3u,1}, {4u,1}, {5u,2}, {6u,2}, {7u,3}, {8u,3} };
    distribute_mesh(bulk, idProcVec);
  }

  std::pair<unsigned, unsigned> hingeCount = stk::tools::impl::get_hinge_count(bulk);
  EXPECT_EQ(1u, hingeCount.first);
  EXPECT_EQ(0u, hingeCount.second);
  output_mesh(bulk);
}

TEST(DetectHinge3D, SingleBlockTwoTetsOneNodeHinge)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 2)
    return;
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_2tet_1node_hinge(bulk);

  two_elements_decomposition(bulk);

  std::pair<unsigned, unsigned> hingeCount = stk::tools::impl::get_hinge_count(bulk);
  EXPECT_EQ(1u, hingeCount.first);
  EXPECT_EQ(0u, hingeCount.second);
  output_mesh(bulk);
}

TEST(DetectHinge3D, SingleBlockTwoHexOneEdgeHinge)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 2)
    return;
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_2hex_1edge_hinge(bulk);

  two_elements_decomposition(bulk);

  std::pair<unsigned, unsigned> hingeCount = stk::tools::impl::get_hinge_count(bulk);
  EXPECT_EQ(0u, hingeCount.first);
  EXPECT_EQ(1u, hingeCount.second);
  output_mesh(bulk);
}

TEST(DetectHinge3D, SingleBlockThreeHexOneEdgeHinge)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 3)
    return;
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_3hex_1edge_hinge(bulk);

  three_elements_decomposition(bulk);

  std::pair<unsigned, unsigned> hingeCount = stk::tools::impl::get_hinge_count(bulk);
  EXPECT_EQ(0u, hingeCount.first);
  EXPECT_EQ(1u, hingeCount.second);
  output_mesh(bulk);
}

TEST(DetectHinge3D, AreNodesPartOfASide)
{
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_SELF);
  setup_mesh_1block_1hex(bulk);
  stk::mesh::EntityVector nodes = get_nodes_from_id_range(bulk, 8);
  stk::mesh::EntityVector elems = get_elements_from_id_range(bulk, 1);
  EXPECT_TRUE(stk::tools::impl::common_nodes_are_part_of_a_side (bulk, elems, stk::mesh::EntityVector{nodes[0],nodes[1],nodes[2],nodes[3]}));
  EXPECT_TRUE(stk::tools::impl::common_nodes_are_part_of_a_side (bulk, elems, stk::mesh::EntityVector{nodes[4],nodes[5],nodes[6],nodes[7]}));
  EXPECT_TRUE(stk::tools::impl::common_nodes_are_part_of_a_side (bulk, elems, stk::mesh::EntityVector{nodes[0],nodes[1],nodes[5],nodes[4]}));
  EXPECT_TRUE(stk::tools::impl::common_nodes_are_part_of_a_side (bulk, elems, stk::mesh::EntityVector{nodes[3],nodes[2],nodes[6],nodes[7]}));
  EXPECT_TRUE(stk::tools::impl::common_nodes_are_part_of_a_side (bulk, elems, stk::mesh::EntityVector{nodes[1],nodes[2],nodes[6],nodes[5]}));
  EXPECT_TRUE(stk::tools::impl::common_nodes_are_part_of_a_side (bulk, elems, stk::mesh::EntityVector{nodes[0],nodes[3],nodes[7],nodes[4]}));

  EXPECT_FALSE(stk::tools::impl::common_nodes_are_part_of_a_side (bulk, elems, stk::mesh::EntityVector{nodes[0],nodes[2],nodes[6],nodes[4]}));
  EXPECT_FALSE(stk::tools::impl::common_nodes_are_part_of_a_side (bulk, elems, stk::mesh::EntityVector{nodes[1],nodes[5],nodes[7],nodes[3]}));
  EXPECT_FALSE(stk::tools::impl::common_nodes_are_part_of_a_side (bulk, elems, stk::mesh::EntityVector{nodes[1],nodes[2],nodes[7],nodes[4]}));
  EXPECT_FALSE(stk::tools::impl::common_nodes_are_part_of_a_side (bulk, elems, stk::mesh::EntityVector{nodes[0],nodes[5],nodes[6],nodes[3]}));
  EXPECT_FALSE(stk::tools::impl::common_nodes_are_part_of_a_side (bulk, elems, stk::mesh::EntityVector{nodes[1],nodes[6],nodes[7],nodes[0]}));
  EXPECT_FALSE(stk::tools::impl::common_nodes_are_part_of_a_side (bulk, elems, stk::mesh::EntityVector{nodes[5],nodes[2],nodes[3],nodes[4]}));
}

TEST(DetectHinge3D, AreNodesPartOfAnEdge)
{
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_SELF);
  setup_mesh_1block_1hex(bulk);

  stk::mesh::EntityVector nodes = get_nodes_from_id_range(bulk, 8);

  EXPECT_TRUE(stk::tools::impl::common_nodes_are_part_of_an_edge(bulk, nodes[0], nodes[1]));
  EXPECT_TRUE(stk::tools::impl::common_nodes_are_part_of_an_edge(bulk, nodes[1], nodes[2]));
  EXPECT_TRUE(stk::tools::impl::common_nodes_are_part_of_an_edge(bulk, nodes[2], nodes[3]));
  EXPECT_TRUE(stk::tools::impl::common_nodes_are_part_of_an_edge(bulk, nodes[3], nodes[0]));
  EXPECT_TRUE(stk::tools::impl::common_nodes_are_part_of_an_edge(bulk, nodes[4], nodes[5]));
  EXPECT_TRUE(stk::tools::impl::common_nodes_are_part_of_an_edge(bulk, nodes[5], nodes[6]));
  EXPECT_TRUE(stk::tools::impl::common_nodes_are_part_of_an_edge(bulk, nodes[6], nodes[7]));
  EXPECT_TRUE(stk::tools::impl::common_nodes_are_part_of_an_edge(bulk, nodes[7], nodes[3]));
  EXPECT_TRUE(stk::tools::impl::common_nodes_are_part_of_an_edge(bulk, nodes[0], nodes[4]));
  EXPECT_TRUE(stk::tools::impl::common_nodes_are_part_of_an_edge(bulk, nodes[1], nodes[5]));
  EXPECT_TRUE(stk::tools::impl::common_nodes_are_part_of_an_edge(bulk, nodes[2], nodes[6]));
  EXPECT_TRUE(stk::tools::impl::common_nodes_are_part_of_an_edge(bulk, nodes[3], nodes[7]));

  EXPECT_FALSE(stk::tools::impl::common_nodes_are_part_of_an_edge(bulk, nodes[0], nodes[2]));
  EXPECT_FALSE(stk::tools::impl::common_nodes_are_part_of_an_edge(bulk, nodes[1], nodes[3]));
  EXPECT_FALSE(stk::tools::impl::common_nodes_are_part_of_an_edge(bulk, nodes[4], nodes[6]));
  EXPECT_FALSE(stk::tools::impl::common_nodes_are_part_of_an_edge(bulk, nodes[5], nodes[7]));
  EXPECT_FALSE(stk::tools::impl::common_nodes_are_part_of_an_edge(bulk, nodes[1], nodes[6]));
  EXPECT_FALSE(stk::tools::impl::common_nodes_are_part_of_an_edge(bulk, nodes[2], nodes[5]));
  EXPECT_FALSE(stk::tools::impl::common_nodes_are_part_of_an_edge(bulk, nodes[0], nodes[7]));
  EXPECT_FALSE(stk::tools::impl::common_nodes_are_part_of_an_edge(bulk, nodes[3], nodes[4]));
  EXPECT_FALSE(stk::tools::impl::common_nodes_are_part_of_an_edge(bulk, nodes[2], nodes[7]));
  EXPECT_FALSE(stk::tools::impl::common_nodes_are_part_of_an_edge(bulk, nodes[3], nodes[6]));
  EXPECT_FALSE(stk::tools::impl::common_nodes_are_part_of_an_edge(bulk, nodes[0], nodes[5]));
  EXPECT_FALSE(stk::tools::impl::common_nodes_are_part_of_an_edge(bulk, nodes[1], nodes[4]));
  EXPECT_FALSE(stk::tools::impl::common_nodes_are_part_of_an_edge(bulk, nodes[3], nodes[5]));
  EXPECT_FALSE(stk::tools::impl::common_nodes_are_part_of_an_edge(bulk, nodes[0], nodes[6]));
  EXPECT_FALSE(stk::tools::impl::common_nodes_are_part_of_an_edge(bulk, nodes[1], nodes[7]));
  EXPECT_FALSE(stk::tools::impl::common_nodes_are_part_of_an_edge(bulk, nodes[2], nodes[4]));
}


TEST(DetectHinge3D, SingleBlockThreeElementsOneNodeHingeOneEdgeHinge)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 3)
    return;
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_3hex_1node_hinge_1edge_hinge(bulk);

  three_elements_decomposition(bulk);

  std::pair<unsigned, unsigned> hingeCount = stk::tools::impl::get_hinge_count(bulk);
  EXPECT_EQ(1u, hingeCount.first);
  EXPECT_EQ(1u, hingeCount.second);
  output_mesh(bulk);
}


TEST(DetectHinge3D, SingleBlockThreeElementsOneNodeHingeOneEdgeHinge2)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 3)
    return;
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_3hex_1node_hinge_1edge_hinge2(bulk);

  three_elements_decomposition(bulk);

  std::pair<unsigned, unsigned> hingeCount = stk::tools::impl::get_hinge_count(bulk);
  EXPECT_EQ(1u, hingeCount.first);
  EXPECT_EQ(1u, hingeCount.second);
  output_mesh(bulk);
}

TEST(DetectHinge3D, SingleBlockThreeElementsOneNodeHingeOneEdgeHinge3)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 3)
    return;
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_3hex_1node_hinge_1edge_hinge3(bulk);

  three_elements_decomposition(bulk);

  std::pair<unsigned, unsigned> hingeCount = stk::tools::impl::get_hinge_count(bulk);
  EXPECT_EQ(1u, hingeCount.first);
  EXPECT_EQ(1u, hingeCount.second);
  output_mesh(bulk);
}


TEST(DetectHinge3D, SingleBlockFourElementsBowtieOneEdgeHinge)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 4)
    return;
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_4hex_bowtie_1edge_hinge(bulk);

  four_elements_decomposition(bulk);

  std::pair<unsigned, unsigned> hingeCount = stk::tools::impl::get_hinge_count(bulk);
  EXPECT_EQ(0u, hingeCount.first);
  EXPECT_EQ(1u, hingeCount.second);
  output_mesh(bulk);
}


TEST(DetectHinge3D, SingleBlockTwoByTwoHexTwoEdgeHinge_Decomp1)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 4)
    return;
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_two_by_two_hex_2edge_hinge(bulk);

  four_elements_decomposition(bulk);

  std::pair<unsigned, unsigned> hingeCount = stk::tools::impl::get_hinge_count(bulk);
  EXPECT_EQ(0u, hingeCount.first);
  EXPECT_EQ(2u, hingeCount.second);
  output_mesh(bulk);
}

TEST(DetectHinge3D, SingleBlockTwoByTwoHexTwoEdgeHinge_Decomp2)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 4)
    return;
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_two_by_two_hex_2edge_hinge(bulk);

  four_elements_decomposition2(bulk);

  std::pair<unsigned, unsigned> hingeCount = stk::tools::impl::get_hinge_count(bulk);
  EXPECT_EQ(0u, hingeCount.first);
  EXPECT_EQ(2u, hingeCount.second);
  output_mesh(bulk);
}

TEST(DetectHinge3D, SingleBlockFourHexOneEdgeOneNodeHinge_Decomp1)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 4)
    return;
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_four_hex_one_edge_one_node_hinge(bulk);

  four_elements_decomposition(bulk);

  std::pair<unsigned, unsigned> hingeCount = stk::tools::impl::get_hinge_count(bulk);
  EXPECT_EQ(1u, hingeCount.first);
  EXPECT_EQ(1u, hingeCount.second);
  output_mesh(bulk);
}

TEST(DetectHinge3D, SingleBlockFourHexOneEdgeOneNodeHinge_Decomp2)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 4)
    return;
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_four_hex_one_edge_one_node_hinge(bulk);

  four_elements_decomposition2(bulk);

  std::pair<unsigned, unsigned> hingeCount = stk::tools::impl::get_hinge_count(bulk);
  EXPECT_EQ(1u, hingeCount.first);
  EXPECT_EQ(1u, hingeCount.second);
  output_mesh(bulk);
}

TEST(DetectHinge3D, SingleBlockFourHexTwoNodeHinges_Decomp1)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 4)
    return;
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_four_hex_2node_hinge(bulk);

  four_elements_decomposition(bulk);

  std::pair<unsigned, unsigned> hingeCount = stk::tools::impl::get_hinge_count(bulk);
  EXPECT_EQ(2u, hingeCount.first);
  EXPECT_EQ(0u, hingeCount.second);
  output_mesh(bulk);
}

TEST(DetectHinge3D, SingleBlockFourHexTwoNodeHinges_Decomp2)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 4)
    return;
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_four_hex_2node_hinge(bulk);

  four_elements_decomposition2(bulk);

  std::pair<unsigned, unsigned> hingeCount = stk::tools::impl::get_hinge_count(bulk);
  EXPECT_EQ(2u, hingeCount.first);
  EXPECT_EQ(0u, hingeCount.second);
  output_mesh(bulk);
}

TEST(DetectHinge3D, inputFile)
{
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  std::string inputFileName = stk::unit_test_util::get_option("--inputFile", "");
  if(!inputFileName.empty()) {
    double startTime = stk::wall_time();
    stk::io::fill_mesh(inputFileName, bulk);
    double meshReadTime = stk::wall_time();
    stk::tools::impl::HingeNodeVector hingeNodes;
    stk::tools::impl::HingeEdgeVector hingeEdges;
    fill_mesh_hinges(bulk, hingeNodes, hingeEdges);
    double detectTime = stk::wall_time();
    print_hinge_info(bulk, hingeNodes, hingeEdges);
    output_mesh(bulk);
    double meshWriteTime = stk::wall_time();

    if (bulk.parallel_rank() == 0) {
      std::cout << " Mesh read time = " << (meshReadTime - startTime) << " s" << std::endl;
      std::cout << "Disconnect time = " << (detectTime - meshReadTime) << " s" << std::endl;
      std::cout << "Mesh write time = " << (meshWriteTime - detectTime) << " s" << std::endl;
    }
  }
}

TEST(DetectHinge3D, GeneratedMesh)
{
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  std::ostringstream os;
  unsigned nproc = stk::parallel_machine_size(MPI_COMM_WORLD);
  os << "generated:" << nproc << "x" << nproc << "x" << nproc;
  stk::io::fill_mesh(os.str(), bulk);
  stk::tools::impl::HingeNodeVector hingeNodes;
  stk::tools::impl::HingeEdgeVector hingeEdges;
  fill_mesh_hinges(bulk, hingeNodes, hingeEdges);
  print_hinge_info(bulk, hingeNodes, hingeEdges);

  EXPECT_EQ(0u, hingeNodes.size());
  EXPECT_EQ(0u, hingeEdges.size());
  output_mesh(bulk);
}

TEST(DetectHinge3D, SingleBlockFourHexTwoNodeHingeOneEdgeHinge)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 3)
    return;
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_four_hex_2node_one_edge_hinge(bulk);

  stk::mesh::EntityIdProcVec idProcVec{ {1u,2}, {2u,2}, {4u,1} };
  distribute_mesh(bulk, idProcVec);

  std::pair<unsigned, unsigned> hingeCount = stk::tools::impl::get_hinge_count(bulk);
  EXPECT_EQ(2u, hingeCount.first);
  EXPECT_EQ(1u, hingeCount.second);
  output_mesh(bulk);
}

TEST(DetectHinge3D, SingleBlockFourHexTwoNodeHingeOneEdgeHinge_Manual)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 3)
    return;
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_four_hex_2node_one_edge_hinge_manual(bulk);

  std::pair<unsigned, unsigned> hingeCount = stk::tools::impl::get_hinge_count(bulk);

  EXPECT_EQ(2u, hingeCount.first);
  EXPECT_EQ(1u, hingeCount.second);
  output_mesh(bulk);
}

TEST(ElementGroups2D, SingleBlockFourQuadOneNodeHinge)
{
  stk::mesh::MetaData meta(2);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_4quad_bowtie_1hinge(bulk);
  stk::mesh::Entity node = bulk.get_entity(stk::topology::NODE_RANK, 4u);
  stk::tools::impl::HingeGroupVector groupings = stk::tools::impl::get_convex_groupings(bulk, node);
  if(!groupings.empty()) {
    EXPECT_EQ(groupings.size(), 4u);
    EXPECT_EQ(groupings[0].size(), 1u);
  }
  output_mesh(bulk);
}

TEST(ElementGroups2D, SingleBlockThreeQuadOneNodeHinge)
{
  stk::mesh::MetaData meta(2);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_3quad_1hinge(bulk);
  stk::mesh::Entity node = bulk.get_entity(stk::topology::NODE_RANK, 4u);
  stk::tools::impl::HingeGroupVector groupings = stk::tools::impl::get_convex_groupings(bulk, node);
  if(!groupings.empty()) {
    EXPECT_EQ(groupings.size(), 2u);
    EXPECT_EQ(groupings[0].size(), 2u);
    EXPECT_EQ(groupings[1].size(), 1u);
  }
  output_mesh(bulk);
}

TEST(ElementGroups3D, EmptyMesh)
{
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  stk::mesh::Entity entity;
  stk::tools::impl::HingeGroupVector groupings = stk::tools::impl::get_convex_groupings(bulk, entity);
  if(!groupings.empty()) {
    EXPECT_EQ(groupings.size(), 0u);
  }
}

TEST(ElementGroups3D, SingleBlockOneHex)
{
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_1hex(bulk);
  stk::mesh::Entity node = bulk.get_entity(stk::topology::NODE_RANK, 1u);
  stk::mesh::Entity elem = bulk.get_entity(stk::topology::ELEMENT_RANK, 1u);
  stk::tools::impl::HingeGroupVector groupings = stk::tools::impl::get_convex_groupings(bulk, node);
  if(!groupings.empty()) {
    EXPECT_EQ(groupings.size(), 1u);
    EXPECT_EQ(groupings[0].size(), 1u);
    EXPECT_EQ(groupings[0][0], elem);
  }
}

TEST(ElementGroups3D, SingleBlockTwoHex)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 2) {
    return;
  }
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_2hex(bulk);
  two_elements_decomposition(bulk);
  stk::mesh::Entity node = bulk.get_entity(stk::topology::NODE_RANK, 1u);
  stk::tools::impl::HingeGroupVector groupings = stk::tools::impl::get_convex_groupings(bulk, node);
  if(!groupings.empty()) {
    EXPECT_EQ(groupings.size(), 1u);
    EXPECT_EQ(groupings[0].size(), 1u);
  }
}

TEST(ElementGroups3D, SingleBlockTwoHexOneNodeHinge)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 2) {
    return;
  }
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_2hex_1node_hinge(bulk);
  two_elements_decomposition(bulk);
  stk::mesh::Entity node = bulk.get_entity(stk::topology::NODE_RANK, 5u);

  test_two_element_one_hinge_grouping(bulk, stk::tools::impl::convert_to_hinge_node(bulk, node));
  output_mesh(bulk);
}

TEST(ElementGroups3D, InsertGroupTest)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 2) {
    return;
  }
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  stk::io::fill_mesh("generated:2x2x2", bulk);

  stk::mesh::Entity elem1 = bulk.get_entity(stk::topology::ELEMENT_RANK, 1u);
  stk::mesh::Entity elem2 = bulk.get_entity(stk::topology::ELEMENT_RANK, 2u);
  stk::mesh::Entity elem3 = bulk.get_entity(stk::topology::ELEMENT_RANK, 3u);
  stk::mesh::Entity elem4 = bulk.get_entity(stk::topology::ELEMENT_RANK, 4u);
  stk::tools::impl::PairwiseSideInfoVector infoVec;
  infoVec.push_back(stk::tools::impl::PairwiseSideInfo(bulk, elem1, elem2));
  infoVec.push_back(stk::tools::impl::PairwiseSideInfo(bulk, elem3, elem4));
  infoVec.push_back(stk::tools::impl::PairwiseSideInfo(bulk, elem1, elem3));

  stk::tools::impl::HingeGroupVector groupings;
  insert_into_group(infoVec, groupings);

  if(!groupings.empty()) {
    EXPECT_EQ(groupings.size(), 1u);
    EXPECT_EQ(groupings[0].size(), 4u);
    EXPECT_EQ(groupings[0][0], elem1);
    EXPECT_EQ(groupings[0][1], elem2);
    EXPECT_EQ(groupings[0][2], elem3);
    EXPECT_EQ(groupings[0][3], elem4);
  }
  output_mesh(bulk);
}

TEST(ElementGroups3D, MergeTest)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 2) {
    return;
  }
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  stk::io::fill_mesh("generated:2x2x2", bulk);

  stk::tools::impl::HingeGroupVector groupings;
  stk::mesh::Entity elem1 = bulk.get_entity(stk::topology::ELEMENT_RANK, 1u);
  stk::mesh::Entity elem2 = bulk.get_entity(stk::topology::ELEMENT_RANK, 2u);
  stk::mesh::Entity elem3 = bulk.get_entity(stk::topology::ELEMENT_RANK, 3u);
  stk::mesh::Entity elem4 = bulk.get_entity(stk::topology::ELEMENT_RANK, 4u);
  groupings.push_back( {elem1,elem2} );
  groupings.push_back( {elem2,elem3} );
  stk::mesh::EntityLess compare(bulk);
  stk::tools::impl::merge_groups(groupings, 5, 6, compare);
  if(!groupings.empty()) {
    EXPECT_EQ(groupings.size(), 2u);
  }

  stk::tools::impl::merge_groups(groupings, 0, 1, compare);
  if(!groupings.empty()) {
    EXPECT_EQ(groupings.size(), 1u);
    EXPECT_EQ(groupings[0].size(), 3u);
    EXPECT_EQ(groupings[0][0], elem1);
    EXPECT_EQ(groupings[0][1], elem2);
    EXPECT_EQ(groupings[0][2], elem3);
  }

  groupings.clear();
  groupings.push_back( {elem1,elem2} );
  groupings.push_back( {elem3,elem4} );

  stk::tools::impl::merge_groups(groupings, 0, 1, compare);
  if(!groupings.empty()) {
    EXPECT_EQ(groupings.size(), 1u);
    EXPECT_EQ(groupings[0].size(), 4u);
    EXPECT_EQ(groupings[0][0], elem1);
    EXPECT_EQ(groupings[0][1], elem2);
    EXPECT_EQ(groupings[0][2], elem3);
    EXPECT_EQ(groupings[0][3], elem4);
  }
}

TEST(ElementGroups2D, TwoBlockFiveElementsOneHinge)
{
  stk::mesh::MetaData meta(2);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_2block_3quad_2tri_1hinge(bulk);
  five_elements_decomposition(bulk);

  stk::tools::impl::HingeNodeVector hingeNodes;
  stk::tools::impl::HingeEdgeVector hingeEdges;
  fill_mesh_hinges(bulk, hingeNodes, hingeEdges);

  int ranTest = 0;

  if(hingeNodes.size() == 1) {
    ranTest = 1;
    stk::tools::impl::HingeGroupVector groupings = stk::tools::impl::get_convex_groupings(bulk, hingeNodes[0].get_node());
    EXPECT_EQ(2u, groupings.size());
    EXPECT_EQ(5u, groupings[0].size() + groupings[1].size() );
    if(groupings[0].size() == 2) {
      EXPECT_EQ(3u, groupings[1].size());
    } else {
      EXPECT_EQ(2u, groupings[1].size());
    }
  }
  verify_test_run(bulk.parallel(), ranTest);
}

TEST(ElementGroups3D, OneBlockFourElementsTwoNodeHinge)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 4) {
    return;
  }
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_four_hex_2node_hinge(bulk);
  four_elements_decomposition2(bulk);
  stk::tools::impl::HingeNodeVector hingeNodes;
  stk::tools::impl::HingeEdgeVector hingeEdges;
  fill_mesh_hinges(bulk, hingeNodes, hingeEdges);

  stk::mesh::Entity elem1 = bulk.get_entity(stk::topology::ELEM_RANK, 1u);
  stk::mesh::Entity elem2 = bulk.get_entity(stk::topology::ELEM_RANK, 2u);

  for(const stk::tools::impl::HingeNode& node : hingeNodes) {
    stk::tools::impl::HingeGroupVector groupings = stk::tools::impl::get_convex_groupings(bulk, node.get_node());
    if(!groupings.empty()) {
      EXPECT_EQ(2u, groupings.size());
      EXPECT_EQ(3u, groupings[0].size() + groupings[1].size());
      if(groupings[0].size() == 1u) {
        stk::mesh::EntityId id = bulk.identifier(groupings[0][0]);
        EXPECT_TRUE( (id == 3u) || (id == 4u) );
        EXPECT_EQ(2u, groupings[1].size());
        bool foundElem1 = std::find(groupings[1].begin(), groupings[1].end(), elem1) != groupings[1].end();
        bool foundElem2 = std::find(groupings[1].begin(), groupings[1].end(), elem2) != groupings[1].end();

        EXPECT_TRUE( (foundElem1 && foundElem2) );
      }
      if(groupings[1].size() == 1u) {
        stk::mesh::EntityId id = bulk.identifier(groupings[1][0]);
        EXPECT_TRUE( (id == 3u) || (id == 4u) );
        EXPECT_EQ(2u, groupings[0].size());
        bool foundElem1 = std::find(groupings[0].begin(), groupings[0].end(), elem1) != groupings[0].end();
        bool foundElem2 = std::find(groupings[0].begin(), groupings[0].end(), elem2) != groupings[0].end();

        EXPECT_TRUE( (foundElem1 && foundElem2) );
      }
    }
  }
}

TEST(ElementGroups3D, OneBlockEightElementsOneNodeHinge)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 4) {
    return;
  }
  stk::mesh::MetaData meta(2);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_eight_tri_1node_hinge(bulk);

  if(bulk.parallel_size() == 2) {
    stk::mesh::EntityIdProcVec idProcVec{ {1u,1}, {3u,1} };
    distribute_mesh(bulk, idProcVec);
  }
  else if(bulk.parallel_size() == 3) {
    stk::mesh::EntityIdProcVec idProcVec{ {5u,1}, {7u,2} };
    distribute_mesh(bulk, idProcVec);
  }
  else if(bulk.parallel_size() == 4) {
    stk::mesh::EntityIdProcVec idProcVec{ {2u,1}, {3u,2}, {4u,3}, {5u,3} };
    distribute_mesh(bulk, idProcVec);
  }

  stk::tools::impl::HingeNodeVector hingeNodes;
  stk::tools::impl::HingeEdgeVector hingeEdges;
  fill_mesh_hinges(bulk, hingeNodes, hingeEdges);

  int ranTest = 0;

  if(hingeNodes.size() > 0) {
    ranTest = 1;
    for(const stk::tools::impl::HingeNode& node : hingeNodes) {
      stk::tools::impl::HingeGroupVector groupings = stk::tools::impl::get_convex_groupings(bulk, node);
      if(!groupings.empty()) {
        EXPECT_EQ(4u, groupings.size());
        EXPECT_EQ(8u, groupings[0].size() + groupings[1].size() + groupings[2].size() + groupings[3].size());
        EXPECT_EQ(2u, groupings[1].size());
        EXPECT_EQ(2u, groupings[2].size());
        EXPECT_EQ(2u, groupings[3].size());
      }
    }
  }

  verify_test_run(bulk.parallel(), ranTest);
}


TEST(ElementGroups3D, SingleBlockTwoHexOneEdgeHinge)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 2)
    return;
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_2hex_1edge_hinge(bulk);

  two_elements_decomposition(bulk);

  stk::tools::impl::HingeNodeVector hingeNodes;
  stk::tools::impl::HingeEdgeVector hingeEdges;
  fill_mesh_hinges(bulk, hingeNodes, hingeEdges);

  EXPECT_EQ(0u, hingeNodes.size());
  EXPECT_EQ(1u, hingeEdges.size());

  test_two_element_one_hinge_grouping(bulk, hingeEdges[0]);
  output_mesh(bulk);
}

TEST(ElementGroups3D, SingleBlockFourElementsBowtieOneEdgeHinge)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 4)
    return;
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_4hex_bowtie_1edge_hinge(bulk);

  four_elements_decomposition2(bulk);

  stk::tools::impl::HingeNodeVector hingeNodes;
  stk::tools::impl::HingeEdgeVector hingeEdges;
  fill_mesh_hinges(bulk, hingeNodes, hingeEdges);
  int ranTest = 0;

  if(hingeEdges.size() > 0) {
    ranTest = 1;
    stk::tools::impl::HingeGroupVector groupings = stk::tools::impl::get_convex_groupings(bulk, hingeEdges[0]);
    if(!groupings.empty()) {
      EXPECT_EQ(4u, groupings.size());
      EXPECT_EQ(4u, groupings[0].size() + groupings[1].size() + groupings[2].size() + groupings[3].size());
      EXPECT_EQ(1u, groupings[0].size());
      EXPECT_EQ(1u, groupings[1].size());
      EXPECT_EQ(1u, groupings[2].size());
      EXPECT_EQ(1u, groupings[3].size());
      EXPECT_TRUE(groupings[0][0] != groupings[1][0]);
      EXPECT_TRUE(groupings[1][0] != groupings[2][0]);
      EXPECT_TRUE(groupings[2][0] != groupings[3][0]);
      EXPECT_TRUE(groupings[3][0] != groupings[0][0]);
    }
  }

  verify_test_run(bulk.parallel(), ranTest);
}

TEST(ElementGroups3D, SingleBlockTwoByTwoHexTwoEdgeHinge_Decomp2)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 4)
    return;
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_two_by_two_hex_2edge_hinge(bulk);

  four_elements_decomposition2(bulk);

  stk::tools::impl::HingeNodeVector hingeNodes;
  stk::tools::impl::HingeEdgeVector hingeEdges;
  fill_mesh_hinges(bulk, hingeNodes, hingeEdges);
  int ranTest = 0;

  if(hingeEdges.size() > 0) {
    ranTest = 1;
    stk::tools::impl::HingeGroupVector groupings = stk::tools::impl::get_convex_groupings(bulk, hingeEdges[0]);
    if(!groupings.empty()) {
      EXPECT_EQ(2u, groupings.size());
      EXPECT_EQ(2u, groupings[0].size() + groupings[1].size());
      EXPECT_EQ(1u, groupings[0].size());
      EXPECT_EQ(1u, groupings[1].size());
      EXPECT_TRUE(groupings[0] != groupings[1]);
    }
    groupings = stk::tools::impl::get_convex_groupings(bulk, hingeEdges[1]);
    if(!groupings.empty()) {
      EXPECT_EQ(2u, groupings.size());
      EXPECT_EQ(2u, groupings[0].size() + groupings[1].size());
      EXPECT_EQ(1u, groupings[0].size());
      EXPECT_EQ(1u, groupings[1].size());
      EXPECT_TRUE(groupings[0] != groupings[1]);
    }
  }

  verify_test_run(bulk.parallel(), ranTest);
  output_mesh(bulk);
}

TEST(ElementGroups3D, SingleBlockThreeElementsOneNodeHingeOneEdgeHinge)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 3)
    return;
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_1block_3hex_1node_hinge_1edge_hinge(bulk);

  three_elements_decomposition(bulk);

  stk::tools::impl::HingeNodeVector hingeNodes;
  stk::tools::impl::HingeEdgeVector hingeEdges;
  fill_mesh_hinges(bulk, hingeNodes, hingeEdges);


  int ranTest = 0;

  if(hingeNodes.size() > 0) {
    ranTest = 1;
    stk::tools::impl::HingeGroupVector groupings = stk::tools::impl::get_convex_groupings(bulk, hingeNodes[0]);
    if(!groupings.empty()) {
      EXPECT_EQ(3u, groupings.size());
      EXPECT_EQ(3u, groupings[0].size() + groupings[1].size() + groupings[2].size());
      EXPECT_EQ(1u, groupings[0].size());
      EXPECT_EQ(1u, groupings[1].size());
      EXPECT_EQ(1u, groupings[2].size());
      EXPECT_TRUE(groupings[0] != groupings[1]);
      EXPECT_TRUE(groupings[1] != groupings[2]);
    }
  }

  verify_test_run(bulk.parallel(), ranTest);

  ranTest = 0;
  if(hingeEdges.size() > 0) {
    ranTest = 1;
    stk::tools::impl::HingeGroupVector groupings = stk::tools::impl::get_convex_groupings(bulk, hingeEdges[0]);
    if(!groupings.empty()) {
      EXPECT_EQ(2u, groupings.size());
      EXPECT_EQ(2u, groupings[0].size() + groupings[1].size());
      EXPECT_EQ(1u, groupings[0].size());
      EXPECT_EQ(1u, groupings[1].size());
      EXPECT_TRUE(groupings[0] != groupings[1]);
    }
  }

  verify_test_run(bulk.parallel(), ranTest);
}

void test_snipping_result(const stk::mesh::BulkData& bulk,
                          unsigned expectedNumElems, unsigned expectedNumNodes,
                          unsigned expectedHingeNodeCount, unsigned expectedHingeEdgeCount)
{
  std::pair<unsigned,unsigned> hingeCounts;
  std::pair<unsigned,unsigned> entityCounts;

  hingeCounts = stk::tools::impl::get_hinge_count(bulk);

  EXPECT_EQ(expectedHingeNodeCount, hingeCounts.first);
  EXPECT_EQ(expectedHingeEdgeCount, hingeCounts.second);

  entityCounts = get_reduced_entity_counts(bulk);

  EXPECT_EQ(expectedNumElems, entityCounts.first);
  EXPECT_EQ(expectedNumNodes, entityCounts.second);
}

void snip_hinges(stk::mesh::BulkData& bulk)
{
  bool debug = is_debug();
  stk::tools::impl::snip_all_hinges_between_blocks(bulk, debug);
}

// Hinge snipping tests
TEST(SnipHinge2D, TwoBlockTwoElementsNoHinge)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 2)
    return;
  stk::mesh::MetaData meta(2);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_2block_2quad(bulk);
  two_elements_decomposition(bulk);

  snip_hinges(bulk);
  test_snipping_result(bulk, 2u, 6u, 0u, 0u);
  output_mesh(bulk);
}

TEST(SnipHinge2D, TwoBlockTwoElementsOneNodeHinge)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 2)
    return;
  stk::mesh::MetaData meta(2);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_2block_2quad_1node_hinge(bulk);
  two_elements_decomposition(bulk);

  snip_hinges(bulk);
  test_snipping_result(bulk, 2u, 8u, 0u, 0u);
  output_mesh(bulk);
}

TEST(SnipHinge2D, TwoBlockTwoElementsTwoNodeHinges)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 2)
    return;
  stk::mesh::MetaData meta(2);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_2block_2quad_2hinge(bulk);
  two_elements_decomposition(bulk);

  snip_hinges(bulk);
  test_snipping_result(bulk, 2u, 8u, 0u, 0u);
  output_mesh(bulk);
}

TEST(SnipHinge2D, ThreeBlockThreeElementsOneHinge_LinearStack)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 3)
    return;
  stk::mesh::MetaData meta(2);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_3block_3quad_1hinge_linear_stack(bulk);
  three_elements_decomposition(bulk);

  snip_hinges(bulk);
  test_snipping_result(bulk, 3u, 10u, 0u, 0u);
  output_mesh(bulk);
}

TEST(SnipHinge2D, FourBlockFourElementsBowtie_Decomp2)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 4)
    return;
  stk::mesh::MetaData meta(2);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_4block_4quad_bowtie_1hinge(bulk);
  four_elements_decomposition(bulk);

  snip_hinges(bulk);
  test_snipping_result(bulk, 4u, 16u, 0u, 0u);
  output_mesh(bulk);
}

TEST(SnipHinge2D, ThreeBlockThreeElementsOneHinge_HorizontalCut)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 3)
    return;
  stk::mesh::MetaData meta(2);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_3block_3quad_1hinge(bulk);
  three_elements_decomposition(bulk);

  snip_hinges(bulk);
  test_snipping_result(bulk, 3u, 10u, 0u, 0u);
  output_mesh(bulk);
}

TEST(SnipHinge2D, FourBlockFourElementsTwoHinge_Decomp2)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 4)
    return;
  stk::mesh::MetaData meta(2);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_4block_4quad_2hinge(bulk);
  four_elements_decomposition(bulk);

  snip_hinges(bulk);
  test_snipping_result(bulk, 4u, 12u, 0u, 0u);
  output_mesh(bulk);
}

TEST(SnipHinge2D, FourBlockFourElementsFourHinge_Decomp2)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 4)
    return;
  stk::mesh::MetaData meta(2);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_4block_4quad_4hinge(bulk);
  four_elements_decomposition2(bulk);

  snip_hinges(bulk);
  test_snipping_result(bulk, 4u, 16u, 0u, 0u);
  output_mesh(bulk);
}

TEST(SnipHinge2D, FourBlockFourElementsNoHingePacman)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 4)
    return;
  stk::mesh::MetaData meta(2);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_4block_4quad_pacman(bulk);
  four_elements_decomposition2(bulk);

  snip_hinges(bulk);
  test_snipping_result(bulk, 4u, 10u, 0u, 0u);
  output_mesh(bulk);
}

TEST(SnipHinge2D, FourBlockFourElementsOneHinge)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 4)
    return;
  stk::mesh::MetaData meta(2);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_4block_4quad_1hinge(bulk);
  four_elements_decomposition2(bulk);

  snip_hinges(bulk);
  test_snipping_result(bulk, 4u, 12u, 0u, 0u);
  output_mesh(bulk);
}

TEST(SnipHinge2D, FiveBlockFiveElementsOneHinge)
{
  stk::mesh::MetaData meta(2);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_5block_3quad_2tri_1hinge(bulk);
  five_elements_decomposition(bulk);

  snip_hinges(bulk);
  test_snipping_result(bulk, 5u, 12u, 0u, 0u);
  output_mesh(bulk);
}

TEST(SnipHinge3D, TwoBlockTwoHexOneEdgeHinge)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 2)
    return;
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_2block_2hex_1edge_hinge(bulk);
  two_elements_decomposition(bulk);

  snip_hinges(bulk);
  test_snipping_result(bulk, 2u, 16u, 0u, 0u);
  output_mesh(bulk);
}

TEST(SnipHinge3D, ThreeBlockThreeHexOneEdgeHinge)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 3)
    return;
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_3block_3hex_1edge_hinge(bulk);
  three_elements_decomposition(bulk);

  snip_hinges(bulk);
  test_snipping_result(bulk, 3u, 20u, 0u, 0u);
  output_mesh(bulk);
}

TEST(SnipHinge3D, ThreeBlockThreeElementsOneNodeHingeOneEdgeHinge)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 3)
    return;
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_3block_3hex_1node_hinge_1edge_hinge(bulk);
  three_elements_decomposition(bulk);

  snip_hinges(bulk);
  test_snipping_result(bulk, 3u, 24u, 0u, 0u);
  output_mesh(bulk);
}

TEST(SnipHinge3D, ThreeBlockThreeElementsOneNodeHingeOneEdgeHinge2)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 3)
    return;
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_3block_3hex_1node_hinge_1edge_hinge2(bulk);
  three_elements_decomposition(bulk);

  snip_hinges(bulk);
  test_snipping_result(bulk, 3u, 24u, 0u, 0u);
  output_mesh(bulk);
}

TEST(SnipHinge3D, ThreeBlockThreeElementsOneNodeHingeOneEdgeHinge3)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 3)
    return;
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_3block_3hex_1node_hinge_1edge_hinge3(bulk);
  three_elements_decomposition(bulk);

  snip_hinges(bulk);
  test_snipping_result(bulk, 3u, 24u, 0u, 0u);
  output_mesh(bulk);
}

TEST(SnipHinge3D, FourBlockFourElementsBowtieOneEdgeHinge)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 4)
    return;
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_4block_4hex_bowtie_1edge_hinge(bulk);
  four_elements_decomposition(bulk);

  snip_hinges(bulk);
  test_snipping_result(bulk, 4u, 32u, 0u, 0u);
  output_mesh(bulk);
}

TEST(SnipHinge3D, FourBlockTwoByTwoHexTwoEdgeHinge_Decomp2)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 4)
    return;
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_4block_two_by_two_hex_2edge_hinge(bulk);
  four_elements_decomposition2(bulk);

  snip_hinges(bulk);
  test_snipping_result(bulk, 4u, 24u, 0u, 0u);
  output_mesh(bulk);
}

TEST(SnipHinge3D, FourBlockFourHexOneEdgeOneNodeHinge_Decomp2)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 4)
    return;
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_4block_four_hex_one_edge_one_node_hinge(bulk);
  four_elements_decomposition2(bulk);

  snip_hinges(bulk);
  test_snipping_result(bulk, 4u, 28u, 0u, 0u);
  output_mesh(bulk);
}

TEST(SnipHinge3D, FourBlockFourHexTwoNodeHinges_Decomp2)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 4)
    return;
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  setup_mesh_4block_four_hex_2node_hinge(bulk);
  four_elements_decomposition2(bulk);

  snip_hinges(bulk);
  test_snipping_result(bulk, 4u, 28u, 0u, 0u);
  output_mesh(bulk);
}

TEST(SnipHinge, inputFile)
{
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  std::string inputFileName = stk::unit_test_util::get_option("--inputFile", "");
  if(!inputFileName.empty()) {
    double startTime = stk::wall_time();
    stk::io::fill_mesh(inputFileName, bulk);
    double meshReadTime = stk::wall_time();
    snip_hinges(bulk);
    double snipTime = stk::wall_time();
    output_mesh(bulk);
    double meshWriteTime = stk::wall_time();

    if (bulk.parallel_rank() == 0) {
      std::cout << " Mesh read time = " << (meshReadTime - startTime) << " s" << std::endl;
      std::cout << " Mesh snip time = " << (snipTime - meshReadTime) << " s" << std::endl;
      std::cout << "Mesh write time = " << (meshWriteTime - snipTime) << " s" << std::endl;
    }
  }
}


