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

#include "DisconnectBlocksMeshConstruction.hpp"
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>
#include "stk_unit_test_utils/getOption.h"
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/util/SortAndUnique.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/WriteMesh.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_tools/mesh_tools/DetectHingesImpl.hpp>
#include <stk_tools/mesh_tools/DisconnectBlocks.hpp>
#include <stk_tools/mesh_tools/DisconnectBlocksImpl.hpp>
#include <stk_tools/mesh_tools/DisconnectUtils.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <string>
#include <algorithm>
#include <fstream>
#include <libgen.h>

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
      STK_ThrowRequire(bulk.is_valid(elem));
      if(idProc.second != bulk.parallel_rank()) {
        procVec.push_back(std::make_pair(elem, idProc.second));
//        std::cout << "p:" << bulk.parallel_rank() << " element: " << bulk.entity_key(elem) << std::endl;
      }

      stk::mesh::Entity const * elemNodes = bulk.begin_nodes(elem);
      for(unsigned i = 0; i < bulk.num_nodes(elem); i++) {
        stk::mesh::Entity node = elemNodes[i];
        STK_ThrowRequire(bulk.is_valid(node));
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
  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

  return {&block1};
}

stk::mesh::PartVector setup_mesh_1block_2quad(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  std::string meshDesc = "0,1,QUAD_4_2D,1,2,4,3,block_1\n"
                         "0,2,QUAD_4_2D,2,5,6,4,block_1";
  std::vector<double> coordinates = { 0,0, 1,0, 0,1, 1,1, 2,0, 2,1 };
  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

  return {&block1};
}

stk::mesh::PartVector setup_mesh_1block_2quad_1node_hinge(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  std::string meshDesc = "0,1,QUAD_4_2D,1,2,4,3,block_1\n"
                         "0,2,QUAD_4_2D,7,5,6,4,block_1";
  std::vector<double> coordinates = { 0,0, (1-EPS),0, 0,1, 1,1, 2,0, 2,1, (1+EPS),0 };
  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

  return {&block1};
}

stk::mesh::PartVector setup_mesh_2block_2quad_1node_hinge(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);
  std::string meshDesc = "0,1,QUAD_4_2D,1,2,4,3,block_1\n"
                         "0,2,QUAD_4_2D,7,5,6,4,block_2";
  std::vector<double> coordinates = { 0,0, (1-EPS),0, 0,1, 1,1, 2,0, 2,1, (1+EPS),0 };
  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

  return {&block1,&block2};
}

stk::mesh::PartVector setup_mesh_1block_2quad_2hinge(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  std::string meshDesc = "0,1,QUAD_4_2D,3,4,1,2,block_1\n"
                         "0,2,QUAD_4_2D,2,6,4,5,block_1";
  std::vector<double> coordinates = { 0,2, 2,1, 1,2, 2,3, 3,2, 4,2 };
  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

  return {&block1};
}

stk::mesh::PartVector setup_mesh_2block_2quad_2hinge(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);
  std::string meshDesc = "0,1,QUAD_4_2D,3,4,1,2,block_1\n"
                         "0,2,QUAD_4_2D,2,6,4,5,block_2";
  std::vector<double> coordinates = { 0,2, 2,1, 1,2, 2,3, 3,2, 4,2 };
  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

  return {&block1,&block2};
}

stk::mesh::PartVector setup_mesh_1block_3quad_1hinge(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  std::string meshDesc = "0,1,QUAD_4_2D,1,2,4,3,block_1\n"
                         "0,2,QUAD_4_2D,2,5,6,4,block_1\n"
                         "0,3,QUAD_4_2D,4,7,8,9,block_1";
  std::vector<double> coordinates = { 0,0, 1,0, 0,1, 1,1, 2,0, 2,(1-EPS), 2,(1+EPS), 2,2, 1,2 };
  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

  return {&block1};
}

stk::mesh::PartVector setup_mesh_1block_3quad_1hinge_linear_stack(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  std::string meshDesc = "0,1,QUAD_4_2D,1,2,4,3,block_1\n"
                         "0,2,QUAD_4_2D,2,5,6,4,block_1\n"
                         "0,3,QUAD_4_2D,7,8,9,6,block_1";
  std::vector<double> coordinates = { 0,0, 1,0, 0,1, 1,1, (2-EPS),0, 2,1, (2+EPS),0, 3,0, 3,1 };
  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

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
  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

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
  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

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
  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

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
  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

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
  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

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
  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

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
  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

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
  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

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
  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

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
  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

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
  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

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
  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

  return {&block1,&block2,&block3,&block4,&block5};
}

stk::mesh::PartVector setup_mesh_1block_1hex(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_1", 1);
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1";
  std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1};

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

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

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

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

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

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

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

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

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

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

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

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

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

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

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

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

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

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

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

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

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

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

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

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

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

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

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

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

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

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

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

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

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

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

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

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

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

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

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

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

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

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

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

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

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

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

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

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

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

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

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

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

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

  return {&block1};
}

stk::mesh::ConnectivityOrdinal destroy_element_node_relation(stk::mesh::BulkData& bulk, stk::mesh::Entity element, stk::mesh::Entity node)
{
  const stk::mesh::Entity* nodes = bulk.begin_nodes(element);
  stk::mesh::ConnectivityOrdinal const *ordinals = bulk.begin_ordinals(element, stk::topology::NODE_RANK);
  unsigned numNodes = bulk.num_nodes(element);

  for(unsigned i = 0; i < numNodes; i++) {
    if(nodes[i] == node) {
      bulk.destroy_relation(element, node, ordinals[i]);
      return ordinals[i];
    }
  }

  return stk::mesh::INVALID_CONNECTIVITY_ORDINAL;
}

void setup_mesh_with_hinge_ring(stk::mesh::BulkData& bulk)
{
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stk::mesh::Part& block2 = create_part(meta, stk::topology::HEX_8, "block_2", 2);
  stk::io::fill_mesh("generated:3x3x1", bulk);

  bulk.modification_begin();

  stk::mesh::Entity node22 = bulk.get_entity(stk::topology::NODE_RANK, 22);
  stk::mesh::Entity node23 = bulk.get_entity(stk::topology::NODE_RANK, 23);
  stk::mesh::Entity node26 = bulk.get_entity(stk::topology::NODE_RANK, 26);
  stk::mesh::Entity node27 = bulk.get_entity(stk::topology::NODE_RANK, 27);

  stk::mesh::Entity newNode122 = bulk.declare_node(122);
  stk::mesh::Entity newNode123 = bulk.declare_node(123);
  stk::mesh::Entity newNode126 = bulk.declare_node(126);
  stk::mesh::Entity newNode127 = bulk.declare_node(127);

  stk::mesh::Entity elem5 = bulk.get_entity(stk::topology::ELEMENT_RANK, 5);
  stk::mesh::ConnectivityOrdinal ordinal;
  ordinal = destroy_element_node_relation(bulk, elem5, node22);
  bulk.declare_relation(elem5, newNode122, ordinal);
  ordinal = destroy_element_node_relation(bulk, elem5, node23);
  bulk.declare_relation(elem5, newNode123, ordinal);
  ordinal = destroy_element_node_relation(bulk, elem5, node26);
  bulk.declare_relation(elem5, newNode126, ordinal);
  ordinal = destroy_element_node_relation(bulk, elem5, node27);
  bulk.declare_relation(elem5, newNode127, ordinal);

  stk::mesh::Part* block1 = meta.get_part("block_1");
  bulk.change_entity_parts(elem5, stk::mesh::PartVector{&block2}, stk::mesh::PartVector{block1});

  bulk.modification_end();

  const stk::mesh::FieldBase* coords = meta.coordinate_field();
  double* coordNode22 = (double*)stk::mesh::field_data(*coords, node22);
  double* coordNode23 = (double*)stk::mesh::field_data(*coords, node23);
  double* coordNode26 = (double*)stk::mesh::field_data(*coords, node26);
  double* coordNode27 = (double*)stk::mesh::field_data(*coords, node27);

  double* coordNode122 = (double*)stk::mesh::field_data(*coords, newNode122);
  double* coordNode123 = (double*)stk::mesh::field_data(*coords, newNode123);
  double* coordNode126 = (double*)stk::mesh::field_data(*coords, newNode126);
  double* coordNode127 = (double*)stk::mesh::field_data(*coords, newNode127);

  coordNode122[0] = coordNode22[0] + EPS;
  coordNode122[1] = coordNode22[1] + EPS;
  coordNode122[2] = coordNode22[2];
  coordNode123[0] = coordNode23[0] - EPS;
  coordNode123[1] = coordNode23[1] + EPS;
  coordNode123[2] = coordNode23[2];
  coordNode126[0] = coordNode26[0] + EPS;
  coordNode126[1] = coordNode26[1] - EPS;
  coordNode126[2] = coordNode26[2];
  coordNode127[0] = coordNode27[0] - EPS;
  coordNode127[1] = coordNode27[1] - EPS;
  coordNode127[2] = coordNode27[2];
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

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

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

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

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
  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

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
  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

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
      os << "\tHinge node id: " << bulk.identifier(node.get_node()) << std::endl;
    }
  }
  if(hingeEdges.size() > 0) {
    os << "PRINTING HINGE EDGES on Proc " << bulk.parallel_rank() << " : " << std::endl;
    for(const stk::tools::impl::HingeEdge& edge : hingeEdges) {
      os << "\tHinge edge ids: " << bulk.identifier(edge.first.get_node())
         << ", " << bulk.identifier(edge.second.get_node()) << std::endl;
    }
  }

  for(int i = 0; i < bulk.parallel_size(); i++) {
    if(i == bulk.parallel_rank()) {
      std::cout << os.str() << std::endl;
    }
    MPI_Barrier(bulk.parallel());
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

stk::tools::BlockPairVector convert_connection_vector_to_pair_vector(const stk::mesh::BulkData& bulk, const BlockConnectionVector& disconnectConnVector)
{
  stk::tools::BlockPairVector pairVector;
  const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stk::tools::impl::PartPairLess comparator;

  for(const BlockConnection& connection : disconnectConnVector) {
    stk::mesh::Part* part1 = meta.get_part(connection.block1);
    stk::mesh::Part* part2 = meta.get_part(connection.block2);

    stk::util::insert_keep_sorted_and_unique(stk::tools::impl::get_block_pair(part1, part2), pairVector, comparator);
  }

  return pairVector;
}

stk::tools::BlockPairVector get_local_reconnect_list(const stk::mesh::BulkData& bulk, const BlockConnectionVector& disconnectList)
{
  stk::tools::BlockPairVector convertedDisconnectList = convert_connection_vector_to_pair_vector(bulk, disconnectList);

  return stk::tools::impl::get_local_reconnect_list(bulk, convertedDisconnectList);
}

void create_sides_between_blocks(stk::mesh::BulkData& bulk,
                                 const std::string& block1Name,
                                 const std::string& block2Name,
                                 const std::string& sidePartName)
{
  stk::mesh::Part& block1 = *bulk.mesh_meta_data().get_part(block1Name);
  stk::mesh::Part& block2 = *bulk.mesh_meta_data().get_part(block2Name);
  stk::mesh::Part& sidePart = *bulk.mesh_meta_data().get_part(sidePartName);

  stk::mesh::Selector blockSelector = block1 | block2;
  stk::mesh::create_interior_block_boundary_sides(bulk, blockSelector, stk::mesh::PartVector{&sidePart});
}

void create_all_boundary_sides(stk::mesh::BulkData& bulk, const std::string& sidePartName)
{
  stk::mesh::Part& sidePart = *bulk.mesh_meta_data().get_part(sidePartName);
  stk::mesh::PartVector allBlocksInMesh;
  stk::tools::impl::get_all_blocks_in_mesh(bulk, allBlocksInMesh);

  stk::mesh::Selector blockSelector = stk::mesh::selectUnion(allBlocksInMesh);
  stk::mesh::create_exposed_block_boundary_sides(bulk, blockSelector, stk::mesh::PartVector{&sidePart});
}

unsigned get_num_surface_nodes(const stk::mesh::BulkData& bulk, const std::vector<std::string>& blockPartNames)
{
  stk::mesh::PartVector blockParts;
  stk::mesh::EntityRank rank = bulk.mesh_meta_data().side_rank();
  for(const std::string& blockName : blockPartNames) {
    stk::mesh::Part* part = bulk.mesh_meta_data().get_part(blockName);
    STK_ThrowRequire(part != nullptr);
    STK_ThrowRequire(part->primary_entity_rank() == rank);
    blockParts.push_back(part);
  }
  stk::mesh::Selector blockSelector =  stk::mesh::selectUnion(blockParts) & bulk.mesh_meta_data().locally_owned_part();
  unsigned localCount = stk::mesh::count_selected_entities(blockSelector, bulk.buckets(stk::topology::NODE_RANK));
  return stk::get_global_sum(MPI_COMM_WORLD, localCount);
}

void create_sideset(stk::mesh::BulkData& bulk,
                    const std::string& surfacePartName,
                    const std::vector<std::string>& blockPartNames)
{
  stk::mesh::ConstPartVector blockParts;
  for(const std::string& blockName : blockPartNames) {
    stk::mesh::Part* part = bulk.mesh_meta_data().get_part(blockName);
    STK_ThrowRequire(part != nullptr);
    blockParts.push_back(part);
  }
  stk::mesh::Part& surfacePart = *bulk.mesh_meta_data().get_part(surfacePartName);

  bulk.mesh_meta_data().set_surface_to_block_mapping(&surfacePart, blockParts);
  bulk.create_sideset(surfacePart);
}

void move_elems_from_block_to_block(stk::mesh::BulkData& bulk,
                                    const std::vector<stk::mesh::EntityId>& elemIDs,
                                    const std::string& fromBlockName,
                                    const std::string& toBlockName)
{
  stk::mesh::Part& fromBlock = *bulk.mesh_meta_data().get_part(fromBlockName);
  stk::mesh::Part& toBlock = *bulk.mesh_meta_data().get_part(toBlockName);

  stk::mesh::EntityVector elems;
  for(stk::mesh::EntityId elemID : elemIDs) {
    stk::mesh::Entity elem = bulk.get_entity(stk::topology::ELEM_RANK, elemID);
    if (bulk.is_valid(elem) && bulk.bucket(elem).owned()) {
      elems.push_back(elem);
    }
  }

  bulk.batch_change_entity_parts(elems, stk::mesh::PartVector{&toBlock}, stk::mesh::PartVector{&fromBlock});
}

unsigned get_num_common_entities(const stk::mesh::BulkData & bulk, const stk::mesh::PartVector & blocks, const stk::mesh::EntityRank rank)
{
  stk::mesh::Selector selector;
  if(!blocks.empty()) {
    for (size_t i = 0; i < blocks.size()-1; ++i) {
      for (size_t j = i+1; j < blocks.size(); ++j) {
        selector |= (*blocks[i] & *blocks[j]);
      }
    }
  }

  selector &= bulk.mesh_meta_data().locally_owned_part();

  unsigned localNumCommonEntities = stk::mesh::count_selected_entities(selector, bulk.buckets(rank));

  return stk::get_global_sum(MPI_COMM_WORLD, localNumCommonEntities);
}

unsigned get_num_intersecting_nodes(const stk::mesh::BulkData & bulk, const stk::mesh::PartVector & blocks)
{
  return get_num_common_entities(bulk, blocks, stk::topology::NODE_RANK);
}

unsigned get_num_common_sides(const stk::mesh::BulkData & bulk, const stk::mesh::PartVector & blocks)
{
  return get_num_common_entities(bulk, blocks, bulk.mesh_meta_data().side_rank());
}

unsigned get_num_total_entities(const stk::mesh::BulkData & bulk, const stk::mesh::EntityRank rank)
{
  unsigned localNumTotalEntities = stk::mesh::count_selected_entities(bulk.mesh_meta_data().locally_owned_part(), bulk.buckets(rank));

  return stk::get_global_sum(MPI_COMM_WORLD, localNumTotalEntities);
}

unsigned get_num_total_nodes(const stk::mesh::BulkData & bulk)
{
  return get_num_total_entities(bulk, stk::topology::NODE_RANK);
}

unsigned get_num_total_sides(const stk::mesh::BulkData & bulk)
{
  return get_num_total_entities(bulk, bulk.mesh_meta_data().side_rank());
}

stk::mesh::EntityVector get_element_side_nodes(stk::mesh::BulkData & bulk,
                                               stk::mesh::Entity elem,
                                               stk::mesh::ConnectivityOrdinal sideOrdinal)
{
  auto elemTopology = bulk.bucket(elem).topology();
  auto sideTopology = elemTopology.side_topology(sideOrdinal);

  std::vector<stk::mesh::ConnectivityOrdinal> elementSideNodeOrdinalVector(sideTopology.num_nodes());
  elemTopology.side_node_ordinals(sideOrdinal, elementSideNodeOrdinalVector.data());

  stk::mesh::EntityVector elementSideNodeVector;
  auto elemNodes = bulk.begin_nodes(elem);
  for(auto nodeIndex : elementSideNodeOrdinalVector) {
    elementSideNodeVector.push_back(elemNodes[nodeIndex]);
  }

  return elementSideNodeVector;
}

bool verify_attached_faces(stk::mesh::BulkData & bulk)
{
  for (stk::mesh::Bucket * bucket : bulk.buckets(bulk.mesh_meta_data().side_rank())) {
    for (stk::mesh::Entity face : *bucket) {
      auto numElems = bulk.num_elements(face);
      auto elems = bulk.begin_elements(face);
      auto ordinals = bulk.begin_ordinals(face, stk::topology::ELEM_RANK);
      stk::mesh::EntityVector faceNodes(bulk.begin_nodes(face), bulk.begin_nodes(face)+bulk.num_nodes(face));
      std::sort(faceNodes.begin(), faceNodes.end());

      for (unsigned i = 0; i < numElems; ++i) {
        auto elem = elems[i];
        auto sideOrdinal = ordinals[i];

        stk::mesh::EntityVector elementSideNodeVector = get_element_side_nodes(bulk, elem, sideOrdinal);
        std::sort(elementSideNodeVector.begin(), elementSideNodeVector.end());

        if(faceNodes != elementSideNodeVector) {
          std::ostringstream oss;

          oss << "P" << bulk.parallel_rank()
              << ": Could not match nodes on face: " << bulk.entity_key(face)
              << " with element: " << bulk.entity_rank(elem)
              << " on ordinal: " << sideOrdinal << std::endl;

          oss << "\tFace nodes\n";
          for(auto node : faceNodes) {
            oss << "\t\t" << bulk.entity_key(node) << std::endl;
          }

          oss << "\n\tElement side nodes\n";
          for(auto node : elementSideNodeVector) {
            oss << "\t\t" << bulk.entity_key(node) << std::endl;
          }

          std::cout << oss.str();
          return false;
        }
      }
    }
  }

  return true;
}

bool check_orphaned_nodes(stk::mesh::BulkData & bulk)
{
  bool foundOrphanedNode = false;
  for (stk::mesh::Bucket * bucket : bulk.buckets(stk::topology::NODE_RANK)) {
    for (stk::mesh::Entity node : *bucket) {
      const unsigned numElems = bulk.num_elements(node);
      if (numElems == 0u) {
        foundOrphanedNode = true;
        std::cout << "[p" << bulk.parallel_rank() << "] Found orphaned node: " << bulk.entity_key(node) << std::endl;
      }
    }
  }
  return foundOrphanedNode;
}

void output_mesh(stk::mesh::BulkData & bulk, const std::string & fileName)
{
  std::string writeOutput = stk::unit_test_util::get_option("--output", "off");
  if (writeOutput == "on") {
    stk::io::write_mesh(fileName, bulk);
  }
}

void output_mesh(stk::mesh::BulkData & bulk)
{
  const std::string fileName = std::string(::testing::UnitTest::GetInstance()->current_test_info()->name()) + ".g";
  output_mesh(bulk, fileName);
}

int get_debug_level()
{
  int level = stk::unit_test_util::get_command_line_option("--debug", 0);
  return std::max(level, 0);
}

stk::mesh::PartVector setup_mesh_2block_1quad(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);

  std::string meshDesc = "0,1,QUAD_4_2D,1,2,4,3,block_1";
  std::vector<double> coordinates = { 0,0, 1,0, 0,1, 1,1 };
  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

  EXPECT_EQ(0u, get_num_intersecting_nodes(bulk, {&block1, &block2}));
  EXPECT_EQ(4u, get_num_total_nodes(bulk));

  return {&block1, &block2};
}

stk::mesh::PartVector setup_mesh_2block_2quad_only_on_proc_0(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);
  std::string meshDesc = "0,1,QUAD_4_2D,1,2,4,3,block_1\n"
                         "0,2,QUAD_4_2D,2,5,6,4,block_2";
  std::vector<double> coordinates = { 0,0, 1,0, 0,1, 1,1, 2,0, 2,1 };
  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

  return {&block1,&block2};
}

stk::mesh::PartVector setup_mesh_2block_2quad(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);

  std::string meshDesc;
  if (bulk.parallel_size() == 1) {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
               "0,2,QUAD_4_2D,2,3,6,5,block_2";
  }
  else {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
               "1,2,QUAD_4_2D,2,3,6,5,block_2";
  }
  std::vector<double> coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1 };
  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

  EXPECT_EQ(2u, get_num_intersecting_nodes(bulk, {&block1, &block2}));
  EXPECT_EQ(6u, get_num_total_nodes(bulk));

  return {&block1, &block2};
}

stk::mesh::PartVector setup_mesh_2block_2quad_reversed(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);

  std::string meshDesc;
  if (bulk.parallel_size() == 1) {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_2\n"
               "0,2,QUAD_4_2D,2,3,6,5,block_1";
  }
  else {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_2\n"
               "1,2,QUAD_4_2D,2,3,6,5,block_1";
  }
  std::vector<double> coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1 };
  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

  EXPECT_EQ(2u, get_num_intersecting_nodes(bulk, {&block1, &block2}));
  EXPECT_EQ(6u, get_num_total_nodes(bulk));

  return {&block1, &block2};
}

stk::mesh::PartVector setup_mesh_2block_4quad_corner(stk::mesh::BulkData& bulk, int decompPattern)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);

  std::string meshDesc;
  if (bulk.parallel_size() == 1) {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
               "0,2,QUAD_4_2D,2,3,6,5,block_2\n"
               "0,3,QUAD_4_2D,4,5,8,7,block_1\n"
               "0,4,QUAD_4_2D,5,6,9,8,block_1";
  }
  else {
    if (decompPattern == 1) {
      // p0 for block_1 and p1 for block_2
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                 "1,2,QUAD_4_2D,2,3,6,5,block_2\n"
                 "0,3,QUAD_4_2D,4,5,8,7,block_1\n"
                 "0,4,QUAD_4_2D,5,6,9,8,block_1";
    }
    else if (decompPattern == 2) {
      // p0 for bottom half and p1 for top half
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                 "0,2,QUAD_4_2D,2,3,6,5,block_2\n"
                 "1,3,QUAD_4_2D,4,5,8,7,block_1\n"
                 "1,4,QUAD_4_2D,5,6,9,8,block_1";
    }
    else if (decompPattern == 3) {
      // p0 for non-face-adjacent block_1 element, p1 for block_2 and face-adjacent block_1 elements
      meshDesc = "1,1,QUAD_4_2D,1,2,5,4,block_1\n"
                 "1,2,QUAD_4_2D,2,3,6,5,block_2\n"
                 "0,3,QUAD_4_2D,4,5,8,7,block_1\n"
                 "1,4,QUAD_4_2D,5,6,9,8,block_1";
    }
    else if (decompPattern == 4) {
      // p0 diagonal, p1 off-diagonal (checkerboard)
      meshDesc = "1,1,QUAD_4_2D,1,2,5,4,block_1\n"
                 "0,2,QUAD_4_2D,2,3,6,5,block_2\n"
                 "0,3,QUAD_4_2D,4,5,8,7,block_1\n"
                 "1,4,QUAD_4_2D,5,6,9,8,block_1";
    }
    else {
      std::cerr << "ERROR: Unexpected decomposition pattern (" << decompPattern << ")!" << std::endl;
      exit(1);
    }
  }
  std::vector<double> coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

  EXPECT_EQ(3u, get_num_intersecting_nodes(bulk, {&block1, &block2}));
  EXPECT_EQ(9u, get_num_total_nodes(bulk));

  return {&block1, &block2};
}

stk::mesh::PartVector setup_mesh_2block_4quad_swappedCorner(stk::mesh::BulkData& bulk, int decompPattern)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);

  std::string meshDesc;
  if (bulk.parallel_size() == 1) {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_2\n"
               "0,2,QUAD_4_2D,2,3,6,5,block_1\n"
               "0,3,QUAD_4_2D,4,5,8,7,block_2\n"
               "0,4,QUAD_4_2D,5,6,9,8,block_2";
  }
  else if (bulk.parallel_size() == 2) {
    if (decompPattern == 1) {
      // p0 for block_1 and p1 for block_2
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_2\n"
                 "1,2,QUAD_4_2D,2,3,6,5,block_1\n"
                 "0,3,QUAD_4_2D,4,5,8,7,block_2\n"
                 "0,4,QUAD_4_2D,5,6,9,8,block_2";
    }
    else if (decompPattern == 2) {
      // p0 for bottom half and p1 for top half
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_2\n"
                 "0,2,QUAD_4_2D,2,3,6,5,block_1\n"
                 "1,3,QUAD_4_2D,4,5,8,7,block_2\n"
                 "1,4,QUAD_4_2D,5,6,9,8,block_2";
    }
    else if (decompPattern == 3) {
      // p0 for non-face-adjacent block_2 element, p1 for block_1 and face-adjacent block_2 elements
      meshDesc = "1,1,QUAD_4_2D,1,2,5,4,block_2\n"
                 "1,2,QUAD_4_2D,2,3,6,5,block_1\n"
                 "0,3,QUAD_4_2D,4,5,8,7,block_2\n"
                 "1,4,QUAD_4_2D,5,6,9,8,block_2";
    }
    else if (decompPattern == 4) {
      // p0 diagonal, p1 off-diagonal (checkerboard)
      meshDesc = "1,1,QUAD_4_2D,1,2,5,4,block_2\n"
                 "0,2,QUAD_4_2D,2,3,6,5,block_1\n"
                 "0,3,QUAD_4_2D,4,5,8,7,block_2\n"
                 "1,4,QUAD_4_2D,5,6,9,8,block_2";
    }
    else {
      std::cerr << "ERROR: Unexpected decomposition pattern (" << decompPattern << ")!" << std::endl;
      exit(1);
    }
  }
  else {
    if (decompPattern == 1) {
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_2\n"
                 "0,2,QUAD_4_2D,2,3,6,5,block_1\n"
                 "1,3,QUAD_4_2D,4,5,8,7,block_2\n"
                 "2,4,QUAD_4_2D,5,6,9,8,block_2";
    }
    else if (decompPattern == 2) {
      meshDesc = "1,1,QUAD_4_2D,1,2,5,4,block_2\n"
                 "0,2,QUAD_4_2D,2,3,6,5,block_1\n"
                 "1,3,QUAD_4_2D,4,5,8,7,block_2\n"
                 "2,4,QUAD_4_2D,5,6,9,8,block_2";
    }
    else if (decompPattern == 3) {
      meshDesc = "1,1,QUAD_4_2D,1,2,5,4,block_2\n"
                 "0,2,QUAD_4_2D,2,3,6,5,block_1\n"
                 "0,3,QUAD_4_2D,4,5,8,7,block_2\n"
                 "2,4,QUAD_4_2D,5,6,9,8,block_2";
    }
    else if (decompPattern == 4) {
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_2\n"
                 "1,2,QUAD_4_2D,2,3,6,5,block_1\n"
                 "2,3,QUAD_4_2D,4,5,8,7,block_2\n"
                 "0,4,QUAD_4_2D,5,6,9,8,block_2";
    }
  }
  std::vector<double> coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

  EXPECT_EQ(3u, get_num_intersecting_nodes(bulk, {&block1, &block2}));
  EXPECT_EQ(9u, get_num_total_nodes(bulk));

  return {&block1, &block2};
}

stk::mesh::PartVector create_3_blocks_order1(stk::mesh::BulkData& bulk)
{
  stk::mesh::PartVector blocks(3);
  blocks[0] = &create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  blocks[1] = &create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);
  blocks[2] = &create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_3", 3);
  return blocks;
}

stk::mesh::PartVector create_3_blocks_order2(stk::mesh::BulkData& bulk)
{
  stk::mesh::PartVector blocks(3);
  blocks[0] = &create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);
  blocks[1] = &create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_3", 3);
  blocks[2] = &create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  return blocks;
}

stk::mesh::PartVector create_3_blocks_order3(stk::mesh::BulkData& bulk)
{
  stk::mesh::PartVector blocks(3);
  blocks[0] = &create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_3", 3);
  blocks[1] = &create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  blocks[2] = &create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);
  return blocks;
}

stk::mesh::PartVector create_3_blocks_order4(stk::mesh::BulkData& bulk)
{
  stk::mesh::PartVector blocks(3);
  blocks[0] = &create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_3", 3);
  blocks[1] = &create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);
  blocks[2] = &create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  return blocks;
}

stk::mesh::PartVector create_3_blocks_order5(stk::mesh::BulkData& bulk)
{
  stk::mesh::PartVector blocks(3);
  blocks[0] = &create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);
  blocks[1] = &create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  blocks[2] = &create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_3", 3);
  return blocks;
}

stk::mesh::PartVector create_3_blocks_order6(stk::mesh::BulkData& bulk)
{
  stk::mesh::PartVector blocks(3);
  blocks[0] = &create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  blocks[1] = &create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_3", 3);
  blocks[2] = &create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);
  return blocks;
}

void setup_mesh_3block_4quad_base(stk::mesh::BulkData& bulk, stk::mesh::PartVector & blocks, unsigned decompPattern)
{
  std::string meshDesc;
  if (bulk.parallel_size() == 1) {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
               "0,2,QUAD_4_2D,2,3,6,5,block_2\n"
               "0,3,QUAD_4_2D,4,5,8,7,block_1\n"
               "0,4,QUAD_4_2D,5,6,9,8,block_3";
  }
  else if (bulk.parallel_size() == 2) {
    if (decompPattern == 1) {
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                 "1,2,QUAD_4_2D,2,3,6,5,block_2\n"
                 "0,3,QUAD_4_2D,4,5,8,7,block_1\n"
                 "1,4,QUAD_4_2D,5,6,9,8,block_3";
    }
    else if (decompPattern == 2) {
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                 "0,2,QUAD_4_2D,2,3,6,5,block_2\n"
                 "1,3,QUAD_4_2D,4,5,8,7,block_1\n"
                 "1,4,QUAD_4_2D,5,6,9,8,block_3";
    }
    else if (decompPattern == 3) {
      meshDesc = "1,1,QUAD_4_2D,1,2,5,4,block_1\n"
                 "0,2,QUAD_4_2D,2,3,6,5,block_2\n"
                 "0,3,QUAD_4_2D,4,5,8,7,block_1\n"
                 "0,4,QUAD_4_2D,5,6,9,8,block_3";
    }
    else {
      std::cerr << "ERROR: Unexpected decomposition pattern (" << decompPattern << ")!" << std::endl;
      exit(1);
    }
  }
  else {
    if (decompPattern == 1) {
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                 "2,2,QUAD_4_2D,2,3,6,5,block_2\n"
                 "0,3,QUAD_4_2D,4,5,8,7,block_1\n"
                 "1,4,QUAD_4_2D,5,6,9,8,block_3";
    }
    else if (decompPattern == 2) {
      meshDesc = "2,1,QUAD_4_2D,1,2,5,4,block_1\n"
                 "0,2,QUAD_4_2D,2,3,6,5,block_2\n"
                 "1,3,QUAD_4_2D,4,5,8,7,block_1\n"
                 "0,4,QUAD_4_2D,5,6,9,8,block_3";
//
//                4
//  1*------------*------------*5
//   |   E1,B1    |   E3,B1    |
//   |    P2      |    P1      |
//  2*-----------5*------------*8
//   |   E2,B2    |   E4,B3    |
//   |    P0      |    P0      |
//  3*------------*------------*9
//                6
//
    }
    else if (decompPattern == 3) {
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                 "0,2,QUAD_4_2D,2,3,6,5,block_2\n"
                 "1,3,QUAD_4_2D,4,5,8,7,block_1\n"
                 "2,4,QUAD_4_2D,5,6,9,8,block_3";
    }
    else {
      std::cerr << "ERROR: Unexpected decomposition pattern (" << decompPattern << ")!" << std::endl;
      exit(1);
    }
  }
  std::vector<double> coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

  EXPECT_EQ(4u, get_num_intersecting_nodes(bulk, blocks));
  EXPECT_EQ(9u, get_num_total_nodes(bulk));
}

stk::mesh::PartVector setup_mesh_3block_4quad(stk::mesh::BulkData& bulk, unsigned blockOrder, unsigned decompPattern)
{
  stk::mesh::PartVector blocks;
  if (blockOrder == 1) {
    blocks = create_3_blocks_order1(bulk);
  } else if (blockOrder == 2) {
    blocks = create_3_blocks_order2(bulk);
  } else if (blockOrder == 3) {
    blocks = create_3_blocks_order3(bulk);
  } else if (blockOrder == 4) {
    blocks = create_3_blocks_order4(bulk);
  } else if (blockOrder == 5) {
    blocks = create_3_blocks_order5(bulk);
  } else if (blockOrder == 6) {
    blocks = create_3_blocks_order6(bulk);
  } else {
    std::cerr << "ERROR: Unexpected part ordinal ordering!!!" << std::endl;
    exit(1);
  }

  setup_mesh_3block_4quad_base(bulk, blocks, decompPattern);
  return blocks;
}

void test_mesh_3block_4quad(stk::mesh::BulkData& bulk, unsigned blockOrder, unsigned decompPattern) {
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(bulk, blockOrder, decompPattern);

  output_mesh(bulk, "disconnect_3block_4quad_blockOrder" + std::to_string(blockOrder) + "_decomp" + std::to_string(decompPattern) + "_init.g");

  stk::tools::disconnect_all_blocks(bulk);

  EXPECT_EQ(0u,  get_num_intersecting_nodes(bulk, blocks));
  EXPECT_EQ(14u, get_num_total_nodes(bulk));
  EXPECT_FALSE(check_orphaned_nodes(bulk));

  output_mesh(bulk, "disconnect_3block_4quad_blockOrder" + std::to_string(blockOrder) + "_decomp" + std::to_string(decompPattern) + ".g");
}

stk::mesh::PartVector setup_mesh_3block_4quad_reverse_ordinal(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & vl = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "vl", 1);
  stk::mesh::Part & radax = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "radax", 2);
  stk::mesh::Part & lateral = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "lateral", 3);

  std::string meshDesc = "0,1,QUAD_4_2D,1,2,5,4,lateral\n"
                         "0,2,QUAD_4_2D,2,3,6,5,vl\n"
                         "0,3,QUAD_4_2D,4,5,8,7,lateral\n"
                         "0,4,QUAD_4_2D,5,6,9,8,radax";

  std::vector<double> coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

  EXPECT_EQ(4u, get_num_intersecting_nodes(bulk, {&vl, &radax, &lateral}));
  EXPECT_EQ(9u, get_num_total_nodes(bulk));

  return {&vl, &radax, &lateral};
}

stk::mesh::PartVector setup_mesh_3block_4quad_keepLowerRight(stk::mesh::BulkData& bulk, unsigned decompPattern)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);
  stk::mesh::Part & block3 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_3", 3);

  std::string meshDesc;
  if (bulk.parallel_size() == 1) {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_3\n"
               "0,2,QUAD_4_2D,2,3,6,5,block_1\n"
               "0,3,QUAD_4_2D,4,5,8,7,block_3\n"
               "0,4,QUAD_4_2D,5,6,9,8,block_2";
  }
  else if (bulk.parallel_size() == 2) {
    if (decompPattern == 1) {
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_3\n"
                 "1,2,QUAD_4_2D,2,3,6,5,block_1\n"
                 "0,3,QUAD_4_2D,4,5,8,7,block_3\n"
                 "1,4,QUAD_4_2D,5,6,9,8,block_2";
    }
    else if (decompPattern == 2) {
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_3\n"
                 "0,2,QUAD_4_2D,2,3,6,5,block_1\n"
                 "1,3,QUAD_4_2D,4,5,8,7,block_3\n"
                 "1,4,QUAD_4_2D,5,6,9,8,block_2";
    }
    else if (decompPattern == 3) {
      meshDesc = "1,1,QUAD_4_2D,1,2,5,4,block_3\n"
                 "0,2,QUAD_4_2D,2,3,6,5,block_1\n"
                 "0,3,QUAD_4_2D,4,5,8,7,block_3\n"
                 "0,4,QUAD_4_2D,5,6,9,8,block_2";
    }
    else {
      std::cerr << "ERROR: Unexpected decomposition pattern (" << decompPattern << ")!" << std::endl;
      exit(1);
    }
  }
  else {
    if (decompPattern == 1) {
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_3\n"
                 "2,2,QUAD_4_2D,2,3,6,5,block_1\n"
                 "0,3,QUAD_4_2D,4,5,8,7,block_3\n"
                 "1,4,QUAD_4_2D,5,6,9,8,block_2";
    }
    else if (decompPattern == 2) {
      meshDesc = "2,1,QUAD_4_2D,1,2,5,4,block_3\n"
                 "0,2,QUAD_4_2D,2,3,6,5,block_1\n"
                 "1,3,QUAD_4_2D,4,5,8,7,block_3\n"
                 "0,4,QUAD_4_2D,5,6,9,8,block_2";
    }
    else if (decompPattern == 3) {
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_3\n"
                 "0,2,QUAD_4_2D,2,3,6,5,block_1\n"
                 "1,3,QUAD_4_2D,4,5,8,7,block_3\n"
                 "2,4,QUAD_4_2D,5,6,9,8,block_2";
    }
    else {
      std::cerr << "ERROR: Unexpected decomposition pattern (" << decompPattern << ")!" << std::endl;
      exit(1);
    }
  }
  std::vector<double> coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

  EXPECT_EQ(4u, get_num_intersecting_nodes(bulk, {&block1, &block2, &block3}));
  EXPECT_EQ(9u, get_num_total_nodes(bulk));

  return {&block1, &block2, &block3};
}

stk::mesh::PartVector setup_mesh_2block_4quad_checkerboard(stk::mesh::BulkData& bulk, unsigned decompPattern)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);

  std::string meshDesc;
  if (bulk.parallel_size() == 1) {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
               "0,2,QUAD_4_2D,2,3,6,5,block_2\n"
               "0,3,QUAD_4_2D,4,5,8,7,block_2\n"
               "0,4,QUAD_4_2D,5,6,9,8,block_1";
  }
  else {
    if (decompPattern == 1) {
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                 "0,2,QUAD_4_2D,2,3,6,5,block_2\n"
                 "1,3,QUAD_4_2D,4,5,8,7,block_2\n"
                 "1,4,QUAD_4_2D,5,6,9,8,block_1";
    }
    else if (decompPattern == 2) {
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                 "1,2,QUAD_4_2D,2,3,6,5,block_2\n"
                 "1,3,QUAD_4_2D,4,5,8,7,block_2\n"
                 "0,4,QUAD_4_2D,5,6,9,8,block_1";
    }
    else if (decompPattern == 3) {
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                 "1,2,QUAD_4_2D,2,3,6,5,block_2\n"
                 "1,3,QUAD_4_2D,4,5,8,7,block_2\n"
                 "1,4,QUAD_4_2D,5,6,9,8,block_1";
    }
    else {
      std::cerr << "ERROR: Unexpected decomposition pattern (" << decompPattern << ")!" << std::endl;
      exit(1);
    }
  }
  std::vector<double> coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

  EXPECT_EQ(5u, get_num_intersecting_nodes(bulk, {&block1, &block2}));
  EXPECT_EQ(9u, get_num_total_nodes(bulk));

  return {&block1, &block2};
}

stk::mesh::PartVector setup_mesh_3block_4quad_checkerboard(stk::mesh::BulkData& bulk, unsigned decompPattern)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);
  stk::mesh::Part & block3 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_3", 3);

  std::string meshDesc;
  if (bulk.parallel_size() == 1) {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
               "0,2,QUAD_4_2D,2,3,6,5,block_3\n"
               "0,3,QUAD_4_2D,4,5,8,7,block_2\n"
               "0,4,QUAD_4_2D,5,6,9,8,block_1";
  }
  else {
    if (decompPattern == 1) {
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                 "0,2,QUAD_4_2D,2,3,6,5,block_3\n"
                 "1,3,QUAD_4_2D,4,5,8,7,block_2\n"
                 "1,4,QUAD_4_2D,5,6,9,8,block_1";
    }
    else if (decompPattern == 2) {
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                 "1,2,QUAD_4_2D,2,3,6,5,block_3\n"
                 "1,3,QUAD_4_2D,4,5,8,7,block_2\n"
                 "0,4,QUAD_4_2D,5,6,9,8,block_1";
    }
    else if (decompPattern == 3) {
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                 "1,2,QUAD_4_2D,2,3,6,5,block_3\n"
                 "1,3,QUAD_4_2D,4,5,8,7,block_2\n"
                 "1,4,QUAD_4_2D,5,6,9,8,block_1";
    }
    else {
      std::cerr << "ERROR: Unexpected decomposition pattern (" << decompPattern << ")!" << std::endl;
      exit(1);
    }
  }
  std::vector<double> coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

  EXPECT_EQ(5u, get_num_intersecting_nodes(bulk, {&block1, &block2, &block3}));
  EXPECT_EQ(9u, get_num_total_nodes(bulk));

  return {&block1, &block2, &block3};
}

stk::mesh::PartVector setup_mesh_2block_2quad_diagonal(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);

  std::string meshDesc;
  if (bulk.parallel_size() == 1) {
    meshDesc = "0,1,QUAD_4_2D,1,2,4,3,block_1\n"
               "0,2,QUAD_4_2D,4,5,7,6,block_2";
  }
  else {
    meshDesc = "0,1,QUAD_4_2D,1,2,4,3,block_1\n"
               "1,2,QUAD_4_2D,4,5,7,6,block_2";
  }
  std::vector<double> coordinates = { 0,0, 1,0, 0,1, 1,1, 2,1, 1,2, 2,2 };

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

  EXPECT_EQ(1u, get_num_intersecting_nodes(bulk, {&block1, &block2}));
  EXPECT_EQ(7u, get_num_total_nodes(bulk));

  return {&block1, &block2};
}

stk::mesh::PartVector setup_mesh_3block_4quad_bowtie(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);
  stk::mesh::Part & block3 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_3", 3);

  std::string meshDesc = "0,1,QUAD_4_2D,1,2, 6, 5,block_1\n"
                         "0,2,QUAD_4_2D,3,4, 7, 6,block_2\n"
                         "0,3,QUAD_4_2D,8,6,11,10,block_3\n"
                         "0,4,QUAD_4_2D,6,9,13,12,block_1";
  std::vector<double> coordinates = { 0,0, 0.9,0, 1.1,0, 2,0, 0,0.9, 1,1, 2,0.9, 0,1.1, 2,1.1, 0,2, 0.9,2, 1.1,2, 2,2 };
  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

  EXPECT_EQ( 1u, get_num_intersecting_nodes(bulk, {&block1, &block2, &block3}));
  EXPECT_EQ(13u, get_num_total_nodes(bulk));

  return {&block1, &block2, &block3};
}

void fill_mesh_description_4block_4quad_np1(stk::mesh::BulkData& bulk, unsigned blockOrder,
                                            std::string& meshDesc, std::vector<double>& coordinates)
{
  STK_ThrowRequire(bulk.parallel_size() == 1);

  if (blockOrder == 1) {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_4\n"
               "0,2,QUAD_4_2D,2,3,6,5,block_3\n"
               "0,3,QUAD_4_2D,4,5,8,7,block_2\n"
               "0,4,QUAD_4_2D,5,6,9,8,block_1";
  }
  else if (blockOrder == 2) {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_4\n"
               "0,2,QUAD_4_2D,2,3,6,5,block_3\n"
               "0,3,QUAD_4_2D,4,5,8,7,block_1\n"
               "0,4,QUAD_4_2D,5,6,9,8,block_2";
  }
  else if (blockOrder == 3) {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
               "0,2,QUAD_4_2D,2,3,6,5,block_2\n"
               "0,3,QUAD_4_2D,4,5,8,7,block_3\n"
               "0,4,QUAD_4_2D,5,6,9,8,block_4";
  }
  coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };
}

void fill_mesh_description_4block_4quad_np2(stk::mesh::BulkData& bulk, unsigned blockOrder,
                                            std::string& meshDesc, std::vector<double>& coordinates)
{
  STK_ThrowRequire(bulk.parallel_size() == 2);
  STK_ThrowRequire(blockOrder <= 3 && blockOrder > 0);

  if (blockOrder == 1) {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_4\n"
               "0,2,QUAD_4_2D,2,3,6,5,block_3\n"
               "1,3,QUAD_4_2D,4,5,8,7,block_2\n"
               "1,4,QUAD_4_2D,5,6,9,8,block_1";
  }
  else if (blockOrder == 2) {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_4\n"
               "0,2,QUAD_4_2D,2,3,6,5,block_3\n"
               "1,3,QUAD_4_2D,4,5,8,7,block_1\n"
               "1,4,QUAD_4_2D,5,6,9,8,block_2";
  }
  else if (blockOrder == 3) {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
               "0,2,QUAD_4_2D,2,3,6,5,block_2\n"
               "1,3,QUAD_4_2D,4,5,8,7,block_3\n"
               "1,4,QUAD_4_2D,5,6,9,8,block_4";
  }
  coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };
}

void fill_mesh_description_4block_4quad_np3(stk::mesh::BulkData& bulk, unsigned blockOrder,
                                            std::string& meshDesc, std::vector<double>& coordinates)
{
  STK_ThrowRequire(bulk.parallel_size() == 3);
  STK_ThrowRequire(blockOrder <= 3 && blockOrder > 0);

  if (blockOrder == 1) {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_4\n"
               "1,2,QUAD_4_2D,2,3,6,5,block_3\n"
               "1,3,QUAD_4_2D,4,5,8,7,block_2\n"
               "2,4,QUAD_4_2D,5,6,9,8,block_1";
  }
  else if (blockOrder == 2) {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_4\n"
               "1,2,QUAD_4_2D,2,3,6,5,block_3\n"
               "1,3,QUAD_4_2D,4,5,8,7,block_1\n"
               "2,4,QUAD_4_2D,5,6,9,8,block_2";
  }
  else if (blockOrder == 3) {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
               "1,2,QUAD_4_2D,2,3,6,5,block_2\n"
               "1,3,QUAD_4_2D,4,5,8,7,block_3\n"
               "2,4,QUAD_4_2D,5,6,9,8,block_4";
  }
  coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };
}

void fill_mesh_description_4block_4quad_np4(stk::mesh::BulkData& bulk, unsigned blockOrder,
                                            std::string& meshDesc, std::vector<double>& coordinates)
{
  STK_ThrowRequire(bulk.parallel_size() == 4);
  STK_ThrowRequire(blockOrder <= 4 && blockOrder > 0);

  if (blockOrder == 1) {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
               "1,2,QUAD_4_2D,2,3,6,5,block_2\n"
               "2,3,QUAD_4_2D,4,5,8,7,block_3\n"
               "3,4,QUAD_4_2D,5,6,9,8,block_4";
  }
  else if (blockOrder == 2) {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_4\n"
               "1,2,QUAD_4_2D,2,3,6,5,block_1\n"
               "2,3,QUAD_4_2D,4,5,8,7,block_2\n"
               "3,4,QUAD_4_2D,5,6,9,8,block_3";
  }
  else if (blockOrder == 3) {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_3\n"
               "1,2,QUAD_4_2D,2,3,6,5,block_4\n"
               "2,3,QUAD_4_2D,4,5,8,7,block_1\n"
               "3,4,QUAD_4_2D,5,6,9,8,block_2";
  }
  else if (blockOrder == 4) {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_2\n"
               "1,2,QUAD_4_2D,2,3,6,5,block_3\n"
               "2,3,QUAD_4_2D,4,5,8,7,block_4\n"
               "3,4,QUAD_4_2D,5,6,9,8,block_1";
  }

  coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };
}

stk::mesh::PartVector setup_mesh_4block_4quad(stk::mesh::BulkData& bulk, unsigned blockOrder)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);
  stk::mesh::Part & block3 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_3", 3);
  stk::mesh::Part & block4 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_4", 4);

  std::string meshDesc;
  std::vector<double> coordinates;
  switch(bulk.parallel_size()) {
  case 1:
    fill_mesh_description_4block_4quad_np1(bulk, blockOrder, meshDesc, coordinates);
    break;
  case 2:
    fill_mesh_description_4block_4quad_np2(bulk, blockOrder, meshDesc, coordinates);
    break;
  case 3:
    fill_mesh_description_4block_4quad_np3(bulk, blockOrder, meshDesc, coordinates);
    break;
  case 4:
    fill_mesh_description_4block_4quad_np4(bulk, blockOrder, meshDesc, coordinates);
    break;
  default:
    STK_ThrowRequireMsg(false, "Unexpected proc count for this test\n");
    break;
  }

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

  EXPECT_EQ(5u, get_num_intersecting_nodes(bulk, {&block1, &block2, &block3, &block4}));
  EXPECT_EQ(9u, get_num_total_nodes(bulk));

  return {&block1, &block2, &block3, &block4};
}

void fill_mesh_description_6block_6quad_np1(stk::mesh::BulkData& bulk, std::string& meshDesc, std::vector<double>& coordinates)
{
  STK_ThrowRequire(bulk.parallel_size() == 1);

  meshDesc = "0,1,QUAD_4_2D,1,2,6,5,block_1\n"
             "0,2,QUAD_4_2D,2,3,7,6,block_2\n"
             "0,3,QUAD_4_2D,3,4,8,7,block_3\n"
             "0,4,QUAD_4_2D,5,6,10,9,block_4\n"
             "0,5,QUAD_4_2D,6,7,11,10,block_5\n"
             "0,6,QUAD_4_2D,7,8,12,11,block_6\n";

  coordinates = { 0,0, 1,0, 2,0, 3,0,
                  0,1, 1,1, 2,1, 3,1,
                  0,2, 1,2, 2,2, 3,2};
}

stk::mesh::PartVector setup_mesh_6block_6quad(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);
  stk::mesh::Part & block3 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_3", 3);
  stk::mesh::Part & block4 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_4", 4);
  stk::mesh::Part & block5 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_5", 5);
  stk::mesh::Part & block6 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_6", 6);

  std::string meshDesc;
  std::vector<double> coordinates;
  switch(bulk.parallel_size()) {
  case 1:
    fill_mesh_description_6block_6quad_np1(bulk, meshDesc, coordinates);
    break;
  default:
    STK_ThrowRequireMsg(false, "Unexpected proc count for this test\n");
    break;
  }

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

  EXPECT_EQ(8u, get_num_intersecting_nodes(bulk, {&block1, &block2, &block3, &block4, &block5, &block6}));
  EXPECT_EQ(12u, get_num_total_nodes(bulk));

  output_mesh(bulk, "initial.g");

  return {&block1, &block2, &block3, &block4, &block5, &block6};
}

void fill_mesh_description_9block_9quad_np1(stk::mesh::BulkData& bulk, std::string& meshDesc, std::vector<double>& coordinates)
{
  STK_ThrowRequire(bulk.parallel_size() == 1);

  meshDesc = "0,1,QUAD_4_2D,1,2,6,5,block_1\n"
             "0,2,QUAD_4_2D,2,3,7,6,block_2\n"
             "0,3,QUAD_4_2D,3,4,8,7,block_3\n"
             "0,4,QUAD_4_2D,5,6,10,9,block_4\n"
             "0,5,QUAD_4_2D,6,7,11,10,block_5\n"
             "0,6,QUAD_4_2D,7,8,12,11,block_6\n"
             "0,7,QUAD_4_2D,9,10,14,13,block_7\n"
             "0,8,QUAD_4_2D,10,11,15,14,block_8\n"
             "0,9,QUAD_4_2D,11,12,16,15,block_9\n";

  coordinates = { 0,0, 1,0, 2,0, 3,0,
                  0,1, 1,1, 2,1, 3,1,
                  0,2, 1,2, 2,2, 3,2,
                  0,3, 1,3, 2,3, 3,3};
}

stk::mesh::PartVector setup_mesh_9block_9quad(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);
  stk::mesh::Part & block3 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_3", 3);
  stk::mesh::Part & block4 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_4", 4);
  stk::mesh::Part & block5 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_5", 5);
  stk::mesh::Part & block6 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_6", 6);
  stk::mesh::Part & block7 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_7", 7);
  stk::mesh::Part & block8 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_8", 8);
  stk::mesh::Part & block9 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_9", 9);

  std::string meshDesc;
  std::vector<double> coordinates;
  switch(bulk.parallel_size()) {
  case 1:
    fill_mesh_description_9block_9quad_np1(bulk, meshDesc, coordinates);
    break;
  default:
    STK_ThrowRequireMsg(false, "Unexpected proc count for this test\n");
    break;
  }

  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

  EXPECT_EQ(12u, get_num_intersecting_nodes(bulk, {&block1, &block2, &block3, &block4, &block5, &block6, &block7, &block8, &block9}));
  EXPECT_EQ(16u, get_num_total_nodes(bulk));

  output_mesh(bulk, "initial.g");

  return {&block1, &block2, &block3, &block4, &block5, &block6, &block7, &block8, &block9};
}

stk::mesh::PartVector setup_mesh_2block_1hex(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part * block2 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_2", 2);
  stk::io::fill_mesh("generated:1x1x1", bulk);
  stk::mesh::Part * block1 = bulk.mesh_meta_data().get_part("block_1");

  EXPECT_EQ(0u, get_num_intersecting_nodes(bulk, {block1, block2}));
  EXPECT_EQ(8u, get_num_total_nodes(bulk));
  EXPECT_FALSE(check_orphaned_nodes(bulk));

  return {block1, block2};
}

stk::mesh::PartVector setup_mesh_2block_2hex(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part * block2 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_2", 2);
  stk::io::fill_mesh("generated:1x1x2", bulk);
  stk::mesh::Part * block1 = bulk.mesh_meta_data().get_part("block_1");

  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{2}, "block_1", "block_2");

  EXPECT_EQ( 4u, get_num_intersecting_nodes(bulk, {block1, block2}));
  EXPECT_EQ(12u, get_num_total_nodes(bulk));
  EXPECT_FALSE(check_orphaned_nodes(bulk));

  return {block1, block2};
}

stk::mesh::PartVector setup_mesh_2block_2hex_with_internal_sides(stk::mesh::BulkData& bulk, bool loadMeshFirst)
{
  stk::mesh::Part * block2 = nullptr;

  if(!loadMeshFirst) {
    block2 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_2", 2);
    create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4, "surface_1", 1);
  }

  stk::io::fill_mesh("generated:1x1x2", bulk);

  if(loadMeshFirst) {
    block2 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_2", 2);
    create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4, "surface_1", 1);
  }

  stk::mesh::Part * block1 = bulk.mesh_meta_data().get_part("block_1");
  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{2}, "block_1", "block_2");

  create_sideset(bulk, "surface_1", {"block_1"});
  create_sides_between_blocks(bulk, "block_1", "block_2", "surface_1");

  EXPECT_EQ( 4u, get_num_intersecting_nodes(bulk, {block1, block2}));
  EXPECT_EQ( 1u, get_num_common_sides(bulk, {block1, block2}));
  EXPECT_EQ(12u, get_num_total_nodes(bulk));
  EXPECT_EQ( 1u, get_num_total_sides(bulk));
  EXPECT_FALSE(check_orphaned_nodes(bulk));

  return {block1, block2};
}

stk::mesh::PartVector setup_mesh_2block_2hex_with_internal_and_external_sides(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part * block2 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_2", 2);
  create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4, "surface_1", 1);
  create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4, "surface_2", 2);

  stk::io::fill_mesh("generated:1x1x2", bulk);
  stk::mesh::Part * block1 = bulk.mesh_meta_data().get_part("block_1");

  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{2}, "block_1", "block_2");

  create_sideset(bulk, "surface_1", {"block_1"});
  create_sides_between_blocks(bulk, "block_1", "block_2", "surface_1");
  create_all_boundary_sides(bulk, "surface_2");

  EXPECT_EQ( 4u, get_num_intersecting_nodes(bulk, {block1, block2}));
  EXPECT_EQ( 1u, get_num_common_sides(bulk, {block1, block2}));
  EXPECT_EQ(12u, get_num_total_nodes(bulk));
  EXPECT_EQ(11u, get_num_total_sides(bulk));
  EXPECT_FALSE(check_orphaned_nodes(bulk));

  return {block1, block2};
}

stk::mesh::PartVector setup_mesh_2block_2hex_with_dual_internal_and_external_sides(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part * block2 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_2", 2);
  create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4, "surface_1", 1);
  create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4, "surface_2", 2);

  stk::io::fill_mesh("generated:1x1x2", bulk);
  stk::mesh::Part * block1 = bulk.mesh_meta_data().get_part("block_1");

  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{2}, "block_1", "block_2");

  create_sideset(bulk, "surface_1", {"block_1", "block_2"});
  create_sides_between_blocks(bulk, "block_1", "block_2", "surface_1");
  create_all_boundary_sides(bulk, "surface_2");

  EXPECT_EQ( 4u, get_num_intersecting_nodes(bulk, {block1, block2}));
  EXPECT_EQ( 1u, get_num_common_sides(bulk, {block1, block2}));
  EXPECT_EQ(12u, get_num_total_nodes(bulk));
  EXPECT_EQ(11u, get_num_total_sides(bulk));
  EXPECT_FALSE(check_orphaned_nodes(bulk));

  return {block1, block2};
}

stk::mesh::PartVector setup_mesh_2block_2hex_with_external_sides(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part * block2 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_2", 2);
  create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4, "surface_1", 1);

  stk::io::fill_mesh("generated:1x1x2", bulk);
  stk::mesh::Part * block1 = bulk.mesh_meta_data().get_part("block_1");

  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{2}, "block_1", "block_2");

  create_sideset(bulk, "surface_1", {"block_1", "block_2"});
  create_all_boundary_sides(bulk, "surface_1");

  EXPECT_EQ( 4u, get_num_intersecting_nodes(bulk, {block1, block2}));
  EXPECT_EQ(12u, get_num_total_nodes(bulk));
  EXPECT_EQ(10u, get_num_total_sides(bulk));
  EXPECT_FALSE(check_orphaned_nodes(bulk));

  return {block1, block2};
}

stk::mesh::PartVector setup_mesh_3block_4hex(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part * block2 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_2", 2);
  stk::mesh::Part * block3 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_3", 3);

  stk::io::fill_mesh("generated:1x2x2", bulk);
  stk::mesh::Part * block1 = bulk.mesh_meta_data().get_part("block_1");

  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{2}, "block_1", "block_2");
  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{4}, "block_1", "block_3");

  EXPECT_EQ( 8u, get_num_intersecting_nodes(bulk, {block1, block2, block3}));
  EXPECT_EQ(18u, get_num_total_nodes(bulk));
  EXPECT_FALSE(check_orphaned_nodes(bulk));

  return {block1, block2, block3};
}

stk::mesh::PartVector setup_mesh_3block_4hex_with_internal_sides(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part * block2 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_2", 2);
  stk::mesh::Part * block3 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_3", 3);
  create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4, "surface_1", 1);
  create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4, "surface_2", 2);

  stk::io::fill_mesh("generated:1x2x2", bulk);
  stk::mesh::Part * block1 = bulk.mesh_meta_data().get_part("block_1");

  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{2}, "block_1", "block_2");
  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{4}, "block_1", "block_3");

  create_sideset(bulk, "surface_1", {"block_1"});
  create_sides_between_blocks(bulk, "block_1", "block_2", "surface_1");
  create_sides_between_blocks(bulk, "block_1", "block_3", "surface_1");
  create_sideset(bulk, "surface_2", {"block_2"});
  create_sides_between_blocks(bulk, "block_2", "block_3", "surface_2");

  EXPECT_EQ( 8u, get_num_intersecting_nodes(bulk, {block1, block2, block3}));
  EXPECT_EQ( 3u, get_num_common_sides(bulk, {block1, block2, block3}));
  EXPECT_EQ(18u, get_num_total_nodes(bulk));
  EXPECT_EQ( 3u, get_num_total_sides(bulk));
  EXPECT_FALSE(check_orphaned_nodes(bulk));

  return {block1, block2, block3};
}

stk::mesh::PartVector setup_mesh_2block_2cubeOfTet(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part * block2 = &create_part(bulk.mesh_meta_data(), stk::topology::TET_4, "block_2", 2);

  stk::io::fill_mesh("generated:1x1x2|tets", bulk);
  stk::mesh::Part * block1 = bulk.mesh_meta_data().get_part("block_1");

  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{7, 8, 9, 10, 11, 12},
                                 "block_1", "block_2");

  EXPECT_EQ( 4u, get_num_intersecting_nodes(bulk, {block1, block2}));
  EXPECT_EQ(12u, get_num_total_nodes(bulk));
  EXPECT_FALSE(check_orphaned_nodes(bulk));

  return {block1, block2};
}

stk::mesh::PartVector setup_mesh_2block_2cubeOfTet_with_internal_sides(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part * block2 = &create_part(bulk.mesh_meta_data(), stk::topology::TET_4, "block_2", 2);
  create_part(bulk.mesh_meta_data(), stk::topology::TRI_3, "surface_1", 1);

  stk::io::fill_mesh("generated:1x1x2|tets", bulk);
  stk::mesh::Part * block1 = bulk.mesh_meta_data().get_part("block_1");

  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{7, 8, 9, 10, 11, 12},
                                 "block_1", "block_2");

  create_sideset(bulk, "surface_1", {"block_1"});
  create_sides_between_blocks(bulk, "block_1", "block_2", "surface_1");

  EXPECT_EQ( 4u, get_num_intersecting_nodes(bulk, {block1, block2}));
  EXPECT_EQ( 2u, get_num_common_sides(bulk, {block1, block2}));
  EXPECT_EQ(12u, get_num_total_nodes(bulk));
  EXPECT_EQ( 2u, get_num_total_sides(bulk));
  EXPECT_FALSE(check_orphaned_nodes(bulk));

  return {block1, block2};
}

stk::mesh::PartVector setup_mesh_3block_4cubeOfTet(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part * block2 = &create_part(bulk.mesh_meta_data(), stk::topology::TET_4, "block_2", 2);
  stk::mesh::Part * block3 = &create_part(bulk.mesh_meta_data(), stk::topology::TET_4, "block_3", 3);

  stk::io::fill_mesh("generated:1x2x2|tets", bulk);
  stk::mesh::Part * block1 = bulk.mesh_meta_data().get_part("block_1");

  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{7, 8, 9, 10, 11, 12}, "block_1", "block_2");
  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{19, 20, 21, 22, 23, 24}, "block_1", "block_3");

  EXPECT_EQ( 8u, get_num_intersecting_nodes(bulk, {block1, block2, block3}));
  EXPECT_EQ(18u, get_num_total_nodes(bulk));
  EXPECT_FALSE(check_orphaned_nodes(bulk));

  return {block1, block2, block3};
}

stk::mesh::PartVector setup_mesh_3block_4cubeOfTet_with_internal_sides(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part * block2 = &create_part(bulk.mesh_meta_data(), stk::topology::TET_4, "block_2", 2);
  stk::mesh::Part * block3 = &create_part(bulk.mesh_meta_data(), stk::topology::TET_4, "block_3", 3);
  create_part(bulk.mesh_meta_data(), stk::topology::TRI_3, "surface_1", 1);
  create_part(bulk.mesh_meta_data(), stk::topology::TRI_3, "surface_2", 2);

  stk::io::fill_mesh("generated:1x2x2|tets", bulk);
  stk::mesh::Part * block1 = bulk.mesh_meta_data().get_part("block_1");

  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{7, 8, 9, 10, 11, 12}, "block_1", "block_2");
  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{19, 20, 21, 22, 23, 24}, "block_1", "block_3");

  create_sideset(bulk, "surface_1", {"block_1"});
  create_sides_between_blocks(bulk, "block_1", "block_2", "surface_1");
  create_sides_between_blocks(bulk, "block_1", "block_3", "surface_1");
  create_sideset(bulk, "surface_2", {"block_2"});
  create_sides_between_blocks(bulk, "block_2", "block_3", "surface_2");

  EXPECT_EQ( 8u, get_num_intersecting_nodes(bulk, {block1, block2, block3}));
  EXPECT_EQ( 6u, get_num_common_sides(bulk, {block1, block2, block3}));
  EXPECT_EQ(18u, get_num_total_nodes(bulk));
  EXPECT_EQ( 6u, get_num_total_sides(bulk));
  EXPECT_FALSE(check_orphaned_nodes(bulk));

  return {block1, block2, block3};
}

stk::mesh::PartVector setup_mesh_4block_4hex(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part * block2 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_2", 2);
  stk::mesh::Part * block3 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_3", 3);
  stk::mesh::Part * block4 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_4", 4);

  stk::io::fill_mesh("generated:1x2x2", bulk);
  stk::mesh::Part * block1 = bulk.mesh_meta_data().get_part("block_1");

  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{2}, "block_1", "block_2");
  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{3}, "block_1", "block_3");
  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{4}, "block_1", "block_4");

  EXPECT_EQ(10u, get_num_intersecting_nodes(bulk, {block1, block2, block3, block4}));
  EXPECT_EQ(18u, get_num_total_nodes(bulk));
  EXPECT_FALSE(check_orphaned_nodes(bulk));

  return {block1, block2, block3, block4};
}

stk::mesh::PartVector setup_mesh_4block_4hex_vertical_stack(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part * block2 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_2", 2);
  stk::mesh::Part * block3 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_3", 3);
  stk::mesh::Part * block4 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_4", 4);

  stk::io::fill_mesh("generated:1x1x4", bulk);
  stk::mesh::Part * block1 = bulk.mesh_meta_data().get_part("block_1");

  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{2}, "block_1", "block_2");
  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{3}, "block_1", "block_3");
  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{4}, "block_1", "block_4");

  EXPECT_EQ(12u, get_num_intersecting_nodes(bulk, {block1, block2, block3, block4}));
  EXPECT_EQ(20u, get_num_total_nodes(bulk));
  EXPECT_FALSE(check_orphaned_nodes(bulk));

  return {block1, block2, block3, block4};
}

stk::mesh::PartVector setup_mesh_4block_8hex_cube(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part * block2 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_2", 2);
  stk::mesh::Part * block3 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_3", 3);
  stk::mesh::Part * block4 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_4", 4);

  stk::io::fill_mesh("generated:2x2x2", bulk);
  stk::mesh::Part * block1 = bulk.mesh_meta_data().get_part("block_1");

  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{3}, "block_1", "block_2");
  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{4}, "block_1", "block_2");
  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{5}, "block_1", "block_3");
  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{6}, "block_1", "block_3");
  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{7}, "block_1", "block_4");
  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{8}, "block_1", "block_4");

  EXPECT_EQ(15u, get_num_intersecting_nodes(bulk, {block1, block2, block3, block4}));
  EXPECT_EQ(27u, get_num_total_nodes(bulk));
  EXPECT_FALSE(check_orphaned_nodes(bulk));

  return {block1, block2, block3, block4};
}

stk::mesh::PartVector setup_mesh_8block_8hex_cube(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part * block2 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_2", 2);
  stk::mesh::Part * block3 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_3", 3);
  stk::mesh::Part * block4 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_4", 4);
  stk::mesh::Part * block5 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_5", 5);
  stk::mesh::Part * block6 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_6", 6);
  stk::mesh::Part * block7 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_7", 7);
  stk::mesh::Part * block8 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_8", 8);

  stk::io::fill_mesh("generated:2x2x2", bulk);
  stk::mesh::Part * block1 = bulk.mesh_meta_data().get_part("block_1");

  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{2}, "block_1", "block_2");
  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{3}, "block_1", "block_3");
  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{4}, "block_1", "block_4");
  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{5}, "block_1", "block_5");
  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{6}, "block_1", "block_6");
  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{7}, "block_1", "block_7");
  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{8}, "block_1", "block_8");

  EXPECT_EQ(19u, get_num_intersecting_nodes(bulk, {block1, block2, block3, block4, block5, block6, block7, block8}));
  EXPECT_EQ(27u, get_num_total_nodes(bulk));
  EXPECT_FALSE(check_orphaned_nodes(bulk));

  return {block1, block2, block3, block4, block5, block6, block7, block8};
}

std::vector<std::string> get_part_names(const stk::mesh::PartVector& parts)
{
  std::vector<std::string> partNames;

  for(stk::mesh::Part* part : parts) {
    partNames.push_back(part->name());
  }

  return partNames;
}

stk::mesh::PartVector setup_mesh_8block_8hex_with_external_sides(stk::mesh::BulkData& bulk)
{
  create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4, "surface_1", 1);

  stk::mesh::PartVector createdParts = setup_mesh_8block_8hex_cube(bulk);

  create_sideset(bulk, "surface_1", get_part_names(createdParts));
  create_all_boundary_sides(bulk, "surface_1");

  EXPECT_EQ(27u, get_num_total_nodes(bulk));
  EXPECT_EQ(24u, get_num_total_sides(bulk));
  EXPECT_FALSE(check_orphaned_nodes(bulk));

  return createdParts;
}
