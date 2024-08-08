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

/*
 * ConstructedMesh.hpp
 *
 *  Created on: Oct 9, 2020
 *      Author: tookusa
 */

#ifndef STK_STK_UNIT_TEST_UTILS_STK_UNIT_TEST_UTILS_CONSTRUCTEDMESH_HPP_
#define STK_STK_UNIT_TEST_UTILS_STK_UNIT_TEST_UTILS_CONSTRUCTEDMESH_HPP_

#include <gtest/gtest.h>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <string>
#include <vector>

namespace stk
{
namespace unit_test_util
{

struct ConstructedElementBlock
{
  ConstructedElementBlock()
    : topology(stk::topology::INVALID_TOPOLOGY),
      name(""),
      id(-1)
  { }

  ConstructedElementBlock(stk::topology topology_, const std::string& name_, int id_)
    : topology(topology_),
      name(name_),
      id(id_)
  { }

  ConstructedElementBlock(stk::topology topology_, const std::string& name_, int id_, const std::vector< std::vector<unsigned> >& connectivityIndex_)
    : topology(topology_),
      name(name_),
      id(id_),
      connectivityIndex(connectivityIndex_)
  { }

  ConstructedElementBlock(const ConstructedElementBlock& block)
    : topology(block.topology),
      name(block.name),
      id(block.id),
      connectivityIndex(block.connectivityIndex)
  { }

  void add_connectivity(const std::vector<unsigned>& connectivity)
  {
    ASSERT_EQ(topology.num_nodes(), connectivity.size());
    connectivityIndex.push_back(connectivity);
  }

  void set_connectivity(const std::vector< std::vector<unsigned> >& connectivities)
  {
    for(const std::vector<unsigned>& connectivity : connectivities) {
      add_connectivity(connectivity);
    }
  }

  stk::topology topology;
  std::string name;
  int id;
  std::vector< std::vector<unsigned> > connectivityIndex;
};

class ConstructedMesh
{
public:
  ConstructedMesh(unsigned spatialDimension)
    : m_spatialDimension(spatialDimension)
  { }

  void set_x_coordinates(const std::vector< double >& xCoordinates)
  {
    m_xCoords = xCoordinates;
  }

  void set_y_coordinates(const std::vector< double >& yCoordinates)
  {
    m_yCoords = yCoordinates;
  }

  void set_z_coordinates(const std::vector< double >& zCoordinates)
  {
    m_zCoords = zCoordinates;
  }

  void set_node_ids(const stk::mesh::EntityIdVector& nodeIds)
  {
    m_nodeIds = nodeIds;
  }

  void set_elem_ids(const stk::mesh::EntityIdVector& elemIds)
  {
    m_elemIds = elemIds;
  }

  void set_elem_distribution(const stk::mesh::EntityIdProcVec& elemDistribution)
  {
    m_elemDistribution = elemDistribution;
  }


  void add_elem_block(const ConstructedElementBlock& block)
  {
    m_elemBlocks.push_back(block);
  }

  void create_mesh(stk::mesh::BulkData& bulk)
  {
    verify_mesh_data(bulk.mesh_meta_data());
    populate_bulk_data(bulk);
  }

private:
  unsigned m_spatialDimension = 0;
  std::vector< double > m_xCoords;
  std::vector< double > m_yCoords;
  std::vector< double > m_zCoords;
  stk::mesh::EntityIdVector m_nodeIds;
  stk::mesh::EntityIdVector m_elemIds;
  stk::mesh::EntityIdProcVec m_elemDistribution;
  std::vector<ConstructedElementBlock> m_elemBlocks;

  ConstructedMesh() {}

  void create_block_elements_and_nodes(stk::mesh::BulkData& bulk, const ConstructedElementBlock& block, const unsigned elemIdOffset);

  void verify_mesh_data(const stk::mesh::MetaData& meta);

  void populate_bulk_data(stk::mesh::BulkData& bulk);
};

namespace simple_fields {

struct ConstructedElementBlock
{
  ConstructedElementBlock()
    : topology(stk::topology::INVALID_TOPOLOGY),
      name(""),
      id(-1)
  { }

  ConstructedElementBlock(stk::topology topology_, const std::string& name_, int id_)
    : topology(topology_),
      name(name_),
      id(id_)
  { }

  ConstructedElementBlock(stk::topology topology_, const std::string& name_, int id_, const std::vector< std::vector<unsigned> >& connectivityIndex_)
    : topology(topology_),
      name(name_),
      id(id_),
      connectivityIndex(connectivityIndex_)
  { }

  ConstructedElementBlock(const ConstructedElementBlock& block)
    : topology(block.topology),
      name(block.name),
      id(block.id),
      connectivityIndex(block.connectivityIndex)
  { }

  void add_connectivity(const std::vector<unsigned>& connectivity)
  {
    ASSERT_EQ(topology.num_nodes(), connectivity.size());
    connectivityIndex.push_back(connectivity);
  }

  void set_connectivity(const std::vector< std::vector<unsigned> >& connectivities)
  {
    for(const std::vector<unsigned>& connectivity : connectivities) {
      add_connectivity(connectivity);
    }
  }

  stk::topology topology;
  std::string name;
  int id;
  std::vector< std::vector<unsigned> > connectivityIndex;
};

class ConstructedMesh
{
public:
  ConstructedMesh(unsigned spatialDimension)
    : m_spatialDimension(spatialDimension)
  { }

  void set_x_coordinates(const std::vector< double >& xCoordinates)
  {
    m_xCoords = xCoordinates;
  }

  void set_y_coordinates(const std::vector< double >& yCoordinates)
  {
    m_yCoords = yCoordinates;
  }

  void set_z_coordinates(const std::vector< double >& zCoordinates)
  {
    m_zCoords = zCoordinates;
  }

  void set_node_ids(const stk::mesh::EntityIdVector& nodeIds)
  {
    m_nodeIds = nodeIds;
  }

  void set_elem_ids(const stk::mesh::EntityIdVector& elemIds)
  {
    m_elemIds = elemIds;
  }

  void set_elem_distribution(const stk::mesh::EntityIdProcVec& elemDistribution)
  {
    m_elemDistribution = elemDistribution;
  }


  void add_elem_block(const ConstructedElementBlock& block)
  {
    m_elemBlocks.push_back(block);
  }

  void create_mesh(stk::mesh::BulkData& bulk)
  {
    verify_mesh_data(bulk.mesh_meta_data());
    populate_bulk_data(bulk);
  }

private:
  unsigned m_spatialDimension = 0;
  std::vector< double > m_xCoords;
  std::vector< double > m_yCoords;
  std::vector< double > m_zCoords;
  stk::mesh::EntityIdVector m_nodeIds;
  stk::mesh::EntityIdVector m_elemIds;
  stk::mesh::EntityIdProcVec m_elemDistribution;
  std::vector<ConstructedElementBlock> m_elemBlocks;

  ConstructedMesh() {}

  void create_block_elements_and_nodes(stk::mesh::BulkData& bulk, const ConstructedElementBlock& block, const unsigned elemIdOffset);

  void verify_mesh_data(const stk::mesh::MetaData& meta);

  void populate_bulk_data(stk::mesh::BulkData& bulk);
};

} // namespace simple_fields

} // namespace unit_test_util
} // namespace stk




#endif /* STK_STK_UNIT_TEST_UTILS_STK_UNIT_TEST_UTILS_CONSTRUCTEDMESH_HPP_ */
