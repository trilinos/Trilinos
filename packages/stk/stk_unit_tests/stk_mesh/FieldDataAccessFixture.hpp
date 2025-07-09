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

#ifndef STK_FIELDDATAACCESSFIXTURE_HPP
#define STK_FIELDDATAACCESSFIXTURE_HPP

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for field_bytes_per_entity, etc
#include "stk_topology/topology.hpp"    // for topology, etc

namespace {

class FieldDataAccessFixture : public stk::unit_test_util::MeshFixture
{
public:
  FieldDataAccessFixture()
    : m_field(nullptr),
      m_leftField(nullptr),
      m_rightField(nullptr)
  {}

  stk::mesh::Entity create_node(stk::mesh::EntityId nodeId)
  {
    get_bulk().modification_begin();
    stk::mesh::Entity node = get_bulk().declare_node(nodeId);
    get_bulk().modification_end();

    return node;
  }

  void create_single_element_mesh()
  {
    stk::io::fill_mesh("generated:1x1x1", get_bulk());
  }

  void build_mesh_with_scalar_field()
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);

    m_field = &get_meta().declare_field<int>(stk::topology::NODE_RANK, "nodeField");
    stk::mesh::put_field_on_mesh(*m_field, get_meta().universal_part(), nullptr);

    stk::io::fill_mesh("generated:2x1x1", get_bulk());
  }

  void build_mesh_with_scalar_left_field()
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);

    m_leftField = &get_meta().declare_field<int, stk::mesh::Layout::Left>(stk::topology::NODE_RANK, "nodeLeftField");
    stk::mesh::put_field_on_mesh(*m_leftField, get_meta().universal_part(), nullptr);

    stk::io::fill_mesh("generated:2x1x1", get_bulk());
  }

  void build_mesh_with_scalar_right_field()
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);

    m_rightField = &get_meta().declare_field<int, stk::mesh::Layout::Right>(stk::topology::NODE_RANK, "nodeRightField");
    stk::mesh::put_field_on_mesh(*m_rightField, get_meta().universal_part(), nullptr);

    stk::io::fill_mesh("generated:2x1x1", get_bulk());
  }

  void build_mesh_with_multi_component_field()
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);

    m_field = &get_meta().declare_field<int>(stk::topology::NODE_RANK, "nodeField");
    stk::mesh::put_field_on_mesh(*m_field, get_meta().universal_part(), 3, nullptr);

    stk::io::fill_mesh("generated:2x1x1", get_bulk());
  }

  void build_mesh_with_multi_component_left_field()
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);

    m_leftField = &get_meta().declare_field<int, stk::mesh::Layout::Left>(stk::topology::NODE_RANK, "nodeLeftField");
    stk::mesh::put_field_on_mesh(*m_leftField, get_meta().universal_part(), 3, nullptr);

    stk::io::fill_mesh("generated:2x1x1", get_bulk());
  }

  void build_mesh_with_multi_component_right_field()
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);

    m_rightField = &get_meta().declare_field<int, stk::mesh::Layout::Right>(stk::topology::NODE_RANK, "nodeRightField");
    stk::mesh::put_field_on_mesh(*m_rightField, get_meta().universal_part(), 3, nullptr);

    stk::io::fill_mesh("generated:2x1x1", get_bulk());
  }

  void build_mesh_with_multi_copy_field()
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);

    m_field = &get_meta().declare_field<int>(stk::topology::NODE_RANK, "nodeField");
    stk::mesh::put_field_on_mesh(*m_field, get_meta().universal_part(), 1, 8, nullptr);

    stk::io::fill_mesh("generated:2x1x1", get_bulk());
  }

  void build_mesh_with_multi_copy_left_field()
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);

    m_leftField = &get_meta().declare_field<int, stk::mesh::Layout::Left>(stk::topology::NODE_RANK, "nodeLeftField");
    stk::mesh::put_field_on_mesh(*m_leftField, get_meta().universal_part(), 1, 8, nullptr);

    stk::io::fill_mesh("generated:2x1x1", get_bulk());
  }

  void build_mesh_with_multi_copy_right_field()
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);

    m_rightField = &get_meta().declare_field<int, stk::mesh::Layout::Right>(stk::topology::NODE_RANK, "nodeRightField");
    stk::mesh::put_field_on_mesh(*m_rightField, get_meta().universal_part(), 1, 8, nullptr);

    stk::io::fill_mesh("generated:2x1x1", get_bulk());
  }

  void build_mesh_with_multi_copy_multi_component_field()
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);

    m_field = &get_meta().declare_field<int>(stk::topology::NODE_RANK, "nodeField");
    stk::mesh::put_field_on_mesh(*m_field, get_meta().universal_part(), 3, 8, nullptr);

    stk::io::fill_mesh("generated:2x1x1", get_bulk());
  }

  void build_mesh_with_multi_copy_multi_component_left_field()
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);

    m_leftField = &get_meta().declare_field<int, stk::mesh::Layout::Left>(stk::topology::NODE_RANK, "nodeLeftField");
    stk::mesh::put_field_on_mesh(*m_leftField, get_meta().universal_part(), 3, 8, nullptr);

    stk::io::fill_mesh("generated:2x1x1", get_bulk());
  }

  void build_mesh_with_multi_copy_multi_component_right_field()
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);

    m_rightField = &get_meta().declare_field<int, stk::mesh::Layout::Right>(stk::topology::NODE_RANK, "nodeRightField");
    stk::mesh::put_field_on_mesh(*m_rightField, get_meta().universal_part(), 3, 8, nullptr);

    stk::io::fill_mesh("generated:2x1x1", get_bulk());
  }

protected:
  stk::mesh::Field<int>* m_field;
  stk::mesh::Field<int, stk::mesh::Layout::Left>* m_leftField;
  stk::mesh::Field<int, stk::mesh::Layout::Right>* m_rightField;
};

}

#endif // STK_FIELDDATAACCESSFIXTURE_HPP
