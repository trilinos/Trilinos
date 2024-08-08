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

#include "GenerateALefRAMesh.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_unit_test_utils/FaceTestingUtils.hpp"
#include "stk_unit_test_utils/ReadWriteSidesetTester.hpp"
#include "stk_unit_test_utils/TextMesh.hpp"
#include <gtest/gtest.h>
#include <stk_io/StkIoUtils.hpp>
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_unit_test_utils/ioUtils.hpp>

namespace stk {
namespace unit_test_util {

stk::mesh::PartVector create_sideset_parts(stk::mesh::MetaData &meta, const std::vector<std::string>&names)
{
  stk::mesh::PartVector parts;

  int id = 1;
  for(const std::string &name : names)
  {
    stk::mesh::Part &part = meta.declare_part_with_topology(name, stk::topology::QUAD_4);
    meta.set_part_id(part, id);
    stk::io::put_io_part_attribute(part);
    parts.push_back(&part);
    ++id;
  }

  return parts;
}

void create_AA_mesh(stk::mesh::BulkData &bulk, ElementOrdering elemOrdering)
{
  std::string meshDesc;
  if (elemOrdering == INCREASING) {
    meshDesc = "0,1,HEX_8,1,2,3,4,5, 6, 7, 8,block_1\n"
               "0,2,HEX_8,5,6,7,8,9,10,11,12,block_1";
  }
  else if (elemOrdering == DECREASING) {
    meshDesc = "0,1,HEX_8,5,6,7,8,9,10,11,12,block_1\n"
               "0,2,HEX_8,1,2,3,4,5, 6, 7, 8,block_1";
  }

  std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0,
                                      0,0,1, 1,0,1, 1,1,1, 0,1,1,
                                      0,0,2, 1,0,2, 1,1,2, 0,1,2 };

  bulk.initialize_face_adjacent_element_graph();
  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));
}

void create_AB_mesh(stk::mesh::BulkData &bulk, ElementOrdering elemOrdering)
{
  std::string meshDesc;
  if (elemOrdering == INCREASING) {
    meshDesc = "0,1,HEX_8,1,2,3,4,5, 6, 7, 8,block_1\n"
               "0,2,HEX_8,5,6,7,8,9,10,11,12,block_2";
  }
  else if (elemOrdering == DECREASING) {
    meshDesc = "0,1,HEX_8,5,6,7,8,9,10,11,12,block_2\n"
               "0,2,HEX_8,1,2,3,4,5, 6, 7, 8,block_1";
  }

  std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0,
                                      0,0,1, 1,0,1, 1,1,1, 0,1,1,
                                      0,0,2, 1,0,2, 1,1,2, 0,1,2 };

  bulk.initialize_face_adjacent_element_graph();
  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));
}

void populate_elem_sides(SidesetDirection direction,
                         ElementOrdering elemOrdering,
                         stk::mesh::EntityIdVector &elem,
                         std::vector<int> &ordinal)
{
  const stk::mesh::EntityId leftElementId  = (elemOrdering == INCREASING) ? 1 : 2;
  const stk::mesh::EntityId rightElementId = (elemOrdering == INCREASING) ? 2 : 1;

  if (direction == LEFT)
  {
    elem.push_back(leftElementId);
    ordinal.push_back(5);
  }
  else if (direction == RIGHT)
  {
    elem.push_back(rightElementId);
    ordinal.push_back(4);
  }
  else
  {
    elem.push_back(leftElementId);
    ordinal.push_back(5);

    elem.push_back(rightElementId);
    ordinal.push_back(4);
  }
}

void populate_sideset_names(SidesetDirection direction,
                            std::vector<std::string> &names)
{
  if(direction == LEFT)
  {
    names = {"surface_1", "surface_block_1_QUAD4_1"};
  }
  else if(direction == RIGHT)
  {
    names = {"surface_1", "surface_block_2_QUAD4_1"};
  }
  else
  {
    names = {"surface_1", "surface_block_1_QUAD4_1", "surface_block_2_QUAD4_1"};
  }
}

void populate_AA_sideset(stk::mesh::BulkData& bulk,
                         SidesetDirection direction,
                         ElementOrdering elemOrdering,
                         const stk::mesh::PartVector& parts)
{
  stk::mesh::MetaData &meta = bulk.mesh_meta_data();
  stk::mesh::EntityIdVector elem;
  std::vector<int> ordinal;
  populate_elem_sides(direction, elemOrdering, elem, ordinal);
  STK_ThrowRequire(elem.size() == ordinal.size());

  stk::mesh::SideSet& sideSet = bulk.create_sideset(*parts[0]);
  sideSet.set_accept_all_internal_non_coincident_entries(false);
  for (unsigned i = 0; i < elem.size(); ++i) {
    sideSet.add( { bulk.get_entity(stk::topology::ELEMENT_RANK, elem[i]),
                   ordinal[i] });
  }

  stk::mesh::Part* block_1 = meta.get_part("block_1");
  EXPECT_TRUE(block_1 != nullptr);

  meta.set_part_id(*block_1, 1);
  std::vector<const stk::mesh::Part*> touchingParts { block_1 };
  meta.set_surface_to_block_mapping(parts[0], touchingParts);

  // tookusa: Order is important for the incremental sideset updater ... surface to block mapping must be set first
  bulk.create_side_entities(sideSet, parts);
}

stk::mesh::Part* create_AA_mesh_with_sideset(stk::mesh::BulkData &bulk,
                                             SidesetDirection direction,
                                             ElementOrdering elemOrdering)
{
  stk::mesh::MetaData &meta = bulk.mesh_meta_data();
  stk::mesh::PartVector parts = create_sideset_parts(meta, std::vector<std::string>{"surface_1"});

  create_AA_mesh(bulk, elemOrdering);

  populate_AA_sideset(bulk, direction, elemOrdering, parts);

  return parts[0];
}

stk::mesh::Part* create_AA_mesh_with_sideset_and_field(stk::mesh::BulkData &bulk,
                                                       SidesetDirection direction,
                                                       ElementOrdering elemOrdering,
                                                       const std::string & fieldName)
{
  stk::mesh::MetaData &meta = bulk.mesh_meta_data();
  stk::mesh::PartVector parts = create_sideset_parts(meta, std::vector<std::string>{"surface_1"});

  const unsigned numberOfStates = 1;
  stk::mesh::Field<double> & field = meta.declare_field<double>(stk::topology::NODE_RANK, fieldName, numberOfStates);
  const double initValue = 123;
  stk::mesh::put_field_on_mesh(field, meta.universal_part(), &initValue);

  create_AA_mesh(bulk, elemOrdering);

  populate_AA_sideset(bulk, direction, elemOrdering, parts);

  return parts[0];
}

void populate_AB_sideset(stk::mesh::BulkData& bulk,
                         SidesetDirection direction,
                         ElementOrdering elemOrdering,
                         const stk::mesh::PartVector& parts)
{
  stk::mesh::MetaData &meta = bulk.mesh_meta_data();
  stk::mesh::EntityIdVector elem;
  std::vector<int> ordinal;

  populate_elem_sides(direction, elemOrdering, elem, ordinal);
  STK_ThrowRequire(elem.size() == ordinal.size());

  stk::mesh::SideSet &sideSet = bulk.create_sideset(*parts[0]);
  sideSet.set_accept_all_internal_non_coincident_entries(false);
  for(unsigned i=0; i<elem.size(); ++i) {
    const stk::mesh::Entity element = bulk.get_entity(stk::topology::ELEMENT_RANK, elem[i]);
    if (bulk.is_valid(element)) {
      sideSet.add({element, ordinal[i]});
    }
  }

  stk::mesh::Part* block_1 = meta.get_part("block_1");
  EXPECT_TRUE(block_1 != nullptr);
  stk::mesh::Part* block_2 = meta.get_part("block_2");
  EXPECT_TRUE(block_2 != nullptr);

  meta.set_part_id(*block_1, 1);
  meta.set_part_id(*block_2, 2);

  std::vector<const stk::mesh::Part*> touchingParts;

  if(direction == LEFT)
  {
    touchingParts = {block_1};
    meta.set_surface_to_block_mapping(parts[1], touchingParts);
  }
  else if(direction == RIGHT)
  {
    touchingParts = {block_2};
    meta.set_surface_to_block_mapping(parts[1], touchingParts);
  }
  else
  {
    touchingParts = {block_1};
    meta.set_surface_to_block_mapping(parts[1], touchingParts);

    touchingParts = {block_2};
    meta.set_surface_to_block_mapping(parts[2], touchingParts);
  }

  // tookusa: Order is important for the incremental sideset updater ... surface to block mapping must be set first
  bulk.create_side_entities(sideSet, parts);
}

stk::mesh::Part* create_AB_mesh_with_sideset(stk::mesh::BulkData &bulk,
                                             SidesetDirection direction,
                                             ElementOrdering elemOrdering)
{
  stk::mesh::MetaData &meta = bulk.mesh_meta_data();
  std::vector<std::string> sideSetNames;
  populate_sideset_names(direction, sideSetNames);

  stk::mesh::PartVector parts = create_sideset_parts(meta, sideSetNames);

  for(unsigned i = 0; i<parts.size(); ++i) {
    meta.set_part_id(*parts[i], 1);
  }

  for(unsigned i=1; i<parts.size(); ++i) {
    meta.declare_part_subset(*parts[0], *parts[i]);
  }

  create_AB_mesh(bulk, elemOrdering);

  populate_AB_sideset(bulk, direction, elemOrdering, parts);

  return parts[0];
}

stk::mesh::Part* create_AB_mesh_with_sideset_and_field(stk::mesh::BulkData &bulk,
                                                       SidesetDirection direction,
                                                       ElementOrdering elemOrdering,
                                                       const std::string & fieldName)
{
  stk::mesh::MetaData &meta = bulk.mesh_meta_data();
  std::vector<std::string> sideSetNames;
  populate_sideset_names(direction, sideSetNames);

  stk::mesh::PartVector parts = create_sideset_parts(meta, sideSetNames);

  for(unsigned i = 0; i<parts.size(); ++i) {
    meta.set_part_id(*parts[i], 1);
  }

  for(unsigned i=1; i<parts.size(); ++i) {
    meta.declare_part_subset(*parts[0], *parts[i]);
  }

  const unsigned numberOfStates = 1;
  stk::mesh::Field<double> & field = meta.declare_field<double>(stk::topology::NODE_RANK, fieldName, numberOfStates);
  const double initValue = 123;
  stk::mesh::put_field_on_mesh(field, meta.universal_part(), &initValue);

  create_AB_mesh(bulk, elemOrdering);

  populate_AB_sideset(bulk, direction, elemOrdering, parts);

  return parts[0];
}

stk::mesh::Part* create_AB_mesh_with_sideset_and_distribution_factors(stk::mesh::BulkData &bulk,
                                                                      SidesetDirection direction,
                                                                      ElementOrdering elemOrdering,
                                                                      const std::string & fieldName,
                                                                      const double initValue)
{
  stk::mesh::MetaData &meta = bulk.mesh_meta_data();
  std::vector<std::string> sideSetNames;
  populate_sideset_names(direction, sideSetNames);

  stk::mesh::PartVector parts = create_sideset_parts(meta, sideSetNames);

  for(unsigned i = 0; i<parts.size(); ++i) {
    meta.set_part_id(*parts[i], 1);
  }

  for(unsigned i=1; i<parts.size(); ++i) {
    meta.declare_part_subset(*parts[0], *parts[i]);
  }

  stk::mesh::Part* sidePart = parts[0];
  STK_ThrowRequire(nullptr != sidePart);
  const unsigned numberOfStates = 1;
  stk::mesh::Field<double> & ssField = meta.declare_field<double>(stk::topology::FACE_RANK, fieldName, numberOfStates);
  for (stk::mesh::Part* part : parts)
  {
    stk::io::set_distribution_factor_field(*part, ssField);
  }
  stk::topology sideTopo = sidePart->topology();
  STK_ThrowRequireMsg(sideTopo != stk::topology::INVALID_TOPOLOGY, "sidePart "<<sidePart->name()<<" has invalid topology.");
  unsigned numNodes = sideTopo.num_nodes();
  std::vector<double> initValVec(numNodes, initValue);
  stk::mesh::put_field_on_mesh(ssField, *sidePart, numNodes, initValVec.data());

  create_AB_mesh(bulk, elemOrdering);

  populate_AB_sideset(bulk, direction, elemOrdering, parts);

  return parts[0];
}

namespace simple_fields {

stk::mesh::PartVector create_sideset_parts(stk::mesh::MetaData &meta, const std::vector<std::string>&names)
{
  stk::mesh::PartVector parts;

  int id = 1;
  for(const std::string &name : names)
  {
    stk::mesh::Part &part = meta.declare_part_with_topology(name, stk::topology::QUAD_4);
    meta.set_part_id(part, id);
    stk::io::put_io_part_attribute(part);
    parts.push_back(&part);
    ++id;
  }

  return parts;
}

void create_AA_mesh(stk::mesh::BulkData &bulk, ElementOrdering elemOrdering)
{
  std::string meshDesc;
  if (elemOrdering == INCREASING) {
    meshDesc = "0,1,HEX_8,1,2,3,4,5, 6, 7, 8,block_1\n"
               "0,2,HEX_8,5,6,7,8,9,10,11,12,block_1";
  }
  else if (elemOrdering == DECREASING) {
    meshDesc = "0,1,HEX_8,5,6,7,8,9,10,11,12,block_1\n"
               "0,2,HEX_8,1,2,3,4,5, 6, 7, 8,block_1";
  }

  std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0,
                                      0,0,1, 1,0,1, 1,1,1, 0,1,1,
                                      0,0,2, 1,0,2, 1,1,2, 0,1,2 };

  bulk.initialize_face_adjacent_element_graph();
  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));
}

void create_AB_mesh(stk::mesh::BulkData &bulk, ElementOrdering elemOrdering)
{
  std::string meshDesc;
  if (elemOrdering == INCREASING) {
    meshDesc = "0,1,HEX_8,1,2,3,4,5, 6, 7, 8,block_1\n"
               "0,2,HEX_8,5,6,7,8,9,10,11,12,block_2";
  }
  else if (elemOrdering == DECREASING) {
    meshDesc = "0,1,HEX_8,5,6,7,8,9,10,11,12,block_2\n"
               "0,2,HEX_8,1,2,3,4,5, 6, 7, 8,block_1";
  }

  std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0,
                                      0,0,1, 1,0,1, 1,1,1, 0,1,1,
                                      0,0,2, 1,0,2, 1,1,2, 0,1,2 };

  bulk.initialize_face_adjacent_element_graph();
  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));
}

void populate_elem_sides(SidesetDirection direction,
                         ElementOrdering elemOrdering,
                         stk::mesh::EntityIdVector &elem,
                         std::vector<int> &ordinal)
{
  const stk::mesh::EntityId leftElementId  = (elemOrdering == INCREASING) ? 1 : 2;
  const stk::mesh::EntityId rightElementId = (elemOrdering == INCREASING) ? 2 : 1;

  if (direction == LEFT)
  {
    elem.push_back(leftElementId);
    ordinal.push_back(5);
  }
  else if (direction == RIGHT)
  {
    elem.push_back(rightElementId);
    ordinal.push_back(4);
  }
  else
  {
    elem.push_back(leftElementId);
    ordinal.push_back(5);

    elem.push_back(rightElementId);
    ordinal.push_back(4);
  }
}

void populate_sideset_names(SidesetDirection direction,
                            std::vector<std::string> &names)
{
  if(direction == LEFT)
  {
    names = {"surface_1", "surface_block_1_QUAD4_1"};
  }
  else if(direction == RIGHT)
  {
    names = {"surface_1", "surface_block_2_QUAD4_1"};
  }
  else
  {
    names = {"surface_1", "surface_block_1_QUAD4_1", "surface_block_2_QUAD4_1"};
  }
}

void populate_AA_sideset(stk::mesh::BulkData& bulk,
                         SidesetDirection direction,
                         ElementOrdering elemOrdering,
                         const stk::mesh::PartVector& parts)
{
  stk::mesh::MetaData &meta = bulk.mesh_meta_data();
  stk::mesh::EntityIdVector elem;
  std::vector<int> ordinal;
  stk::unit_test_util::populate_elem_sides(direction, elemOrdering, elem, ordinal);
  STK_ThrowRequire(elem.size() == ordinal.size());

  stk::mesh::SideSet& sideSet = bulk.create_sideset(*parts[0]);
  sideSet.set_accept_all_internal_non_coincident_entries(false);
  for (unsigned i = 0; i < elem.size(); ++i) {
    sideSet.add( { bulk.get_entity(stk::topology::ELEMENT_RANK, elem[i]),
                   ordinal[i] });
  }

  stk::mesh::Part* block_1 = meta.get_part("block_1");
  EXPECT_TRUE(block_1 != nullptr);

  meta.set_part_id(*block_1, 1);
  std::vector<const stk::mesh::Part*> touchingParts { block_1 };
  meta.set_surface_to_block_mapping(parts[0], touchingParts);

  // tookusa: Order is important for the incremental sideset updater ... surface to block mapping must be set first
  bulk.create_side_entities(sideSet, parts);
}

stk::mesh::Part* create_AA_mesh_with_sideset(stk::mesh::BulkData &bulk,
                                             SidesetDirection direction,
                                             ElementOrdering elemOrdering)
{
  stk::mesh::MetaData &meta = bulk.mesh_meta_data();
  stk::mesh::PartVector parts = create_sideset_parts(meta, std::vector<std::string>{"surface_1"});

  stk::unit_test_util::create_AA_mesh(bulk, elemOrdering);

  stk::unit_test_util::populate_AA_sideset(bulk, direction, elemOrdering, parts);

  return parts[0];
}

stk::mesh::Part* create_AA_mesh_with_sideset_and_field(stk::mesh::BulkData &bulk,
                                                       SidesetDirection direction,
                                                       ElementOrdering elemOrdering,
                                                       const std::string & fieldName)
{
  stk::mesh::MetaData &meta = bulk.mesh_meta_data();
  stk::mesh::PartVector parts = create_sideset_parts(meta, std::vector<std::string>{"surface_1"});

  const unsigned numberOfStates = 1;
  stk::mesh::Field<double> & field = meta.declare_field<double>(stk::topology::NODE_RANK, fieldName, numberOfStates);
  const double initValue = 123;
  stk::mesh::put_field_on_mesh(field, meta.universal_part(), &initValue);

  stk::unit_test_util::create_AA_mesh(bulk, elemOrdering);

  stk::unit_test_util::populate_AA_sideset(bulk, direction, elemOrdering, parts);

  return parts[0];
}

void populate_AB_sideset(stk::mesh::BulkData& bulk,
                         SidesetDirection direction,
                         ElementOrdering elemOrdering,
                         const stk::mesh::PartVector& parts)
{
  stk::mesh::MetaData &meta = bulk.mesh_meta_data();
  stk::mesh::EntityIdVector elem;
  std::vector<int> ordinal;

  stk::unit_test_util::populate_elem_sides(direction, elemOrdering, elem, ordinal);
  STK_ThrowRequire(elem.size() == ordinal.size());

  stk::mesh::SideSet &sideSet = bulk.create_sideset(*parts[0]);
  sideSet.set_accept_all_internal_non_coincident_entries(false);
  for(unsigned i=0; i<elem.size(); ++i) {
    const stk::mesh::Entity element = bulk.get_entity(stk::topology::ELEMENT_RANK, elem[i]);
    if (bulk.is_valid(element)) {
      sideSet.add({element, ordinal[i]});
    }
  }

  stk::mesh::Part* block_1 = meta.get_part("block_1");
  EXPECT_TRUE(block_1 != nullptr);
  stk::mesh::Part* block_2 = meta.get_part("block_2");
  EXPECT_TRUE(block_2 != nullptr);

  meta.set_part_id(*block_1, 1);
  meta.set_part_id(*block_2, 2);

  std::vector<const stk::mesh::Part*> touchingParts;

  if(direction == LEFT)
  {
    touchingParts = {block_1};
    meta.set_surface_to_block_mapping(parts[1], touchingParts);
  }
  else if(direction == RIGHT)
  {
    touchingParts = {block_2};
    meta.set_surface_to_block_mapping(parts[1], touchingParts);
  }
  else
  {
    touchingParts = {block_1};
    meta.set_surface_to_block_mapping(parts[1], touchingParts);

    touchingParts = {block_2};
    meta.set_surface_to_block_mapping(parts[2], touchingParts);
  }

  // tookusa: Order is important for the incremental sideset updater ... surface to block mapping must be set first
  bulk.create_side_entities(sideSet, parts);
}

stk::mesh::Part* create_AB_mesh_with_sideset(stk::mesh::BulkData &bulk,
                                             SidesetDirection direction,
                                             ElementOrdering elemOrdering)
{
  stk::mesh::MetaData &meta = bulk.mesh_meta_data();
  std::vector<std::string> sideSetNames;
  stk::unit_test_util::populate_sideset_names(direction, sideSetNames);

  stk::mesh::PartVector parts = create_sideset_parts(meta, sideSetNames);

  for(unsigned i = 0; i<parts.size(); ++i) {
    meta.set_part_id(*parts[i], 1);
  }

  for(unsigned i=1; i<parts.size(); ++i) {
    meta.declare_part_subset(*parts[0], *parts[i]);
  }

  stk::unit_test_util::create_AB_mesh(bulk, elemOrdering);

  stk::unit_test_util::populate_AB_sideset(bulk, direction, elemOrdering, parts);

  return parts[0];
}

stk::mesh::Part* create_AB_mesh_with_sideset_and_field(stk::mesh::BulkData &bulk,
                                                       SidesetDirection direction,
                                                       ElementOrdering elemOrdering,
                                                       const std::string & fieldName)
{
  stk::mesh::MetaData &meta = bulk.mesh_meta_data();
  std::vector<std::string> sideSetNames;
  stk::unit_test_util::populate_sideset_names(direction, sideSetNames);

  stk::mesh::PartVector parts = create_sideset_parts(meta, sideSetNames);

  for(unsigned i = 0; i<parts.size(); ++i) {
    meta.set_part_id(*parts[i], 1);
  }

  for(unsigned i=1; i<parts.size(); ++i) {
    meta.declare_part_subset(*parts[0], *parts[i]);
  }

  const unsigned numberOfStates = 1;
  stk::mesh::Field<double> & field = meta.declare_field<double>(stk::topology::NODE_RANK, fieldName, numberOfStates);
  const double initValue = 123;
  stk::mesh::put_field_on_mesh(field, meta.universal_part(), &initValue);

  stk::unit_test_util::create_AB_mesh(bulk, elemOrdering);

  stk::unit_test_util::populate_AB_sideset(bulk, direction, elemOrdering, parts);

  return parts[0];
}

stk::mesh::Part* create_AB_mesh_with_sideset_and_distribution_factors(stk::mesh::BulkData &bulk,
                                                                      SidesetDirection direction,
                                                                      ElementOrdering elemOrdering,
                                                                      const std::string & fieldName,
                                                                      const double initValue)
{
  stk::mesh::MetaData &meta = bulk.mesh_meta_data();
  std::vector<std::string> sideSetNames;
  stk::unit_test_util::populate_sideset_names(direction, sideSetNames);

  stk::mesh::PartVector parts = create_sideset_parts(meta, sideSetNames);

  for(unsigned i = 0; i<parts.size(); ++i) {
    meta.set_part_id(*parts[i], 1);
  }

  for(unsigned i=1; i<parts.size(); ++i) {
    meta.declare_part_subset(*parts[0], *parts[i]);
  }

  stk::mesh::Part* sidePart = parts[0];
  STK_ThrowRequire(nullptr != sidePart);
  const unsigned numberOfStates = 1;
  stk::mesh::Field<double> & ssField = meta.declare_field<double>(stk::topology::FACE_RANK, fieldName, numberOfStates);
  for (stk::mesh::Part* part : parts)
  {
    stk::io::set_distribution_factor_field(*part, ssField);
  }
  stk::topology sideTopo = sidePart->topology();
  STK_ThrowRequireMsg(sideTopo != stk::topology::INVALID_TOPOLOGY, "sidePart "<<sidePart->name()<<" has invalid topology.");
  unsigned numNodes = sideTopo.num_nodes();
  std::vector<double> initValVec(numNodes, initValue);
  stk::mesh::put_field_on_mesh(ssField, *sidePart, numNodes, initValVec.data());

  stk::unit_test_util::create_AB_mesh(bulk, elemOrdering);

  stk::unit_test_util::populate_AB_sideset(bulk, direction, elemOrdering, parts);

  return parts[0];
}

} // namespace simple_fields

}
}

