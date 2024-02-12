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

#include "stk_performance_tests/stk_mesh/multi_block.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/SkinBoundary.hpp"
#include "stk_mesh/base/ExodusTranslator.hpp"
#include "stk_io/IossBridge.hpp"

namespace stk {
namespace performance_tests {

void setup_multiple_blocks(stk::mesh::MetaData& meta, unsigned numBlocks)
{
  std::string block1Name("block_1");
  stk::mesh::Part& block1Part = meta.declare_part_with_topology(block1Name, stk::topology::HEX_8);
  meta.set_part_id(block1Part, 1);
  stk::io::put_io_part_attribute(block1Part);

  unsigned partId = 10;
  for (unsigned i = 1; i < numBlocks; i++) {
    std::string blockName = "block_" + std::to_string(partId);
    stk::mesh::Part& part = meta.declare_part_with_topology(blockName, stk::topology::HEX_8);
    meta.set_part_id(part, partId);
    stk::io::put_io_part_attribute(part);
    ++partId;
  }
}

void setup_elem_fields_on_blocks(stk::mesh::MetaData& meta, unsigned numFields, bool allFieldsSameSize)
{
  {
    std::vector<stk::mesh::FieldBase*> fields(numFields);
    for(unsigned f=0; f<numFields; ++f) {
      fields[f] = &meta.declare_field<double>(stk::topology::ELEM_RANK, std::string("elemField"+std::to_string(f)));
    }
  }

  unsigned scalarsPerEntity = 1;
  stk::mesh::PartVector elemBlocks;
  stk::mesh::fill_element_block_parts(meta, stk::topology::HEX_8, elemBlocks);
  for(stk::mesh::Part* elemBlock : elemBlocks) {
    for(stk::mesh::FieldBase* field : meta.get_fields(stk::topology::ELEM_RANK)) {
      stk::mesh::put_field_on_mesh(*field, *elemBlock, scalarsPerEntity, nullptr);
    }
    if(!allFieldsSameSize) {
      ++scalarsPerEntity;
    }
  }
}

std::string sideset_name_between_blocks(unsigned leftBlock, unsigned rightBlock)
{
  return "surfaceBetween_" + std::to_string(leftBlock) + "_and_" + std::to_string(rightBlock);
}

void setup_sidesets_between_blocks(stk::mesh::MetaData& meta)
{
  stk::mesh::PartVector elemBlocks;
  stk::mesh::fill_element_block_parts(meta, stk::topology::HEX_8, elemBlocks);
  const unsigned numBlocks = elemBlocks.size();
  unsigned partId = 10;
  for(unsigned i = 1; i < numBlocks; i++) {
    unsigned leftBlockId = elemBlocks[i-1]->id();
    unsigned rightBlockId = elemBlocks[i]->id();
    std::string sidesetName = sideset_name_between_blocks(leftBlockId, rightBlockId)+"_ss1";
    stk::mesh::Part& part = meta.declare_part(sidesetName, meta.side_rank());
    meta.set_part_id(part, partId++);

    std::string sidesetName2 = sideset_name_between_blocks(leftBlockId, rightBlockId)+"_ss2";
    stk::mesh::Part& part2 = meta.declare_part(sidesetName2, meta.side_rank());
    meta.set_part_id(part2, partId++);

    stk::mesh::Part& leftBlockPart = *elemBlocks[i-1];
    stk::mesh::Part& rightBlockPart = *elemBlocks[i];

    meta.set_surface_to_block_mapping(&part, {&leftBlockPart, &rightBlockPart});
    meta.set_surface_to_block_mapping(&part2, {&leftBlockPart, &rightBlockPart});

    stk::io::put_io_part_attribute(part);
    stk::io::put_io_part_attribute(part2);
  }
}

void setup_sidesets_for_blocks(stk::mesh::MetaData& meta)
{
  stk::mesh::PartVector elemBlocks;
  stk::mesh::fill_element_block_parts(meta, stk::topology::HEX_8, elemBlocks);
  const unsigned numBlocks = elemBlocks.size();
  for(unsigned i = 0; i < numBlocks; i++) {
    stk::mesh::Part& blockPart = *elemBlocks[i];
    unsigned partId = blockPart.id();

    std::string sidesetName = "surface_" + std::to_string(partId);
    stk::mesh::Part& part = meta.declare_part(sidesetName, meta.side_rank());
    meta.set_part_id(part, partId);

    meta.set_surface_to_block_mapping(&part, {&blockPart});

    stk::io::put_io_part_attribute(part);
  }
}

void move_elements_to_other_blocks(stk::mesh::BulkData& bulk, unsigned numElemsInDimX)
{
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stk::mesh::PartVector elemBlocks;
  stk::mesh::fill_element_block_parts(meta, stk::topology::HEX_8, elemBlocks);
  const unsigned numBlocks = elemBlocks.size();

  STK_ThrowRequireMsg(((numElemsInDimX % numBlocks) == 0),
                  "Number of blocks (" << numBlocks << ") must divide evenly into numElemsInDimX (" << numElemsInDimX << ")");

  stk::mesh::EntityVector elems;
  stk::mesh::Selector ownedHexes = meta.get_topology_root_part(stk::topology::HEX_8) &
                                   meta.locally_owned_part();
  stk::mesh::get_selected_entities(ownedHexes, bulk.buckets(stk::topology::ELEMENT_RANK), elems);
  const stk::mesh::Part* block1Part = elemBlocks[0];

  std::vector<stk::mesh::EntityVector> elemsToMove(elemBlocks.size());

  unsigned blockIdx = 0;
  for(stk::mesh::Entity elem : elems) {
    elemsToMove[blockIdx].push_back(elem);
    ++blockIdx;
    if (blockIdx == numBlocks) {
      blockIdx = 0;
    }
  }

  bulk.modification_begin();

  for(unsigned i = 1; i < numBlocks; i++) {
    stk::mesh::Part* newBlock = elemBlocks[i];

    bulk.change_entity_parts(elemsToMove[i], stk::mesh::ConstPartVector{newBlock}, stk::mesh::ConstPartVector{block1Part});
  }

  bulk.modification_end();
}

void move_elements_to_other_contiguous_blocks(stk::mesh::BulkData& bulk, unsigned numElemsInDimX)
{
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stk::mesh::PartVector elemBlocks;
  stk::mesh::fill_element_block_parts(meta, stk::topology::HEX_8, elemBlocks);
  const unsigned numBlocks = elemBlocks.size();

  STK_ThrowRequireMsg(((numElemsInDimX % numBlocks) == 0),
                  "Number of blocks (" << numBlocks << ") must divide evenly into numElemsInDimX (" << numElemsInDimX << ")");

  stk::mesh::EntityVector elems;
  stk::mesh::Selector ownedHexes = meta.get_topology_root_part(stk::topology::HEX_8) &
                                   meta.locally_owned_part();
  stk::mesh::get_selected_entities(ownedHexes, bulk.buckets(stk::topology::ELEMENT_RANK), elems);
  const stk::mesh::Part* block1Part = elemBlocks[0];

  std::vector<stk::mesh::EntityVector> elemsToMove(elemBlocks.size());
  unsigned numElemsInDimXPerBlock = numElemsInDimX / numBlocks;

  for(stk::mesh::Entity elem : elems) {
    stk::mesh::EntityId id = bulk.identifier(elem);
    unsigned blockIdx = ((id - 1)%numElemsInDimX)/numElemsInDimXPerBlock;

    STK_ThrowRequire(blockIdx < numBlocks);
    elemsToMove[blockIdx].push_back(elem);
  }

  bulk.modification_begin();

  for(unsigned i = 1; i < numBlocks; i++) {
    stk::mesh::Part* newBlock = elemBlocks[i];

    bulk.change_entity_parts(elemsToMove[i], stk::mesh::ConstPartVector{newBlock}, stk::mesh::ConstPartVector{block1Part});
  }

  bulk.modification_end();
}

void fill_sideset(stk::mesh::BulkData& bulk,
                  stk::mesh::Part& partToPutSidesInto,
                  const stk::mesh::Selector& blockSelector)
{
  stk::mesh::create_interior_block_boundary_sides(bulk, blockSelector,
                                                 stk::mesh::PartVector{&partToPutSidesInto});
  stk::mesh::SideSet* sideset = nullptr;
  if (bulk.does_sideset_exist(partToPutSidesInto)) {
    sideset = &bulk.get_sideset(partToPutSidesInto);
  }
  else {
    sideset = &bulk.create_sideset(partToPutSidesInto);
  }

  stk::mesh::EntityVector sides;
  stk::mesh::get_selected_entities(partToPutSidesInto, bulk.buckets(stk::topology::FACE_RANK), sides);

  for(stk::mesh::Entity side : sides) {
    const unsigned numElems = bulk.num_elements(side);
    const stk::mesh::Entity* elems = bulk.begin_elements(side);
    const stk::mesh::ConnectivityOrdinal* ords = bulk.begin_element_ordinals(side);
    for(unsigned i=0; i<numElems; ++i) {
      sideset->add({elems[i], ords[i]});
    } 
  }
}

void fill_sidesets_between_blocks(stk::mesh::BulkData& bulk)
{
  const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stk::mesh::PartVector elemBlocks;
  stk::mesh::fill_element_block_parts(meta, stk::topology::HEX_8, elemBlocks);
  const unsigned numBlocks = elemBlocks.size();
  for(unsigned i = 1; i < numBlocks; i++) {
    stk::mesh::Part& leftBlockPart = *elemBlocks[i-1];
    stk::mesh::Part& rightBlockPart = *elemBlocks[i];
    std::string sidesetName = sideset_name_between_blocks(leftBlockPart.id(), rightBlockPart.id())+"_ss1";
    stk::mesh::Part& sidesetPart = *meta.get_part(sidesetName);

    fill_sideset(bulk, sidesetPart, stk::mesh::Selector(leftBlockPart|rightBlockPart));

    std::string sidesetName2 = sideset_name_between_blocks(leftBlockPart.id(), rightBlockPart.id())+"_ss2";
    stk::mesh::Part& sidesetPart2 = *meta.get_part(sidesetName2);

    fill_sideset(bulk, sidesetPart2, stk::mesh::Selector(leftBlockPart|rightBlockPart));
  }
}

}}
