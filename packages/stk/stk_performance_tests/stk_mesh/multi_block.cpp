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
#include "stk_io/IossBridge.hpp"

namespace stk {
namespace performance_tests {

void setup_multiple_blocks(stk::mesh::MetaData& meta, unsigned numBlocks)
{
  for(unsigned i = 1; i < numBlocks; i++) {
    std::string blockName = "block_" + std::to_string(i+1);
    stk::mesh::Part& part = meta.declare_part_with_topology(blockName, stk::topology::HEX_8);
    stk::io::put_io_part_attribute(part);
  }
}

void move_elements_to_other_blocks(stk::mesh::BulkData& bulk, unsigned numElemsPerDim, unsigned numBlocks)
{
  ThrowRequireMsg(numBlocks <= numElemsPerDim,
                  "Number of blocks (" << numBlocks << ") cannot be greater than numElemsPerDim (" << numElemsPerDim << ")");
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stk::mesh::EntityVector elems;
  stk::mesh::EntityVector elemsToMove;
  unsigned numElemsPerBlock = (numElemsPerDim * numElemsPerDim) / bulk.parallel_size();
  const stk::mesh::Part* block1Part = meta.get_part("block_1");

  stk::mesh::get_selected_entities(*block1Part, bulk.buckets(stk::topology::ELEMENT_RANK), elems);

  bulk.modification_begin();

  for(unsigned i = 1; i < numBlocks; i++) {
    std::string blockName = "block_" + std::to_string(i+1);

    stk::mesh::Part* newBlock = meta.get_part(blockName);

    elemsToMove.clear();

    unsigned startingIndex = i*numElemsPerBlock;

    for(unsigned j = startingIndex; j < (startingIndex+numElemsPerBlock); j++) {
      elemsToMove.push_back(elems[j]);
    }
    bulk.change_entity_parts(elemsToMove, stk::mesh::ConstPartVector{newBlock}, stk::mesh::ConstPartVector{block1Part});
  }
  bulk.modification_end();
}

}}
