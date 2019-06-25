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
#include <gtest/gtest.h>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_tools/block_extractor/ExtractBlocks.hpp>
#include "stk_mesh/base/Entity.hpp"

namespace
{

class MeshWithTwoBlocks : public stk::unit_test_util::MeshFixture
{
protected:
    void switch_half_mesh_to_part(stk::mesh::Part &addPart, stk::mesh::Part &removePart)
    {
        get_bulk().modification_begin();
        switch_elem_to_part(1, addPart, removePart);
        switch_elem_to_part(2, addPart, removePart);
        get_bulk().modification_end();
    }

    void switch_elem_to_part(stk::mesh::EntityId elemId, stk::mesh::Part &addPart, stk::mesh::Part &removePart)
    {
        stk::mesh::Entity elem = get_bulk().get_entity(stk::topology::ELEM_RANK, elemId);
        if(get_bulk().is_valid(elem) && get_bulk().bucket(elem).owned())
            get_bulk().change_entity_parts(elem, stk::mesh::ConstPartVector{&addPart}, stk::mesh::ConstPartVector{&removePart});
    }

    void expect_num_elems_in_part(stk::mesh::BulkData &bulk, size_t expectedNum, stk::mesh::Part& block)
    {
        stk::mesh::Selector selector(block);
        std::vector<size_t> counts;
        stk::mesh::comm_mesh_counts(bulk, counts, &selector);
        EXPECT_EQ(expectedNum, counts[stk::topology::ELEM_RANK]);
    }
};

TEST_F(MeshWithTwoBlocks, extractBlock2_onlyHaveBlock2)
{
    stk::mesh::Part &block2 = get_meta().declare_part_with_topology("block_2", stk::topology::HEX_8);
    setup_mesh("generated:1x1x4", stk::mesh::BulkData::AUTO_AURA);
    stk::mesh::Part &block1 = *get_meta().get_part("block_1");
    switch_half_mesh_to_part(block2, block1);

    expect_num_elems_in_part(get_bulk(), 2, block1);
    expect_num_elems_in_part(get_bulk(), 2, block2);

    stk::mesh::MetaData newMeta;
    stk::mesh::BulkData newBulk(newMeta, get_comm());
    stk::tools::extract_blocks(get_bulk(), newBulk, {"block_2"});

    expect_num_elems_in_part(newBulk, 0, block1);
    expect_num_elems_in_part(newBulk, 2, block2);
}

class MeshWithOneBlock : public stk::unit_test_util::MeshFixture
{

};

TEST_F(MeshWithOneBlock, extractBlock2_throws)
{
    setup_mesh("generated:1x1x4", stk::mesh::BulkData::AUTO_AURA);
    stk::mesh::MetaData newMeta;
    stk::mesh::BulkData newBulk(newMeta, get_comm());
    EXPECT_THROW(stk::tools::extract_blocks(get_bulk(), newBulk, {"block_2"}), std::exception);
}

}
