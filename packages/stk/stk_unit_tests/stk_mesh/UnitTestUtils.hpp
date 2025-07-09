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
#ifndef _UnitTestUtils_hpp_
#define _UnitTestUtils_hpp_

#include <gtest/gtest.h>
#include <stddef.h>  // for size_t
#include <stdlib.h>  // for exit

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <iostream>
#include <vector>

namespace stk { namespace mesh { namespace unit_test {

inline
void proc0_ghost_node_to_proc1_and_proc2(stk::mesh::BulkData& bulk,
                                         stk::mesh::EntityId nodeID,
                                         const std::string& ghostingName = "myCustomGhosting")
{
  bulk.modification_begin();

  stk::mesh::Ghosting& myGhosting = bulk.create_ghosting(ghostingName);
  std::vector<stk::mesh::EntityProc> nodesToGhost;

  stk::mesh::Entity node1 = bulk.get_entity(stk::topology::NODE_RANK, nodeID);
  if (bulk.parallel_rank() == 0) {
    EXPECT_TRUE(bulk.is_valid(node1));
    EXPECT_TRUE(bulk.bucket(node1).owned());
    nodesToGhost.push_back(stk::mesh::EntityProc(node1, 1));
    nodesToGhost.push_back(stk::mesh::EntityProc(node1, 2));
  }

  bulk.change_ghosting(myGhosting, nodesToGhost);

  bulk.modification_end();
}

inline
void proc0_ghost_node5_to_proc2(stk::mesh::BulkData& bulk,
                                const std::string& ghostingName = "myCustomGhosting")
{
  bulk.modification_begin();

  stk::mesh::Ghosting& myGhosting = bulk.create_ghosting(ghostingName);
  std::vector<stk::mesh::EntityProc> nodesToGhost;

  stk::mesh::Entity node5 = bulk.get_entity(stk::topology::NODE_RANK, 5);
  if (bulk.parallel_rank() == 0) {
    EXPECT_TRUE(bulk.is_valid(node5));
    EXPECT_TRUE(bulk.bucket(node5).owned());
    nodesToGhost.push_back(stk::mesh::EntityProc(node5, 2));
  }

  bulk.change_ghosting(myGhosting, nodesToGhost);

  bulk.modification_end();
}

} } } // namespace stk mesh unit_test

#endif
