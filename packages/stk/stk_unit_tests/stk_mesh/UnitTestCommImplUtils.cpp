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
#include <stk_mesh/baseImpl/CommImplUtils.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stddef.h>
#include <vector>

TEST ( CommImplUtils, fill_procs_shared_then_ghosted__empty )
{
  std::vector<stk::mesh::EntityCommInfo> info = {};
  stk::mesh::PairIterEntityComm pec(info.data(), info.data()+info.size());

  std::vector<int> commProcs;
  unsigned numShared = stk::mesh::impl::fill_procs_shared_then_ghosted(pec, commProcs);
  EXPECT_EQ(0u, numShared);
  EXPECT_TRUE(commProcs.empty());
}

TEST ( CommImplUtils, fill_procs_shared_then_ghosted__only_shared )
{
  std::vector<stk::mesh::EntityCommInfo> info = {
    stk::mesh::EntityCommInfo(stk::mesh::BulkData::SHARED, 0),
    stk::mesh::EntityCommInfo(stk::mesh::BulkData::SHARED, 1),
  };
  stk::mesh::PairIterEntityComm pec(info.data(), info.data()+info.size());

  std::vector<int> commProcs;
  unsigned numShared = stk::mesh::impl::fill_procs_shared_then_ghosted(pec, commProcs);
  EXPECT_EQ(2u, numShared);
  EXPECT_EQ(2u, commProcs.size());
}

TEST ( CommImplUtils, fill_procs_shared_then_ghosted__both_shared_and_ghosted )
{
  std::vector<stk::mesh::EntityCommInfo> info = {
    stk::mesh::EntityCommInfo(stk::mesh::BulkData::SHARED, 0),
    stk::mesh::EntityCommInfo(stk::mesh::BulkData::SHARED, 2),
    stk::mesh::EntityCommInfo(stk::mesh::BulkData::AURA, 1),
    stk::mesh::EntityCommInfo(stk::mesh::BulkData::AURA, 3),
  };
  stk::mesh::PairIterEntityComm pec(info.data(), info.data()+info.size());

  std::vector<int> commProcs;
  unsigned numShared = stk::mesh::impl::fill_procs_shared_then_ghosted(pec, commProcs);
  EXPECT_EQ(2u, numShared);
  EXPECT_EQ(4u, commProcs.size());
}

TEST ( CommImplUtils, fill_procs_shared_then_ghosted__same_procs_shared_and_ghosted )
{
  std::vector<stk::mesh::EntityCommInfo> info = {
    stk::mesh::EntityCommInfo(stk::mesh::BulkData::SHARED, 0),
    stk::mesh::EntityCommInfo(stk::mesh::BulkData::SHARED, 2),
    stk::mesh::EntityCommInfo(stk::mesh::BulkData::AURA, 0),
    stk::mesh::EntityCommInfo(stk::mesh::BulkData::AURA, 2),
  };
  stk::mesh::PairIterEntityComm pec(info.data(), info.data()+info.size());

  std::vector<int> commProcs;
  unsigned numShared = stk::mesh::impl::fill_procs_shared_then_ghosted(pec, commProcs);
  EXPECT_EQ(2u, numShared);
  EXPECT_EQ(2u, commProcs.size());
}

TEST ( CommImplUtils, fill_procs_shared_then_ghosted__2_ghostings_same_proc )
{
  std::vector<stk::mesh::EntityCommInfo> info = {
    stk::mesh::EntityCommInfo(stk::mesh::BulkData::SHARED, 0),
    stk::mesh::EntityCommInfo(stk::mesh::BulkData::SHARED, 2),
    stk::mesh::EntityCommInfo(stk::mesh::BulkData::AURA, 1),
    stk::mesh::EntityCommInfo(                        2, 1),
  };
  stk::mesh::PairIterEntityComm pec(info.data(), info.data()+info.size());

  std::vector<int> commProcs;
  unsigned numShared = stk::mesh::impl::fill_procs_shared_then_ghosted(pec, commProcs);
  EXPECT_EQ(2u, numShared);
  EXPECT_EQ(3u, commProcs.size());
  EXPECT_EQ(1, commProcs[2]);
}

TEST ( CommImplUtils, fill_procs_shared_then_ghosted__both_shared_and_ghosted_procs_unsorted)
{
  std::vector<stk::mesh::EntityCommInfo> info = {
    stk::mesh::EntityCommInfo(stk::mesh::BulkData::SHARED, 4),
    stk::mesh::EntityCommInfo(stk::mesh::BulkData::SHARED, 3),
    stk::mesh::EntityCommInfo(stk::mesh::BulkData::AURA, 1),
    stk::mesh::EntityCommInfo(stk::mesh::BulkData::AURA, 2),
  };
  stk::mesh::PairIterEntityComm pec(info.data(), info.data()+info.size());

  std::vector<int> commProcs;
  unsigned numShared = stk::mesh::impl::fill_procs_shared_then_ghosted(pec, commProcs);
  EXPECT_EQ(2u, numShared);
  EXPECT_EQ(4u, commProcs.size());
  EXPECT_EQ(1, commProcs[2]);
  EXPECT_EQ(2, commProcs[3]);
}

TEST ( CommImplUtils, fill_procs_shared_then_ghosted__proc_shared_and_ghost )
{
  std::vector<stk::mesh::EntityCommInfo> info = {
    stk::mesh::EntityCommInfo(stk::mesh::BulkData::SHARED, 0),
    stk::mesh::EntityCommInfo(stk::mesh::BulkData::SHARED, 1),
    stk::mesh::EntityCommInfo(stk::mesh::BulkData::SHARED, 2),
    stk::mesh::EntityCommInfo(stk::mesh::BulkData::AURA, 1),
    stk::mesh::EntityCommInfo(stk::mesh::BulkData::AURA, 3),
  };
  stk::mesh::PairIterEntityComm pec(info.data(), info.data()+info.size());

  std::vector<int> commProcs;
  unsigned numShared = stk::mesh::impl::fill_procs_shared_then_ghosted(pec, commProcs);
  EXPECT_EQ(3u, numShared);
  EXPECT_EQ(4u, commProcs.size());
  EXPECT_EQ(3, commProcs[3]);
}
