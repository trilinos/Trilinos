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
#include "gtest/gtest.h"
#include "stk_util/parallel/MPITagManager.hpp"

namespace {

class TagManagerTester : public ::testing::Test
{
  protected:
    TagManagerTester() :
      tag_manager(stk::get_mpi_tag_manager())
    {
      comm1 = MPI_COMM_WORLD;
      MPI_Comm_dup(comm1, &comm2);
    }

    MPI_Comm comm1;
    MPI_Comm comm2;
    stk::MPITagManager& tag_manager;
};

}  // namespace

TEST_F(TagManagerTester, TagUniqueness)
{
  auto tag1 = tag_manager.get_tag(comm1, 1);
  auto tag2 = tag_manager.get_tag(comm1, 1);
  auto tag3 = tag_manager.get_tag(comm1, 1);
  auto tag4 = tag_manager.get_tag(comm1, 1);

  EXPECT_NE(tag1, tag2);
  EXPECT_NE(tag1, tag3);
  EXPECT_NE(tag1, tag4);

  EXPECT_NE(tag2, tag1);
  EXPECT_NE(tag2, tag3);
  EXPECT_NE(tag2, tag4);

  EXPECT_NE(tag3, tag1);
  EXPECT_NE(tag3, tag2);
  EXPECT_NE(tag3, tag4);

  EXPECT_NE(tag3, tag1);
  EXPECT_NE(tag3, tag2);
  EXPECT_NE(tag3, tag4);
}


TEST_F(TagManagerTester, TagUniquenessWithHoles)
{
  auto tag1 = tag_manager.get_tag(comm1, 1);
  auto tag2 = tag_manager.get_tag(comm1, 4);
  auto tag3 = tag_manager.get_tag(comm1, 1);
  auto tag4 = tag_manager.get_tag(comm1, 1);
  auto tag5 = tag_manager.get_tag(comm1, 5);

  std::vector<stk::MPITag> tags{tag1, tag2, tag3, tag4, tag5};

  for (unsigned int i=0; i < tags.size(); ++i)
    for (unsigned int j=0; j < tags.size(); ++j) {
      if (i == j) {
        EXPECT_EQ(tags[i], tags[j]);
      } else {
        EXPECT_NE(tags[i], tags[j]);
      }
    }
}


TEST_F(TagManagerTester, TagReuse)
{
  {
    auto tag1 = tag_manager.get_tag(comm1, 1);
    auto tag2 = tag_manager.get_tag(comm1, 2);
    EXPECT_EQ(static_cast<int>(tag1), 1);
    EXPECT_EQ(static_cast<int>(tag2), 2);
  }

  // the tags should have be freed at the end of the block
  {
    auto tag1 = tag_manager.get_tag(comm1, 1);
    auto tag2 = tag_manager.get_tag(comm1, 2);
    EXPECT_EQ(static_cast<int>(tag1), 1);
    EXPECT_EQ(static_cast<int>(tag2), 2);
  }
}


TEST_F(TagManagerTester, TagCopying)
{
  {
    auto tag1 = tag_manager.get_tag(comm1, 1);
    auto tag2 = tag1;
    auto tag3 = tag2;
    EXPECT_EQ(static_cast<int>(tag1), 1);
    EXPECT_EQ(static_cast<int>(tag2), 1);
    EXPECT_EQ(static_cast<int>(tag2), 1);
  } // if a tag was double freed, an error will be thrown here
}


TEST_F(TagManagerTester, CommFree)
{
  {
    auto tag1 = tag_manager.get_tag(comm2, 1);
    auto tag2 = tag_manager.get_tag(comm2, 2);
    auto tag3 = tag_manager.get_tag(comm2, 3);

    MPI_Comm_free(&comm2);
  } // if a tag was double freed, an error will be thrown here
}


TEST_F(TagManagerTester, CommFree2)
{
  {
    auto tag1 = tag_manager.get_tag(comm2, 1);
    auto tag2 = tag_manager.get_tag(comm2, 2);
    auto tag3 = tag_manager.get_tag(comm2, 3);

    MPI_Comm_free(&comm2);
    MPI_Comm_dup(comm1, &comm2);
    
    auto tag4 = tag_manager.get_tag(comm2, 1);
    EXPECT_EQ(static_cast<int>(tag4), 1);
  } // if a tag was double freed, an error will be thrown here
}

TEST_F(TagManagerTester, CommFreeMultiple)
{
  MPI_Comm comm3, comm4, comm5;

  MPI_Comm_dup(comm1, &comm3);
  MPI_Comm_dup(comm1, &comm4);
  MPI_Comm_dup(comm1, &comm5);

  auto tag11 = tag_manager.get_tag(comm1);
  auto tag21 = tag_manager.get_tag(comm2);
  auto tag31 = tag_manager.get_tag(comm3);
  auto tag41 = tag_manager.get_tag(comm4);
  auto tag51 = tag_manager.get_tag(comm5);

  MPI_Comm_free(&comm3);
  MPI_Comm_free(&comm4);

  MPI_Comm_dup(comm1, &comm3);
  MPI_Comm_dup(comm1, &comm4);

  auto tag32 = tag_manager.get_tag(comm3);
  auto tag42 = tag_manager.get_tag(comm4);

  std::vector<stk::MPITag> tags{tag11, tag21, tag32, tag42, tag51};

  for (unsigned int i=0; i < tags.size(); ++i)
    for (unsigned int j=0; j < tags.size(); ++j) {
      if (i == j) {
        EXPECT_EQ(tags[i], tags[j]);
      } else {
        EXPECT_NE(tags[i], tags[j]);
      }
    }

  EXPECT_NE(tag31, tag31);
  EXPECT_NE(tag41, tag41);
  for (unsigned int i=0; i < tags.size(); ++i) {
    EXPECT_NE(tag31, tags[i]);
    EXPECT_NE(tag41, tags[i]);
  }
}