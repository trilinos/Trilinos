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

#include "stk_util/stk_config.h"            // for STK_HAS_MPI
#ifdef STK_HAS_MPI

#include "gtest/gtest.h"
#include "stk_util/parallel/MPITagManager.hpp"
#include <thread>

namespace {

MPI_Comm copy_comm(MPI_Comm comm)
{
  MPI_Comm newcomm;
  MPI_Comm_dup(comm, &newcomm);
  return newcomm;
}


class MPITagManagerForTesting : public stk::MPITagManager
{
  public:
    MPITagManagerForTesting(int deletionGroupSize, int delayCount) :
      MPITagManager(deletionGroupSize, delayCount)
    {}
};


std::shared_ptr<stk::MPITagManager> make_tag_manager_for_testing(int deletionGroupSize, int delayCount)
{
  return std::make_shared<MPITagManagerForTesting>(deletionGroupSize, delayCount);
}


class TagManagerTester : public ::testing::Test
{
  protected:
    TagManagerTester() :
      tag_manager_ptr(make_tag_manager_for_testing(1, 0)),
      tag_manager(*tag_manager_ptr)

    {
      comm1 = copy_comm(MPI_COMM_WORLD);
      comm2 = copy_comm(MPI_COMM_WORLD);
    }

    ~TagManagerTester()
    {
      MPI_Comm_free(&comm1);
      MPI_Comm_free(&comm2);
    }

    MPI_Comm comm1;
    MPI_Comm comm2;
    std::shared_ptr<stk::MPITagManager> tag_manager_ptr;
    stk::MPITagManager& tag_manager;
};

void sleep(int milliseconds)
{
  std::this_thread::sleep_for(std::chrono::milliseconds(milliseconds));
}

}  // namespace

TEST_F(TagManagerTester, TagComm)
{
  stk::MPITag tag1 = tag_manager.get_tag(comm1);
  stk::MPITag tag2 = tag_manager.get_tag(comm2);
  stk::MPITag tag3 = tag_manager.get_tag(comm1);
  stk::MPITag tag4 = tag_manager.get_tag(comm2);

  EXPECT_EQ(tag1.get_comm(), comm1);
  EXPECT_EQ(tag2.get_comm(), comm2);
  EXPECT_EQ(tag3.get_comm(), comm1);
  EXPECT_EQ(tag4.get_comm(), comm2);

  std::cout << "about to exit function" << std::endl;
}


TEST_F(TagManagerTester, TagUniqueness)
{
  std::set<int> unique_tags;
  stk::MPITag tag1 = tag_manager.get_tag(comm1, 1); unique_tags.insert(tag1);
  stk::MPITag tag2 = tag_manager.get_tag(comm1, 1); unique_tags.insert(tag2);
  stk::MPITag tag3 = tag_manager.get_tag(comm1, 1); unique_tags.insert(tag3);
  stk::MPITag tag4 = tag_manager.get_tag(comm1, 1); unique_tags.insert(tag4);

  EXPECT_EQ(unique_tags.size(), 4u);
}


TEST_F(TagManagerTester, TagUniquenessWithHoles)
{
  stk::MPITag tag1 = tag_manager.get_tag(comm1, 1);
  stk::MPITag tag2 = tag_manager.get_tag(comm1, 4);
  stk::MPITag tag3 = tag_manager.get_tag(comm1, 1);
  stk::MPITag tag4 = tag_manager.get_tag(comm1, 1);
  stk::MPITag tag5 = tag_manager.get_tag(comm1, 5);

  std::vector<stk::MPITag> tags{tag1, tag2, tag3, tag4, tag5};

  for (unsigned int i=0; i < tags.size(); ++i) {
    for (unsigned int j=0; j < tags.size(); ++j) {
      if (i == j) {
        EXPECT_EQ(tags[i], tags[j]);
      } else {
        EXPECT_NE(tags[i], tags[j]);
      }
    }
  }
}


TEST_F(TagManagerTester, TagReuse)
{
  {
    stk::MPITag tag1 = tag_manager.get_tag(comm1, 1);
    stk::MPITag tag2 = tag_manager.get_tag(comm1, 2);

    EXPECT_EQ(static_cast<int>(tag1), 1);
    EXPECT_EQ(static_cast<int>(tag2), 2);
  }  // tags get freed

  {
    stk::MPITag tag1 = tag_manager.get_tag(comm1, 1);
    stk::MPITag tag2 = tag_manager.get_tag(comm1, 2);

    EXPECT_EQ(static_cast<int>(tag1), 1);
    EXPECT_EQ(static_cast<int>(tag2), 2);
  }
}

TEST_F(TagManagerTester, TagReuseMultipleComms)
{
  {
    stk::MPITag tag11 = tag_manager.get_tag(comm1, 1);
    stk::MPITag tag12 = tag_manager.get_tag(comm1, 2);

    stk::MPITag tag21 = tag_manager.get_tag(comm2, 1);
    stk::MPITag tag22 = tag_manager.get_tag(comm2, 2);

    EXPECT_EQ(static_cast<int>(tag11), 1);
    EXPECT_EQ(static_cast<int>(tag12), 2);
    EXPECT_EQ(static_cast<int>(tag21), 1);
    EXPECT_EQ(static_cast<int>(tag22), 2);
  }  // tags get freed

  {
    stk::MPITag tag11 = tag_manager.get_tag(comm1, 1);
    stk::MPITag tag12 = tag_manager.get_tag(comm1, 2);

    stk::MPITag tag21 = tag_manager.get_tag(comm2, 1);
    stk::MPITag tag22 = tag_manager.get_tag(comm2, 2);

    EXPECT_EQ(static_cast<int>(tag11), 1);
    EXPECT_EQ(static_cast<int>(tag12), 2);
    EXPECT_EQ(static_cast<int>(tag21), 1);
    EXPECT_EQ(static_cast<int>(tag22), 2);
  }

}


TEST_F(TagManagerTester, TagCopying)
{
  {
    stk::MPITag tag1 = tag_manager.get_tag(comm1, 1);
    stk::MPITag tag2 = tag1;
    stk::MPITag tag3 = tag2;
    EXPECT_EQ(static_cast<int>(tag1), 1);
    EXPECT_EQ(static_cast<int>(tag2), 1);
    EXPECT_EQ(static_cast<int>(tag3), 1);
  }
}


TEST_F(TagManagerTester, CommFree)
{
  {
    MPI_Comm  comm3 = copy_comm(comm1);
    stk::MPITag tag1 = tag_manager.get_tag(comm3, 1);
    stk::MPITag tag2 = tag_manager.get_tag(comm3, 2);
    stk::MPITag tag3 = tag_manager.get_tag(comm3, 3);

    MPI_Comm_free(&comm3);
  }
}


TEST_F(TagManagerTester, CommFree_ResetTagRange)
{
  {
    MPI_Comm comm3 = copy_comm(comm1);
    stk::MPITag tag1 = tag_manager.get_tag(comm3, 1);
    stk::MPITag tag2 = tag_manager.get_tag(comm3, 2);
    stk::MPITag tag3 = tag_manager.get_tag(comm3, 3);

    MPI_Comm_free(&comm3);
    MPI_Comm_dup(comm1, &comm3);

    stk::MPITag tag4 = tag_manager.get_tag(comm3, 1);
    EXPECT_EQ(static_cast<int>(tag4), 1);
    MPI_Comm_free(&comm3);
  }
}

TEST_F(TagManagerTester, CommFree_ResetTagRangeNewScope)
{
  {
    MPI_Comm comm3 = copy_comm(comm1);
    stk::MPITag tag1 = tag_manager.get_tag(comm3, 1);
    stk::MPITag tag2 = tag_manager.get_tag(comm3, 2);
    stk::MPITag tag3 = tag_manager.get_tag(comm3, 3);

    MPI_Comm_free(&comm3);
  }

  {
    MPI_Comm comm3 = copy_comm(comm1);
    stk::MPITag tag4 = tag_manager.get_tag(comm3, 1);
    EXPECT_EQ(static_cast<int>(tag4), 1);
    MPI_Comm_free(&comm3);
  }
}

TEST_F(TagManagerTester, CommFreeMultiple)
{
  MPI_Comm comm3 = copy_comm(comm1);
  MPI_Comm comm4 = copy_comm(comm1);
  MPI_Comm comm5 = copy_comm(comm1);
  stk::MPITag tag11 = tag_manager.get_tag(comm1);
  stk::MPITag tag21 = tag_manager.get_tag(comm2);
  stk::MPITag tag31 = tag_manager.get_tag(comm3);
  stk::MPITag tag41 = tag_manager.get_tag(comm4);
  stk::MPITag tag51 = tag_manager.get_tag(comm5);

  MPI_Comm_free(&comm3);
  MPI_Comm_free(&comm4);

  MPI_Comm_dup(comm1, &comm3);
  MPI_Comm_dup(comm1, &comm4);

  stk::MPITag tag32 = tag_manager.get_tag(comm3);
  stk::MPITag tag42 = tag_manager.get_tag(comm4);

  std::vector<stk::MPITag> tags{tag11, tag21, tag32, tag42, tag51};

  for (unsigned int i=0; i < tags.size(); ++i) {
    for (unsigned int j=0; j < tags.size(); ++j) {
      if (i == j) {
        EXPECT_EQ(tags[i], tags[j]);
      } else {
        EXPECT_NE(tags[i], tags[j]);
      }
    }
  }
  EXPECT_NE(tag31, tag31);
  EXPECT_NE(tag41, tag41);

  for (unsigned int i=0; i < tags.size(); ++i) {
    EXPECT_NE(tag31, tags[i]);
    EXPECT_NE(tag41, tags[i]);
  }

  MPI_Comm_free(&comm3);
  MPI_Comm_free(&comm4);
  MPI_Comm_free(&comm5);
}


TEST(TagManager, NonSortedCommOrder_DoesntHang)
{
  MPI_Comm comm1=MPI_COMM_WORLD, comm2, comm3, comm4;
  int myRank;
  MPI_Comm_rank(comm1, &myRank);

  int color = myRank == 0 ? MPI_UNDEFINED : 1;
  MPI_Comm_split(comm1, color, 0, &comm2);
  MPI_Comm_dup(comm1, &comm3);
  MPI_Comm_dup(comm1, &comm4);

  {
    int deletionGroupSize = 1, delayCount = 100;
    auto tagManager = make_tag_manager_for_testing(deletionGroupSize, delayCount);
    tagManager->get_tag(comm1);
    if (color != MPI_UNDEFINED) {
      tagManager->get_tag(comm2);
    }

    tagManager->get_tag(comm3);
    if (color != MPI_UNDEFINED) {
      MPI_Comm_free(&comm2);
    }

    tagManager->get_tag(comm4);
  }
  MPI_Comm_free(&comm3);
  MPI_Comm_free(&comm4);
}

namespace {
void check_tags_same_on_all_procs(MPI_Comm comm, const std::vector<int>& vals)
{
  int myRank, commSize;
  MPI_Comm_rank(comm, &myRank);
  MPI_Comm_size(comm, &commSize);

  int root = 0;
  int niters = vals.size();
  std::vector<int> vals_all(niters * commSize);
  MPI_Gather(vals.data(), niters, MPI_INT, vals_all.data(), niters, MPI_INT, root, comm);
  if (myRank == root) {
    for (int i=0; i < niters; ++i) {
      for (int rank=0; rank < commSize; ++rank) {
        EXPECT_EQ(vals[i], vals_all[rank * niters + i]);
      }
    }
  }
}
}

TEST(TagManager, MisalignedProcs)
{
  MPI_Comm comm1 = MPI_COMM_WORLD;
  int myRank, commSize;
  MPI_Comm_rank(comm1, &myRank);
  MPI_Comm_size(comm1, &commSize);

  int deletionGroupSize=16, delayCount=8;
  auto tagManager = make_tag_manager_for_testing(deletionGroupSize, delayCount);
  std::vector<int> tag_vals;

  int niters = 50;
  for (int i=0; i < niters; ++i) {
    tag_vals.push_back(tagManager->get_tag(comm1));
    sleep(10*myRank);
  }

  check_tags_same_on_all_procs(comm1, tag_vals);
}

TEST(TagManager, MisalignedProcsMultipleComms)
{
  MPI_Comm comm1 = copy_comm(MPI_COMM_WORLD);
  int myRank1, commSize1;
  MPI_Comm_rank(comm1, &myRank1);
  MPI_Comm_size(comm1, &commSize1);

  MPI_Comm comm2;
  int color = 0;
  if (commSize1 >= 2) {
    color = myRank1 / (commSize1 / 2);
  }
  MPI_Comm_split(comm1, color, 0, &comm2);

  int deletionGroupSize=16, delayCount=8;
  auto tagManager = make_tag_manager_for_testing(deletionGroupSize, delayCount);
  std::vector<int> tag_vals1, tag_vals2;

  int niters = 50;
  for (int i=0; i < niters; ++i) {
    tag_vals1.push_back(tagManager->get_tag(comm1));
    tag_vals2.push_back(tagManager->get_tag(comm2));
    sleep(10*myRank1);
  }

  check_tags_same_on_all_procs(comm1, tag_vals1);
  check_tags_same_on_all_procs(comm2, tag_vals2);

  MPI_Comm_free(&comm1);
  MPI_Comm_free(&comm2);
}


TEST(TagManager, MisalignedProcsFreePeriodically)
{
  MPI_Comm comm1 = copy_comm(MPI_COMM_WORLD);
  int myRank, commSize;
  MPI_Comm_rank(comm1, &myRank);
  MPI_Comm_size(comm1, &commSize);

  int deletionGroupSize=16, delayCount=8;
  auto tagManager = make_tag_manager_for_testing(deletionGroupSize, delayCount);
  std::vector<stk::MPITag> tags;
  std::vector<std::vector<int>> tag_vals;

  int niters = 50, free_every = 22;
  for (int i=0; i < niters; ++i) {
    tags.push_back(tagManager->get_tag(comm1));
    sleep(10*myRank);

    if (i % free_every == 0) {
      tag_vals.emplace_back(tags.begin(), tags.end());
      tags.clear();
    }
  }

  for (auto& tag_vals_group : tag_vals) {
    check_tags_same_on_all_procs(comm1, tag_vals_group);
  }

  MPI_Comm_free(&comm1);
}



TEST(TagManager, StaticInstance)
{
  stk::get_mpi_tag_manager().get_tag(MPI_COMM_WORLD);

  // if this doesn't segfault when MPI_Finalize is called, then the test passes
}

#endif
