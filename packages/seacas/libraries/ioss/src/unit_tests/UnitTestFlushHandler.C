// Copyright(C) 2024, 2025 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "Ioss_CodeTypes.h"
#include <memory>
#include <string>
#include <vector>

#ifdef SEACAS_HAVE_MPI
#include "mpi.h"
#endif

#include "Ionit_Initializer.h"
#include "Ioss_FlushHandler.h"
#include "Ioss_IOFactory.h"
#include "exodus/Ioex_BaseDatabaseIO.h"
#include <chrono>
#include <gtest/gtest.h>
#include <thread>

class FlushTest : public ::testing::Test
{
protected:
  Ioss::ParallelUtils util_;
  Ioss::FlushHandler  fh{util_};
  int                 rank;
  int                 size;

  void verifyFlush(int state)
  {
    fh.resetTimeStepBegin();
    ASSERT_TRUE(fh.doFlush(state));
  }

  void verifyNoFlush(int state)
  {
    fh.resetTimeStepBegin();
    ASSERT_FALSE(fh.doFlush(state));
  }

#ifdef SEACAS_HAVE_MPI
  MPI_Comm get_split_comm()
  {
    int color = rank % 2;
    int key   = rank;

    MPI_Comm sub_comm;
    MPI_Comm_split(MPI_COMM_WORLD, color, key, &sub_comm);
    return sub_comm;
  }

  void compare_mpi_comms(MPI_Comm comm1, MPI_Comm comm2)
  {
    int cmp = MPI_UNEQUAL;
    MPI_Comm_compare(comm1, comm2, &cmp);
    EXPECT_TRUE(cmp == MPI_IDENT || cmp == MPI_CONGRUENT);
  }
#endif

  void SetUp() override
  {
    rank = 0;
    size = 1;
#ifdef SEACAS_HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
  }

  void TearDown() override {}
};

TEST_F(FlushTest, EmptyTest) {}

TEST_F(FlushTest, DefaultFlushHandler)
{
  ASSERT_EQ(-1, fh.getFlushInterval());
  ASSERT_TRUE(fh.getIsParallel());
  ASSERT_FALSE(fh.getFlushOnFirstOutput());
  ASSERT_EQ(10, fh.getTimeLastFlushInterval());
  ASSERT_EQ(10, fh.getTimeStepBeginFlushInterval());
}

TEST_F(FlushTest, DefaultFlushHandlerDoFlush)
{
  fh.setTimeLastFlushInterval(1);
  std::this_thread::sleep_for(std::chrono::seconds(2));
  verifyFlush(0);
  verifyNoFlush(20);
}

TEST_F(FlushTest, FlushIntervalOneDoFlush)
{
  fh.setFlushInterval(1);
  verifyFlush(0);
  verifyFlush(10);
}

TEST_F(FlushTest, FlushIntervalZeroDoFlush)
{
  fh.setFlushInterval(0);
  verifyNoFlush(0);
  verifyNoFlush(10);
}

TEST_F(FlushTest, FlushIntervalStepDoFlush)
{
  fh.setFlushInterval(20);
  verifyFlush(0);
  verifyNoFlush(10);
  verifyFlush(20);
  verifyNoFlush(35);
  verifyFlush(100);
}

TEST_F(FlushTest, FlushIntervalStepTimeBeginStepDoFlush)
{
  fh.setFlushInterval(10);
  fh.setTimeStepBeginFlushInterval(1);
  fh.resetTimeStepBegin();
  std::this_thread::sleep_for(std::chrono::seconds(2));
  ASSERT_TRUE(fh.doFlush(15));
}

TEST_F(FlushTest, DefaultFlushHandlerTimeBeginStepDoFlush)
{
  fh.setTimeStepBeginFlushInterval(1);
  fh.resetTimeStepBegin();
  std::this_thread::sleep_for(std::chrono::seconds(2));
  ASSERT_TRUE(fh.doFlush(0));
}

TEST_F(FlushTest, DefaultFlushHandlerFlushOnFirstOutput)
{
  fh.setTimeStepBeginFlushInterval(1);
  fh.setFlushOnFirstOutput(true);
  fh.resetTimeStepBegin();
  ASSERT_TRUE(fh.doFlush(0));
  ASSERT_FALSE(fh.doFlush(10));
  std::this_thread::sleep_for(std::chrono::seconds(2));
  ASSERT_TRUE(fh.doFlush(20));
}

TEST_F(FlushTest, DefaultFlushHandlerMultipleRanks)
{
  if (size < 2) {
    GTEST_SKIP() << "Skipping parallel test\n";
  }

  fh.setTimeLastFlushInterval(1);
  if (rank != 0) {
    std::this_thread::sleep_for(std::chrono::seconds(2));
  }
  verifyFlush(10);
}

TEST_F(FlushTest, FlushHandlerSplitMPICommunicator)
{
  if (size < 2) {
    GTEST_SKIP() << "Skipping parallel test\n";
  }

#ifdef SEACAS_HAVE_MPI
  MPI_Comm sub_comm = get_split_comm();

  Ioss::ParallelUtils util_(sub_comm);
  fh = Ioss::FlushHandler(util_);

  compare_mpi_comms(fh.util().communicator(), sub_comm);
#endif
}

TEST_F(FlushTest, FlushHandlerSplitMPICommunicatorDatabase)
{
  if (size < 2) {
    GTEST_SKIP() << "Skipping parallel test\n";
  }

#ifdef SEACAS_HAVE_MPI
  MPI_Comm sub_comm = get_split_comm();

  Ioss::Init::Initializer init_db;
  Ioss::PropertyManager   properties;
  Ioex::BaseDatabaseIO   *dbi = dynamic_cast<Ioex::BaseDatabaseIO *>(
      Ioss::IOFactory::create("exodus", "mesh.exo", Ioss::WRITE_RESTART, sub_comm, properties));
  ASSERT_NE(dbi, nullptr);

  compare_mpi_comms(dbi->get_flush_handler().util().communicator(), sub_comm);
  delete dbi;
#endif
}
