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
#include <stk_util/parallel/Parallel.hpp>
#include <stk_coupling/Utils.hpp>
#include <stk_coupling/SplitComms.hpp>
#include <stk_coupling/SyncInfo.hpp>
#include <stk_coupling/SplitCommsSingleton.hpp>

//BEGINsplit_comms_example
TEST(StkCouplingDocTest, split_comms)
{
  auto commWorld = MPI_COMM_WORLD;
  auto rank = stk::parallel_machine_rank(commWorld);
  auto commSize = stk::parallel_machine_size(commWorld);

  if (commSize < 2) GTEST_SKIP();

  auto color = rank % 2;

  stk::coupling::SplitComms splitComms(commWorld, color);
  splitComms.set_free_comms_in_destructor(true);

  auto subComm = splitComms.get_split_comm();
  std::vector<int> otherColors = splitComms.get_other_colors();
  EXPECT_EQ(1u, otherColors.size());

  for (auto otherColor : otherColors) {
    auto otherComm = splitComms.get_pairwise_comm(otherColor);

    int result;
    MPI_Comm_compare(subComm, otherComm, &result);
    if (color != otherColor) {
      EXPECT_NE(MPI_IDENT, result);
    } else {
      EXPECT_EQ(MPI_IDENT, result);
    }

    EXPECT_EQ(splitComms.get_parent_comm(), commWorld);
  }
}
//ENDsplit_comms_example

//BEGINsplit_comms_singleton_example
TEST(StkCouplingDocTest, split_comms_singleton)
{
  auto commWorld = MPI_COMM_WORLD;
  auto rank = stk::parallel_machine_rank(commWorld);
  auto color = rank % 2;

  stk::coupling::SplitComms splitComms(commWorld, color);
  splitComms.set_free_comms_in_destructor(true);
  stk::coupling::set_split_comms_singleton(splitComms);

  auto singletonComms = stk::coupling::get_split_comms_singleton();
  EXPECT_TRUE(singletonComms.is_initialized());

  int result;
  MPI_Comm_compare(splitComms.get_split_comm(), singletonComms.get_split_comm(), &result);
  EXPECT_EQ(MPI_IDENT, result);
}
//ENDsplit_comms_singleton_example

TEST(StkCouplingDocTest, sync_info_set_get_values)
{
  using stk::coupling::SyncInfo;
  using stk::coupling::SyncMode;

  double value = -2.0;

  SyncInfo syncInfo("sync_info");
  syncInfo.set_value("SyncModeValue", SyncMode::Minimum);
  syncInfo.set_value("someDoubleValue", value);
  EXPECT_TRUE(syncInfo.has_value<double>("someDoubleValue"));
  EXPECT_DOUBLE_EQ(value, syncInfo.get_value<double>("someDoubleValue"));
}

//BEGINsync_info_exchange_two_colors
TEST(StkCouplingDocTest, sync_info_exchange_two_colors)
{
  using stk::coupling::SplitComms;
  using stk::coupling::SyncInfo;

  auto commWorld = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(commWorld) != 4) GTEST_SKIP();

  auto rank = stk::parallel_machine_rank(commWorld);
  auto color = rank % 2;

  SplitComms splitComms(commWorld, color);
  SyncInfo syncInfo("exchange_info");

  std::string stringValue("DataFrom" + std::to_string(color));
  std::vector<int> intVector = (color == 0) ? std::vector<int>{1, 3, 5} : std::vector<int>{2, 4, 6, 8};
  std::vector<std::pair<std::string, double>> color0_vectorPairStringDouble = {{"one", 1.0}, {"two", 2.0}};
  std::vector<std::pair<std::string, double>> color1_vectorPairStringDouble = {{"three", 3.0}};
  std::vector<std::pair<std::string, double>> vectorPairStringDouble =
      (color == 0) ? color0_vectorPairStringDouble : color1_vectorPairStringDouble;

  syncInfo.set_value("stringToExchange", stringValue);
  syncInfo.set_value("vectorOfIntToExchange", intVector);
  syncInfo.set_value("vectorOfPairToExchange", vectorPairStringDouble);

  auto otherColors = splitComms.get_other_colors();
  SyncInfo exchangeInfo = syncInfo.exchange(splitComms, otherColors[0]);

  std::string exepctedStringValue("DataFrom" + std::to_string(otherColors[0]));
  std::vector<int> expectedIntVector = (color == 1) ? std::vector<int>{1, 3, 5} : std::vector<int>{2, 4, 6, 8};
  std::vector<std::pair<std::string, double>> expectedVectorPairStringDouble =
      (color == 1) ? color0_vectorPairStringDouble : color1_vectorPairStringDouble;

  auto recvString = exchangeInfo.get_value<std::string>("stringToExchange");
  auto recvVectorOfInt = exchangeInfo.get_value<std::vector<int>>("vectorOfIntToExchange");
  auto recvVectorOfPair = exchangeInfo.get_value<std::vector<std::pair<std::string, double>>>("vectorOfPairToExchange");

  EXPECT_EQ(exepctedStringValue, recvString);
  EXPECT_EQ(expectedIntVector, recvVectorOfInt);
  EXPECT_EQ(expectedVectorPairStringDouble, recvVectorOfPair);
}
//ENDsync_info_exchange_two_colors

//BEGINsync_info_exchange_multi_colors
TEST(StkCouplingDocTest, sync_info_exchange_multi_colors)
{
  using stk::coupling::SplitComms;
  using stk::coupling::SyncInfo;

  auto commWorld = MPI_COMM_WORLD;
  auto commSize = stk::parallel_machine_size(commWorld);

  if (commSize % 3 != 0) GTEST_SKIP();

  auto rank = stk::parallel_machine_rank(commWorld);
  auto color = rank;
  auto intToExchange = color * (commSize / 3);

  SyncInfo syncInfo("exchange_info");
  SplitComms splitComms(commWorld, color);
  splitComms.set_free_comms_in_destructor(true);

  syncInfo.set_value<int>("value", intToExchange);

  SyncInfo::ColorToSyncInfoMap otherInfos = syncInfo.exchange(splitComms);

  std::for_each(otherInfos.begin(), otherInfos.end(), [&](const auto &mapElement) {
    [[maybe_unused]] int otherColor = mapElement.first;
    SyncInfo otherSyncInfo = mapElement.second;

    EXPECT_NE(intToExchange, otherSyncInfo.get_value<int>("value"));
    EXPECT_TRUE(otherSyncInfo.has_value<int>("value"));
  });
}
//ENDsync_info_exchange_multi_colors

//BEGINsync_info_choose_values
TEST(StkCouplingDocTest, sync_info_choose_values)
{
  using stk::coupling::SplitComms;
  using stk::coupling::SyncInfo;
  using stk::coupling::SyncMode;

  SyncInfo syncInfo("sync_info");
  SyncInfo otherInfo("other_sync_info");
  const std::string parameterName = "time_step";

  syncInfo.set_value(parameterName, 1.0);
  otherInfo.set_value(parameterName, 2.0);

  EXPECT_DOUBLE_EQ(stk::coupling::choose_value(syncInfo, otherInfo, parameterName, SyncMode::Send), 1.0);
  EXPECT_DOUBLE_EQ(stk::coupling::choose_value(syncInfo, otherInfo, parameterName, SyncMode::Receive), 2.0);
  EXPECT_DOUBLE_EQ(stk::coupling::choose_value(syncInfo, otherInfo, parameterName, SyncMode::Minimum), 1.0);
  EXPECT_DOUBLE_EQ(stk::coupling::choose_value(syncInfo, otherInfo, parameterName, SyncMode::Any), 1.0);
}
//ENDsync_info_choose_values

TEST(StkCouplingDocTest, sync_info_check_sync_mode)
{ 
  using stk::coupling::SplitComms;
  using stk::coupling::SyncInfo;
  using stk::coupling::SyncMode;

  auto commWorld = MPI_COMM_WORLD;
  auto commSize = stk::parallel_machine_size(commWorld);

  if (commSize != 2) GTEST_SKIP();

  auto rank = stk::parallel_machine_rank(commWorld);
  auto color = rank;
  SyncMode modeToCheck = (color == 0) ? stk::coupling::Send : stk::coupling::Receive;

  SplitComms splitComms(commWorld, color);
  splitComms.set_free_comms_in_destructor(true);
  SyncInfo syncInfo("exchange_info");
  syncInfo.set_value(stk::coupling::TimeSyncMode, modeToCheck);

  SyncInfo exchangeInfo = syncInfo.exchange(splitComms, splitComms.get_other_colors()[0]);
  EXPECT_NO_THROW(stk::coupling::check_sync_mode_consistency(syncInfo, exchangeInfo));

  SyncMode modeFail = stk::coupling::Send;
  syncInfo.set_value(stk::coupling::TimeSyncMode, modeFail);
  exchangeInfo = syncInfo.exchange(splitComms, splitComms.get_other_colors()[0]);
  EXPECT_ANY_THROW(stk::coupling::check_sync_mode_consistency(syncInfo, exchangeInfo));
}
