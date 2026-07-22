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

#include <gtest/gtest.h>                // for AssertHelper, ASSERT_TRUE, etc
#include <stk_util/parallel/Parallel.hpp>
#include <stk_coupling/Utils.hpp>
#include <stk_coupling/SplitComms.hpp>
#include <stk_coupling/SyncInfo.hpp>
#include <stdexcept>

namespace {

TEST(UnitTestSyncInfo, string_to_sync_mode)
{
  EXPECT_EQ(stk::coupling::Minimum, stk::coupling::string_to_sync_mode("Minimum"));
  EXPECT_EQ(stk::coupling::Minimum, stk::coupling::string_to_sync_mode("minIMum"));
  EXPECT_EQ(stk::coupling::Receive, stk::coupling::string_to_sync_mode("Receive"));
  EXPECT_EQ(stk::coupling::Send, stk::coupling::string_to_sync_mode("Send"));
  EXPECT_EQ(stk::coupling::Any, stk::coupling::string_to_sync_mode("Any"));
  EXPECT_ANY_THROW(stk::coupling::string_to_sync_mode("BadInput"));
}

TEST(UnitTestSyncInfo, defaultConstructorStillWorks)
{
  stk::coupling::SyncInfo info;
}

TEST(UnitTestSyncInfo, get_and_set)
{
  using stk::coupling::SyncMode;

  stk::coupling::SyncInfo info("get_and_set");

  bool boolValue = true;
  int intValue = 3;
  SyncMode syncModeValue = SyncMode::Minimum;
  double doubleValue = 3.14;
  std::string stringValue("foo");
  std::vector<int> vectorIntValue = {6, 5, 4, 3};

  info.set_value("boolValue", boolValue);
  info.set_value("intValue", intValue);
  info.set_value("SyncModeValue", syncModeValue);
  info.set_value("doubleValue", doubleValue);
  info.set_value("stringValue", stringValue);
  info.set_value("vectorIntValue", vectorIntValue);

  EXPECT_FALSE(info.has_value<bool>("garbage"));
  EXPECT_FALSE(info.has_value<int>("garbage"));
  EXPECT_FALSE(info.has_value<SyncMode>("garbage"));
  EXPECT_FALSE(info.has_value<double>("garbage"));
  EXPECT_FALSE(info.has_value<std::string>("garbage"));
  EXPECT_FALSE(info.has_value<std::vector<int>>("garbage"));

  EXPECT_TRUE(info.has_value<bool>("boolValue"));
  EXPECT_TRUE(info.has_value<int>("intValue"));
  EXPECT_TRUE(info.has_value<SyncMode>("SyncModeValue"));
  EXPECT_TRUE(info.has_value<double>("doubleValue"));
  EXPECT_TRUE(info.has_value<std::string>("stringValue"));
  EXPECT_TRUE(info.has_value<std::vector<int>>("vectorIntValue"));

  EXPECT_ANY_THROW(info.get_value<bool>("garbage"));
  EXPECT_ANY_THROW(info.get_value<int>("garbage"));
  EXPECT_ANY_THROW(info.get_value<SyncMode>("garbage"));
  EXPECT_ANY_THROW(info.get_value<double>("garbage"));
  EXPECT_ANY_THROW(info.get_value<std::string>("garbage"));
  EXPECT_ANY_THROW(info.get_value<std::vector<int>>("garbage"));

  EXPECT_EQ(boolValue, info.get_value<bool>("boolValue"));
  EXPECT_EQ(intValue, info.get_value<int>("intValue"));
  EXPECT_EQ(syncModeValue, info.get_value<SyncMode>("SyncModeValue"));
  EXPECT_NEAR(doubleValue, info.get_value<double>("doubleValue"), 1.e-12);
  EXPECT_EQ(stringValue, info.get_value<std::string>("stringValue"));
  EXPECT_EQ(vectorIntValue, info.get_value<std::vector<int>>("vectorIntValue"));
}

TEST(UnitTestSyncInfo, getDefaultValue)
{
  using stk::coupling::SyncMode;

  stk::coupling::SyncInfo info("getDefaultValue");

  bool boolValue = true;
  int intValue = 3;
  SyncMode syncModeValue = SyncMode::Minimum;
  double doubleValue = 3.14;
  std::string stringValue("foo");
  std::vector<int> vectorIntValue = {3, 4, 7, 6};

  info.set_value("boolValue", boolValue);
  info.set_value("intValue", intValue);
  info.set_value("SyncModeValue", syncModeValue);
  info.set_value("doubleValue", doubleValue);
  info.set_value("stringValue", stringValue);
  info.set_value("vectorIntValue", vectorIntValue);

  EXPECT_EQ(boolValue, info.get_value<bool>("boolValue", false));
  EXPECT_EQ(intValue, info.get_value<int>("intValue", 1));
  EXPECT_EQ(syncModeValue, info.get_value<SyncMode>("SyncModeValue", SyncMode::Receive));
  EXPECT_NEAR(doubleValue, info.get_value<double>("doubleValue", 0.1), 1.e-12);
  EXPECT_EQ(stringValue, info.get_value<std::string>("stringValue", "default"));
  EXPECT_EQ(vectorIntValue, info.get_value<std::vector<int>>("vectorIntValue", std::vector<int>{7,8,9}));

  EXPECT_EQ(false, info.get_value<bool>("garbage", false));
  EXPECT_EQ(1, info.get_value<int>("garbage", 1));
  EXPECT_EQ(SyncMode::Receive, info.get_value<SyncMode>("garbage", SyncMode::Receive));
  EXPECT_NEAR(0.1, info.get_value<double>("garbage", 0.1), 1.e-12);
  EXPECT_EQ("default", info.get_value<std::string>("garbage", "default"));
  EXPECT_EQ((std::vector<int>{7,8,9}), (info.get_value<std::vector<int>>("garbage", std::vector<int>{7,8,9})));
}

TEST(UnitTestSyncInfo, exchange_2way)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 4) { return; }

  using stk::coupling::SyncMode;

  int myRank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  int color = myRank < 2 ? 0 : 1;
  int otherColor = color==0 ? 1 : 0;
  stk::coupling::SplitComms splitComms(MPI_COMM_WORLD, color);

  stk::coupling::SyncInfo myInfo("exchange");

  bool boolValue = color == 0;
  int intValue = color == 0 ? 3 : 99;
  SyncMode syncModeValue = color == 0 ? SyncMode::Send : SyncMode::Receive;
  double doubleValue = color == 0 ? 3.14 : 99.99;
  std::string stringValue("foo"+std::to_string(color));
  std::vector<int> vectorIntValue = (color == 0) ? std::vector<int>{1, 3, 5} : std::vector<int>{2, 4, 6, 8};
  std::vector<std::pair<std::string,int>> color0_vectorPairStringInt = {{"one",1}, {"two",2}};
  std::vector<std::pair<std::string,int>> color1_vectorPairStringInt = {{"this is a string that is longer than 32 chars",3}};
  std::vector<std::pair<std::string,int>> vectorPairStringIntValue =
     (color == 0) ? color0_vectorPairStringInt : color1_vectorPairStringInt;
  std::vector<std::pair<std::string,double>> color0_vectorPairStringDouble = {{"one",1.0}, {"two",2.0}};
  std::vector<std::pair<std::string,double>> color1_vectorPairStringDouble = {{"this is a string that is longer than 32 chars",3.0}};
  std::vector<std::pair<std::string,double>> vectorPairStringDoubleValue =
     (color == 0) ? color0_vectorPairStringDouble : color1_vectorPairStringDouble;

  myInfo.set_value("boolValue", boolValue);
  myInfo.set_value("intValue", intValue);
  myInfo.set_value("SyncModeValue", syncModeValue);
  myInfo.set_value("doubleValue", doubleValue);
  myInfo.set_value("stringValue", stringValue);
  myInfo.set_value("vectorIntValue", vectorIntValue);
  myInfo.set_value("vectorPairStringIntValue", vectorPairStringIntValue);
  myInfo.set_value("vectorPairStringDoubleValue", vectorPairStringDoubleValue);

  stk::coupling::SyncInfo otherInfo = myInfo.exchange(splitComms, otherColor);

  bool expectedBoolValue = otherColor == 0;
  int expectedIntValue = color == 0 ? 99 : 3;
  SyncMode expectedSyncModeValue = color == 0 ? SyncMode::Receive : SyncMode::Send;
  double expectedDoubleValue = color == 0 ? 99.99 : 3.14;
  std::string expectedStringValue("foo"+std::to_string(otherColor));
  std::vector<int> expectedVectorIntValue = (color == 0) ? std::vector<int>{2, 4, 6, 8} : std::vector<int>{1, 3, 5};
  std::vector<std::pair<std::string,int>> expectedVectorPairStringInt =
    (color == 0) ? color1_vectorPairStringInt : color0_vectorPairStringInt;
  std::vector<std::pair<std::string,double>> expectedVectorPairStringDouble =
    (color == 0) ? color1_vectorPairStringDouble : color0_vectorPairStringDouble;

  EXPECT_EQ(expectedBoolValue, otherInfo.get_value<bool>("boolValue"));
  EXPECT_EQ(expectedIntValue, otherInfo.get_value<int>("intValue"));
  EXPECT_EQ(expectedSyncModeValue, otherInfo.get_value<SyncMode>("SyncModeValue"));
  EXPECT_NEAR(expectedDoubleValue, otherInfo.get_value<double>("doubleValue"), 1.e-12);
  EXPECT_EQ(expectedStringValue, otherInfo.get_value<std::string>("stringValue"));
  EXPECT_EQ(expectedVectorIntValue, otherInfo.get_value<std::vector<int>>("vectorIntValue"));
  EXPECT_EQ(expectedVectorPairStringInt, (otherInfo.get_value<std::vector<std::pair<std::string,int>>>("vectorPairStringIntValue")));
  EXPECT_EQ(expectedVectorPairStringDouble, (otherInfo.get_value<std::vector<std::pair<std::string,double>>>("vectorPairStringDoubleValue")));
}

TEST(UnitTestSyncInfo, exchange_2way_SplitComms_API)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 4) { GTEST_SKIP(); }

  using stk::coupling::SyncMode;

  int myRank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  int color = myRank < 2 ? 0 : 1;
  int otherColor = color == 0 ? 1 : 0;
  stk::coupling::SplitComms splitComms(MPI_COMM_WORLD, color);

  stk::coupling::SyncInfo myInfo("exchange");

  bool boolValue = color == 0;
  int intValue = color == 0 ? 3 : 99;
  SyncMode syncModeValue = color == 0 ? SyncMode::Send : SyncMode::Receive;
  double doubleValue = color == 0 ? 3.14 : 99.99;
  std::string stringValue("foo"+std::to_string(color));
  std::vector<int> vectorIntValue = (color == 0) ? std::vector<int>{1, 3, 5} : std::vector<int>{2, 4, 6, 8};
  std::vector<std::pair<std::string,int>> color0_vectorPairStringInt = {{"one",1}, {"two",2}};
  std::vector<std::pair<std::string,int>> color1_vectorPairStringInt = {{"this is a string that is longer than 32 chars",3}};
  std::vector<std::pair<std::string,int>> vectorPairStringIntValue =
     (color == 0) ? color0_vectorPairStringInt : color1_vectorPairStringInt;

  myInfo.set_value("boolValue", boolValue);
  myInfo.set_value("intValue", intValue);
  myInfo.set_value("SyncModeValue", syncModeValue);
  myInfo.set_value("doubleValue", doubleValue);
  myInfo.set_value("stringValue", stringValue);
  myInfo.set_value("vectorIntValue", vectorIntValue);
  myInfo.set_value("vectorPairStringIntValue", vectorPairStringIntValue);

  stk::coupling::SyncInfo otherInfo = myInfo.exchange(splitComms, otherColor);

  bool expectedBoolValue = otherColor == 0;
  int expectedIntValue = color == 0 ? 99 : 3;
  SyncMode expectedSyncModeValue = color == 0 ? SyncMode::Receive : SyncMode::Send;
  double expectedDoubleValue = color == 0 ? 99.99 : 3.14;
  std::string expectedStringValue("foo"+std::to_string(otherColor));
  std::vector<int> expectedVectorIntValue = (color == 0) ? std::vector<int>{2, 4, 6, 8} : std::vector<int>{1, 3, 5};
  std::vector<std::pair<std::string,int>> expectedVectorPairStringInt =
    (color == 0) ? color1_vectorPairStringInt : color0_vectorPairStringInt;

  EXPECT_EQ(expectedBoolValue, otherInfo.get_value<bool>("boolValue"));
  EXPECT_EQ(expectedIntValue, otherInfo.get_value<int>("intValue"));
  EXPECT_EQ(expectedSyncModeValue, otherInfo.get_value<SyncMode>("SyncModeValue"));
  EXPECT_NEAR(expectedDoubleValue, otherInfo.get_value<double>("doubleValue"), 1.e-12);
  EXPECT_EQ(expectedStringValue, otherInfo.get_value<std::string>("stringValue"));
  EXPECT_EQ(expectedVectorIntValue, otherInfo.get_value<std::vector<int>>("vectorIntValue"));
  EXPECT_EQ(expectedVectorPairStringInt, (otherInfo.get_value<std::vector<std::pair<std::string,int>>>("vectorPairStringIntValue")));
}

TEST(UnitTestSyncInfo, exchangeAsymmetric_2way)
{
  stk::ParallelMachine global = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(global) != 4) return;

  int globalRank = stk::parallel_machine_rank(global);
  const int color = (globalRank >= 2) ? 1 : 0;
  const int otherColor = color==0 ? 1 : 0;
  stk::coupling::SplitComms splitComms(MPI_COMM_WORLD, color);

  stk::coupling::SyncInfo config;
  if (color == 0)
  {
    config.set_value<int>("rank", 2);
    config.set_value<double>("score", 12.75);
  }
  else
  {
    config.set_value<int>("rank", 7);
    config.set_value<std::string>("name", "codeB");
    config.set_value<double>("meta", -3.1);
  }

  stk::coupling::SyncInfo otherConfig = config.exchange(splitComms, otherColor);

  if (color == 0)
  {
    EXPECT_EQ("codeB", otherConfig.get_value<std::string>("name"));
    EXPECT_NEAR(-3.1, otherConfig.get_value<double>("meta"), 1e-12);
    EXPECT_EQ(7, otherConfig.get_value<int>("rank"));
  }
  else
  {
    ASSERT_ANY_THROW(otherConfig.get_value<std::string>("name"));
    EXPECT_NEAR(12.75, otherConfig.get_value<double>("score"), 1e-12);
    EXPECT_EQ(2, otherConfig.get_value<int>("rank"));
  }
}

class TestSyncInfo : public ::testing::Test
{
public:
  TestSyncInfo()
  : m_syncInfo(),
    rankConstant(2),
    nameConstant("code"),
    metaConstant(3.14159)
  {}

  void setup_sync_info(int color)
  {
    m_syncInfo.set_value<int>("rank", rankConstant * color);
    m_syncInfo.set_value<std::string>("name", nameConstant + std::to_string(color));
    m_syncInfo.set_value<double>("meta", metaConstant * color);

    if(color == 0) {
      m_syncInfo.set_value<double>("meta2", 1.5 * metaConstant);
    }
  }

  void check_sync_info(stk::coupling::SyncInfo::ColorToSyncInfoMap syncInfoMap) const
  {
    std::for_each(syncInfoMap.begin(), syncInfoMap.end(), [this](const auto &mapElement) {
      const int color = mapElement.first;
      stk::coupling::SyncInfo syncInfo = mapElement.second;
      EXPECT_EQ(color * rankConstant, syncInfo.get_value<int>("rank"));
      EXPECT_EQ(nameConstant + std::to_string(color), syncInfo.get_value<std::string>("name"));
      EXPECT_NEAR(color * metaConstant, syncInfo.get_value<double>("meta"), 1.e-12);

      if(color == 0) {
        EXPECT_NEAR(1.5 * metaConstant, syncInfo.get_value<double>("meta2"), 1.e-12);
      }
      else {
        EXPECT_FALSE(syncInfo.has_value<double>("meta2"));
      }
    });
  }

  stk::coupling::SyncInfo m_syncInfo;

  int rankConstant;
  std::string nameConstant;
  double metaConstant;

};

TEST_F(TestSyncInfo, exchange_multi_way)
{
  stk::ParallelMachine commWorld = MPI_COMM_WORLD;

  int numColors = stk::parallel_machine_size(commWorld);
  int myColor = stk::parallel_machine_rank(commWorld);
  setup_sync_info(myColor);

  stk::coupling::SplitComms splitComms(commWorld, myColor);
  splitComms.set_free_comms_in_destructor(true);

  using ColorToSyncInfoMap = stk::coupling::SyncInfo::ColorToSyncInfoMap;

  ColorToSyncInfoMap otherInfos = m_syncInfo.exchange(splitComms);
  size_t expectedNumInfos = numColors - 1;
  EXPECT_EQ(expectedNumInfos, otherInfos.size());
  check_sync_info(otherInfos);
}

TEST(UnitTestSyncInfo, choose_value_mine_smaller)
{
  using stk::coupling::SyncMode;

  stk::coupling::SyncInfo myInfo;
  stk::coupling::SyncInfo otherInfo;
  const std::string parameterName = "time step";

  myInfo.set_value(parameterName, 0.01);
  otherInfo.set_value(parameterName, 0.1);

  EXPECT_DOUBLE_EQ(stk::coupling::choose_value(myInfo, otherInfo, parameterName, SyncMode::Send), 0.01);
  EXPECT_DOUBLE_EQ(stk::coupling::choose_value(myInfo, otherInfo, parameterName, SyncMode::Receive), 0.1);
  EXPECT_DOUBLE_EQ(stk::coupling::choose_value(myInfo, otherInfo, parameterName, SyncMode::Minimum), 0.01);
  EXPECT_DOUBLE_EQ(stk::coupling::choose_value(myInfo, otherInfo, parameterName, SyncMode::Any), 0.01);
}

TEST(UnitTestSyncInfo, choose_value_mine_larger)
{
  using stk::coupling::SyncMode;

  stk::coupling::SyncInfo myInfo;
  stk::coupling::SyncInfo otherInfo;
  const std::string parameterName = "time step";

  myInfo.set_value(parameterName, 0.1);
  otherInfo.set_value(parameterName, 0.01);

  EXPECT_DOUBLE_EQ(stk::coupling::choose_value(myInfo, otherInfo, parameterName, SyncMode::Send), 0.1);
  EXPECT_DOUBLE_EQ(stk::coupling::choose_value(myInfo, otherInfo, parameterName, SyncMode::Receive), 0.01);
  EXPECT_DOUBLE_EQ(stk::coupling::choose_value(myInfo, otherInfo, parameterName, SyncMode::Minimum), 0.01);
  EXPECT_DOUBLE_EQ(stk::coupling::choose_value(myInfo, otherInfo, parameterName, SyncMode::Any), 0.1);
}

TEST(UnitTestSyncInfo, choose_value_other_specifies)
{
  using stk::coupling::SyncMode;

  stk::coupling::SyncInfo myInfo;
  stk::coupling::SyncInfo otherInfo;
  const std::string parameterName = "time step";

  myInfo.set_value(parameterName, 0.1);
  otherInfo.set_value(parameterName, 0.01);

  otherInfo.set_value(stk::coupling::TimeSyncMode, stk::coupling::Send);
  EXPECT_DOUBLE_EQ(stk::coupling::choose_value(myInfo, otherInfo, parameterName, SyncMode::Any), 0.01);
  otherInfo.set_value(stk::coupling::TimeSyncMode, stk::coupling::Receive);
  EXPECT_DOUBLE_EQ(stk::coupling::choose_value(myInfo, otherInfo, parameterName, SyncMode::Any), 0.1);
  otherInfo.set_value(stk::coupling::TimeSyncMode, stk::coupling::Minimum);
  EXPECT_DOUBLE_EQ(stk::coupling::choose_value(myInfo, otherInfo, parameterName, SyncMode::Any), 0.01);
}

TEST(UnitTestSyncInfo, choose_value_other_is_Any)
{
  using stk::coupling::SyncMode;

  stk::coupling::SyncInfo myInfo;
  stk::coupling::SyncInfo otherInfo;
  const std::string parameterName = "time step";

  myInfo.set_value(parameterName, 0.1);
  otherInfo.set_value(parameterName, 0.01);
  otherInfo.set_value(stk::coupling::TimeSyncMode, stk::coupling::Any);

  EXPECT_DOUBLE_EQ(stk::coupling::choose_value(myInfo, otherInfo, parameterName, SyncMode::Send), 0.1);
  EXPECT_DOUBLE_EQ(stk::coupling::choose_value(myInfo, otherInfo, parameterName, SyncMode::Receive), 0.01);
  EXPECT_DOUBLE_EQ(stk::coupling::choose_value(myInfo, otherInfo, parameterName, SyncMode::Minimum), 0.01);
  EXPECT_ANY_THROW(stk::coupling::choose_value(myInfo, otherInfo, parameterName, SyncMode::Any));
}

TEST(UnitTestSyncInfo, choose_value_other_missing)
{
  using stk::coupling::SyncMode;

  stk::coupling::SyncInfo myInfo;
  stk::coupling::SyncInfo otherInfo;
  const std::string parameterName = "time step";

  myInfo.set_value(parameterName, 0.1);

  EXPECT_DOUBLE_EQ(stk::coupling::choose_value(myInfo, otherInfo, parameterName, SyncMode::Send), 0.1);
  EXPECT_ANY_THROW(stk::coupling::choose_value(myInfo, otherInfo, parameterName, SyncMode::Receive));
  EXPECT_ANY_THROW(stk::coupling::choose_value(myInfo, otherInfo, parameterName, SyncMode::Minimum));
  EXPECT_DOUBLE_EQ(stk::coupling::choose_value(myInfo, otherInfo, parameterName, SyncMode::Any), 0.1);
}

TEST(UnitTestSyncInfo, choose_value_mine_missing)
{
  using stk::coupling::SyncMode;

  stk::coupling::SyncInfo myInfo;
  stk::coupling::SyncInfo otherInfo;
  const std::string parameterName = "time step";

  otherInfo.set_value(parameterName, 0.1);

  EXPECT_ANY_THROW(stk::coupling::choose_value(myInfo, otherInfo, parameterName, SyncMode::Send));
  EXPECT_DOUBLE_EQ(stk::coupling::choose_value(myInfo, otherInfo, parameterName, SyncMode::Receive), 0.1);
  EXPECT_ANY_THROW(stk::coupling::choose_value(myInfo, otherInfo, parameterName, SyncMode::Minimum));
  EXPECT_ANY_THROW(stk::coupling::choose_value(myInfo, otherInfo, parameterName, SyncMode::Any));
}

TEST(UnitTestSyncInfo, choose_value_both_missing)
{
  using stk::coupling::SyncMode;

  stk::coupling::SyncInfo myInfo;
  stk::coupling::SyncInfo otherInfo;
  const std::string parameterName = "time step";

  EXPECT_ANY_THROW(stk::coupling::choose_value(myInfo, otherInfo, parameterName, SyncMode::Send));
  EXPECT_ANY_THROW(stk::coupling::choose_value(myInfo, otherInfo, parameterName, SyncMode::Receive));
  EXPECT_ANY_THROW(stk::coupling::choose_value(myInfo, otherInfo, parameterName, SyncMode::Minimum));
  EXPECT_ANY_THROW(stk::coupling::choose_value(myInfo, otherInfo, parameterName, SyncMode::Any));
}

class TestSyncMode : public ::testing::Test
{
public:
  TestSyncMode()
  : myColor(stk::parallel_machine_rank(MPI_COMM_WORLD)),
    otherColor(1 - myColor),
    splitComms(MPI_COMM_WORLD, myColor)
  {
  }

  void setup(stk::coupling::SyncMode myMode)
  {
    EXPECT_EQ(1u, splitComms.get_other_colors().size());
    EXPECT_EQ(otherColor, splitComms.get_other_colors()[0]);
    myInfo.set_value(stk::coupling::TimeSyncMode, myMode);
    otherInfo = myInfo.exchange(splitComms, otherColor);
  }

  void setup_color0_no_mode(stk::coupling::SyncMode color1Mode)
  {
    EXPECT_EQ(1u, splitComms.get_other_colors().size());
    EXPECT_EQ(otherColor, splitComms.get_other_colors()[0]);
    if (myColor == 1) {
      myInfo.set_value(stk::coupling::TimeSyncMode, color1Mode);
    }
    otherInfo = myInfo.exchange(splitComms, otherColor);
  }

  int myColor;
  int otherColor;
  stk::coupling::SplitComms splitComms;
  stk::coupling::SyncInfo myInfo;
  stk::coupling::SyncInfo otherInfo;
};

TEST_F(TestSyncMode, check_sync_mode_consistency_Send_Receive)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { GTEST_SKIP(); }
  setup(myColor==0 ? stk::coupling::Send : stk::coupling::Receive);
  EXPECT_NO_THROW(stk::coupling::check_sync_mode_consistency(myInfo, otherInfo));
}

TEST_F(TestSyncMode, check_sync_mode_consistency_Send_Send)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { GTEST_SKIP(); }
  setup(stk::coupling::Send);
  EXPECT_ANY_THROW(stk::coupling::check_sync_mode_consistency(myInfo, otherInfo));
}

TEST_F(TestSyncMode, check_sync_mode_consistency_Receive_Receive)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { GTEST_SKIP(); }
  setup(stk::coupling::Receive);
  EXPECT_ANY_THROW(stk::coupling::check_sync_mode_consistency(myInfo, otherInfo));
}

TEST_F(TestSyncMode, check_sync_mode_consistency_Minimum)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { GTEST_SKIP(); }
  setup(stk::coupling::Minimum);
  EXPECT_NO_THROW(stk::coupling::check_sync_mode_consistency(myInfo, otherInfo));
}

TEST_F(TestSyncMode, check_sync_mode_consistency_Send_Minimum)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { GTEST_SKIP(); }
  setup(myColor==0 ? stk::coupling::Send : stk::coupling::Minimum);
  EXPECT_ANY_THROW(stk::coupling::check_sync_mode_consistency(myInfo, otherInfo));
}

TEST_F(TestSyncMode, check_sync_mode_consistency_Minimum_Receive)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { GTEST_SKIP(); }
  setup(myColor==0 ? stk::coupling::Minimum : stk::coupling::Receive);
  EXPECT_ANY_THROW(stk::coupling::check_sync_mode_consistency(myInfo, otherInfo));
}

TEST_F(TestSyncMode, check_sync_mode_consistency_Send_Any)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { GTEST_SKIP(); }
  setup(myColor==0 ? stk::coupling::Send : stk::coupling::Any);
  EXPECT_NO_THROW(stk::coupling::check_sync_mode_consistency(myInfo, otherInfo));
}

TEST_F(TestSyncMode, check_sync_mode_consistency_Receive_Any)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { GTEST_SKIP(); }
  setup(myColor==0 ? stk::coupling::Receive : stk::coupling::Any);
  EXPECT_NO_THROW(stk::coupling::check_sync_mode_consistency(myInfo, otherInfo));
}

TEST_F(TestSyncMode, check_sync_mode_consistency_Minimum_Any)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { GTEST_SKIP(); }
  setup(myColor==0 ? stk::coupling::Minimum : stk::coupling::Any);
  EXPECT_NO_THROW(stk::coupling::check_sync_mode_consistency(myInfo, otherInfo));
}

TEST_F(TestSyncMode, check_sync_mode_consistency_none_Any)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { GTEST_SKIP(); }
  setup_color0_no_mode(stk::coupling::Any);
  EXPECT_NO_THROW(stk::coupling::check_sync_mode_consistency(myInfo, otherInfo));
}

TEST_F(TestSyncMode, check_sync_mode_consistency_Any_Any)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { GTEST_SKIP(); }
  setup(stk::coupling::Any);
  EXPECT_ANY_THROW(stk::coupling::check_sync_mode_consistency(myInfo, otherInfo));
}

TEST_F(TestSyncMode, check_sync_mode_consistency_none_none)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { GTEST_SKIP(); }
  EXPECT_ANY_THROW(stk::coupling::check_sync_mode_consistency(myInfo, otherInfo));
}

}
