#include <Kokkos_Core.hpp>

#include <gtest/gtest.h>

#include <cstddef>
#include <iostream>

namespace Test {

class ctest_environment : public ::testing::Test {
protected:
  void SetUp();
};

void ctest_environment::SetUp()
{
  setenv("CTEST_KOKKOS_DEVICE_TYPE", "gpus", 1);
  setenv("CTEST_PROCESS_COUNT", "10", 1);
  unsetenv("CTEST_PROCESS_0");
  setenv("CTEST_PROCESS_1", "threads", 1);
  setenv("CTEST_PROCESS_2", "threads,cores", 1);

  setenv("CTEST_PROCESS_3", "gpus", 1);
  unsetenv("CTEST_PROCESS_3_GPUS");

  setenv("CTEST_PROCESS_4", "gpus", 1);
  setenv("CTEST_PROCESS_4_GPUS", "id:2", 1);

  setenv("CTEST_PROCESS_5", "gpus", 1);
  setenv("CTEST_PROCESS_5_GPUS", "slots:1,id:2", 1);

  setenv("CTEST_PROCESS_6", "gpus", 1);
  setenv("CTEST_PROCESS_6_GPUS", "id:2,slots:1", 1);

  setenv("CTEST_PROCESS_7", "threads,gpus", 1);
  setenv("CTEST_PROCESS_7_GPUS", "id:3,slots:1", 1);

  setenv("CTEST_PROCESS_8", "gpus,threads", 1);
  setenv("CTEST_PROCESS_8_GPUS", "id:1,slots:1", 1);

  setenv("CTEST_PROCESS_9", "cores,gpus,threads", 1);
  setenv("CTEST_PROCESS_9_GPUS", "id:4,slots:1", 1);
}

TEST_F(ctest_environment, no_device_type)
{
  unsetenv("CTEST_KOKKOS_DEVICE_TYPE");
  EXPECT_EQ(Kokkos::get_ctest_gpu("0"), 0);
}

TEST_F(ctest_environment, no_process_count)
{
  unsetenv("CTEST_PROCESS_COUNT");
  EXPECT_EQ(Kokkos::get_ctest_gpu("0"), 0);
}

TEST_F(ctest_environment, invalid_rank)
{
  EXPECT_EQ(Kokkos::get_ctest_gpu("10"), 0);
}

TEST_F(ctest_environment, no_type_str)
{
  EXPECT_EQ(Kokkos::get_ctest_gpu("0"), 0);
}

TEST_F(ctest_environment, missing_type)
{
  EXPECT_EQ(Kokkos::get_ctest_gpu("1"), 0);
  EXPECT_EQ(Kokkos::get_ctest_gpu("2"), 0);
}

TEST_F(ctest_environment, no_id_str)
{
  EXPECT_EQ(Kokkos::get_ctest_gpu("3"), 0);
}

TEST_F(ctest_environment, invalid_id_str)
{
  EXPECT_EQ(Kokkos::get_ctest_gpu("4"), 0);
  EXPECT_EQ(Kokkos::get_ctest_gpu("5"), 0);
}

TEST_F(ctest_environment, good)
{
  EXPECT_EQ(Kokkos::get_ctest_gpu("6"), 2);
  EXPECT_EQ(Kokkos::get_ctest_gpu("7"), 3);
  EXPECT_EQ(Kokkos::get_ctest_gpu("8"), 1);
  EXPECT_EQ(Kokkos::get_ctest_gpu("9"), 4);
}

} // namespace Test
