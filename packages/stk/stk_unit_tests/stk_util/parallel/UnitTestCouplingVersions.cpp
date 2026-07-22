#include "gtest/gtest.h"
#include "stk_util/parallel/CouplingVersions.hpp"
#include "stk_util/parallel/CouplingVersions_impl.hpp"
#include "stk_util/parallel/Parallel.hpp"

#ifdef STK_HAS_MPI
namespace {

class CouplingVersionsTester : public ::testing::Test
{
  protected:

    CouplingVersionsTester()
    {
      stk::util::impl::set_error_on_reset(false);
    }
    ~CouplingVersionsTester()
    {
      stk::util::impl::set_coupling_version(MPI_COMM_WORLD, STK_MAX_COUPLING_VERSION);
      stk::util::impl::reset_global_max_coupling_version();
    }

};
}

TEST_F(CouplingVersionsTester, CompatibileRangeGetter)
{
  EXPECT_EQ(stk::util::get_local_max_coupling_version(), STK_MAX_COUPLING_VERSION);
  EXPECT_EQ(stk::util::get_local_min_coupling_version(), STK_MIN_COUPLING_VERSION);
}

TEST_F(CouplingVersionsTester, DefaultVersion)
{
  EXPECT_EQ(stk::util::get_common_coupling_version(), STK_MAX_COUPLING_VERSION);
  EXPECT_EQ(stk::util::get_global_max_coupling_version(), STK_MAX_COUPLING_VERSION);
}

TEST_F(CouplingVersionsTester, NewVersion)
{
  stk::util::impl::set_coupling_version(MPI_COMM_WORLD, STK_MAX_COUPLING_VERSION);
  EXPECT_EQ(stk::util::get_common_coupling_version(), STK_MAX_COUPLING_VERSION);
  EXPECT_EQ(stk::util::get_global_max_coupling_version(), STK_MAX_COUPLING_VERSION);
  EXPECT_FALSE(stk::util::is_local_stk_coupling_deprecated());

}

TEST_F(CouplingVersionsTester, OldVersion)
{
  stk::util::impl::set_coupling_version(MPI_COMM_WORLD, STK_MIN_COUPLING_VERSION);
  EXPECT_EQ(stk::util::get_common_coupling_version(), STK_MIN_COUPLING_VERSION);
  EXPECT_EQ(stk::util::get_global_max_coupling_version(), STK_MAX_COUPLING_VERSION);
}

TEST_F(CouplingVersionsTester, MixedVersion)
{
  int version = stk::parallel_machine_rank(MPI_COMM_WORLD) > 0 ? STK_MIN_COUPLING_VERSION : STK_MAX_COUPLING_VERSION;
  stk::util::impl::set_coupling_version(MPI_COMM_WORLD, version);

  int expectedVal = stk::parallel_machine_size(MPI_COMM_WORLD) == 1 ? STK_MAX_COUPLING_VERSION : STK_MIN_COUPLING_VERSION;
  EXPECT_EQ(stk::util::get_common_coupling_version(), expectedVal);
  EXPECT_EQ(stk::util::get_global_max_coupling_version(), STK_MAX_COUPLING_VERSION);


}

TEST_F(CouplingVersionsTester, DeprecatedVersionCheck)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
    GTEST_SKIP();
  }

  int version = stk::parallel_machine_rank(MPI_COMM_WORLD) == 0 ? STK_MAX_COUPLING_VERSION : STK_MAX_COUPLING_VERSION + 1;
  stk::util::impl::set_coupling_version(MPI_COMM_WORLD, version);

  if (stk::parallel_machine_rank(MPI_COMM_WORLD) == 0) {
    EXPECT_TRUE(stk::util::is_local_stk_coupling_deprecated());
  }
}


TEST_F(CouplingVersionsTester, NewVersionComm)
{
  stk::util::set_coupling_version(MPI_COMM_WORLD);
  EXPECT_EQ(stk::util::get_common_coupling_version(), STK_MAX_COUPLING_VERSION);
  EXPECT_EQ(stk::util::get_global_max_coupling_version(), STK_MAX_COUPLING_VERSION);

}

TEST_F(CouplingVersionsTester, NonincreasingVersion)
{
  stk::util::impl::set_coupling_version(MPI_COMM_WORLD, STK_MIN_COUPLING_VERSION);
  EXPECT_EQ(stk::util::get_common_coupling_version(), STK_MIN_COUPLING_VERSION);
  EXPECT_EQ(stk::util::get_global_max_coupling_version(), STK_MAX_COUPLING_VERSION);


  stk::util::set_coupling_version(MPI_COMM_WORLD);
  EXPECT_EQ(stk::util::get_common_coupling_version(), STK_MIN_COUPLING_VERSION);
  EXPECT_EQ(stk::util::get_global_max_coupling_version(), STK_MAX_COUPLING_VERSION);
}


TEST_F(CouplingVersionsTester, DeprecationDate)
{
  for (int i=STK_MIN_COUPLING_VERSION; i <= STK_MAX_COUPLING_VERSION; ++i)
  {
    std::string date = stk::util::get_deprecation_date(i);
    if (i == STK_MAX_COUPLING_VERSION) {
      EXPECT_EQ(date.size(), 0u);
    } else
    {
      EXPECT_GE(date.size(), 0u);
    }
  }
}

TEST(CouplingVersions, VersionWarningOutput)
{
  testing::internal::CaptureStderr();
  stk::util::print_unsupported_version_warning(-1, __LINE__, __FILE__);
  std::string stderrString = testing::internal::GetCapturedStderr();
  EXPECT_GE(stderrString.size(), 0u);
}

TEST(CouplingVersions, VersionWarningNoOutput)
{
  testing::internal::CaptureStderr();
  stk::util::print_unsupported_version_warning(STK_MIN_COUPLING_VERSION, __LINE__, __FILE__);
  std::string stderrString = testing::internal::GetCapturedStderr();
  EXPECT_EQ(stderrString.size(), 0u);
}
#endif