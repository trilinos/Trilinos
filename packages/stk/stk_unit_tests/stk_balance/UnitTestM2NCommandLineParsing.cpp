#include "gtest/gtest.h"
#include "stk_unit_test_utils/MeshFixture.hpp"
#include "stk_balance/setup/M2NParser.hpp"
#include "stk_balance/balanceUtils.hpp"

namespace {

class M2NBalanceCommandLine : public stk::unit_test_util::MeshFixture
{
protected:
  M2NBalanceCommandLine()
    : m_execName("stk_balance_m2n"),
      m_parser(get_comm())
  { }

  const stk::balance::M2NBalanceSettings & get_balance_settings(const std::vector<std::string>& options)
  {
    std::vector<const char*> args = assemble_args(options);
    m_parser.parse_command_line_options(args.size(), args.data(), m_balanceSettings);
    return m_balanceSettings;
  }

  std::vector<const char*> assemble_args(const std::vector<std::string>& options) const
  {
    std::vector<const char*> args = {m_execName.c_str()};

    for (const std::string& option : options) {
      args.push_back(option.c_str());
    }

    return args;
  }

private:
  std::string m_execName;
  stk::balance::M2NParser m_parser;
  stk::balance::M2NBalanceSettings m_balanceSettings;
};

TEST_F(M2NBalanceCommandLine, missingAllRequiredArguments)
{
  testing::internal::CaptureStderr();
  EXPECT_THROW(get_balance_settings({}), std::logic_error);
  testing::internal::GetCapturedStderr();
}

TEST_F(M2NBalanceCommandLine, missingInfilePositionalArgument)
{
  testing::internal::CaptureStderr();
  EXPECT_THROW(get_balance_settings({"16"}), std::logic_error);
  testing::internal::GetCapturedStderr();
}

TEST_F(M2NBalanceCommandLine, missingNprocsPositionalArgument)
{
  testing::internal::CaptureStderr();
  EXPECT_THROW(get_balance_settings({"mesh.g"}), std::logic_error);
  testing::internal::GetCapturedStderr();
}

TEST_F(M2NBalanceCommandLine, missingInfileArgument)
{
  testing::internal::CaptureStderr();
  EXPECT_THROW(get_balance_settings({"--nprocs=16"}), std::logic_error);
  testing::internal::GetCapturedStderr();
}

TEST_F(M2NBalanceCommandLine, missingNprocsArgument)
{
  testing::internal::CaptureStderr();
  EXPECT_THROW(get_balance_settings({"--infile=mesh.g"}), std::logic_error);
  testing::internal::GetCapturedStderr();
}

TEST_F(M2NBalanceCommandLine, misspelledInfileArgument)
{
  testing::internal::CaptureStderr();
  EXPECT_THROW(get_balance_settings({"--infle=mesh.g --nprocs=16"}), std::logic_error);
  testing::internal::GetCapturedStderr();
}

TEST_F(M2NBalanceCommandLine, misspelledNprocsArgument)
{
  testing::internal::CaptureStderr();
  EXPECT_THROW(get_balance_settings({"--infile=mesh.g --procs=16"}), std::logic_error);
  testing::internal::GetCapturedStderr();
}

TEST_F(M2NBalanceCommandLine, normalPostionalArguments)
{
  const stk::balance::M2NBalanceSettings & balanceSettings = get_balance_settings({"mesh.g", "16"});

  const int initialNumProcs = stk::parallel_machine_size(MPI_COMM_WORLD);

  EXPECT_EQ(balanceSettings.get_input_filename(), "mesh.g");
  EXPECT_EQ(balanceSettings.get_log_filename(), "mesh." + std::to_string(initialNumProcs) + "_to_16.log");
  EXPECT_EQ(balanceSettings.get_num_output_processors(), 16u);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(), false);
}

TEST_F(M2NBalanceCommandLine, normalArguments)
{
  const stk::balance::M2NBalanceSettings & balanceSettings = get_balance_settings({"--infile=mesh.g", "--nprocs=32"});

  const int initialNumProcs = stk::parallel_machine_size(MPI_COMM_WORLD);
  
  EXPECT_EQ(balanceSettings.get_input_filename(), "mesh.g");
  EXPECT_EQ(balanceSettings.get_log_filename(), "mesh." + std::to_string(initialNumProcs) + "_to_32.log");
  EXPECT_EQ(balanceSettings.get_num_output_processors(), 32u);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(), false);
}

TEST_F(M2NBalanceCommandLine, customLogFile)
{
  const stk::balance::M2NBalanceSettings & balanceSettings = get_balance_settings({"mesh.g", "16", "--logfile=custom.log"});

  EXPECT_EQ(balanceSettings.get_input_filename(), "mesh.g");
  EXPECT_EQ(balanceSettings.get_log_filename(), "custom.log");
  EXPECT_EQ(balanceSettings.get_num_output_processors(), 16u);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(), false);
}

TEST_F(M2NBalanceCommandLine, shortCustomLogFile)
{
  const stk::balance::M2NBalanceSettings & balanceSettings = get_balance_settings({"mesh.g", "16", "-l", "custom.log"});

  EXPECT_EQ(balanceSettings.get_input_filename(), "mesh.g");
  EXPECT_EQ(balanceSettings.get_log_filename(), "custom.log");
  EXPECT_EQ(balanceSettings.get_num_output_processors(), 16u);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(), false);
}

TEST_F(M2NBalanceCommandLine, coutLogFile)
{
  const stk::balance::M2NBalanceSettings & balanceSettings = get_balance_settings({"mesh.g", "16", "-l", "cout"});

  EXPECT_EQ(balanceSettings.get_input_filename(), "mesh.g");
  EXPECT_EQ(balanceSettings.get_log_filename(), "cout");
  EXPECT_EQ(balanceSettings.get_num_output_processors(), 16u);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(), false);
}

TEST_F(M2NBalanceCommandLine, useNestedDecomp_initialOneProc)
{
  if (get_parallel_size() != 1) return;
  const stk::balance::M2NBalanceSettings & balanceSettings = get_balance_settings({"--infile=mesh.g", "--nprocs=2", "--use-nested-decomp"});

  EXPECT_EQ(balanceSettings.get_input_filename(), "mesh.g");
  EXPECT_EQ(balanceSettings.get_num_output_processors(), 2u);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(), true);
}

TEST_F(M2NBalanceCommandLine, useNestedDecomp_sameFinalNumProcs)
{
  if (get_parallel_size() != 2) return;
  const stk::balance::M2NBalanceSettings & balanceSettings = get_balance_settings({"--infile=mesh.g", "--nprocs=2", "--use-nested-decomp"});

  EXPECT_EQ(balanceSettings.get_input_filename(), "mesh.g");
  EXPECT_EQ(balanceSettings.get_num_output_processors(), 2u);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(), true);
}

TEST_F(M2NBalanceCommandLine, useNestedDecomp_largerFinalNumProcs)
{
  if (get_parallel_size() != 2) return;
  const stk::balance::M2NBalanceSettings & balanceSettings = get_balance_settings({"--infile=mesh.g", "--nprocs=8", "--use-nested-decomp"});

  EXPECT_EQ(balanceSettings.get_input_filename(), "mesh.g");
  EXPECT_EQ(balanceSettings.get_num_output_processors(), 8u);
  EXPECT_EQ(balanceSettings.get_use_nested_decomp(), true);
}

TEST_F(M2NBalanceCommandLine, useNestedDecomp_lowerFinalNumProcs)
{
  if (get_parallel_size() != 2) return;
  EXPECT_THROW(get_balance_settings({"--infile=mesh.g", "--nprocs=1", "--use-nested-decomp"}), std::logic_error);
}

TEST_F(M2NBalanceCommandLine, useNestedDecomp_largerInvalidFinalNumProcs)
{
  if (get_parallel_size() != 2) return;
  EXPECT_THROW(get_balance_settings({"--infile=mesh.g", "--nprocs=7", "--use-nested-decomp"}), std::logic_error);
}

}
