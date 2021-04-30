#include "gtest/gtest.h"
#include "stk_unit_test_utils/MeshFixture.hpp"
#include "stk_balance/setup/M2NParser.hpp"

namespace {

class M2NBalanceCommandLine : public stk::unit_test_util::MeshFixture
{
protected:
  M2NBalanceCommandLine()
    : m_execName("stk_balance_m2n"),
      m_parser(get_comm())
  {
  }

  const stk::balance::M2NParsedOptions & get_parsed_options(const std::vector<std::string>& options)
  {
    std::vector<const char*> args = assemble_args(options);
    m_parser.parse_command_line_options(args.size(), args.data(), m_parsedOptions);
    return m_parsedOptions;
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
  stk::balance::M2NParsedOptions m_parsedOptions;
};

TEST_F(M2NBalanceCommandLine, missingAllRequiredArguments)
{
  testing::internal::CaptureStderr();
  EXPECT_THROW(get_parsed_options({}), std::logic_error);
  testing::internal::GetCapturedStderr();
}

TEST_F(M2NBalanceCommandLine, missingInfilePositionalArgument)
{
  testing::internal::CaptureStderr();
  EXPECT_THROW(get_parsed_options({"16"}), std::logic_error);
  testing::internal::GetCapturedStderr();
}

TEST_F(M2NBalanceCommandLine, missingNprocsPositionalArgument)
{
  testing::internal::CaptureStderr();
  EXPECT_THROW(get_parsed_options({"mesh.g"}), std::logic_error);
  testing::internal::GetCapturedStderr();
}

TEST_F(M2NBalanceCommandLine, missingInfileArgument)
{
  testing::internal::CaptureStderr();
  EXPECT_THROW(get_parsed_options({"--nprocs=16"}), std::logic_error);
  testing::internal::GetCapturedStderr();
}

TEST_F(M2NBalanceCommandLine, missingNprocsArgument)
{
  testing::internal::CaptureStderr();
  EXPECT_THROW(get_parsed_options({"--infile=mesh.g"}), std::logic_error);
  testing::internal::GetCapturedStderr();
}

TEST_F(M2NBalanceCommandLine, misspelledInfileArgument)
{
  testing::internal::CaptureStderr();
  EXPECT_THROW(get_parsed_options({"--infle=mesh.g --nprocs=16"}), std::logic_error);
  testing::internal::GetCapturedStderr();
}

TEST_F(M2NBalanceCommandLine, misspelledNprocsArgument)
{
  testing::internal::CaptureStderr();
  EXPECT_THROW(get_parsed_options({"--infile=mesh.g --procs=16"}), std::logic_error);
  testing::internal::GetCapturedStderr();
}

TEST_F(M2NBalanceCommandLine, normalPostionalArguments)
{
  const stk::balance::M2NParsedOptions & parsedOptions = get_parsed_options({"mesh.g", "16"});

  EXPECT_EQ(parsedOptions.inFile, "mesh.g");
  EXPECT_EQ(parsedOptions.targetNumProcs, 16);
  EXPECT_EQ(parsedOptions.useNestedDecomp, false);
}

TEST_F(M2NBalanceCommandLine, normalArguments)
{
  const stk::balance::M2NParsedOptions & parsedOptions = get_parsed_options({"--infile=mesh.g", "--nprocs=32"});

  EXPECT_EQ(parsedOptions.inFile, "mesh.g");
  EXPECT_EQ(parsedOptions.targetNumProcs, 32);
  EXPECT_EQ(parsedOptions.useNestedDecomp, false);
}

TEST_F(M2NBalanceCommandLine, useNestedDecomp_initialOneProc)
{
  if (get_parallel_size() != 1) return;
  const stk::balance::M2NParsedOptions & parsedOptions = get_parsed_options({"--infile=mesh.g", "--nprocs=2", "--use-nested-decomp"});

  EXPECT_EQ(parsedOptions.inFile, "mesh.g");
  EXPECT_EQ(parsedOptions.targetNumProcs, 2);
  EXPECT_EQ(parsedOptions.useNestedDecomp, true);
}

TEST_F(M2NBalanceCommandLine, useNestedDecomp_sameFinalNumProcs)
{
  if (get_parallel_size() != 2) return;
  const stk::balance::M2NParsedOptions & parsedOptions = get_parsed_options({"--infile=mesh.g", "--nprocs=2", "--use-nested-decomp"});

  EXPECT_EQ(parsedOptions.inFile, "mesh.g");
  EXPECT_EQ(parsedOptions.targetNumProcs, 2);
  EXPECT_EQ(parsedOptions.useNestedDecomp, true);
}

TEST_F(M2NBalanceCommandLine, useNestedDecomp_largerFinalNumProcs)
{
  if (get_parallel_size() != 2) return;
  const stk::balance::M2NParsedOptions & parsedOptions = get_parsed_options({"--infile=mesh.g", "--nprocs=8", "--use-nested-decomp"});

  EXPECT_EQ(parsedOptions.inFile, "mesh.g");
  EXPECT_EQ(parsedOptions.targetNumProcs, 8);
  EXPECT_EQ(parsedOptions.useNestedDecomp, true);
}

TEST_F(M2NBalanceCommandLine, useNestedDecomp_lowerFinalNumProcs)
{
  if (get_parallel_size() != 2) return;
  EXPECT_THROW(get_parsed_options({"--infile=mesh.g", "--nprocs=1", "--use-nested-decomp"}), std::logic_error);
}

TEST_F(M2NBalanceCommandLine, useNestedDecomp_largerInvalidFinalNumProcs)
{
  if (get_parallel_size() != 2) return;
  EXPECT_THROW(get_parsed_options({"--infile=mesh.g", "--nprocs=7", "--use-nested-decomp"}), std::logic_error);
}

}
