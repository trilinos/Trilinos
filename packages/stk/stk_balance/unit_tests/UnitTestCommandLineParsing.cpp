#include <iosfwd>

#include <gtest/gtest.h>
#include <stk_util/parallel/Parallel.hpp>

#include <stk_balance/internal/Inputs.hpp>
#include <stk_util/command_line/CommandLineParserUtils.hpp>
#include <fstream>

namespace
{

void touch_new_file(const std::string& filename)
{
    std::ofstream out(filename, std::ios_base::app);
    out << "hi";
    out.close();
}

stk::balance::Inputs get_inputs_for_no_options()
{
    const int numInputs = 1;
    const char* inputs[numInputs] = { "stk_balance" };
    stk::balance::Inputs parsedInput(numInputs, inputs);
    return parsedInput;
}

stk::balance::Inputs get_inputs_for_filename_only()
{
    const int numInputs = 2;
    const char* inputs[numInputs] = { "stk_balance", "junk.exo"};
    stk::balance::Inputs parsedInput(numInputs, inputs);
    return parsedInput;
}

stk::balance::Inputs get_inputs_for_valid_inputs()
{
    const int numInputs = 3;
    const char* inputs[numInputs] = { "stk_balance", "junk.exo", "alpha" };
    stk::balance::Inputs parsedInput(numInputs, inputs);
    return parsedInput;
}

TEST(CommandLineParse, testNoOptions)
{
    stk::balance::Inputs parsedInput = get_inputs_for_no_options();
    EXPECT_EQ("stk_balance", parsedInput.get_executable_name());
    EXPECT_EQ("", parsedInput.get_exodus_filename());
    EXPECT_EQ(".", parsedInput.get_output_directory());
}

TEST(CommandLineParse, testWithOneOption)
{
    stk::balance::Inputs parsedInput= get_inputs_for_filename_only();
    EXPECT_EQ("stk_balance", parsedInput.get_executable_name());
    EXPECT_EQ("junk.exo", parsedInput.get_exodus_filename());
    EXPECT_EQ(".", parsedInput.get_output_directory());
}

TEST(CommandLineParse, testOptions)
{
    stk::balance::Inputs parsedInput= get_inputs_for_valid_inputs();
    EXPECT_EQ("stk_balance", parsedInput.get_executable_name());
    EXPECT_EQ("junk.exo", parsedInput.get_exodus_filename());
    EXPECT_EQ("alpha", parsedInput.get_output_directory());
}

TEST(ParsedInput, withNoInputsWriteUsageInfo)
{
    stk::balance::Inputs parsedInput = get_inputs_for_no_options();
    EXPECT_TRUE( stk::balance::should_write_usage_info( parsedInput.get_exodus_filename() ) );
}

TEST(ParsedInput, checkIfFilenameExistsWhenFileExists)
{
    if(stk::parallel_machine_size(MPI_COMM_WORLD)==1)
    {
        stk::balance::Inputs parsedInput= get_inputs_for_filename_only();
        touch_new_file(parsedInput.get_exodus_filename());
        EXPECT_TRUE( stk::parallel::does_file_exist( parsedInput.get_exodus_filename() ));
        unlink(parsedInput.get_exodus_filename().c_str());
    }
}

TEST(ParsedInput, checkIfFilenameExistsWhenFileDoesNotExist)
{
    stk::balance::Inputs parsedInput= get_inputs_for_filename_only();
    EXPECT_FALSE( stk::parallel::does_file_exist( parsedInput.get_exodus_filename() ));
}


TEST(ParsedInput, outputDirectoryExistsAfterCallToCreatePath)
{
    if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
    {
        stk::balance::Inputs parsedInput = get_inputs_for_valid_inputs();
        EXPECT_TRUE(stk::balance::create_path(parsedInput.get_output_directory()));
    }
}


}
