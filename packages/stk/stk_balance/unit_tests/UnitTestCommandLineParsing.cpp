#include <iosfwd>

#include <gtest/gtest.h>
#include <stk_balance/internal/balanceCommandLine.hpp>
#include <stk_util/command_line/CommandLineParserParallel.hpp>
#include <stk_util/parallel/Parallel.hpp>

#include <stk_util/command_line/CommandLineParserUtils.hpp>
#include <fstream>

namespace
{

class BalanceCommandLine : public ::testing::Test
{
protected:
    std::vector<std::string> argv;
    BalanceCommandLine() :
        argv({"exe", "infile"})
    {

    }
    virtual stk::balance::ParsedOptions test_command_line(stk::balance::AppTypeDefaults expectedAppDefaults)
    {
        const std::string quickExample = "";

        MPI_Comm comm(MPI_COMM_WORLD);

        std::vector<const char *> argvCstr = get_argv_c_strings();
        stk::balance::ParsedOptions balanceOptions = stk::balance::parse_balance_command_line(argvCstr.size(), argvCstr.data(), argv[0], comm);

        EXPECT_EQ(expectedAppDefaults, balanceOptions.appTypeDefaults);
        EXPECT_EQ(argv[1], balanceOptions.inFile);

        return balanceOptions;
    }
    std::vector<const char *> get_argv_c_strings()
    {
        std::vector<const char *> argvCstr(argv.size());
        for(size_t i=0; i<argv.size(); i++)
            argvCstr[i] = argv[i].c_str();
        return argvCstr;
    }
    void add_defaults_flag(const std::string &defaultsName)
    {
        std::string defaultsFlag = stk::dash_it(defaultsName);
        argv.push_back(defaultsFlag.c_str());
    }
};

TEST_F(BalanceCommandLine, parse_getValues)
{
    test_command_line(stk::balance::NO_DEFAULTS);
}

TEST_F(BalanceCommandLine, parseSmDefaultOption_smDefaultSet)
{
    add_defaults_flag(stk::balance::CommandLineOptions().smDefaults.name);
    test_command_line(stk::balance::SM_DEFAULTS);
}

TEST_F(BalanceCommandLine, parseSmDefaultOption_sdDefaultSet)
{
    add_defaults_flag(stk::balance::CommandLineOptions().sdDefaults.name);
    test_command_line(stk::balance::SD_DEFAULTS);
}

TEST_F(BalanceCommandLine, parseSmAndSdDefaultOptions_throws)
{
    add_defaults_flag(stk::balance::CommandLineOptions().smDefaults.name);
    add_defaults_flag(stk::balance::CommandLineOptions().sdDefaults.name);
    EXPECT_THROW(test_command_line(stk::balance::SD_DEFAULTS), std::exception);
}

class BalanceCommandLineWithOutputDir : public BalanceCommandLine
{
protected:
    BalanceCommandLineWithOutputDir()
    {
        argv.push_back("outputDir");
    }
    virtual stk::balance::ParsedOptions test_command_line(stk::balance::AppTypeDefaults expectedAppDefaults)
    {
        stk::balance::ParsedOptions balanceOptions = BalanceCommandLine::test_command_line(expectedAppDefaults);
        EXPECT_EQ(argv[2], balanceOptions.outputDirectory);
        return balanceOptions;
    }
};

TEST_F(BalanceCommandLineWithOutputDir, parseSmDefaultOption_smDefaultSet)
{
    add_defaults_flag(stk::balance::CommandLineOptions().smDefaults.name);
    test_command_line(stk::balance::SM_DEFAULTS);
}

}
