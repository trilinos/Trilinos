#include <gtest/gtest.h>
#include <stk_util/command_line/OptionsSpecification.hpp>
#include <stk_util/command_line/ParsedOptions.hpp>
#include <stk_util/command_line/ParseCommandLineArgs.hpp>

namespace {

class ParseCommandLineOnlyExe : public ::testing::Test
{
protected:
  static constexpr int argc = 1;
  const char *argv[argc] = {"exeName"};
};

TEST_F(ParseCommandLineOnlyExe, parseCmdLine_varMapEmpty)
{
    stk::OptionsSpecification optionsSpec;
    stk::ParsedOptions parsedOptions;

    stk::parse_command_line_args(argc, argv, optionsSpec, parsedOptions);

    EXPECT_TRUE(parsedOptions.empty());
}

class ParseCommandLineOneFlag : public ::testing::Test
{
protected:
  static constexpr int argc = 2;
  const char *argv[argc] = {"exeName", "--flag"};
};

TEST_F(ParseCommandLineOneFlag, parseCmdLine_varMapHasFlag)
{
    stk::OptionsSpecification optionsSpec;
    stk::ParsedOptions parsedOptions;

    optionsSpec.add_options()("flag,f","this is a non-required flag");

    stk::parse_command_line_args(argc, argv, optionsSpec, parsedOptions);

    EXPECT_FALSE(parsedOptions.empty());
    EXPECT_EQ(1u, parsedOptions.count("flag"));
    EXPECT_TRUE(parsedOptions["flag"].empty());
}

class ParseCommandLineFlagAndOption : public ::testing::Test
{
protected:
  static constexpr int argc = 4;
  const char *argv[argc] = {"exeName", "--flag", "--option", "42"};
};

TEST_F(ParseCommandLineFlagAndOption, parseCmdLine_varMapHasFlagAndOption)
{
    stk::OptionsSpecification optionsSpec;
    stk::ParsedOptions parsedOptions;

    optionsSpec.add_options()
        ("flag,f","this is a non-required flag")
        ("option,o", "this is an option", static_cast<int>(24));

    stk::parse_command_line_args(argc, argv, optionsSpec, parsedOptions);

    EXPECT_FALSE(parsedOptions.empty());
    EXPECT_EQ(1u, parsedOptions.count("flag"));
    EXPECT_EQ(1u, parsedOptions.count("option"));
    int expectedValue = 42;
    EXPECT_EQ(expectedValue, parsedOptions["option"].as<int>());
}

class ParseCommandLineFlagAndAbbrevOption : public ::testing::Test
{
protected:
  static constexpr int argc = 4;
  const char *argv[argc] = {"exeName", "--flag", "-o", "42"};
};

TEST_F(ParseCommandLineFlagAndAbbrevOption, parseCmdLine_varMapHasFlagAndOption)
{
    stk::OptionsSpecification optionsSpec;
    stk::ParsedOptions parsedOptions;

    optionsSpec.add_options()
        ("flag,f","this is a non-required flag")
        ("option,o", "this is an option", static_cast<int>(24));

    stk::parse_command_line_args(argc, argv, optionsSpec, parsedOptions);

    EXPECT_FALSE(parsedOptions.empty());
    EXPECT_EQ(1u, parsedOptions.count("flag"));
    EXPECT_TRUE(parsedOptions["flag"].empty());
    EXPECT_EQ(1u, parsedOptions.count("option"));
    int expectedValue = 42;
    EXPECT_EQ(expectedValue, parsedOptions["option"].as<int>());
}

TEST_F(ParseCommandLineFlagAndAbbrevOption, parseCmdLine_varMapHasFlagAndOptionAndDefaultValue)
{
    stk::OptionsSpecification optionsSpec;
    stk::ParsedOptions parsedOptions;

    optionsSpec.add_options()
        ("flag,f","this is a non-required flag")
        ("option,o", "this is an option", static_cast<int>(24))
        ("defaultOption,d", "this is an option with a default value", static_cast<double>(99.9));

    stk::parse_command_line_args(argc, argv, optionsSpec, parsedOptions);

    EXPECT_FALSE(parsedOptions.empty());
    EXPECT_EQ(1u, parsedOptions.count("flag"));
    EXPECT_TRUE(parsedOptions["flag"].empty());
    EXPECT_EQ(1u, parsedOptions.count("option"));
    int expectedValue = 42;
    EXPECT_EQ(expectedValue, parsedOptions["option"].as<int>());
    double expectedDoubleValue = 99.9;
    EXPECT_EQ(1u, parsedOptions.count("defaultOption"));
    EXPECT_EQ(expectedDoubleValue, parsedOptions["defaultOption"].as<double>());
    EXPECT_EQ(expectedDoubleValue, parsedOptions["d"].as<double>());
}

}
