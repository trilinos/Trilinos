#include "gtest/gtest.h"
#include "stk_util/environment/OptionsSpecification.hpp"  // for OptionsSpecification, DefaultValue
#include "stk_util/environment/ParseCommandLineArgs.hpp"  // for parse_command_line_args
#include "stk_util/environment/ParsedOptions.hpp"         // for ParsedOptions, VariableType
#include <string>                                         // for string, basic_string
#include <vector>                                         // for vector

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
        ("option,o", "this is an option", stk::DefaultValue<int>(24));

    stk::parse_command_line_args(argc, argv, optionsSpec, parsedOptions);

    EXPECT_FALSE(parsedOptions.empty());
    EXPECT_EQ(1u, parsedOptions.count("flag"));
    EXPECT_EQ(1u, parsedOptions.count("option"));
    int expectedValue = 42;
    EXPECT_EQ(expectedValue, parsedOptions["option"].as<int>());
}

TEST_F(ParseCommandLineFlagAndOption, parseCmdLine_varMapHasFlagAndOption_storeOptionInVariable)
{
    stk::OptionsSpecification optionsSpec;
    stk::ParsedOptions parsedOptions;

    int parsedOptionValue = 0;

    optionsSpec.add_options()
        ("flag,f","this is a non-required flag")
        ("option,o", "this is an option", stk::TargetPointer<int>(&parsedOptionValue), stk::DefaultValue<int>(24));

    stk::parse_command_line_args(argc, argv, optionsSpec, parsedOptions);

    EXPECT_FALSE(parsedOptions.empty());
    EXPECT_EQ(1u, parsedOptions.count("flag"));
    EXPECT_EQ(1u, parsedOptions.count("option"));
    int expectedValue = 42;
    EXPECT_EQ(expectedValue, parsedOptions["option"].as<int>());
    EXPECT_EQ(expectedValue, parsedOptionValue);
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
        ("option,o", "this is an option", stk::DefaultValue<int>(24));

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
        ("option,o", "this is an option", stk::DefaultValue<int>(24))
        ("defaultOption,d", "this is an option with a default value", stk::DefaultValue<double>(99.9));

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

class ParseCommandLineImplicitOptionAndFlag : public ::testing::Test
{
protected:
  static constexpr int argc = 3;
  const char *argv[argc] = {"exeName", "--aprepro", "--flag"};
};

TEST_F(ParseCommandLineImplicitOptionAndFlag, parseCmdLine_varMapHasFlagAndImplicitOptionOn)
{
    stk::OptionsSpecification optionsSpec;
    stk::ParsedOptions parsedOptions;

    optionsSpec.add_options()
        ("flag,f","this is a non-required flag")
        ("aprepro,a", "this is an option", stk::ImplicitValue<std::string>("on"));

    stk::parse_command_line_args(argc, argv, optionsSpec, parsedOptions);

    EXPECT_FALSE(parsedOptions.empty());
    EXPECT_EQ(1u, parsedOptions.count("flag"));
    EXPECT_TRUE(parsedOptions["flag"].empty());
    EXPECT_EQ(1u, parsedOptions.count("aprepro"));
    std::string expectedValue = "on";
    EXPECT_EQ(expectedValue, parsedOptions["aprepro"].as<std::string>());
}

class ParseCommandLineImplicitOptionOffAndFlag : public ::testing::Test
{
protected:
  static constexpr int argc = 4;
  const char *argv[argc] = {"exeName", "--aprepro", "off", "--flag"};
};

TEST_F(ParseCommandLineImplicitOptionOffAndFlag, parseCmdLine_varMapHasFlagAndImplicitOptionOff)
{
    stk::OptionsSpecification optionsSpec;
    stk::ParsedOptions parsedOptions;

    optionsSpec.add_options()
        ("flag,f","this is a non-required flag")
        ("aprepro,a", "this is an option", stk::ImplicitValue<std::string>("on"));

    stk::parse_command_line_args(argc, argv, optionsSpec, parsedOptions);

    EXPECT_FALSE(parsedOptions.empty());
    EXPECT_EQ(1u, parsedOptions.count("flag"));
    EXPECT_TRUE(parsedOptions["flag"].empty());
    EXPECT_EQ(1u, parsedOptions.count("aprepro"));
    std::string expectedValue = "off";
    EXPECT_EQ(expectedValue, parsedOptions["aprepro"].as<std::string>());
}

class ParseCommandLineImplicitOptionNotPassedAndFlag : public ::testing::Test
{
protected:
  static constexpr int argc = 2;
  const char *argv[argc] = {"exeName", "--flag"};
};

TEST_F(ParseCommandLineImplicitOptionNotPassedAndFlag, parseCmdLine_varMapHasFlagAndNotImplicitOption)
{
    stk::OptionsSpecification optionsSpec;
    stk::ParsedOptions parsedOptions;

    optionsSpec.add_options()
        ("flag,f","this is a non-required flag")
        ("aprepro,a", "this is an option", stk::ImplicitValue<std::string>("on"));

    stk::parse_command_line_args(argc, argv, optionsSpec, parsedOptions);

    EXPECT_FALSE(parsedOptions.empty());
    EXPECT_EQ(1u, parsedOptions.count("flag"));
    EXPECT_TRUE(parsedOptions["flag"].empty());
    EXPECT_EQ(0u, parsedOptions.count("aprepro"));
}

class ParseCommandLineImplicitOptionEqualsOffAndFlag : public ::testing::Test
{
protected:
  static constexpr int argc = 3;
  const char *argv[argc] = {"exeName", "--aprepro=off", "--flag"};
};

TEST_F(ParseCommandLineImplicitOptionEqualsOffAndFlag, parseCmdLine_varMapHasFlagAndImplicitOptionOff)
{
    stk::OptionsSpecification optionsSpec;
    stk::ParsedOptions parsedOptions;

    optionsSpec.add_options()
        ("flag,f","this is a non-required flag")
        ("aprepro,a", "this is an option", stk::ImplicitValue<std::string>("on"));

    stk::parse_command_line_args(argc, argv, optionsSpec, parsedOptions);

    EXPECT_FALSE(parsedOptions.empty());
    EXPECT_EQ(1u, parsedOptions.count("flag"));
    EXPECT_TRUE(parsedOptions["flag"].empty());
    EXPECT_EQ(1u, parsedOptions.count("aprepro"));
    std::string expectedValue = "off";
    EXPECT_EQ(expectedValue, parsedOptions["aprepro"].as<std::string>());
}

}
