#include <gtest/gtest.h>
#include <stk_util/command_line/CommandLineParser.hpp>

namespace {

class EmptyCommandLine : public ::testing::Test
{
protected:
    static constexpr int argc = 1;
    const char *argv[argc] = {"exeName"};
};

TEST_F(EmptyCommandLine, parseNoOptions_noArgs)
{
    stk::CommandLineParser parser;
    parser.parse(argc, argv);

    EXPECT_TRUE(parser.is_empty());
}

TEST_F(EmptyCommandLine, askForOption_throws)
{
    stk::CommandLineParser parser;
    parser.parse(argc, argv);

    EXPECT_THROW(parser.get_option_value<std::string>("oneOpt"), std::exception);
}

void add_one_option(stk::CommandLineParser &parser)
{
    parser.add_required<std::string>({"oneOpt", "o", "one option"});
}

void add_flag(stk::CommandLineParser &parser)
{
    parser.add_flag("flag,f", "a simple flag");
}

void add_positional_argument(stk::CommandLineParser &parser)
{
    parser.add_required_positional<std::string>({"positional", "p", "a positional argument"});
}

void add_positional_argument_with_default(stk::CommandLineParser &parser)
{
    parser.add_optional_positional<std::string>({"positional", "p", "a positional argument"}, "def");
}

TEST_F(EmptyCommandLine, parseOneOption_noArgs)
{
    stk::CommandLineParser parser;
    add_one_option(parser);
    parser.parse(argc, argv);

    EXPECT_TRUE(parser.is_empty());
}

TEST_F(EmptyCommandLine, parseOneOptionDefaultValue_getDefault)
{
    stk::CommandLineParser parser;
    parser.add_optional<std::string>({"oneOpt" ,"o", "one option"}, "default");
    EXPECT_EQ(stk::CommandLineParser::ParseComplete, parser.parse(argc, argv));

    EXPECT_TRUE(!parser.is_empty());
    EXPECT_EQ("default", parser.get_option_value<std::string>("oneOpt"));
}

TEST_F(EmptyCommandLine, parseOneOptionWithoutAbbreviation)
{
    stk::CommandLineParser parser;
    parser.add_optional<std::string>("oneOpt", "one option", "default");
    EXPECT_EQ(stk::CommandLineParser::ParseComplete, parser.parse(argc, argv));

    EXPECT_TRUE(!parser.is_empty());
    EXPECT_EQ("default", parser.get_option_value<std::string>("oneOpt"));
}

TEST_F(EmptyCommandLine, queryIfFlagProvided_notProvided)
{
    stk::CommandLineParser parser;
    add_flag(parser);
    EXPECT_EQ(stk::CommandLineParser::ParseComplete, parser.parse(argc, argv));

    EXPECT_TRUE(!parser.is_option_provided("flag"));
}

TEST_F(EmptyCommandLine, getUsage_nonEmptyString)
{
    stk::CommandLineParser parser;
    std::string usage = parser.get_usage();
    EXPECT_TRUE(!usage.empty());
}

TEST_F(EmptyCommandLine, getUsageSpecified_usageStartsWithSpecified)
{
    const std::string word("myProgram");
    stk::CommandLineParser parser(word);
    std::string usage = parser.get_usage();
    EXPECT_TRUE(usage.find(word) != std::string::npos);
}

TEST_F(EmptyCommandLine, requiredOneOption_parseError)
{
    stk::CommandLineParser parser;
    add_one_option(parser);
    EXPECT_EQ(stk::CommandLineParser::ParseError, parser.parse(argc, argv));
}

class HelpOnlyCommandLine : public ::testing::Test
{
protected:
    static constexpr int argc = 2;
    const char *argv[argc] = {"exeName", "-h"};
};

TEST_F(HelpOnlyCommandLine, parsing_parsesHelp)
{
    stk::CommandLineParser parser;
    EXPECT_EQ(stk::CommandLineParser::ParseHelpOnly, parser.parse(argc, argv));
}

TEST_F(HelpOnlyCommandLine, requiredArgumentNotGiven_parsesHelp)
{
    stk::CommandLineParser parser;
    add_one_option(parser);
    EXPECT_EQ(stk::CommandLineParser::ParseHelpOnly, parser.parse(argc, argv));
}

class VersionOnlyCommandLine : public ::testing::Test
{
protected:
    static constexpr int argc = 2;
    const char *argv[argc] = {"exeName", "-v"};
};

TEST_F(VersionOnlyCommandLine, parsing_parsesVersion)
{
    stk::CommandLineParser parser;
    EXPECT_EQ(stk::CommandLineParser::ParseVersionOnly, parser.parse(argc, argv));
}

TEST_F(VersionOnlyCommandLine, requiredArgumentNotGiven_parsesVersion)
{
    stk::CommandLineParser parser;
    add_one_option(parser);
    EXPECT_EQ(stk::CommandLineParser::ParseVersionOnly, parser.parse(argc, argv));
}

class FlagCommandLine : public ::testing::Test
{
protected:
    static constexpr int argc = 2;
    const char *argv[argc] = {"exeName", "-f"};
};

TEST_F(FlagCommandLine, queryIfFlagProvided_Provided)
{
    stk::CommandLineParser parser;
    add_flag(parser);
    EXPECT_EQ(stk::CommandLineParser::ParseComplete, parser.parse(argc, argv));
    EXPECT_TRUE(parser.is_option_provided("flag"));
}

class OnePositionalArgumentOnCommandLine : public ::testing::Test
{
protected:
    static constexpr int argc = 2;
    const char *argv[argc] = {"exeName", "value"};
};

TEST_F(OnePositionalArgumentOnCommandLine, parsePositional_parseCompletes)
{
    stk::CommandLineParser parser;
    add_positional_argument(parser);
    EXPECT_EQ(stk::CommandLineParser::ParseComplete, parser.parse(argc, argv));
    EXPECT_EQ(argv[1], parser.get_option_value<std::string>("positional"));
}

TEST_F(EmptyCommandLine, parsePositional_parseError)
{
    stk::CommandLineParser parser;
    add_positional_argument(parser);
    EXPECT_EQ(stk::CommandLineParser::ParseError, parser.parse(argc, argv));
}

TEST_F(OnePositionalArgumentOnCommandLine, parsePositionalWithDefault_parseCompletes)
{
    stk::CommandLineParser parser;
    add_positional_argument_with_default(parser);
    EXPECT_EQ(stk::CommandLineParser::ParseComplete, parser.parse(argc, argv));
    EXPECT_EQ(argv[1], parser.get_option_value<std::string>("positional"));
}

TEST_F(EmptyCommandLine, parsePositionalWithDefault_parseCompletesAndGetDefaultValue)
{
    stk::CommandLineParser parser;
    add_positional_argument_with_default(parser);
    EXPECT_EQ(stk::CommandLineParser::ParseComplete, parser.parse(argc, argv));
    EXPECT_EQ("def", parser.get_option_value<std::string>("positional"));
}

class OneOptionOnCommandLine : public ::testing::Test
{
protected:
    static constexpr int argc = 3;
    const char *argv[argc] = {"exeName", "-o", "value"};
};

TEST_F(OneOptionOnCommandLine, parseOneOption_oneArg)
{
    stk::CommandLineParser parser;
    add_one_option(parser);
    parser.parse(argc, argv);

    ASSERT_TRUE(!parser.is_empty());
    EXPECT_TRUE(parser.is_option_provided("oneOpt"));
    EXPECT_EQ("value", parser.get_option_value<std::string>("oneOpt"));
}

TEST_F(OneOptionOnCommandLine, noOptionsSpecified_throws)
{
    stk::CommandLineParser parser;
    parser.add_required<std::string>({"otherOpt", "p", "other option"});
    EXPECT_EQ(stk::CommandLineParser::ParseError, parser.parse(argc, argv));
}

class TwoOptionsOnCommandLine : public ::testing::Test
{
protected:
    static constexpr int argc = 4;
    const char *argv[argc] = {"exeName", "-o", "value", "--twoOpt=2"};
};

void add_two_option(stk::CommandLineParser &parser)
{
    parser.add_required<int>({"twoOpt", "t", "two option"});
}

TEST_F(TwoOptionsOnCommandLine, parseTwoOption_twoArgs)
{
    stk::CommandLineParser parser;
    add_one_option(parser);
    add_two_option(parser);
    parser.parse(argc, argv);

    ASSERT_TRUE(!parser.is_empty());
    EXPECT_EQ("value", parser.get_option_value<std::string>("oneOpt"));
    EXPECT_EQ(2, parser.get_option_value<int>("twoOpt"));
}

}
