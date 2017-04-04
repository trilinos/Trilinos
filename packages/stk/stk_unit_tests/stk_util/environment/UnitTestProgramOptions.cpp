#include <gtest/gtest.h>
#include <stk_util/environment/CommandLineParser.hpp>

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
    parser.add_option<std::string>("oneOpt,o", "one option");
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
    parser.add_option_with_default<std::string>("oneOpt,o", "one option", "default");
    parser.parse(argc, argv);

    EXPECT_TRUE(!parser.is_empty());
    EXPECT_EQ("default", parser.get_option_value<std::string>("oneOpt"));
}

TEST_F(EmptyCommandLine, queryIfOptionProvided_notProvided)
{
    stk::CommandLineParser parser;
    add_one_option(parser);
    parser.parse(argc, argv);

    EXPECT_TRUE(!parser.is_option_provided("oneOpt"));
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
    parser.add_option<std::string>("otherOpt,p", "other option");
    EXPECT_THROW(parser.parse(argc, argv), std::exception);
}

class TwoOptionsOnCommandLine : public ::testing::Test
{
protected:
    static constexpr int argc = 4;
    const char *argv[argc] = {"exeName", "-o", "value", "--twoOpt=2"};
};

void add_two_option(stk::CommandLineParser &parser)
{
    parser.add_option<int>("twoOpt,t", "two option");
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
