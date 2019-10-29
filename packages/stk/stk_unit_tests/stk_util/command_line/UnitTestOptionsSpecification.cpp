#include <gtest/gtest.h>
#include <stk_util/command_line/OptionsSpecification.hpp>

namespace {

TEST(OptionsSpecification, construct_empty)
{
    stk::OptionsSpecification optionsSpec;

    EXPECT_TRUE(optionsSpec.empty());
}

TEST(OptionsSpecification, insert_flag_and_print)
{
    stk::OptionsSpecification optionsSpec;

    optionsSpec.add_options()("flag,f", "this is a flag");

    EXPECT_EQ(1u, optionsSpec.size());
    EXPECT_FALSE(optionsSpec.empty());

    std::ostringstream os;
    os << optionsSpec;
    std::ostringstream oss;
    oss<<"--flag,-f   this is a flag"<<std::endl;
    std::string expected = oss.str();
    EXPECT_EQ(expected, os.str());
}

TEST(OptionsSpecification, insert_required_option_and_print)
{
    stk::OptionsSpecification optionsSpec;

    const bool isFlag = true;
    optionsSpec.add_options()("option,o", isFlag, true, "this is a required option");

    EXPECT_EQ(1u, optionsSpec.size());
    EXPECT_FALSE(optionsSpec.empty());

    std::ostringstream os;
    os << optionsSpec;
    std::ostringstream oss;
    oss<<"--option,-o   this is a required option (required)"<<std::endl;
    std::string expected = oss.str();
    EXPECT_EQ(expected, os.str());
}

TEST(OptionsSpecification, insert_option_with_default_value_and_print)
{
    stk::OptionsSpecification optionsSpec;

    optionsSpec.add_options()("option,o", "this is an option", 99.9);

    EXPECT_EQ(1u, optionsSpec.size());
    EXPECT_FALSE(optionsSpec.empty());

    std::ostringstream os;
    os << optionsSpec;
    std::ostringstream oss;
    oss<<"--option,-o   this is an option default: 99.9"<<std::endl;
    std::string expected = oss.str();
    EXPECT_EQ(expected, os.str());
}

TEST(OptionsSpecification, insert_flag_and_option_with_default_value_and_print)
{
    stk::OptionsSpecification optionsSpec;

    optionsSpec.add_options()("flag,f", "this is a flag");
    optionsSpec.add_options()("option,o", "this is an option", 99.9);

    EXPECT_EQ(2u, optionsSpec.size());
    EXPECT_FALSE(optionsSpec.empty());

    std::ostringstream os;
    os << optionsSpec;
    std::ostringstream oss;
    oss<<"--flag,-f     this is a flag"<<std::endl;
    oss<<"--option,-o   this is an option default: 99.9"<<std::endl;
    std::string expected = oss.str();
    EXPECT_EQ(expected, os.str());
}

}
