#include "gtest/gtest.h"
#include "stk_util/environment/ParsedOptions.hpp"  // for ParsedOptions, VariableType
#include <string>                                  // for string

namespace {

TEST(ParsedOptions, construct_empty)
{
    stk::ParsedOptions vm;

    EXPECT_TRUE(vm.empty());
}

TEST(ParsedOptions, insert_value_and_count)
{
    stk::ParsedOptions vm;

    vm.insert("foo", "bar");

    EXPECT_EQ(1u, vm.count("foo"));
    EXPECT_FALSE(vm.empty());

    EXPECT_EQ(0u, vm.count("garbage"));
}

TEST(ParsedOptions, get_value_as_type)
{
    stk::ParsedOptions vm;

    vm.insert("int", "99");

    int val = vm["int"].as<int>();
    EXPECT_EQ(99, val);

    std::string strval = vm["int"];
    std::string expectedVal("99");
    EXPECT_EQ(expectedVal, strval);
}

}
