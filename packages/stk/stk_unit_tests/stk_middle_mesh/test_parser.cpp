#include "gtest/gtest.h"

#include "stk_middle_mesh/parser.hpp"

namespace stk {
namespace middle_mesh {
namespace impl {

TEST(Parser, splitSemicolons)
{
  {
    auto s = utils::impl::split_semicolons("abc");
    EXPECT_EQ(s.size(), 1u);
    EXPECT_EQ(s[0], std::string("abc"));
  }

  {
    auto s = utils::impl::split_semicolons("abc;");
    EXPECT_EQ(s.size(), 1u);
    EXPECT_EQ(s[0], std::string("abc"));
  }

  {
    auto s = utils::impl::split_semicolons("abc;def");
    EXPECT_EQ(s.size(), 2u);
    EXPECT_EQ(s[0], std::string("abc"));
    EXPECT_EQ(s[1], std::string("def"));
  }

  {
    auto s = utils::impl::split_semicolons("abc;def;");
    EXPECT_EQ(s.size(), 2u);
    EXPECT_EQ(s[0], std::string("abc"));
    EXPECT_EQ(s[1], std::string("def"));
  }

  {
    auto s = utils::impl::split_semicolons("abc12;def;");
    EXPECT_EQ(s.size(), 2u);
    EXPECT_EQ(s[0], std::string("abc12"));
    EXPECT_EQ(s[1], std::string("def"));
  }

  {
    auto s = utils::impl::split_semicolons("abc12;def");
    EXPECT_EQ(s.size(), 2u);
    EXPECT_EQ(s[0], std::string("abc12"));
    EXPECT_EQ(s[1], std::string("def"));
  }

  {
    auto s = utils::impl::split_semicolons("abc;def12;");
    EXPECT_EQ(s.size(), 2u);
    EXPECT_EQ(s[0], std::string("abc"));
    EXPECT_EQ(s[1], std::string("def12"));
  }

  {
    auto s = utils::impl::split_semicolons("abc;def12");
    EXPECT_EQ(s.size(), 2u);
    EXPECT_EQ(s[0], std::string("abc"));
    EXPECT_EQ(s[1], std::string("def12"));
  }
  // TODO: test errors
}

TEST(Parser, splitSemicolonsErrors)
{
  EXPECT_ANY_THROW(utils::impl::split_semicolons(""));
  EXPECT_ANY_THROW(utils::impl::split_semicolons(";"));
  EXPECT_ANY_THROW(utils::impl::split_semicolons(";;"));
  EXPECT_ANY_THROW(utils::impl::split_semicolons(";;;"));
}

TEST(Parser, trimWhitespace)
{
  EXPECT_EQ(utils::impl::trim_whitespace("abc"), std::string("abc"));
  EXPECT_EQ(utils::impl::trim_whitespace(" abc"), std::string("abc"));
  EXPECT_EQ(utils::impl::trim_whitespace("abc "), std::string("abc"));
  EXPECT_EQ(utils::impl::trim_whitespace(" abc "), std::string("abc"));
  EXPECT_EQ(utils::impl::trim_whitespace("  abc"), std::string("abc"));
  EXPECT_EQ(utils::impl::trim_whitespace("abc  "), std::string("abc"));
  EXPECT_EQ(utils::impl::trim_whitespace("  abc  "), std::string("abc"));
}

TEST(Parser, trimWhitespaceErrors)
{
  EXPECT_ANY_THROW(utils::impl::trim_whitespace(" "));
  EXPECT_ANY_THROW(utils::impl::trim_whitespace("  "));
}

TEST(Parser, splitColons)
{
  {
    std::vector<std::string> strs{"abc:def", "abc:def"};
    auto strPairs = utils::impl::split_colons(strs);
    EXPECT_EQ(strPairs.size(), 2u);
    EXPECT_EQ(strPairs[0].first, "abc");
    EXPECT_EQ(strPairs[0].second, "def");
    EXPECT_EQ(strPairs[1].first, "abc");
    EXPECT_EQ(strPairs[1].second, "def");
  }

  {
    std::vector<std::string> strs{"  abc  :  def12  ", "abc:def12"};
    auto strPairs = utils::impl::split_colons(strs);
    EXPECT_EQ(strPairs.size(), 2u);
    EXPECT_EQ(strPairs[0].first, "abc");
    EXPECT_EQ(strPairs[0].second, "def12");
    EXPECT_EQ(strPairs[1].first, "abc");
    EXPECT_EQ(strPairs[1].second, "def12");
  }
}

TEST(Parser, splitColonsErrors)
{
  {
    std::vector<std::string> strs{"abc::def", "abc:def"};
    EXPECT_ANY_THROW(utils::impl::split_colons(strs));
  }

  {
    std::vector<std::string> strs{"abcdef", "abc:def"};
    EXPECT_ANY_THROW(utils::impl::split_colons(strs));
  }
}

TEST(Parser, parseSideSetNames)
{
  auto names = utils::impl::parse_sideset_names("foo:bar");
  EXPECT_EQ(names.size(), 1u);
  EXPECT_EQ(names[0].first, "foo");
  EXPECT_EQ(names[0].second, "bar");

  names = utils::impl::parse_sideset_names(" foo : bar");
  EXPECT_EQ(names.size(), 1u);
  EXPECT_EQ(names[0].first, "foo");
  EXPECT_EQ(names[0].second, "bar");

  names = utils::impl::parse_sideset_names("foo:bar;foo2:bar2");
  EXPECT_EQ(names.size(), 2u);
  EXPECT_EQ(names[0].first, "foo");
  EXPECT_EQ(names[0].second, "bar");
  EXPECT_EQ(names[1].first, "foo2");
  EXPECT_EQ(names[1].second, "bar2");

  names = utils::impl::parse_sideset_names(" foo : bar ; foo2 : bar2 ");
  EXPECT_EQ(names.size(), 2u);
  EXPECT_EQ(names[0].first, "foo");
  EXPECT_EQ(names[0].second, "bar");
  EXPECT_EQ(names[1].first, "foo2");
  EXPECT_EQ(names[1].second, "bar2");

  names = utils::impl::parse_sideset_names("foo:bar;");
  EXPECT_EQ(names.size(), 1u);
  EXPECT_EQ(names[0].first, "foo");
  EXPECT_EQ(names[0].second, "bar");

  names = utils::impl::parse_sideset_names("foo:bar;foo2:bar2;");
  EXPECT_EQ(names.size(), 2u);
  EXPECT_EQ(names[0].first, "foo");
  EXPECT_EQ(names[0].second, "bar");
  EXPECT_EQ(names[1].first, "foo2");
  EXPECT_EQ(names[1].second, "bar2");

  EXPECT_ANY_THROW(utils::impl::parse_sideset_names("foo::bar"));
  EXPECT_ANY_THROW(utils::impl::parse_sideset_names("foo::bar;;"));
  EXPECT_ANY_THROW(utils::impl::parse_sideset_names("foo:bar;;"));
  EXPECT_ANY_THROW(utils::impl::parse_sideset_names("foo:bar;foo2:bar2;;"));
  EXPECT_ANY_THROW(utils::impl::parse_sideset_names(";foo:bar"));
  EXPECT_ANY_THROW(utils::impl::parse_sideset_names("foo:bar:"));
  EXPECT_ANY_THROW(utils::impl::parse_sideset_names("foo;bar:"));
  EXPECT_ANY_THROW(utils::impl::parse_sideset_names("foo:bar::"));
  EXPECT_ANY_THROW(utils::impl::parse_sideset_names("foo:bar; ;"));
  EXPECT_ANY_THROW(utils::impl::parse_sideset_names("foo:bar; : ;"));
  EXPECT_ANY_THROW(utils::impl::parse_sideset_names(""));
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
