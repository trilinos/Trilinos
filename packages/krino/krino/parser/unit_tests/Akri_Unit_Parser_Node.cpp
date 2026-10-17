// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Parser.hpp>

#include <stk_util/environment/RuntimeDoomed.hpp>

#include <gtest/gtest.h>
#include <Akri_YAML.hpp>

#include <string>
#include <vector>

namespace krino {
namespace {

#ifdef KRINO_HAVE_YAML

// Unit tests for the Parser::Node abstraction over YAML::Node.  These are the
// lowest-level building block used by every krino YAML parser, and they need no
// MetaData / Region / mesh / MPI setup, so they exercise the node predicates,
// the get_*_if_present lookups (including their RuntimeDoomed error paths), the
// bool conversion, and iteration in complete isolation.
//
// RuntimeDoomedAdHoc does not throw; it increments a monotonic process-global
// counter.  Following the convention in Akri_Unit_Surface_Parser.cpp, error
// paths capture stk::get_doomed_count() before the call and assert the count
// increased, while happy paths assert the count is unchanged.
class ParserNodeTest : public ::testing::Test
{
protected:
  Parser::Node load_node(const std::string & yamlText)
  {
    return Parser::Node(YAML::Load(yamlText));
  }
};

// -----------------------------------------------------------------------------
// Type predicates and operator bool
// -----------------------------------------------------------------------------

TEST_F(ParserNodeTest, DefaultConstructedNodeIsFalsy)
{
  Parser::Node node;
  EXPECT_FALSE(bool(node));
}

TEST_F(ParserNodeTest, MapNodeIsMap)
{
  Parser::Node node = load_node(R"yaml(
alpha: 1
beta: 2
)yaml");
  EXPECT_TRUE(bool(node));
  EXPECT_TRUE(node.is_map());
  EXPECT_FALSE(node.is_sequence());
  EXPECT_FALSE(node.is_scalar());
  EXPECT_FALSE(node.is_null());
  EXPECT_EQ(2u, node.size());
}

TEST_F(ParserNodeTest, SequenceNodeIsSequence)
{
  Parser::Node node = load_node(R"yaml(
- 1
- 2
- 3
)yaml");
  EXPECT_TRUE(node.is_sequence());
  EXPECT_FALSE(node.is_map());
  EXPECT_EQ(3u, node.size());
}

TEST_F(ParserNodeTest, ScalarNodeIsScalar)
{
  Parser::Node node = load_node("hello");
  EXPECT_TRUE(node.is_scalar());
  EXPECT_FALSE(node.is_map());
  EXPECT_FALSE(node.is_sequence());
}

TEST_F(ParserNodeTest, NullValueNodeIsNull)
{
  Parser::Node map = load_node(R"yaml(
empty_value:
)yaml");
  Parser::Node child = map.get_if_present("empty_value");
  EXPECT_TRUE(bool(child));
  EXPECT_TRUE(child.is_null());
}

// -----------------------------------------------------------------------------
// is_empty_scalar
// -----------------------------------------------------------------------------

TEST_F(ParserNodeTest, IsEmptyScalarTrueForEmptyQuotedString)
{
  Parser::Node map = load_node(R"yaml(
blank: ""
)yaml");
  Parser::Node child = map.get_if_present("blank");
  EXPECT_TRUE(child.is_scalar());
  EXPECT_TRUE(child.is_empty_scalar());
}

TEST_F(ParserNodeTest, IsEmptyScalarFalseForNonEmptyScalar)
{
  Parser::Node node = load_node("value");
  EXPECT_FALSE(node.is_empty_scalar());
}

TEST_F(ParserNodeTest, IsEmptyScalarFalseForNonScalar)
{
  Parser::Node node = load_node(R"yaml(
- 1
- 2
)yaml");
  EXPECT_FALSE(node.is_empty_scalar());
}

// -----------------------------------------------------------------------------
// get_if_present(key) -> Node
// -----------------------------------------------------------------------------

TEST_F(ParserNodeTest, GetIfPresentReturnsChildForExistingKey)
{
  Parser::Node node = load_node(R"yaml(
alpha: 1
)yaml");
  Parser::Node child = node.get_if_present("alpha");
  EXPECT_TRUE(bool(child));
  EXPECT_EQ(1, child.as<int>());
}

TEST_F(ParserNodeTest, GetIfPresentReturnsEmptyForAbsentKeyWithoutDooming)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  Parser::Node node = load_node(R"yaml(
alpha: 1
)yaml");
  Parser::Node child = node.get_if_present("does_not_exist");
  EXPECT_FALSE(bool(child));
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(ParserNodeTest, GetIfPresentReturnsEmptyForFalsyNodeWithoutDooming)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  Parser::Node node; // default-constructed / undefined
  Parser::Node child = node.get_if_present("anything");
  EXPECT_FALSE(bool(child));
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(ParserNodeTest, GetIfPresentOnScalarNodeDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  Parser::Node node = load_node("just_a_scalar");
  Parser::Node child = node.get_if_present("key");
  EXPECT_FALSE(bool(child));
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

// -----------------------------------------------------------------------------
// get_map/sequence/scalar/null_if_present -> type-checked lookups
// -----------------------------------------------------------------------------

TEST_F(ParserNodeTest, GetMapIfPresentReturnsMap)
{
  Parser::Node node = load_node(R"yaml(
child:
  inner: 1
)yaml");
  Parser::Node child = node.get_map_if_present("child");
  EXPECT_TRUE(bool(child));
  EXPECT_TRUE(child.is_map());
}

TEST_F(ParserNodeTest, GetSequenceIfPresentReturnsSequence)
{
  Parser::Node node = load_node(R"yaml(
items:
  - 1
  - 2
)yaml");
  Parser::Node seq = node.get_sequence_if_present("items");
  EXPECT_TRUE(bool(seq));
  EXPECT_TRUE(seq.is_sequence());
  EXPECT_EQ(2u, seq.size());
}

TEST_F(ParserNodeTest, GetScalarIfPresentReturnsScalar)
{
  Parser::Node node = load_node(R"yaml(
name: value
)yaml");
  Parser::Node scalar = node.get_scalar_if_present("name");
  EXPECT_TRUE(bool(scalar));
  EXPECT_TRUE(scalar.is_scalar());
}

TEST_F(ParserNodeTest, GetTypeIfPresentReturnsEmptyForAbsentKeyWithoutDooming)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  Parser::Node node = load_node(R"yaml(
alpha: 1
)yaml");
  Parser::Node child = node.get_map_if_present("missing");
  EXPECT_FALSE(bool(child));
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(ParserNodeTest, GetMapIfPresentOnSequenceValueDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  Parser::Node node = load_node(R"yaml(
child:
  - 1
  - 2
)yaml");
  Parser::Node child = node.get_map_if_present("child");
  EXPECT_FALSE(bool(child));
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(ParserNodeTest, GetSequenceIfPresentOnScalarValueDooms)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  Parser::Node node = load_node(R"yaml(
child: 5
)yaml");
  Parser::Node child = node.get_sequence_if_present("child");
  EXPECT_FALSE(bool(child));
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

// -----------------------------------------------------------------------------
// get_options_if_present -> options nested under a keyword or beside it
// -----------------------------------------------------------------------------

TEST_F(ParserNodeTest, GetOptionsIfPresentReturnsNestedMap)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  Parser::Node node = load_node(R"yaml(
plane:
  offset: 1
)yaml");
  Parser::Node options = node.get_options_if_present("plane");
  ASSERT_TRUE(bool(options));
  EXPECT_TRUE(options.is_map());
  EXPECT_FALSE(bool(options.get_if_present("plane")));
  int offset = 0;
  EXPECT_TRUE(options.get_if_present("offset", offset));
  EXPECT_EQ(1, offset);
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(ParserNodeTest, GetOptionsIfPresentReturnsEnclosingNodeForNullValue)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  Parser::Node node = load_node(R"yaml(
plane:
offset: 1
)yaml");
  Parser::Node options = node.get_options_if_present("plane");
  ASSERT_TRUE(bool(options));
  EXPECT_TRUE(bool(options.get_if_present("plane")));
  int offset = 0;
  EXPECT_TRUE(options.get_if_present("offset", offset));
  EXPECT_EQ(1, offset);
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(ParserNodeTest, GetOptionsIfPresentReturnsEnclosingNodeForEmptyScalarValue)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  Parser::Node node = load_node(R"yaml(
plane: ""
offset: 1
)yaml");
  Parser::Node options = node.get_options_if_present("plane");
  ASSERT_TRUE(bool(options));
  int offset = 0;
  EXPECT_TRUE(options.get_if_present("offset", offset));
  EXPECT_EQ(1, offset);
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(ParserNodeTest, GetOptionsIfPresentReturnsEmptyForAbsentKeyWithoutDooming)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  Parser::Node node = load_node(R"yaml(
offset: 1
)yaml");
  EXPECT_FALSE(bool(node.get_options_if_present("plane")));
  EXPECT_EQ(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(ParserNodeTest, GetOptionsIfPresentOnScalarOrSequenceValueDooms)
{
  for (const std::string yamlText : {"plane: yes\noffset: 1\n", "plane: [1, 2]\noffset: 1\n"})
  {
    const unsigned oldDoomedCount = stk::get_doomed_count();
    Parser::Node node = load_node(yamlText);
    EXPECT_FALSE(bool(node.get_options_if_present("plane"))) << yamlText;
    EXPECT_GT(stk::get_doomed_count(), oldDoomedCount) << yamlText;
  }
}

// -----------------------------------------------------------------------------
// Templated get_if_present(key, T& result)
// -----------------------------------------------------------------------------

TEST_F(ParserNodeTest, GetTypedIntWhenPresent)
{
  Parser::Node node = load_node(R"yaml(
count: 42
)yaml");
  int value = 0;
  EXPECT_TRUE(node.get_if_present("count", value));
  EXPECT_EQ(42, value);
}

TEST_F(ParserNodeTest, GetTypedDoubleWhenPresent)
{
  Parser::Node node = load_node(R"yaml(
scale: 1.5
)yaml");
  double value = 0.0;
  EXPECT_TRUE(node.get_if_present("scale", value));
  EXPECT_DOUBLE_EQ(1.5, value);
}

TEST_F(ParserNodeTest, GetTypedStringWhenPresent)
{
  Parser::Node node = load_node(R"yaml(
label: hello
)yaml");
  std::string value;
  EXPECT_TRUE(node.get_if_present("label", value));
  EXPECT_EQ("hello", value);
}

TEST_F(ParserNodeTest, GetTypedVectorWhenPresent)
{
  Parser::Node node = load_node(R"yaml(
coords: [1.0, 2.0, 3.0]
)yaml");
  std::vector<double> value;
  EXPECT_TRUE(node.get_if_present("coords", value));
  ASSERT_EQ(3u, value.size());
  EXPECT_DOUBLE_EQ(1.0, value[0]);
  EXPECT_DOUBLE_EQ(2.0, value[1]);
  EXPECT_DOUBLE_EQ(3.0, value[2]);
}

TEST_F(ParserNodeTest, GetTypedReturnsFalseAndLeavesResultWhenAbsent)
{
  Parser::Node node = load_node(R"yaml(
present: 1
)yaml");
  int value = 7;
  EXPECT_FALSE(node.get_if_present("absent", value));
  EXPECT_EQ(7, value); // untouched
}

TEST_F(ParserNodeTest, GetTypedBoolWhenPresent)
{
  Parser::Node node = load_node(R"yaml(
enabled: true
)yaml");
  bool value = false;
  EXPECT_TRUE(node.get_if_present("enabled", value));
  EXPECT_TRUE(value);
}

TEST_F(ParserNodeTest, GetTypedNonConvertibleScalarDoomsAndLeavesResult)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  Parser::Node node = load_node(R"yaml(
count: not_an_int
)yaml");
  int value = 7;
  // A present-but-unconvertible scalar must doom (not throw) and return false,
  // leaving the caller's result untouched so parse() can collect and report it.
  EXPECT_FALSE(node.get_if_present("count", value));
  EXPECT_EQ(7, value); // untouched
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

TEST_F(ParserNodeTest, GetTypedNonConvertibleDoubleDoomsAndLeavesResult)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  Parser::Node node = load_node(R"yaml(
scale: abc
)yaml");
  double value = 1.5;
  EXPECT_FALSE(node.get_if_present("scale", value));
  EXPECT_DOUBLE_EQ(1.5, value); // untouched
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

// -----------------------------------------------------------------------------
// as<bool>() truth table + failure path
// -----------------------------------------------------------------------------

TEST_F(ParserNodeTest, AsBoolTrueVariants)
{
  EXPECT_TRUE(load_node("true").as<bool>());
  EXPECT_TRUE(load_node("True").as<bool>());
  EXPECT_TRUE(load_node("YES").as<bool>());
  EXPECT_TRUE(load_node("yes").as<bool>());
  EXPECT_TRUE(load_node("1").as<bool>());
}

TEST_F(ParserNodeTest, AsBoolFalseVariants)
{
  EXPECT_FALSE(load_node("false").as<bool>());
  EXPECT_FALSE(load_node("False").as<bool>());
  EXPECT_FALSE(load_node("NO").as<bool>());
  EXPECT_FALSE(load_node("no").as<bool>());
  EXPECT_FALSE(load_node("0").as<bool>());
}

TEST_F(ParserNodeTest, AsBoolUnrecognizedDoomsAndReturnsFalse)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  const bool result = load_node("maybe").as<bool>();
  EXPECT_FALSE(result);
  EXPECT_GT(stk::get_doomed_count(), oldDoomedCount);
}

// -----------------------------------------------------------------------------
// Iteration
// -----------------------------------------------------------------------------

TEST_F(ParserNodeTest, IterateSequenceYieldsEachElement)
{
  Parser::Node node = load_node(R"yaml(
- 10
- 20
- 30
)yaml");
  std::vector<int> collected;
  for (auto && entry : node)
    collected.push_back(entry.as<int>());
  ASSERT_EQ(3u, collected.size());
  EXPECT_EQ(10, collected[0]);
  EXPECT_EQ(20, collected[1]);
  EXPECT_EQ(30, collected[2]);
}

// -----------------------------------------------------------------------------
// info()
// -----------------------------------------------------------------------------

TEST_F(ParserNodeTest, InfoReturnsNonEmptyDescription)
{
  Parser::Node node = load_node(R"yaml(
alpha: 1
)yaml");
  EXPECT_FALSE(node.info().empty());
}

// -----------------------------------------------------------------------------
// yaml_parse_error_message
// -----------------------------------------------------------------------------

std::string get_yaml_parse_error_message(const std::string & yamlText)
{
  try
  {
    YAML::Load(yamlText);
  }
  catch (const YAML::ParserException & e)
  {
    return Parser::yaml_parse_error_message(e);
  }
  return std::string();
}

std::string get_yaml_parse_exception_what(const std::string & yamlText)
{
  try
  {
    YAML::Load(yamlText);
  }
  catch (const YAML::ParserException & e)
  {
    return e.what();
  }
  return std::string();
}

constexpr const char * quoteHint = "put it in quotes";

TEST(YamlParseErrorMessage, UnquotedLeadingStarPattern_HasYamlErrorAndQuoteHint)
{
  const std::string yamlText = "parts:\n- name: *_dragon\n";
  const std::string message = get_yaml_parse_error_message(yamlText);

  EXPECT_EQ(0u, message.find(get_yaml_parse_exception_what(yamlText))); // keeps the yaml-cpp message and location
  EXPECT_NE(std::string::npos, message.find("line 2"));
  EXPECT_NE(std::string::npos, message.find(quoteHint));
}

TEST(YamlParseErrorMessage, UnquotedBareStar_HasQuoteHint)
{
  const std::string message = get_yaml_parse_error_message("subdomains:\n- *\n");
  EXPECT_NE(std::string::npos, message.find(quoteHint));
}

TEST(YamlParseErrorMessage, OtherYamlError_IsUnchanged)
{
  const std::string yamlText = "subdomains: [a, b\n";
  const std::string message = get_yaml_parse_error_message(yamlText);

  ASSERT_FALSE(message.empty());
  EXPECT_EQ(get_yaml_parse_exception_what(yamlText), message);
  EXPECT_EQ(std::string::npos, message.find(quoteHint));
}

// -----------------------------------------------------------------------------
// input_excerpt
// -----------------------------------------------------------------------------

const std::string tenLines = "l1\nl2\nl3\nl4\nl5\nl6\nl7\nl8\nl9\nl10\n";

TEST(InputExcerpt, ErrorInMiddle_ShowsContextNumberedWithMarkerAndCaret)
{
  // 0-based line 4 is "l5"; with 3 lines of context, lines 2..8 (1-based) are shown.
  const std::string excerpt = Parser::input_excerpt(tenLines, 4, 1, 3);
  const std::string expected =
      "  2 | l2\n"
      "  3 | l3\n"
      "  4 | l4\n"
      "> 5 | l5\n"
      "       ^\n"
      "  6 | l6\n"
      "  7 | l7\n"
      "  8 | l8\n";
  EXPECT_EQ(expected, excerpt);
}

TEST(InputExcerpt, LineNumbersAreRightAlignedToTheWidestShownNumber)
{
  // 0-based line 8 is "l9"; with 1 line of context, lines 8..10 are shown, so numbers are 2 wide.
  const std::string excerpt = Parser::input_excerpt(tenLines, 8, 0, 1);
  const std::string expected =
      "   8 | l8\n"
      ">  9 | l9\n"
      "       ^\n"
      "  10 | l10\n";
  EXPECT_EQ(expected, excerpt);
}

TEST(InputExcerpt, ErrorOnFirstLine_ContextIsClamped)
{
  const std::string excerpt = Parser::input_excerpt(tenLines, 0, 0, 2);
  const std::string expected =
      "> 1 | l1\n"
      "      ^\n"
      "  2 | l2\n"
      "  3 | l3\n";
  EXPECT_EQ(expected, excerpt);
}

TEST(InputExcerpt, ErrorOnLastLineWithoutTrailingNewline_ContextIsClamped)
{
  const std::string excerpt = Parser::input_excerpt("a\nb\nc", 2, 0, 3);
  const std::string expected =
      "  1 | a\n"
      "  2 | b\n"
      "> 3 | c\n"
      "      ^\n";
  EXPECT_EQ(expected, excerpt);
}

TEST(InputExcerpt, LineBeyondEndOfText_MarksLastLineWithoutCaret)
{
  const std::string excerpt = Parser::input_excerpt("a\nb\n", 7, 0, 1);
  const std::string expected =
      "  1 | a\n"
      "> 2 | b\n";
  EXPECT_EQ(expected, excerpt);
}

TEST(InputExcerpt, EmptyText_IsEmpty)
{
  EXPECT_TRUE(Parser::input_excerpt("", 0, 0).empty());
}

TEST(InputExcerpt, UnquotedLeadingStarPattern_ExcerptPointsAtStar)
{
  const std::string yamlText = "parts:\n- name: block_1\n- name: *_dragon\n  subdomains:\n  - a\n";
  try
  {
    YAML::Load(yamlText);
    FAIL() << "Expected a YAML parse error";
  }
  catch (const YAML::ParserException & e)
  {
    const std::string excerpt = Parser::input_excerpt(yamlText, e.mark.line, e.mark.column);
    const std::size_t errorLinePos = excerpt.find("> 3 | - name: *_dragon\n");
    ASSERT_NE(std::string::npos, errorLinePos) << excerpt;

    // The caret line follows the error line, with '^' at the position of the '*'.
    const std::size_t caretLineBegin = excerpt.find('\n', errorLinePos) + 1;
    const std::size_t caretLineEnd = excerpt.find('\n', caretLineBegin);
    const std::string caretLine = excerpt.substr(caretLineBegin, caretLineEnd-caretLineBegin);
    const std::size_t starPos = excerpt.find('*', errorLinePos) - errorLinePos;
    EXPECT_EQ(starPos, caretLine.find('^')) << excerpt;
  }
}

TEST(YamlParseErrorMessage, QuotedStarPatterns_Parse)
{
  EXPECT_EQ("*_dragon", YAML::Load("parts:\n- name: '*_dragon'\n")["parts"][0]["name"].as<std::string>());
  EXPECT_EQ("*", YAML::Load("subdomains:\n- \"*\"\n")["subdomains"][0].as<std::string>());
}

#else

TEST(ParserNodeTest, YamlSupportNotEnabled)
{
  GTEST_SKIP() << "KRINO_HAVE_YAML is not enabled; skipping YAML parser tests.";
}

#endif

} // namespace
} // namespace krino
