#include "gtest/gtest.h"
#include "stk_util/command_line/CommandLineParser.hpp"  // for CommandLineParser, CommandLinePar...
#include "stk_util/command_line/CommandLineParserUtils.hpp"
#include <exception>                                    // for exception
#include <string>                                       // for string
#include <vector>                                       // for vector

namespace {

class Args
{
public:
  Args(const std::vector<std::string> & strArgs)
    : m_stringArgs(strArgs),
      m_argc(m_stringArgs.size()),
      m_argv(strArgs.empty() ? nullptr : new char*[m_argc])
  {
    for (int i = 0; i < m_argc; ++i) {
      m_argv[i] = const_cast<char*>(m_stringArgs[i].c_str());
    }
  }

  ~Args()
  {
    delete [] m_argv;
  }

  int argc() { return m_argc; }
  char** argv() { return m_argv; }

private:
  const std::vector<std::string> m_stringArgs;
  int m_argc;
  char** m_argv;
};

TEST(UnitTestGetOption, get_command_line_option_null)
{
  Args args({});
  int defaultValue = -1;
  int result = stk::get_command_line_option(args.argc(), args.argv(), "foo", defaultValue);
  EXPECT_EQ(defaultValue, result);
}

TEST(UnitTestGetOption, get_command_line_option_bad_arg)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { return; }

  int myRank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  Args args({"exe", "--garbage-color", std::to_string(myRank)});
  int defaultValue = 0;

  testing::internal::CaptureStderr();
  int result = stk::get_command_line_option(args.argc(), args.argv(), "app-color", defaultValue);
  testing::internal::GetCapturedStderr();
  EXPECT_EQ(defaultValue, result);
}

TEST(UnitTestGetOption, get_command_line_option_no_value)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { return; }

  Args args({"exe", "--app-color"});
  int defaultValue = -1;

  testing::internal::CaptureStderr();
  EXPECT_THROW(stk::get_command_line_option(args.argc(), args.argv(), "app-color", defaultValue), std::runtime_error);
  testing::internal::GetCapturedStderr();
}

TEST(UnitTestGetOption, get_command_line_option_non_int_value)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { return; }

  Args args({"exe", "--app-color", "foo"});
  int defaultValue = -1;

  EXPECT_THROW(stk::get_command_line_option(args.argc(), args.argv(), "app-color", defaultValue), std::logic_error);
}

TEST(UnitTestGetOption, get_command_line_option)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { return; }

  int myRank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  Args args({"exe", "--my-color", std::to_string(myRank)});
  int defaultColor = -1;

  int expectedResult = myRank;
  int result = stk::get_command_line_option(args.argc(), args.argv(), "my-color", defaultColor);
  EXPECT_EQ(expectedResult, result);
}


bool messageContains(const std::string & errorMsg, const std::string & containedString) {
  return (errorMsg.find(containedString) != std::string::npos);
}

bool parse_command_line_without_error(stk::CommandLineParser & parser, Args & args)
{
  return parser.parse(args.argc(), const_cast<const char**>(args.argv())) == stk::CommandLineParser::ParseComplete;
}

bool parse_command_line_with_help(stk::CommandLineParser & parser, Args & args)
{
  testing::internal::CaptureStderr();

  const bool requestedHelpDuringParse = parser.parse(args.argc(), const_cast<const char**>(args.argv())) == stk::CommandLineParser::ParseHelpOnly;
  EXPECT_TRUE(requestedHelpDuringParse);

  const std::string output = testing::internal::GetCapturedStderr();
  const bool printedNoOutput = output.empty();
  EXPECT_EQ(printedNoOutput, true) << "Actual output: " << output;

  return requestedHelpDuringParse && printedNoOutput;
}

bool parse_command_line_with_version(stk::CommandLineParser & parser, Args & args)
{
  testing::internal::CaptureStderr();

  const bool requestedVersionDuringParse = parser.parse(args.argc(), const_cast<const char**>(args.argv())) == stk::CommandLineParser::ParseVersionOnly;
  EXPECT_TRUE(requestedVersionDuringParse);

  const std::string output = testing::internal::GetCapturedStderr();
  const bool printedNoOutput = output.empty();
  EXPECT_EQ(printedNoOutput, true) << "Actual output: " << output;

  return requestedVersionDuringParse && printedNoOutput;
}

bool parse_command_line_with_error(stk::CommandLineParser & parser, Args & args, const std::string & expectedErrorText)
{
  testing::internal::CaptureStderr();

  const bool hadErrorDuringParse = parser.parse(args.argc(), const_cast<const char**>(args.argv())) == stk::CommandLineParser::ParseError;
  EXPECT_TRUE(hadErrorDuringParse);

  const std::string errorMessage = testing::internal::GetCapturedStderr();
  const bool hadExpectedErrorMessage = messageContains(errorMessage, expectedErrorText);
  EXPECT_EQ(hadExpectedErrorMessage, true) << "Actual error message: " << errorMessage;

  return hadErrorDuringParse && hadExpectedErrorMessage;
}

//==============================================================================
TEST(CommandLineParser, oneFlag_notProvided_querySaysNotProvided)
{
  stk::CommandLineParser parser;
  parser.add_flag("flag,f", "One flag description");

  Args args({"exe"});
  EXPECT_TRUE(parse_command_line_without_error(parser, args));

  EXPECT_TRUE(parser.is_empty());
  EXPECT_FALSE(parser.is_option_parsed("flag"));
  EXPECT_FALSE(parser.is_option_provided("flag"));
}

TEST(CommandLineParser, oneFlag_provided_querySaysProvided)
{
  stk::CommandLineParser parser;
  parser.add_flag("flag,f", "One flag description");

  Args args({"exe", "--flag"});
  EXPECT_TRUE(parse_command_line_without_error(parser, args));

  EXPECT_FALSE(parser.is_empty());
  EXPECT_TRUE(parser.is_option_parsed("flag"));
  EXPECT_TRUE(parser.is_option_provided("flag"));
}

TEST(CommandLineParser, oneFlag_shortProvided_querySaysProvided)
{
  stk::CommandLineParser parser;
  parser.add_flag("flag,f", "One flag description");

  Args args({"exe", "-f"});
  EXPECT_TRUE(parse_command_line_without_error(parser, args));

  EXPECT_FALSE(parser.is_empty());
  EXPECT_TRUE(parser.is_option_parsed("flag"));
  EXPECT_TRUE(parser.is_option_provided("flag"));
}


//==============================================================================
TEST(CommandLineParser, oneRequiredPositional_notProvided_printsErrorDuringParse_getValueThrows)
{
  stk::CommandLineParser parser;
  parser.add_required_positional<std::string>({"positionalOpt", "p", "One required positional description"});

  Args args({"exe"});
  EXPECT_TRUE(parse_command_line_with_error(parser, args, ""));

  EXPECT_TRUE(parser.is_empty());
  EXPECT_FALSE(parser.is_option_parsed("positionalOpt"));
  EXPECT_FALSE(parser.is_option_provided("positionalOpt"));

  EXPECT_THROW(parser.get_option_value<std::string>("positionalOpt"), std::logic_error);
}

TEST(CommandLineParser, oneRequiredPositional_valueProvided_getValueReturnsProvided)
{
  stk::CommandLineParser parser;
  parser.add_required_positional<std::string>({"positionalOpt", "p", "One required positional description"});

  Args args({"exe", "positionalValue"});
  EXPECT_TRUE(parse_command_line_without_error(parser, args));

  EXPECT_FALSE(parser.is_empty());
  EXPECT_FALSE(parser.is_option_parsed("positionalOpt"));
  EXPECT_TRUE(parser.is_option_provided("positionalOpt"));

  EXPECT_EQ(parser.get_option_value<std::string>("positionalOpt"), "positionalValue");
}

TEST(CommandLineParser, oneRequiredPositional_optionAndValueProvided_getValueReturnsProvided)
{
  stk::CommandLineParser parser;
  parser.add_required_positional<std::string>({"positionalOpt", "p", "One required positional description"});

  Args args({"exe", "--positionalOpt=positionalValue"});
  EXPECT_TRUE(parse_command_line_without_error(parser, args));

  EXPECT_FALSE(parser.is_empty());
  EXPECT_TRUE(parser.is_option_parsed("positionalOpt"));
  EXPECT_TRUE(parser.is_option_provided("positionalOpt"));

  EXPECT_EQ(parser.get_option_value<std::string>("positionalOpt"), "positionalValue");
}


//==============================================================================
TEST(CommandLineParser, oneOptionalPositional_notProvided_getValueReturnsDefault)
{
  stk::CommandLineParser parser;
  parser.add_optional_positional<std::string>({"positionalOpt", "p", "One required positional description"}, "defaultValue");

  Args args({"exe"});
  EXPECT_TRUE(parse_command_line_without_error(parser, args));

  EXPECT_FALSE(parser.is_empty());
  EXPECT_FALSE(parser.is_option_parsed("positionalOpt"));
  EXPECT_TRUE(parser.is_option_provided("positionalOpt"));

  EXPECT_EQ(parser.get_option_value<std::string>("positionalOpt"), "defaultValue");
}

TEST(CommandLineParser, oneOptionalPositional_valueProvided_getValueReturnsProvided)
{
  stk::CommandLineParser parser;
  parser.add_optional_positional<std::string>({"positionalOpt", "p", "One required positional description"}, "defaultValue");

  Args args({"exe", "providedValue"});
  EXPECT_TRUE(parse_command_line_without_error(parser, args));

  EXPECT_FALSE(parser.is_empty());
  EXPECT_FALSE(parser.is_option_parsed("positionalOpt"));
  EXPECT_TRUE(parser.is_option_provided("positionalOpt"));

  EXPECT_EQ(parser.get_option_value<std::string>("positionalOpt"), "providedValue");
}

TEST(CommandLineParser, oneOptionalPositional_optionAndValueProvided_getValueReturnsProvided)
{
  stk::CommandLineParser parser;
  parser.add_optional_positional<std::string>({"positionalOpt", "p", "One required positional description"}, "defaultValue");

  Args args({"exe", "--positionalOpt=providedValue"});
  EXPECT_TRUE(parse_command_line_without_error(parser, args));

  EXPECT_FALSE(parser.is_empty());
  EXPECT_TRUE(parser.is_option_parsed("positionalOpt"));
  EXPECT_TRUE(parser.is_option_provided("positionalOpt"));

  EXPECT_EQ(parser.get_option_value<std::string>("positionalOpt"), "providedValue");
}

TEST(CommandLineParser, oneOptionalPositional_disallowUnrecognized_unrecognizedOptionInPlaceOfPositional_errorDuringParse)
{
  stk::CommandLineParser parser;
  parser.add_optional_positional<std::string>({"positionalOpt", "p", "One positional description"}, "defaultValue");
  parser.disallow_unrecognized();

  Args args({"exe", "--unrecognizedOpt"});
  EXPECT_TRUE(parse_command_line_with_error(parser, args, "Unrecognized option: '--unrecognizedOpt'"));

  EXPECT_FALSE(parser.is_empty());
  EXPECT_FALSE(parser.is_option_parsed("positionalOpt"));
  EXPECT_TRUE(parser.is_option_provided("positionalOpt"));

  EXPECT_EQ(parser.get_option_value<std::string>("positionalOpt"), "defaultValue");
}


//==============================================================================
TEST(CommandLineParser, oneRequired_notProvided_printErrorDuringParse_getValueThrows)
{
  stk::CommandLineParser parser;
  parser.add_required<std::string>({"requiredOpt", "r", "One required option description"});

  Args args({"exe"});
  EXPECT_TRUE(parse_command_line_with_error(parser, args, "Required option '--requiredOpt' not found"));

  EXPECT_TRUE(parser.is_empty());
  EXPECT_FALSE(parser.is_option_parsed("requiredOpt"));
  EXPECT_FALSE(parser.is_option_provided("requiredOpt"));

  EXPECT_THROW(parser.get_option_value<std::string>("requiredOpt"), std::exception);
}

TEST(CommandLineParser, oneRequired_shortOptionAndValueProvided_getValueReturnsProvided)
{
  stk::CommandLineParser parser;
  parser.add_required<std::string>({"requiredOpt", "r", "One required description"});

  Args args({"exe", "-r", "providedValue"});
  EXPECT_TRUE(parse_command_line_without_error(parser, args));

  EXPECT_FALSE(parser.is_empty());
  EXPECT_TRUE(parser.is_option_parsed("requiredOpt"));
  EXPECT_TRUE(parser.is_option_provided("requiredOpt"));

  EXPECT_EQ(parser.get_option_value<std::string>("requiredOpt"), "providedValue");
}

TEST(CommandLineParser, oneRequired_optionAndValueProvided_getValueReturnsProvided)
{
  stk::CommandLineParser parser;
  parser.add_required<std::string>({"requiredOpt", "r", "One required description"});

  Args args({"exe", "--requiredOpt", "providedValue"});
  EXPECT_TRUE(parse_command_line_without_error(parser, args));

  EXPECT_FALSE(parser.is_empty());
  EXPECT_TRUE(parser.is_option_parsed("requiredOpt"));
  EXPECT_TRUE(parser.is_option_provided("requiredOpt"));

  EXPECT_EQ(parser.get_option_value<std::string>("requiredOpt"), "providedValue");
}

TEST(CommandLineParser, oneRequired_optionAndNoValueProvided_errorDuringParse)
{
  stk::CommandLineParser parser;
  parser.add_required<std::string>({"requiredOpt", "r", "One required description"});

  Args args({"exe", "--requiredOpt"});
  EXPECT_TRUE(parse_command_line_with_error(parser, args, "Missing value for option --requiredOpt"));

  EXPECT_TRUE(parser.is_empty());
  EXPECT_FALSE(parser.is_option_parsed("requiredOpt"));
  EXPECT_FALSE(parser.is_option_provided("requiredOpt"));

  EXPECT_THROW(parser.get_option_value<std::string>("requiredOpt"), std::logic_error);
}


//==============================================================================
TEST(CommandLineParser, oneOptional_notProvided_getValueReturnsDefault)
{
  stk::CommandLineParser parser;
  parser.add_optional<std::string>({"optionalOpt" ,"o", "One optional option description"}, "defaultValue");

  Args args({"exe"});
  EXPECT_TRUE(parse_command_line_without_error(parser, args));

  EXPECT_FALSE(parser.is_empty());
  EXPECT_FALSE(parser.is_option_parsed("optionalOpt"));
  EXPECT_TRUE(parser.is_option_provided("optionalOpt"));

  EXPECT_EQ(parser.get_option_value<std::string>("optionalOpt"), "defaultValue");
}

TEST(CommandLineParser, oneOptionalWithoutAbbreviation_notProvided_getValueReturnsDefault)
{
  stk::CommandLineParser parser;
  parser.add_optional<std::string>("optionalOpt", "One optional option description", "defaultValue");

  Args args({"exe"});
  EXPECT_TRUE(parse_command_line_without_error(parser, args));

  EXPECT_FALSE(parser.is_empty());
  EXPECT_FALSE(parser.is_option_parsed("optionalOpt"));
  EXPECT_TRUE(parser.is_option_provided("optionalOpt"));

  EXPECT_EQ(parser.get_option_value<std::string>("optionalOpt"), "defaultValue");
}

TEST(CommandLineParser, oneOptional_notProvided_throwIfGettingWrongTypeFromDefault)
{
  stk::CommandLineParser parser;
  parser.add_optional<std::string>({"optionalOpt", "r", "One required description"}, "defaultValue");

  Args args({"exe"});
  EXPECT_TRUE(parse_command_line_without_error(parser, args));

  EXPECT_FALSE(parser.is_empty());
  EXPECT_FALSE(parser.is_option_parsed("optionalOpt"));
  EXPECT_TRUE(parser.is_option_provided("optionalOpt"));

  EXPECT_THROW(parser.get_option_value<double>("optionalOpt"), std::logic_error);
}

TEST(CommandLineParser, oneOptional_povided_getValueReturnsProvided)
{
  stk::CommandLineParser parser;
  parser.add_optional<std::string>({"optionalOpt", "o", "One optional option description"}, "defaultValue");

  Args args({"exe", "--optionalOpt", "providedValue"});
  EXPECT_TRUE(parse_command_line_without_error(parser, args));

  EXPECT_FALSE(parser.is_empty());
  EXPECT_TRUE(parser.is_option_parsed("optionalOpt"));
  EXPECT_TRUE(parser.is_option_provided("optionalOpt"));

  EXPECT_EQ(parser.get_option_value<std::string>("optionalOpt"), "providedValue");
}

TEST(CommandLineParser, oneOptional_shortPovided_getValueReturnsProvided)
{
  stk::CommandLineParser parser;
  parser.add_optional<std::string>({"optionalOpt", "o", "One optional option description"}, "defaultValue");

  Args args({"exe", "-o", "providedValue"});
  EXPECT_TRUE(parse_command_line_without_error(parser, args));

  EXPECT_FALSE(parser.is_empty());
  EXPECT_TRUE(parser.is_option_parsed("optionalOpt"));
  EXPECT_TRUE(parser.is_option_provided("optionalOpt"));

  EXPECT_EQ(parser.get_option_value<std::string>("optionalOpt"), "providedValue");
}

TEST(CommandLineParser, oneOptional_optionWithEqualsAndNoValueProvided_errorDuringParse)
{
  stk::CommandLineParser parser;
  parser.add_optional<std::string>({"optionalOpt", "r", "One required description"}, "defaultValue");

  Args args({"exe", "--optionalOpt="});
  EXPECT_TRUE(parse_command_line_with_error(parser, args, "Missing value for option --optionalOpt"));

  EXPECT_FALSE(parser.is_empty());
  EXPECT_FALSE(parser.is_option_parsed("optionalOpt"));
  EXPECT_TRUE(parser.is_option_provided("optionalOpt"));

  EXPECT_EQ(parser.get_option_value<std::string>("optionalOpt"), "defaultValue");
}

TEST(CommandLineParser, oneOptional_optionAndNoValueProvided_errorDuringParse)
{
  stk::CommandLineParser parser;
  parser.add_optional<std::string>({"optionalOpt", "r", "One required description"}, "defaultValue");

  Args args({"exe", "--optionalOpt"});
  EXPECT_TRUE(parse_command_line_with_error(parser, args, "Missing value for option --optionalOpt"));

  EXPECT_FALSE(parser.is_empty());
  EXPECT_FALSE(parser.is_option_parsed("optionalOpt"));
  EXPECT_TRUE(parser.is_option_provided("optionalOpt"));

  EXPECT_EQ(parser.get_option_value<std::string>("optionalOpt"), "defaultValue");
}

TEST(CommandLineParser, oneOptional_noDefault_notProvided_getValueThrows)
{
  stk::CommandLineParser parser;
  parser.add_optional<std::string>({"optionalOpt" ,"o", "One optional option description"});

  Args args({"exe"});
  EXPECT_TRUE(parse_command_line_without_error(parser, args));

  EXPECT_TRUE(parser.is_empty());
  EXPECT_FALSE(parser.is_option_parsed("optionalOpt"));
  EXPECT_FALSE(parser.is_option_provided("optionalOpt"));

  EXPECT_THROW(parser.get_option_value<std::string>("optionalOpt"), std::logic_error);
}

TEST(CommandLineParser, oneOptional_noDefault_optionProvided_getValueReturnsProvided)
{
  stk::CommandLineParser parser;
  parser.add_optional<std::string>({"optionalOpt" ,"o", "One optional option description"});

  Args args({"exe", "--optionalOpt", "providedValue"});
  EXPECT_TRUE(parse_command_line_without_error(parser, args));

  EXPECT_FALSE(parser.is_empty());
  EXPECT_TRUE(parser.is_option_parsed("optionalOpt"));
  EXPECT_TRUE(parser.is_option_provided("optionalOpt"));

  EXPECT_EQ(parser.get_option_value<std::string>("optionalOpt"), "providedValue");
}

//==============================================================================
TEST(CommandLineParser, oneOptionalImplicit_provided_getValueReturnsProvided)
{
  stk::CommandLineParser parser;
  parser.add_optional_implicit<std::string>({"implicitOpt", "i", "One optional implicit description"}, "defaultValue");

  Args args({"exe", "--implicitOpt", "providedValue"});
  EXPECT_TRUE(parse_command_line_without_error(parser, args));

  EXPECT_FALSE(parser.is_empty());
  EXPECT_TRUE(parser.is_option_parsed("implicitOpt"));
  EXPECT_TRUE(parser.is_option_provided("implicitOpt"));

  EXPECT_EQ(parser.get_option_value<std::string>("implicitOpt"), "providedValue");
}

TEST(CommandLineParser, oneOptionalImplicit_providedWithoutValue_getValueReturnsDefault)
{
  stk::CommandLineParser parser;
  parser.add_optional_implicit<std::string>({"implicitOpt", "i", "One optional implicit description"}, "defaultValue");

  Args args({"exe", "--implicitOpt"});
  EXPECT_TRUE(parse_command_line_without_error(parser, args));

  EXPECT_FALSE(parser.is_empty());
  EXPECT_TRUE(parser.is_option_parsed("implicitOpt"));
  EXPECT_TRUE(parser.is_option_provided("implicitOpt"));

  EXPECT_EQ(parser.get_option_value<std::string>("implicitOpt"), "defaultValue");
}

TEST(CommandLineParser, oneOptionalImplicit_notProvided_getValueThrows)
{
  stk::CommandLineParser parser;
  parser.add_optional_implicit<std::string>({"implicitOpt", "i", "One optional implicit description"}, "defaultValue");

  Args args({"exe"});
  EXPECT_TRUE(parse_command_line_without_error(parser, args));

  EXPECT_TRUE(parser.is_empty());
  EXPECT_FALSE(parser.is_option_parsed("implicitOpt"));
  EXPECT_FALSE(parser.is_option_provided("implicitOpt"));

  EXPECT_THROW(parser.get_option_value<std::string>("implicitOpt"), std::logic_error);
}


//==============================================================================
TEST(CommandLineParser, twoRequired_bothOptionsAndValuesProvided_getValueReturnsProvided)
{
  stk::CommandLineParser parser;
  parser.add_required<std::string>({"requiredOpt1", "r", "First required description"});
  parser.add_required<std::string>({"requiredOpt2", "o", "Second required description"});

  Args args({"exe", "-r", "firstProvidedValue", "--requiredOpt2=secondProvidedValue"});
  EXPECT_TRUE(parse_command_line_without_error(parser, args));

  EXPECT_FALSE(parser.is_empty());
  EXPECT_TRUE(parser.is_option_parsed("requiredOpt1"));
  EXPECT_TRUE(parser.is_option_parsed("requiredOpt2"));
  EXPECT_TRUE(parser.is_option_provided("requiredOpt1"));
  EXPECT_TRUE(parser.is_option_provided("requiredOpt2"));

  EXPECT_EQ(parser.get_option_value<std::string>("requiredOpt1"), "firstProvidedValue");
  EXPECT_EQ(parser.get_option_value<std::string>("requiredOpt2"), "secondProvidedValue");
}


TEST(CommandLineParser, oneRequiredPositional_disallowUnrecognized_unrecognizedOptionInPlaceOfPositional_errorDuringParse)
{
  stk::CommandLineParser parser;
  parser.add_required_positional<std::string>({"positionalOpt", "p", "One positional description"});
  parser.disallow_unrecognized();

  Args args({"exe", "--unrecognizedOpt"});
  EXPECT_TRUE(parse_command_line_with_error(parser, args, "Required option '--positionalOpt' not found"));

  EXPECT_TRUE(parser.is_empty());
  EXPECT_FALSE(parser.is_option_parsed("positionalOpt"));
  EXPECT_FALSE(parser.is_option_provided("positionalOpt"));

  EXPECT_THROW(parser.get_option_value<std::string>("positionalOpt"), std::exception);
}

TEST(CommandLineParser, oneRequiredPositionalAndOneOptional_disallowUnrecognized_bothOptionsProvidedWithUnrecognized_errorDuringParse)
{
  stk::CommandLineParser parser;
  parser.add_required_positional<std::string>({"positionalOpt", "p", "One positional description"});
  parser.add_optional<std::string>({"optionalOpt", "o", "One optional description"}, "defaultValue");
  parser.disallow_unrecognized();

  Args args({"exe", "positionalValue", "--unrecognizedOpt", "unrecognizedValue", "--optionalOpt", "optionalValue"});
  EXPECT_TRUE(parse_command_line_with_error(parser, args, "Unrecognized option: '--unrecognizedOpt'"));

  EXPECT_FALSE(parser.is_empty());
  EXPECT_FALSE(parser.is_option_parsed("positionalOpt"));
  EXPECT_TRUE(parser.is_option_parsed("optionalOpt"));
  EXPECT_TRUE(parser.is_option_provided("positionalOpt"));
  EXPECT_TRUE(parser.is_option_provided("optionalOpt"));

  EXPECT_EQ(parser.get_option_value<std::string>("positionalOpt"), "positionalValue");
  EXPECT_EQ(parser.get_option_value<std::string>("optionalOpt"), "optionalValue");
}

//==============================================================================
TEST(CommandLineParser, noDefinedArgs_getUsage_displaysDefaultArgs)
{
  stk::CommandLineParser parser;
  std::string usage = parser.get_usage();

  EXPECT_TRUE(messageContains(usage, "Options"));
  EXPECT_TRUE(messageContains(usage, "--help,-h"));
  EXPECT_TRUE(messageContains(usage, "--version,-v"));
}

TEST(CommandLineParser, noDefinedArgs_specifiedPreamble_getUsage_displaysDefaultArgsWithhPreamble)
{
  stk::CommandLineParser parser("myWackyProgram usage:");
  std::string usage = parser.get_usage();

  EXPECT_FALSE(messageContains(usage, "Options"));
  EXPECT_TRUE(messageContains(usage, "myWackyProgram usage:"));
  EXPECT_TRUE(messageContains(usage, "--help,-h"));
  EXPECT_TRUE(messageContains(usage, "--version,-v"));
}

TEST(CommandLineParser, emptyCommandLine_parseNothing)
{
  stk::CommandLineParser parser;

  Args args({"exe"});
  EXPECT_TRUE(parse_command_line_without_error(parser, args));

  EXPECT_TRUE(parser.is_empty());
}

TEST(CommandLineParser, emptyCommandLine_parseNothing_getValueThrows)
{
  stk::CommandLineParser parser;

  Args args({"exe"});
  EXPECT_TRUE(parse_command_line_without_error(parser, args));

  EXPECT_TRUE(parser.is_empty());
  EXPECT_THROW(parser.get_option_value<std::string>("undefinedOpt"), std::exception);
}

TEST(CommandLineParser, noArgsDefined_askHelp_returnsHelpStatus)
{
  stk::CommandLineParser parser;

  Args args({"exe", "--help"});
  EXPECT_TRUE(parse_command_line_with_help(parser, args));
}

TEST(CommandLineParser, noArgsDefined_askShortHelp_returnsHelpStatus)
{
  stk::CommandLineParser parser;

  Args args({"exe", "-h"});
  EXPECT_TRUE(parse_command_line_with_help(parser, args));
}

TEST(CommandLineParser, oneRequired_notProvided_askHelp_returnsHelpStatus)
{
  stk::CommandLineParser parser;
  parser.add_required<std::string>({"requiredOpt", "r", "One required option description"});

  Args args({"exe", "--help"});
  EXPECT_TRUE(parse_command_line_with_help(parser, args));
}

TEST(CommandLineParser, oneRequired_notProvided_askShortHelp_returnsHelpStatus)
{
  stk::CommandLineParser parser;
  parser.add_required<std::string>({"requiredOpt", "r", "One required option description"});

  Args args({"exe", "-h"});
  EXPECT_TRUE(parse_command_line_with_help(parser, args));
}

TEST(CommandLineParser, noArgsDefined_askVersion_returnsVersionStatus)
{
  stk::CommandLineParser parser;

  Args args({"exe", "--version"});
  EXPECT_TRUE(parse_command_line_with_version(parser, args));
}

TEST(CommandLineParser, noArgsDefined_askShortVersion_returnsVersionStatus)
{
  stk::CommandLineParser parser;

  Args args({"exe", "-v"});
  EXPECT_TRUE(parse_command_line_with_version(parser, args));
}

TEST(CommandLineParser, oneRequired_askVersion_returnsVersionStatus)
{
  stk::CommandLineParser parser;
  parser.add_required<std::string>({"requiredOpt", "r", "One required option description"});

  Args args({"exe", "--version"});
  EXPECT_TRUE(parse_command_line_with_version(parser, args));
}

TEST(CommandLineParser, oneRequired_askShortVersion_returnsVersionStatus)
{
  stk::CommandLineParser parser;
  parser.add_required<std::string>({"requiredOpt", "r", "One required option description"});

  Args args({"exe", "-v"});
  EXPECT_TRUE(parse_command_line_with_version(parser, args));
}


//==============================================================================
}
