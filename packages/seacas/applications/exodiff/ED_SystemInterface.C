// Copyright(C) 1999-2024 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "ED_SystemInterface.h"

#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <vector>

#include <climits>
#include <cstdlib>
#include <cstring>

#include "ED_Version.h"
#include "copy_string_cpp.h"
#include "copyright.h"
#include "fmt/ostream.h"
#include "stringx.h"
#include <SL_tokenize.h>

namespace {
  [[noreturn]] void Parse_Die(std::string &sline)
  {
    chop_whitespace(sline);
    Error(fmt::format("parsing input file, currently at \"{}\".\n", sline));
  }

  std::string Parse_Variables(std::string xline, std::ifstream &cmd_file, bool &all_flag,
                              Tolerance &def_tol, std::vector<std::string> &names,
                              std::vector<Tolerance> &toler);

  bool str_equal(const std::string &s1, const std::string &s2)
  {
    return (s1.size() == s2.size()) &&
           std::equal(s1.begin(), s1.end(), s2.begin(),
                      [](char a, char b) { return std::tolower(a) == std::tolower(b); });
  }

  void file_help();
  void tolerance_help();

  double To_Double(const std::string &str_val)
  {
    SMART_ASSERT(!str_val.empty());

    double val = 0;
    try {
      val = std::stod(str_val);
    }
    catch (...) {
      Error(fmt::format(" Problem converting the string '{}'"
                        " to a double value while parsing tolerance.  Aborting...\n",
                        str_val));
    }

    if (val < 0.0) {
      Error(fmt::format(" Parsed a negative value \"{}\".  Aborting...\n", val));
    }
    return val;
  }

  int File_Exists(const std::string &fname)
  {
    if (fname.empty()) {
      return 0;
    }
    std::ifstream file_check(fname, std::ios::in);
    if (file_check.fail()) {
      return 0;
    }
    file_check.close();
    return 1;
  }

  void Parse_Steps_Option(const std::string &option, int &start, int &stop, int &increment)
  {
    //: The defined formats for the count attribute are:<br>
    //:  <ul>
    //:    <li><missing> -- default -- 1 <= count <= oo  (all steps)</li>
    //:    <li>"X"                  -- X <= count <= X  (just step X)</li>
    //:    <li>"X:Y"                -- X to Y by 1</li>
    //:    <li>"X:"                 -- X to oo by 1</li>
    //:    <li>":Y"                 -- 1 to Y by 1</li>
    //:    <li>"::Z"                -- 1 to oo by Z</li>
    //:  </ul>
    //: The count and step must always be >= 0

    // Break into tokens separated by ":"

    // Default is given in constructor above...

    if (option == "last" || option == "LAST") {
      start = -1;
      return;
    }
    const char *tokens = option.c_str();
    if (strchr(tokens, ':') != nullptr) {
      // The string contains a separator

      int vals[3];
      vals[0] = start;
      vals[1] = stop;
      vals[2] = increment;

      int j = 0;
      for (auto &val : vals) {
        // Parse 'i'th field
        char tmp_str[128];
        ;
        int k = 0;

        while (tokens[j] != '\0' && tokens[j] != ':') {
          tmp_str[k++] = tokens[j++];
        }

        tmp_str[k] = '\0';
        if (strlen(tmp_str) > 0) {
          val = strtol(tmp_str, nullptr, 0);
        }

        if (tokens[j++] == '\0') {
          break; // Reached end of string
        }
      }
      start     = vals[0];
      stop      = vals[1];
      increment = vals[2];
    }
    else {
      // Does not contain a separator, min == max
      start = stop = strtol(tokens, nullptr, 0);
    }
  }

  void Check_Parsed_Names(const std::vector<std::string> &names, bool &all_flag)
  {
    int num_include = 0;
    int num_exclude = 0;
    for (const auto &name : names) {
      SMART_ASSERT(!name.empty());
      if (name[0] == '!') {
        ++num_exclude;
      }
      else {
        ++num_include;
      }
    }
    if (!all_flag && num_include > 0 && num_exclude > 0) {
      Error(fmt::format("Parsing error: Cannot specify both "
                        "variables to include and exclude without using the "
                        "'(all)' specifier.  Aborting...\n"));
    }
    if (num_include == 0 && num_exclude > 0) {
      all_flag = true;
    }
  }

  void parseExcludeTimes(std::string exclude_arg, std::vector<int> &exclude_steps)
  {
    std::string arg_copy = exclude_arg;

    int num_excluded_steps = 0;

    // first pass just counts the number of excluded time steps:

    std::string tok = extract_token(exclude_arg, ",");
    while (!tok.empty()) {
      std::string subtok = extract_token(tok, "-");
      SMART_ASSERT(!subtok.empty());

      errno     = 0;
      int ival1 = std::stoi(subtok);
      SMART_ASSERT(errno == 0);

      if (ival1 < 1) {
        Error(fmt::format("parsing exclusion times from command "
                          "line .. value was less than 1\n"));
      }

      ++num_excluded_steps;

      subtok = extract_token(tok, "-");
      if (!subtok.empty()) {
        errno     = 0;
        int ival2 = std::stoi(subtok);
        SMART_ASSERT(errno == 0);

        if (ival2 < 1) {
          Error(fmt::format("parsing exclusion times from command "
                            "line .. value was less than 1\n"));
        }

        if (ival1 < ival2) {
          for (int i = ival1 + 1; i <= ival2; ++i) {
            ++num_excluded_steps;
          }
        }
        else if (ival1 > ival2) {
          Error(fmt::format("parsing exclusion times from command "
                            "line .. first value in a range was greater than the "
                            "second.\n"));
        }
      }

      tok = extract_token(exclude_arg, ",");
    }

    if (num_excluded_steps > 0) {
      exclude_steps.resize(num_excluded_steps);

      // second pass collects the excluded time steps

      exclude_arg        = std::move(arg_copy);
      num_excluded_steps = 0;

      tok = extract_token(exclude_arg, ",");
      while (!tok.empty()) {
        std::string subtok = extract_token(tok, "-");
        SMART_ASSERT(!subtok.empty());

        errno     = 0;
        int ival1 = std::stoi(subtok);
        SMART_ASSERT(errno == 0);

        exclude_steps[num_excluded_steps++] = ival1;

        subtok = extract_token(tok, "-");
        if (!subtok.empty()) {
          errno     = 0;
          int ival2 = std::stoi(subtok);
          SMART_ASSERT(errno == 0);

          for (int i = ival1 + 1; i <= ival2; ++i) {
            exclude_steps[num_excluded_steps++] = i;
          }
        }

        tok = extract_token(exclude_arg, ",");
      }
    }
  }
} // namespace

SystemInterface::SystemInterface() { enroll_options(); }

void SystemInterface::show_version()
{
  fmt::print("EXODIFF\t(Version: {}) Modified: {}\n", version, verdate);
}

void SystemInterface::enroll_options()
{
  options_.usage("[options] file1.exo file2.exo [diffile.exo]                   \n"
                 "\tor:  exodiff -summary <file.exo> (create variable summary)  \n"
                 "\tor:  exodiff [-help] [tolerance|file] (usage)               \n"
                 "\tor:  exodiff [-version]                                     \n");

  options_.enroll(
      "help", GetLongOption::OptionalValue,
      "Print this summary and exit.\n"
      "\t\tEnter \"-help file\" for the syntax of the command file\n"
      "\t\t      \"-help tolerance\" for information on the supported tolerance options.",
      nullptr, "usage");

  options_.enroll("Help", GetLongOption::NoValue, "Print this summary and exit.", nullptr);
  options_.enroll("file", GetLongOption::MandatoryValue,
                  "Use the given file to specify the variables to be considered and to\n"
                  "\t\twhat tolerances. Enter \"-help file\" for the syntax of the command file",
                  nullptr);
  options_.enroll("summary", GetLongOption::NoValue,
                  "Produce a summary in exodiff input format.\n"
                  "\t\tThis will create output with max/min statistics on the data in the format\n"
                  "\t\tof an exodiff input file.",
                  nullptr, nullptr, true);

  // Tolerance options...
  options_.enroll("tolerance", GetLongOption::MandatoryValue,
                  "Overrides the default tolerance of 1.0E-6.", "1.0E-6");

  options_.enroll("Floor", GetLongOption::MandatoryValue,
                  "Overrides the default floor tolerance of 0.0.", "0.0");

  options_.enroll("absolute", GetLongOption::NoValue,
                  "Default tolerance is absolute difference. |a-b| > tolerance", nullptr);
  options_.enroll("relative", GetLongOption::NoValue,
                  "Default tolerance is relative difference. |a-b| > max(|a|,|b|)*tolerance",
                  nullptr);
  options_.enroll("combined", GetLongOption::NoValue,
                  "Default tolerance is combined difference. (-help tolerance for info)", nullptr);
  options_.enroll("ulps_float", GetLongOption::NoValue,
                  "Default tolerance if number of ulps (units last position) of difference\n"
                  "\t\twhen values converted to floats.",
                  nullptr);
  options_.enroll("ulps_double", GetLongOption::NoValue,
                  "Default tolerance is number of ulps (units last position) of difference.",
                  nullptr);
  options_.enroll("eigen_absolute", GetLongOption::NoValue,
                  "Default tolerance is absolute differences of the absolute value of the values.",
                  nullptr);
  options_.enroll("eigen_relative", GetLongOption::NoValue,
                  "Default tolerance is relative differences of the absolute value of the values.",
                  nullptr);
  options_.enroll("eigen_combined", GetLongOption::NoValue,
                  "Default tolerance is combined differences of the absolute value of the values.",
                  nullptr);
  options_.enroll("ignore", GetLongOption::NoValue,
                  "Default tolerance is ignored (turn off all checking by default).", nullptr);
  options_.enroll("coordinate_tolerance", GetLongOption::MandatoryValue,
                  "Overrides the default coordinate comparison tolerance of 1.0E-6.", "1.0E-6",
                  nullptr, true);

  options_.enroll("pedantic", GetLongOption::NoValue, "Be more picky about what is a difference.",
                  nullptr);
  options_.enroll("quiet", GetLongOption::NoValue,
                  "Quiet.  Only errors will be sent to stdout.  Comparison mode will echo\n"
                  "\t\t\"exodiff: Files are the same.\" or \"exodiff: Files are different.\"",
                  nullptr);
  options_.enroll("show_all_diffs", GetLongOption::NoValue,
                  "Show all differences for all variables, not just the maximum.", nullptr, nullptr,
                  true);

  options_.enroll("ignore_steps", GetLongOption::NoValue,
                  "Don't compare any transient data; compare mesh only.", nullptr);
  options_.enroll(
      "x", GetLongOption::MandatoryValue,
      "Exclude time steps.  Does not consider the time steps given in the list of integers.\n"
      "\t\tThe format is comma-separated and ranged integers (with no spaces), such as "
      "\"1,5-9,28\".\n"
      "\t\tThe first time step is the number '1'.",
      nullptr);
  options_.enroll(
      "exclude", GetLongOption::MandatoryValue,
      "Exclude time steps.  Does not consider the time steps given in the list of integers.\n"
      "\t\tThe format is comma-separated and ranged integers (with no spaces), such as "
      "\"1,5-9,28\".\n"
      "\t\tThe first time step is the number '1'.",
      nullptr);
  options_.enroll("steps", GetLongOption::MandatoryValue,
                  "Specify subset of steps to consider. Syntax is beg:end:increment,\n"
                  "\t\tEnter '-steps last' for just the last step. If only beg set, end=beg",
                  nullptr);
  options_.enroll("explicit", GetLongOption::MandatoryValue,
                  "Specify an explicit match of a step on database 1 with a step on database 2.\n"
                  "\t\tSyntax is '-explicit db1_step:db2_step' where 'db*_step' is either\n"
                  "\t\tthe 1-based step number or 'last' for the last step on the database.\n"
                  "\t\tExample: '-explicit 42:last' to match step 42 on database 1 with last step "
                  "on database 2",
                  nullptr);
  options_.enroll("TimeStepOffset", GetLongOption::MandatoryValue,
                  "Timestep 'x+offset' in first file matches timestep 'x' in second file.",
                  nullptr);
  options_.enroll("TA", GetLongOption::NoValue,
                  "Automatic determination of timestep offset -- end at same step.", nullptr);
  options_.enroll(
      "TM", GetLongOption::NoValue,
      "Automatic determination of timestep offset -- closest match to first step on file2.",
      nullptr);
  options_.enroll("interpolate", GetLongOption::NoValue,
                  "Interpolate times on file2 to match times on file1.", nullptr);
  options_.enroll(
      "final_time_tolerance", GetLongOption::MandatoryValue,
      "Tolerance on matching of final times on database when interpolate option specified.\n"
      "\t\tIf final times do not match within this tolerance, files are different.",
      nullptr, nullptr, true);

  options_.enroll("time_scale", GetLongOption::MandatoryValue,
                  "Scale the time values on the input database by the specified value.", nullptr);
  options_.enroll(
      "time_offset", GetLongOption::MandatoryValue,
      "Offset the (possibly scaled) time values on the input database by the specified value.",
      nullptr, nullptr, true);

  options_.enroll("map", GetLongOption::NoValue,
                  "Invokes a matching algorithm to create a mapping between the\n"
                  "\t\tnodes and elements of the two files.  The topology must still be\n"
                  "\t\tthe same (within tolerance), but can be ordered differently.",
                  nullptr);
  options_.enroll("dumpmap", GetLongOption::NoValue,
                  "If the -map switch used, print the resulting node and element maps.", nullptr);
  options_.enroll("partial", GetLongOption::NoValue,
                  "Invokes a matching algorithm similar to the -m option.  However \n"
                  "\t\tthis option ignores unmatched nodes and elements.  This allows \n"
                  "\t\tcomparison of files that only partially overlap.",
                  nullptr);
  options_.enroll("show_unmatched", GetLongOption::NoValue,
                  "If the -partial switch used, print the elements that did not match.", nullptr);
  options_.enroll("ignore_dups", GetLongOption::NoValue,
                  "If two elements/nodes are in the same location in map or partial\n"
                  "                  map case, just return first match instead of aborting.",
                  nullptr);
  options_.enroll("match_ids", GetLongOption::NoValue,
                  "Invokes a matching algorithm using the node and element global id\n"
                  "\t\tmaps in the two files.",
                  nullptr);
  options_.enroll("match_file_order", GetLongOption::NoValue,
                  "Verifies that node and element ids match and are in same order\n"
                  "\t\tin the two files.",
                  nullptr);
  options_.enroll("match_by_name", GetLongOption::NoValue,
                  "Match element blocks, nodesets, and sidesets by name instead of by id.",
                  nullptr);
  options_.enroll("nsmap", GetLongOption::NoValue,
                  "Creates a map between the nodeset nodes in the two files\n"
                  "\t\tif they include the same nodes, but are in different order.",
                  nullptr);
  options_.enroll("ssmap", GetLongOption::NoValue,
                  "Creates a map between the sideset faces in the two files\n"
                  "\t\tif they include the same sides, but are in different order.",
                  nullptr);
  options_.enroll("no_nsmap", GetLongOption::NoValue,
                  "Compare nodeset nodes based on file order only", nullptr);
  options_.enroll("no_ssmap", GetLongOption::NoValue,
                  "Compare sideset faces based on file order only", nullptr, nullptr, true);

  options_.enroll("short", GetLongOption::NoValue,
                  "Short block type compare.  Forces element block type strings to\n"
                  "\t\tbe compared only up to the shortest string length.  For example,\n"
                  "\t\t\"HEX\" and \"HEX8\" will be considered the same. (default)",
                  nullptr);
  options_.enroll("no_short", GetLongOption::NoValue,
                  "Do not do short block type compare.  Forces element block\n"
                  "\t\ttype strings to fully match. For example, \"HEX\" and \"HEX8\"\n"
                  "\t\twill be considered different.",
                  nullptr);
  options_.enroll("ignore_case", GetLongOption::NoValue,
                  "Ignore case.  Variable names are compared case in-sensitive (default).",
                  nullptr);
  options_.enroll("case_sensitive", GetLongOption::NoValue,
                  "Variable names are compared case sensitive.", nullptr);
  options_.enroll("nosymmetric_name_check", GetLongOption::NoValue,
                  "No symmetric variable name checking.  By default, a warning will\n"
                  "\t\tbe produced if a name that is not to be excluded is contained\n"
                  "\t\tin the second file given on the command line but not the first.\n"
                  "\t\tThis \"symmetric\" check can be turned off with this option.",
                  nullptr);
  options_.enroll("allow_name_mismatch", GetLongOption::NoValue,
                  "Allow a variable name that is in the first database to not be in the\n"
                  "\t\tsecond database",
                  nullptr, nullptr, true);

  options_.enroll("ignore_maps", GetLongOption::NoValue,
                  "Output node and element diff summaries using file local implicit ids\n"
                  "\t\tinstead of global ids.",
                  nullptr);
  options_.enroll("ignore_nans", GetLongOption::NoValue, "Don't check data for NaNs", nullptr);
  options_.enroll("ignore_attributes", GetLongOption::NoValue,
                  "Don't compare element attribute values.", nullptr);
  options_.enroll("ignore_sideset_df", GetLongOption::NoValue,
                  "Don't compare sideset distribution factors.", nullptr, nullptr, true);

  options_.enroll("norms", GetLongOption::NoValue,
                  "Calculate L1 and L2 norms of variable differences and output if > 0.0", nullptr);
  options_.enroll("l2norms", GetLongOption::NoValue,
                  "Calculate L2 norm of variable differences and output if > 0.0", nullptr);
  options_.enroll("l1norms", GetLongOption::NoValue,
                  "Calculate L1 norm of variable differences and output if > 0.0", nullptr);

  options_.enroll("status", GetLongOption::NoValue,
                  "Return exit status of 2 if the files are different. (default).", nullptr);
  options_.enroll("ignore_status", GetLongOption::NoValue,
                  "The exit status is always zero unless an error occurs.", nullptr);

  options_.enroll(
      "max_warnings", GetLongOption::MandatoryValue,
      "Maximum number of warnings to output during element/node matching process.  Default 100.",
      "100");
  options_.enroll("use_old_floor", GetLongOption::NoValue,
                  "use the older definition of the floor tolerance.\n"
                  "\t\tOLD: ignore if |a-b| < floor.\n"
                  "\t\tNEW: ignore if |a| < floor && |b| < floor.",
                  nullptr);
  options_.enroll("64-bit", GetLongOption::NoValue,
                  "True if forcing the use of 64-bit integers for the output file in summary mode",
                  nullptr);
  options_.enroll("min_coordinate_separation", GetLongOption::NoValue,
                  "In summary mode, calculate the minimum distance between any two nodes", nullptr);
  options_.enroll("copyright", GetLongOption::NoValue, "Output copyright and license information.",
                  nullptr);
  options_.enroll("version", GetLongOption::NoValue, "Output code version", nullptr);

  options_.enroll("maxnames", GetLongOption::MandatoryValue, "[deprecated -- no longer needed]",
                  "1000");
  options_.enroll("t", GetLongOption::MandatoryValue, "Backward-compatible option for -tolerance",
                  "1.0E-6");
  options_.enroll("m", GetLongOption::NoValue, "Backward-compatible option for -map", nullptr);
  options_.enroll("p", GetLongOption::NoValue, "Backward-compatible option for -partial.", nullptr);
  options_.enroll("s", GetLongOption::NoValue, "Backward-compatible option for -short", nullptr);
  options_.enroll("i", GetLongOption::NoValue, "Backward-compatible option for -ignore_case.",
                  nullptr);
  options_.enroll("f", GetLongOption::MandatoryValue, "Backward-compatible option for -file",
                  nullptr);
  options_.enroll("T", GetLongOption::MandatoryValue,
                  "Backward-compatible option for -TimeStepOffset", nullptr);
}

bool SystemInterface::parse_options(int argc, char **argv)
{
  int option_index = options_.parse(argc, argv);
  if (option_index < 1) {
    return false;
  }

  {
    const char *temp = options_.retrieve("help");
    if (temp != nullptr) {
      if ((str_equal("usage", temp)) || (str_equal("all", temp))) {
        options_.usage();
      }
      if ((str_equal("file", temp)) || (str_equal("all", temp))) {
        file_help();
      }
      if ((str_equal("tolerance", temp)) || (str_equal("all", temp))) {
        tolerance_help();
      }
      fmt::print("\n\t\tCan also set options via EXODIFF_OPTIONS environment variable.\n");
      fmt::print("\n\t\tDocumentation: "
                 "https://sandialabs.github.io/seacas-docs/sphinx/html/index.html#exodiff\n");
      fmt::print("\t\t->->-> Send email to gdsjaar@sandia.gov for exodiff support.<-<-<-\n");
      exit(EXIT_SUCCESS);
    }
  }

  if (options_.retrieve("Help") != nullptr) {
    options_.usage();
    fmt::print("\n\t\tCan also set options via EXODIFF_OPTIONS environment variable.\n");
    fmt::print("\n\t\tDocumentation: "
               "https://sandialabs.github.io/seacas-docs/sphinx/html/index.html#exodiff\n");
    fmt::print("\t\t->->-> Send email to gdsjaar@sandia.gov for exodiff support.<-<-<-\n");
    exit(EXIT_SUCCESS);
  }

  if (options_.retrieve("version") != nullptr) {
    show_version();
    exit(EXIT_SUCCESS);
  }

  if (options_.retrieve("copyright") != nullptr) {
    fmt::print("{}", copyright("2008-2021"));
    exit(EXIT_SUCCESS);
  }

  // Parse remaining options as filenames
  if (option_index < argc) {
    file1 = argv[option_index++];
    if (option_index < argc) {
      file2 = argv[option_index++];
    }
    if (option_index < argc) {
      if (option_index + 1 == argc) {
        diff_file = argv[option_index++];
      }
      else {
        // Check for additional unknown arguments...
        std::ostringstream out;
        fmt::print(
            out,
            "\nexodiff: ERROR: Too many file arguments specified.\n"
            "         Probably caused by options following filenames which is no longer allowed.\n"
            "         Unknown options are: ");
        while (option_index < argc) {
          fmt::print(out, "'{}' ", argv[option_index++]);
        }
        fmt::print(out, "\n\n");
        ERR_OUT(out);
        return false;
      }
    }
  }
  else {
    Error("no files specified\n\n");
  }

  // Get options from environment variable also...
  char *options = getenv("EXODIFF_OPTIONS");
  if (options != nullptr) {
    fmt::print(
        "\nThe following options were specified via the EXODIFF_OPTIONS environment variable:\n"
        "\t\t{}\n\n",
        options);
    options_.parse(options, GetLongOption::basename(*argv));
  }

  if (options_.retrieve("summary") != nullptr) {
    summary_flag = true;
  }

  if (options_.retrieve("min_coordinate_separation") != nullptr) {
    coord_sep = true;
  }

  {
    const char *temp = options_.retrieve("exclude");
    if (temp != nullptr) {
      parseExcludeTimes(temp, exclude_steps);
    }
  }

  {
    const char *temp = options_.retrieve("x");
    if (temp != nullptr) {
      parseExcludeTimes(temp, exclude_steps);
    }
  }

  {
    auto t1 = options_.get_option_value("t", default_tol.value);
    auto t2 = options_.get_option_value("tolerance", default_tol.value);
    if (t1 != default_tol.value) {
      default_tol.value = t1;
    }
    else if (t2 != default_tol.value) {
      default_tol.value = t2;
    }
  }

  coord_tol.value      = options_.get_option_value("coordinate_tolerance", coord_tol.value);
  default_tol.floor    = options_.get_option_value("Floor", default_tol.floor);
  final_time_tol.value = options_.get_option_value("final_time_tolerance", final_time_tol.value);

  time_value_offset = options_.get_option_value("time_offset", time_value_offset);
  time_value_scale  = options_.get_option_value("time_scale", time_value_scale);

  {
    const char *temp = options_.retrieve("TimeStepOffset");
    if (temp != nullptr) {
      errno            = 0;
      time_step_offset = strtol(temp, nullptr, 10);
      SMART_ASSERT(errno == 0);
    }
    else {
      const char *temp2 = options_.retrieve("T");
      if (temp2 != nullptr) {
        errno            = 0;
        time_step_offset = strtol(temp2, nullptr, 10);
        SMART_ASSERT(errno == 0);
      }
    }
  }

  if (options_.retrieve("TA") != nullptr) {
    time_step_offset = -1; // Signifies automatic offset calculation.
  }

  if (options_.retrieve("TM") != nullptr) {
    time_step_offset = -2; // Signifies automatic offset calculation -- closest match
  }

  {
    const char *temp = options_.retrieve("steps");
    if (temp != nullptr) {
      Parse_Steps_Option(temp, time_step_start, time_step_stop, time_step_increment);
    }
  }

  {
    const char *temp = options_.retrieve("explicit");
    if (temp != nullptr) {
      // temp should be of the form <ts1>:<ts2>  where ts# is either a timestep number
      // (1-based) or 'last'
      std::vector<std::string> tokens = SLIB::tokenize(temp, ":");
      if (tokens.size() == 2) {
        if (str_equal(tokens[0], "last")) {
          explicit_steps.first = -1;
        }
        else {
          // Try to convert to integer...
          explicit_steps.first = std::stoi(tokens[0]);
        }

        if (str_equal(tokens[1], "last")) {
          explicit_steps.second = -1;
        }
        else {
          // Try to convert to integer...
          explicit_steps.second = std::stoi(tokens[1]);
        }
      }
      else {
        Error(fmt::format("parse error for -explicit keyword. "
                          "Expected '<int|last>:<int|last>', found '{}' Aborting...\n",
                          temp));
      }
    }
  }

  if (options_.retrieve("quiet") != nullptr) {
    quiet_flag = true;
  }

  if (options_.retrieve("show_all_diffs") != nullptr) {
    show_all_diffs = true;
  }

  if ((options_.retrieve("partial") != nullptr) || (options_.retrieve("p") != nullptr)) {
    map_flag = MapType::PARTIAL;
  }

  if (options_.retrieve("match_ids") != nullptr) {
    map_flag = MapType::USE_FILE_IDS;
  }

  if (options_.retrieve("match_file_order") != nullptr) {
    map_flag = MapType::FILE_ORDER;
  }

  if (options_.retrieve("match_by_name") != nullptr) {
    by_name = true;
  }

  if ((options_.retrieve("map") != nullptr) || (options_.retrieve("m") != nullptr)) {
    map_flag = MapType::DISTANCE;
  }
  if (options_.retrieve("nsmap") != nullptr) {
    nsmap_flag = true;
  }
  if (options_.retrieve("no_nsmap") != nullptr) {
    nsmap_flag = false;
  }
  if (options_.retrieve("ssmap") != nullptr) {
    ssmap_flag = true;
  }
  if (options_.retrieve("no_ssmap") != nullptr) {
    ssmap_flag = false;
  }
  if ((options_.retrieve("short") != nullptr) || (options_.retrieve("s") != nullptr)) {
    short_block_check = true;
  }
  if (options_.retrieve("no_short") != nullptr) {
    short_block_check = false;
  }
  if (options_.retrieve("nosymmetric_name_check") != nullptr) {
    noSymmetricNameCheck = true;
  }
  if (options_.retrieve("norms") != nullptr) {
    doL1Norm = true;
    doL2Norm = true;
  }
  if (options_.retrieve("l2norms") != nullptr) {
    doL2Norm = true;
  }
  if (options_.retrieve("l1norms") != nullptr) {
    doL1Norm = true;
  }
  if (options_.retrieve("pedantic") != nullptr) {
    pedantic = true;
  }
  if (options_.retrieve("interpolate") != nullptr) {
    interpolating = true;
  }

  if (options_.retrieve("allow_name_mismatch") != nullptr) {
    allowNameMismatch = true;
  }
  if ((options_.retrieve("ignore_case") != nullptr) || (options_.retrieve("i") != nullptr)) {
    nocase_var_names = true;
  }
  if (options_.retrieve("case_sensitive") != nullptr) {
    nocase_var_names = false;
  }
  if (options_.retrieve("ignore_maps") != nullptr) {
    ignore_maps = true;
  }
  if (options_.retrieve("ignore_nans") != nullptr) {
    ignore_nans = true;
  }
  if (options_.retrieve("ignore_dups") != nullptr) {
    ignore_dups = true;
  }
  if (options_.retrieve("ignore_steps") != nullptr) {
    ignore_steps = true;
  }
  if (options_.retrieve("64-bit") != nullptr) {
    ints_64_bits = true;
  }
  if (options_.retrieve("ignore_attributes") != nullptr) {
    ignore_attributes = true;
  }
  if (options_.retrieve("ignore_sideset_df") != nullptr) {
    ignore_sideset_df = true;
  }
  if (options_.retrieve("relative") != nullptr) {
    output_type      = ToleranceMode::RELATIVE_; // Change type to relative.
    default_tol.type = ToleranceMode::RELATIVE_;
  }
  if (options_.retrieve("ignore") != nullptr) {
    output_type      = ToleranceMode::IGNORE_; // Change type to ignored
    default_tol.type = ToleranceMode::IGNORE_;
  }
  if (options_.retrieve("absolute") != nullptr) {
    output_type      = ToleranceMode::ABSOLUTE_; // Change type to absolute
    default_tol.type = ToleranceMode::ABSOLUTE_;
  }
  if (options_.retrieve("combined") != nullptr) {
    output_type      = ToleranceMode::COMBINED_; // Change type to combine
    default_tol.type = ToleranceMode::COMBINED_;
  }
  if (options_.retrieve("ulps_float") != nullptr) {
    output_type      = ToleranceMode::ULPS_FLOAT_;
    default_tol.type = ToleranceMode::ULPS_FLOAT_;
  }
  if (options_.retrieve("ulps_double") != nullptr) {
    output_type      = ToleranceMode::ULPS_DOUBLE_;
    default_tol.type = ToleranceMode::ULPS_DOUBLE_;
  }
  if (options_.retrieve("eigen_relative") != nullptr) {
    output_type      = ToleranceMode::EIGEN_REL_; // Change type to relative.
    default_tol.type = ToleranceMode::EIGEN_REL_;
  }
  if (options_.retrieve("eigen_absolute") != nullptr) {
    output_type      = ToleranceMode::EIGEN_ABS_; // Change type to absolute
    default_tol.type = ToleranceMode::EIGEN_ABS_;
  }
  if (options_.retrieve("eigen_combined") != nullptr) {
    output_type      = ToleranceMode::EIGEN_COM_; // Change type to combine
    default_tol.type = ToleranceMode::EIGEN_COM_;
  }
  if (options_.retrieve("dumpmap") != nullptr) {
    dump_mapping = true;
  }
  if (options_.retrieve("show_unmatched") != nullptr) {
    show_unmatched = true;
  }

  {
    const char *temp = options_.retrieve("max_warnings");
    if (temp != nullptr) {
      errno        = 0;
      max_warnings = strtol(temp, nullptr, 10);
      SMART_ASSERT(errno == 0);
    }
  }

  if (options_.retrieve("status") != nullptr) {
    exit_status_switch = true;
  }

  if (options_.retrieve("ignore_status") != nullptr) {
    exit_status_switch = false;
  }

  if (options_.retrieve("use_old_floor") != nullptr) {
    Tolerance::use_old_floor = true; // Change type to relative.
  }

  {
    // Reset default tolerances in case the -t flag was given.
    time_tol         = default_tol;
    glob_var_default = default_tol;
    node_var_default = default_tol;
    elmt_var_default = default_tol;
    elmt_att_default = default_tol;
    ns_var_default   = default_tol;
    ss_var_default   = default_tol;
    eb_var_default   = default_tol;
    fb_var_default   = default_tol;

    const char *temp = options_.retrieve("file");
    if (temp != nullptr) {
      command_file = temp;
      if (!summary_flag && (File_Exists(command_file) == 0)) {
        Error(fmt::format("Can't open file \"{}\".\n", command_file));
      }

      // Command file exists, parse contents...
      Parse_Command_File();
    }
    else {
      const char *t2 = options_.retrieve("f");
      if (t2 != nullptr) {
        command_file = t2;
        if (!summary_flag && (File_Exists(command_file) == 0)) {
          Error(fmt::format("Can't open file \"{}\".\n", command_file));
        }

        // Command file exists, parse contents...
        Parse_Command_File();
      }
      else {
        glob_var_do_all_flag = true;
        node_var_do_all_flag = true;
        elmt_var_do_all_flag = true;
        elmt_att_do_all_flag = true;
        ns_var_do_all_flag   = true;
        ss_var_do_all_flag   = true;
        eb_var_do_all_flag   = true;
        fb_var_do_all_flag   = true;
      }
    }
  }
  return true;
}

void SystemInterface::Parse_Command_File()
{
  int default_tol_specified = 0;

  std::ifstream cmd_file(command_file, std::ios::in);
  SMART_ASSERT(cmd_file.good());

  std::string line;
  std::string xline, tok2, tok3;
  std::getline(cmd_file, line);
  xline = line;
  while (!cmd_file.eof()) {
    std::string tok1;
    // Skip blank lines and comment lines.
    if (count_tokens(xline, " \t") > 0 && (tok1 = extract_token(xline, " \t"))[0] != '#') {
      to_lower(tok1); // Make case insensitive.
      tok2 = extract_token(xline, " \t");
      to_lower(tok2);

      if (abbreviation(tok1, "default", 3) && abbreviation(tok2, "tolerance", 3)) {
        std::string tok = extract_token(xline, " \n\t=,");
        to_lower(tok);
        if (tok.empty()) {
          Parse_Die(line);
        }

        if (abbreviation(tok, "relative", 3)) {
          default_tol.type = ToleranceMode::RELATIVE_;
          tok              = extract_token(xline, " \n\t=,");
        }
        else if (abbreviation(tok, "absolute", 3)) {
          default_tol.type = ToleranceMode::ABSOLUTE_;
          tok              = extract_token(xline, " \n\t=,");
        }
        else if (abbreviation(tok, "combine", 3)) {
          default_tol.type = ToleranceMode::COMBINED_;
          tok              = extract_token(xline, " \n\t=,");
        }
        else if (abbreviation(tok, "eigen_relative", 7)) {
          default_tol.type = ToleranceMode::EIGEN_REL_;
          tok              = extract_token(xline, " \n\t=,");
        }
        else if (abbreviation(tok, "eigen_absolute", 7)) {
          default_tol.type = ToleranceMode::EIGEN_ABS_;
          tok              = extract_token(xline, " \n\t=,");
        }
        else if (abbreviation(tok, "eigen_combine", 7)) {
          default_tol.type = ToleranceMode::EIGEN_COM_;
          tok              = extract_token(xline, " \n\t=,");
        }
        else if (abbreviation(tok, "ignore", 3)) {
          default_tol.type = ToleranceMode::IGNORE_;
          tok              = extract_token(xline, " \n\t=,");
        }
        if (tok.empty()) {
          Parse_Die(line);
        }

        default_tol.value = To_Double(tok);

        tok = extract_token(xline, " \n\t=,");
        to_lower(tok);
        if (abbreviation(tok, "floor", 3)) {
          tok = extract_token(xline, " \n\t=,");
          if (tok.empty()) {
            Parse_Die(line);
          }
          default_tol.floor = To_Double(tok);
        }
        default_tol_specified = 1;
      }
      else if (abbreviation(tok1, "max", 3) && abbreviation(tok2, "names", 3)) {
        ; // Ignored -- no longer needed.
      }
      else if (abbreviation(tok1, "final", 3) && abbreviation(tok2, "time", 3)) {
        tok3 = extract_token(xline, " \t");
        to_lower(tok3);
        if (!abbreviation(tok3, "tolerance", 3)) {
          Error(fmt::format(" expected \"TOLERANCE\" after the \"FINAL TIME\" keyword. "
                            "Found \"{}\" instead. Aborting...\n",
                            tok3));
        }
        std::string tok = extract_token(xline, " \n\t=,");
        if (tok.empty()) {
          Parse_Die(line);
        }
        final_time_tol.value = To_Double(tok);
      }
      else if (abbreviation(tok1, "return", 3) && abbreviation(tok2, "status", 3)) {
        exit_status_switch = true;
      }
      else if (abbreviation(tok1, "ignore", 3) && abbreviation(tok2, "status", 3)) {
        exit_status_switch = false;
      }
      else if (abbreviation(tok1, "exclude", 3) && abbreviation(tok2, "times", 3)) {
        std::string tok = extract_token(xline, " \n\t=");
        if (!tok.empty() && tok[0] != '#') {
          parseExcludeTimes(tok, exclude_steps);
        }
      }
      else if (abbreviation(tok1, "apply", 3) && abbreviation(tok2, "matching", 3)) {
        map_flag = MapType::DISTANCE;
      }
      else if (abbreviation(tok1, "calculate", 3) && abbreviation(tok2, "norms", 3)) {
        doL2Norm = true;
        doL1Norm = true;
      }
      else if (abbreviation(tok1, "calculate", 3) && abbreviation(tok2, "l2norms", 3)) {
        doL2Norm = true;
      }
      else if (abbreviation(tok1, "calculate", 3) && abbreviation(tok2, "l1norms", 3)) {
        doL1Norm = true;
      }
      else if (tok1 == "nodeset" && abbreviation(tok2, "match", 3)) {
        nsmap_flag = true;
      }
      else if (tok1 == "pedantic") {
        pedantic = true;
      }
      else if (tok1 == "interpolate") {
        interpolating = true;
      }
      else if (tok1 == "sideset" && abbreviation(tok2, "match", 3)) {
        ssmap_flag = true;
      }
      else if (abbreviation(tok1, "short", 3) && abbreviation(tok2, "blocks", 3)) {
        short_block_check = true;
      }
      else if (tok1 == "no" && abbreviation(tok2, "short", 3)) {
        short_block_check = false;
      }
      else if (abbreviation(tok1, "ignore", 3) && abbreviation(tok2, "case", 3)) {
        nocase_var_names = true;
      }
      else if (abbreviation(tok1, "case", 3) && abbreviation(tok2, "sensitive", 3)) {
        nocase_var_names = false;
      }
      else if (abbreviation(tok1, "ignore", 3) && abbreviation(tok2, "maps", 3)) {
        ignore_maps = true;
      }
      else if (abbreviation(tok1, "ignore", 3) && abbreviation(tok2, "nans", 3)) {
        ignore_nans = true;
      }
      else if (abbreviation(tok1, "ignore", 3) && abbreviation(tok2, "dups", 3)) {
        ignore_dups = true;
      }
      else if (abbreviation(tok1, "ignore", 3) && abbreviation(tok2, "attributes", 3)) {
        ignore_attributes = true;
      }
      else if (abbreviation(tok1, "ignore", 3) && abbreviation(tok2, "sideset", 3)) {
        tok3 = extract_token(xline, " \t");
        to_lower(tok3);
        if (abbreviation(tok3, "distribution", 3)) {
          ignore_sideset_df = true;
        }
      }
      else if (tok1 == "step" && tok2 == "offset") {
        std::string tok = extract_token(xline, " \n\t=");
        if (abbreviation(tok, "automatic", 4)) {
          time_step_offset = -1;
        }
        else if (abbreviation(tok, "match", 4)) {
          time_step_offset = -2;
        }
        else {
          errno            = 0;
          time_step_offset = std::stoi(tok);
          SMART_ASSERT(errno == 0);
        }
      }
      else if (abbreviation(tok1, "coordinates", 4)) {
        if (default_tol_specified != 0) {
          coord_tol = default_tol;
        }
        else {
          coord_tol.type  = ToleranceMode::ABSOLUTE_; // These should correspond to
          coord_tol.value = 1.e-6;                    // the defaults at the top of
          coord_tol.floor = 0.0;                      // this file.
        }

        if (!tok2.empty() && tok2[0] != '#') {
          // If rel or abs is specified, then the tolerance must
          // be specified.
          if (abbreviation(tok2, "relative", 3)) {
            coord_tol.type = ToleranceMode::RELATIVE_;
            tok2           = extract_token(xline, " \n\t=");
            if (tok2.empty()) {
              Parse_Die(line);
            }
            coord_tol.value = To_Double(tok2);
          }
          else if (abbreviation(tok2, "absolute", 3)) {
            coord_tol.type = ToleranceMode::ABSOLUTE_;
            tok2           = extract_token(xline, " \n\t=");
            if (tok2.empty()) {
              Parse_Die(line);
            }
            coord_tol.value = To_Double(tok2);
          }
          else if (abbreviation(tok2, "combine", 3)) {
            coord_tol.type = ToleranceMode::COMBINED_;
            tok2           = extract_token(xline, " \n\t=");
            if (tok2.empty()) {
              Parse_Die(line);
            }
            coord_tol.value = To_Double(tok2);
          }
          else if (abbreviation(tok2, "eigen_relative", 7)) {
            coord_tol.type = ToleranceMode::EIGEN_REL_;
            tok2           = extract_token(xline, " \n\t=");
            if (tok2.empty()) {
              Parse_Die(line);
            }
            coord_tol.value = To_Double(tok2);
          }
          else if (abbreviation(tok2, "eigen_absolute", 7)) {
            coord_tol.type = ToleranceMode::EIGEN_ABS_;
            tok2           = extract_token(xline, " \n\t=");
            if (tok2.empty()) {
              Parse_Die(line);
            }
            coord_tol.value = To_Double(tok2);
          }
          else if (abbreviation(tok2, "eigen_combine", 7)) {
            coord_tol.type = ToleranceMode::EIGEN_COM_;
            tok2           = extract_token(xline, " \n\t=");
            if (tok2.empty()) {
              Parse_Die(line);
            }
            coord_tol.value = To_Double(tok2);
          }
          else if (abbreviation(tok2, "ignore", 3)) {
            coord_tol.type  = ToleranceMode::IGNORE_;
            coord_tol.value = 0.0;
          }
          else if (abbreviation(tok2, "floor", 3)) {
            tok2 = extract_token(xline, " \n\t=");
            if (tok2.empty()) {
              Parse_Die(line);
            }
            coord_tol.floor = To_Double(tok2);
          }

          tok2 = extract_token(xline, " \n\t=,");
          to_lower(tok2);
          if (abbreviation(tok2, "floor", 3)) {
            tok2 = extract_token(xline, " \n\t=,");
            if (tok2.empty()) {
              Parse_Die(line);
            }
            coord_tol.floor = To_Double(tok2);
          }
        }
      }
      else if (tok1 == "time" && abbreviation(tok2, "steps", 4)) {
        time_tol = default_tol;

        std::string tok = extract_token(xline, " \n\t=");
        to_lower(tok);
        if (!tok.empty() && tok[0] != '#') {
          // If rel or abs is specified, then the tolerance
          // must be specified.
          if (abbreviation(tok, "relative", 3)) {
            time_tol.type = ToleranceMode::RELATIVE_;
            tok           = extract_token(xline, " \n\t=");
            if (tok.empty()) {
              Parse_Die(line);
            }
            time_tol.value = To_Double(tok);
          }
          else if (abbreviation(tok, "absolute", 3)) {
            time_tol.type = ToleranceMode::ABSOLUTE_;
            tok           = extract_token(xline, " \n\t=");
            if (tok.empty()) {
              Parse_Die(line);
            }
            time_tol.value = To_Double(tok);
          }
          else if (abbreviation(tok, "combine", 3)) {
            time_tol.type = ToleranceMode::COMBINED_;
            tok           = extract_token(xline, " \n\t=");
            if (tok.empty()) {
              Parse_Die(line);
            }
            time_tol.value = To_Double(tok);
          }
          else if (abbreviation(tok, "ignore", 3)) {
            time_tol.type  = ToleranceMode::IGNORE_;
            time_tol.value = 0.0;
          }
          else if (abbreviation(tok, "floor", 3)) {
            tok = extract_token(xline, " \n\t=");
            if (tok.empty()) {
              Parse_Die(line);
            }
            time_tol.floor = To_Double(tok);
          }

          tok2 = extract_token(xline, " \n\t=,");
          to_lower(tok2);
          if (abbreviation(tok2, "floor", 3)) {
            tok2 = extract_token(xline, " \n\t=,");
            if (tok2.empty()) {
              Parse_Die(line);
            }
            time_tol.floor = To_Double(tok2);
          }
        }
      }
      else if (abbreviation(tok1, "global", 4) && abbreviation(tok2, "variables", 3)) {
        glob_var_default = default_tol;
        xline            = Parse_Variables(xline, cmd_file, glob_var_do_all_flag, glob_var_default,
                                           glob_var_names, glob_var);

        Check_Parsed_Names(glob_var_names, glob_var_do_all_flag);

        line = xline;
        continue;
      }
      else if (abbreviation(tok1, "nodal", 4) && abbreviation(tok2, "variables", 3)) {
        node_var_default = default_tol;
        xline            = Parse_Variables(xline, cmd_file, node_var_do_all_flag, node_var_default,
                                           node_var_names, node_var);

        Check_Parsed_Names(node_var_names, node_var_do_all_flag);

        line = xline;
        continue;
      }
      else if (abbreviation(tok1, "element", 4) && abbreviation(tok2, "variables", 3)) {
        elmt_var_default = default_tol;
        xline            = Parse_Variables(xline, cmd_file, elmt_var_do_all_flag, elmt_var_default,
                                           elmt_var_names, elmt_var);

        Check_Parsed_Names(elmt_var_names, elmt_var_do_all_flag);

        line = xline;
        continue;
      }
      else if (tok1 == "nodeset" && abbreviation(tok2, "variables", 3)) {
        ns_var_default = default_tol;
        xline = Parse_Variables(xline, cmd_file, ns_var_do_all_flag, ns_var_default, ns_var_names,
                                ns_var);

        Check_Parsed_Names(ns_var_names, ns_var_do_all_flag);

        line = xline;
        continue;
      }
      else if (abbreviation(tok1, "sideset", 4) && abbreviation(tok2, "variables", 3)) {
        ss_var_default = default_tol;
        xline = Parse_Variables(xline, cmd_file, ss_var_do_all_flag, ss_var_default, ss_var_names,
                                ss_var);

        Check_Parsed_Names(ss_var_names, ss_var_do_all_flag);

        line = xline;
        continue;
      }
      else if (abbreviation(tok1, "sideset", 4) && abbreviation(tok2, "distribution", 4)) {
        if (default_tol_specified != 0) {
          ss_df_tol = default_tol;
        }
        else {
          ss_df_tol.type  = ToleranceMode::ABSOLUTE_; // These should correspond to
          ss_df_tol.value = 1.e-6;                    // the defaults at the top of
          ss_df_tol.floor = 0.0;                      // this file.
        }

        if (!tok2.empty() && tok2[0] != '#') {
          // If rel or abs is specified, then the tolerance must
          // be specified.
          if (abbreviation(tok2, "relative", 3)) {
            ss_df_tol.type = ToleranceMode::RELATIVE_;
            tok2           = extract_token(xline, " \n\t=");
            if (tok2.empty()) {
              Parse_Die(line);
            }
            ss_df_tol.value = To_Double(tok2);
          }
          else if (abbreviation(tok2, "absolute", 3)) {
            ss_df_tol.type = ToleranceMode::ABSOLUTE_;
            tok2           = extract_token(xline, " \n\t=");
            if (tok2.empty()) {
              Parse_Die(line);
            }
            ss_df_tol.value = To_Double(tok2);
          }
          else if (abbreviation(tok2, "combine", 3)) {
            ss_df_tol.type = ToleranceMode::COMBINED_;
            tok2           = extract_token(xline, " \n\t=");
            if (tok2.empty()) {
              Parse_Die(line);
            }
            ss_df_tol.value = To_Double(tok2);
          }
          else if (abbreviation(tok2, "eigen_relative", 7)) {
            ss_df_tol.type = ToleranceMode::EIGEN_REL_;
            tok2           = extract_token(xline, " \n\t=");
            if (tok2.empty()) {
              Parse_Die(line);
            }
            ss_df_tol.value = To_Double(tok2);
          }
          else if (abbreviation(tok2, "eigen_absolute", 7)) {
            ss_df_tol.type = ToleranceMode::EIGEN_ABS_;
            tok2           = extract_token(xline, " \n\t=");
            if (tok2.empty()) {
              Parse_Die(line);
            }
            ss_df_tol.value = To_Double(tok2);
          }
          else if (abbreviation(tok2, "eigen_combine", 7)) {
            ss_df_tol.type = ToleranceMode::EIGEN_COM_;
            tok2           = extract_token(xline, " \n\t=");
            if (tok2.empty()) {
              Parse_Die(line);
            }
            ss_df_tol.value = To_Double(tok2);
          }
          else if (abbreviation(tok2, "ignore", 3)) {
            ss_df_tol.type  = ToleranceMode::IGNORE_;
            ss_df_tol.value = 0.0;
          }
          else if (abbreviation(tok2, "floor", 3)) {
            tok2 = extract_token(xline, " \n\t=");
            if (tok2.empty()) {
              Parse_Die(line);
            }
            ss_df_tol.floor = To_Double(tok2);
          }

          tok2 = extract_token(xline, " \n\t=,");
          to_lower(tok2);
          if (abbreviation(tok2, "floor", 3)) {
            tok2 = extract_token(xline, " \n\t=,");
            if (tok2.empty()) {
              Parse_Die(line);
            }
            ss_df_tol.floor = To_Double(tok2);
          }
        }
      }
      else if (abbreviation(tok1, "edgeblock", 4) && abbreviation(tok2, "variables", 3)) {
        eb_var_default = default_tol;
        xline = Parse_Variables(xline, cmd_file, eb_var_do_all_flag, eb_var_default, eb_var_names,
                                eb_var);

        Check_Parsed_Names(eb_var_names, eb_var_do_all_flag);

        line = xline;
        continue;
      }
      else if (abbreviation(tok1, "faceblock", 4) && abbreviation(tok2, "variables", 3)) {
        fb_var_default = default_tol;
        xline = Parse_Variables(xline, cmd_file, fb_var_do_all_flag, fb_var_default, fb_var_names,
                                fb_var);

        Check_Parsed_Names(fb_var_names, fb_var_do_all_flag);

        line = xline;
        continue;
      }
      else if (abbreviation(tok1, "element", 4) && abbreviation(tok2, "attributes", 3)) {
        elmt_att_default = default_tol;
        xline            = Parse_Variables(xline, cmd_file, elmt_att_do_all_flag, elmt_att_default,
                                           elmt_att_names, elmt_att);

        Check_Parsed_Names(elmt_att_names, elmt_att_do_all_flag);

        line = xline;
        continue;
      }
      else {
        Parse_Die(line);
      }
    }

    std::getline(cmd_file, line);
    xline = line;
  }
}

namespace {
  std::string Parse_Variables(std::string xline, std::ifstream &cmd_file, bool &all_flag,
                              Tolerance &def_tol, std::vector<std::string> &names,
                              std::vector<Tolerance> &toler)
  {
    toler.clear();
    names.clear();

    std::string tok = extract_token(xline, " \n\t=,");
    to_lower(tok);
    if (!tok.empty()) {
      if (tok != "(all)" && tok != "all" && !abbreviation(tok, "relative", 3) &&
          !abbreviation(tok, "absolute", 3) && !abbreviation(tok, "combine", 3) &&
          !abbreviation(tok, "ulps_float", 6) && !abbreviation(tok, "ulps_double", 6) &&
          !abbreviation(tok, "eigen_relative", 7) && !abbreviation(tok, "eigen_absolute", 7) &&
          !abbreviation(tok, "eigen_combine", 7) && !abbreviation(tok, "ignore", 3) &&
          !abbreviation(tok, "floor", 3)) {
        Error(fmt::format("in parsing command file: unrecognized keyword \"{}\"\n", tok));
      }

      if (tok == "(all)" || tok == "all") {
        all_flag = true;
        tok      = extract_token(xline, " \n\t=,");
      }

      // If rel or abs is specified, then the tolerance must be specified.
      if (abbreviation(tok, "relative", 3)) {
        def_tol.type = ToleranceMode::RELATIVE_;
        tok          = extract_token(xline, " \n\t=,");
        if (tok == "floor" || tok.empty()) {
          Error(" Input file specifies a tolerance type "
                "but no tolerance\n");
        }
        def_tol.value = To_Double(tok);
        tok           = extract_token(xline, " \n\t=,");
        to_lower(tok);
      }
      else if (abbreviation(tok, "absolute", 3)) {
        def_tol.type = ToleranceMode::ABSOLUTE_;
        tok          = extract_token(xline, " \n\t=,");
        if (tok == "floor" || tok.empty()) {
          Error("Input file specifies a tolerance type "
                "but no tolerance\n");
        }
        def_tol.value = To_Double(tok);
        tok           = extract_token(xline, " \n\t=,");
        to_lower(tok);
      }
      else if (abbreviation(tok, "combine", 3)) {
        def_tol.type = ToleranceMode::COMBINED_;
        tok          = extract_token(xline, " \n\t=,");
        if (tok == "floor" || tok.empty()) {
          Error("Input file specifies a tolerance type "
                "but no tolerance\n");
        }
        def_tol.value = To_Double(tok);
        tok           = extract_token(xline, " \n\t=,");
        to_lower(tok);
      }
      else if (abbreviation(tok, "ulps_float", 6)) {
        def_tol.type = ToleranceMode::ULPS_FLOAT_;
        tok          = extract_token(xline, " \n\t=,");
        if (tok == "floor" || tok.empty()) {
          Error("Input file specifies a tolerance type "
                "but no tolerance\n");
        }
        def_tol.value = To_Double(tok);
        tok           = extract_token(xline, " \n\t=,");
        to_lower(tok);
      }
      else if (abbreviation(tok, "ulps_double", 6)) {
        def_tol.type = ToleranceMode::ULPS_DOUBLE_;
        tok          = extract_token(xline, " \n\t=,");
        if (tok == "floor" || tok.empty()) {
          Error("Input file specifies a tolerance type "
                "but no tolerance\n");
        }
        def_tol.value = To_Double(tok);
        tok           = extract_token(xline, " \n\t=,");
        to_lower(tok);
      }
      else if (abbreviation(tok, "eigen_relative", 7)) {
        def_tol.type = ToleranceMode::EIGEN_REL_;
        tok          = extract_token(xline, " \n\t=,");
        if (tok == "floor" || tok.empty()) {
          Error("Input file specifies a tolerance type "
                "but no tolerance\n");
        }
        def_tol.value = To_Double(tok);
        tok           = extract_token(xline, " \n\t=,");
        to_lower(tok);
      }
      else if (abbreviation(tok, "eigen_absolute", 7)) {
        def_tol.type = ToleranceMode::EIGEN_ABS_;
        tok          = extract_token(xline, " \n\t=,");
        if (tok == "floor" || tok.empty()) {
          Error("Input file specifies a tolerance type "
                "but no tolerance\n");
        }
        def_tol.value = To_Double(tok);
        tok           = extract_token(xline, " \n\t=,");
        to_lower(tok);
      }
      else if (abbreviation(tok, "eigen_combine", 7)) {
        def_tol.type = ToleranceMode::EIGEN_COM_;
        tok          = extract_token(xline, " \n\t=,");
        if (tok == "floor" || tok.empty()) {
          Error("Input file specifies a tolerance type "
                "but no tolerance\n");
        }
        def_tol.value = To_Double(tok);
        tok           = extract_token(xline, " \n\t=,");
        to_lower(tok);
      }
      else if (abbreviation(tok, "ignore", 3)) {
        def_tol.type  = ToleranceMode::IGNORE_;
        def_tol.value = 0.0;
        tok           = extract_token(xline, " \n\t=,");
        to_lower(tok);
      }

      if (abbreviation(tok, "floor", 3)) {
        tok = extract_token(xline, " \n\t=,");
        if (tok.empty() || tok[0] == '#') {
          Error("Floor specified but couldn't find value\n");
        }
        def_tol.floor = To_Double(tok);
      }
    }

    std::string line{};
    std::getline(cmd_file, line);
    xline = line;
    while (!cmd_file.eof()) {
      if (xline.empty() ||
          ((xline[0] != '\t' && xline[0] != ' ') && first_character(xline) != '#')) {
        break;
      }

      if (first_character(xline) != '#') {
        tok = extract_token(xline);
        chop_whitespace(tok);
        if (tok.empty()) {
          continue; // Found tab but no name given.
        }

        if (tok[0] == '!') {
          // A "!" in front of a name means to exclude the name so no
          // need to look for difference type and tolerance.
          std::string tmp = tok;
          if (!extract_token(tmp, "!").empty()) {
            names.push_back(tok);
            toler.push_back(def_tol);
          }
          std::getline(cmd_file, line);
          xline = line;
          continue;
        }

        int idx = names.size();
        names.push_back(tok);
        toler.push_back(def_tol);

        tok = extract_token(xline);
        to_lower(tok);

        if (!tok.empty() && tok[0] != '#') {
          if (abbreviation(tok, "relative", 3)) {
            toler[idx].type = ToleranceMode::RELATIVE_;
            tok             = extract_token(xline, " \n\t=,");
          }
          else if (abbreviation(tok, "absolute", 3)) {
            toler[idx].type = ToleranceMode::ABSOLUTE_;
            tok             = extract_token(xline, " \n\t=,");
          }
          else if (abbreviation(tok, "combine", 3)) {
            toler[idx].type = ToleranceMode::COMBINED_;
            tok             = extract_token(xline, " \n\t=,");
          }
          else if (abbreviation(tok, "eigen_relative", 7)) {
            toler[idx].type = ToleranceMode::EIGEN_REL_;
            tok             = extract_token(xline, " \n\t=,");
          }
          else if (abbreviation(tok, "eigen_absolute", 7)) {
            toler[idx].type = ToleranceMode::EIGEN_ABS_;
            tok             = extract_token(xline, " \n\t=,");
          }

          else if (abbreviation(tok, "eigen_com", 7)) {
            toler[idx].type = ToleranceMode::EIGEN_COM_;
            tok             = extract_token(xline, " \n\t=,");
          }

          if (abbreviation(tok, "floor", 3)) {
            toler[idx].value = def_tol.value;

            tok = extract_token(xline, " \n\t=,");
            if (tok.empty()) {
              Parse_Die(line);
            }
            toler[idx].floor = To_Double(tok);
          }
          else {
            if (tok.empty()) {
              Parse_Die(line);
            }
            toler[idx].value = To_Double(tok);

            tok = extract_token(xline, " \n\t=,");
            to_lower(tok);
            if (abbreviation(tok, "floor", 3)) {
              tok = extract_token(xline, " \n\t=,");
              if (tok.empty()) {
                Parse_Die(line);
              }
              toler[idx].floor = To_Double(tok);
            }
            else {
              toler[idx].floor = def_tol.floor;
            }
          }
        }
        else {
          toler[idx] = def_tol;
        }
      }

      std::getline(cmd_file, line);
      xline = line;
    }

    if (names.empty()) {
      all_flag = true;
    }

    return xline;
  }

  void tolerance_help()
  {
    fmt::print(
        "\n Tolerance Help:\n"
        "\n"
        "\t Relative difference  |val1 - val2| / max(|val1|, |val2|)\n"
        "\t Absolute difference  |val1 - val2|\n"
        "\t Combined difference  |val1 - val2| / max(tol, tol * max(|val1|, |val2|))\n"
        "\t Eigen_relative difference  ||val1| - |val2|| / max(|val1|,|val2|)\n"
        "\t Eigen_absolute difference  ||val1| - |val2||\n"
        "\t Eigen_combined difference  ||val1| - |val2|| / max(tol, tol * max(|val1|, |val2|))\n"
        "\t Ulps_float difference  -- Calculate number of representable floats between the two "
        "values\n"
        "\t Ulps_double difference  -- Calculate number of representable doubles between the "
        "two values\n"
        "\n"
        "\t Values are considered equal if |val1| <= floor and |val2| <= floor;\n"
        "\t where floor is a user-specified value (-Floor option). Otherwise the difference is\n"
        "\t computed using one of the above formulas and compared to a tolerance.\n"
        "\t If the difference is greater than the tolerance, then the databases\n"
        "\t are different.  At the end of execution, a summary of the differences\n"
        "\t found is output.\n"
        "\t \n"
        "\t By default:\n"
        "\t * All results variables and attributes are compared using a relative difference\n"
        "\t   of 10^{{-6}} (about 6 significant digits) and a floor of 0.0.\n"
        "\t * Nodal locations are compared using {{absolute difference}} with\n"
        "\t   a tolerance of 10^{{-6}} and a floor of 0.0.\n"
        "\t * Time step values are compared using relative difference tolerance of 10^{{-6}}\n"
        "\t   and a floor of 10^{{-15}}.\n"
        "\n\n");
  }

  void file_help()
  {
    fmt::print(
        "\n  Command file syntax:\n"
        "\n"
        "                # Anything following a # is a comment.\n"
        "                DEFAULT TOLERANCE relative 1.E-8 floor 1.E-14\n"
        "                COORDINATES absolute 1.E-12\n"
        "                TIME STEPS absolute 1.E-14\n"
        "                GLOBAL VARIABLES relative 1.E-4 floor 1.E-12\n"
        "                NODAL VARIABLES absolute 1.E-8\n"
        "                <tab> DISPLX\n"
        "                <tab> VELX absolute 1.E-6\n"
        "                <tab> VELY relative 1.E-6 floor 1.e-10\n"
        "                ELEMENT VARIABLES\n"
        "                <tab> !SIGYY\n"
        "                <tab> !SIGZZ\n"
        "\n"
        "         - The variable names are case insensitive (unless or CASE SENSITIVE "
        "specified),\n"
        "           All other comparisons are also case insensitive. Abbreviations can be used. "
        "\n"
        "         - All comparisons use the default of relative 1.e-6 for\n"
        "           variables and absolute 1.e-6 for coordinates.  This is overridden\n"
        "           with the DEFAULT TOLERANCE line.  The DEFAULT TOLERANCE values\n"
        "           are overridden by the values given on the VARIABLES line and apply\n"
        "           only to those variables.  Each variable can override all values\n"
        "           by following its name with a value.\n"
        "         - A variable name must start with a tab character.  Only those\n"
        "           variables listed will be considered (unless \"(all)\") specified.\n"
        "           The NOT symbol \"!\" means do not include this variable.\n"
        "           Mixing non-! and ! is not allowed without the \"(all)\" specifier, e.g.:\n\n"
        "                NODAL VARIABLES (all) absolute 1.E-8\n"
        "                <tab> DISPLX\n"
        "                <tab> !VELX\n"
        "                <tab> VELY relative 1.E-6 floor 1.e-10\n\n"
        "           In this case, all variables are considered that are not prepended\n"
        "           with a \"!\" symbol.\n"
        "         - If a variable type (e.g. NODAL VARIABLES) is not specified, no\n"
        "           variables of that type will be considered.\n"
        "\n"
        "  Other Keywords:\n"
        "         - EXCLUDE TIMES <list>: \"-exclude\" Specify the time step exclusion option.\n"
        "               <list> has the same format as in the command line option.\n"
        "         - APPLY MATCHING: \"-map\" matching algorithm.\n"
        "         - NODESET MATCH: \"-nsmap\", nodeset mapping algorithm.\n"
        "         - SIDESET MATCH: \"-ssmap\", sideset mapping algorithm.\n"
        "         - SHORT BLOCKS: \"-short\", short block type compare.\n"
        "         - NO SHORT BLOCKS: \"-no_short\", do not do short block type compare.\n"
        "         - CASE_SENSITIVE: \"-case_sensitive\", variable names are compared case "
        "sensitive.\n"
        "         - IGNORE CASE: \"-ignore_case\", variable names compared case insensitive "
        "[default].\n"
        "         - IGNORE MAPS: \"-ignore_maps\", use local implicit instead of global ids.\n"
        "         - IGNORE NANS: \"-ignore_nans\", do not check data for NaNs\n"
        "         - IGNORE DUPLICATES: \"-ignore_dups\", ignore two elements/nodes in same "
        "position...\n"
        "         - STEP OFFSET <val>: \"-TimeStepOffset\", timestep 'val+offset' in file1 matches "
        "'val' in file2.\n"
        "         - STEP OFFSET AUTOMATIC: \"-TA\", automatic determination of timestep offset, "
        "end at same step.\n"
        "         - STEP OFFSET MATCH: \"-TM\", automatic determination of timestep offset, start "
        "at same step.\n"
        "         - INTERPOLATE: \"-interpolate\", interpolate times on file2 to match times on "
        "file1.\n"
        "         - FINAL TIME TOLERANCE <tol>: \"-final_time_tolerance <tol>\", used with "
        "interpolation.\n"
        "         - CALCULATE NORMS: \"-norms\", calculate L1 and L2 norms of variable differences "
        "and output if > 0.0.\n"
        "         - RETURN STATUS: \"-stat\", return exit status of 2 if the files are different. "
        "(default).\n"
        "         - IGNORE STATUS: \"-ignore_status\", exit status is always zero unless an error "
        "occurs.\n"
        "         - PEDANTIC: \"-pedantic\", be more picky about what is a difference.\n\n");
  }
} // namespace
