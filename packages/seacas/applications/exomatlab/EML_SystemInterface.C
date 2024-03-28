/*
 * Copyright(C) 1999-2023 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "EML_CodeTypes.h" // for StringIdVector
#include "EML_SystemInterface.h"
#include "EML_Version.h"   // for qainfo
#include "GetLongOpt.h"    // for GetLongOption, etc
#include "SL_tokenize.h"   // for tokenize
#include <Ioss_FileInfo.h> // for FileInfo
#include <algorithm>       // for sort, transform
#include <cctype>          // for tolower
#include <copyright.h>
#include <cstddef> // for size_t
#include <cstdlib> // for exit, strtod, EXIT_SUCCESS, etc
#include <cstring> // for strcmp
#include <fmt/format.h>
#include <iosfwd>  // for ostream
#include <utility> // for pair, make_pair
#include <vector>  // for vector

namespace {
  void parse_variable_names(const char *tokens, StringIdVector *variable_list);
} // namespace

SystemInterface::SystemInterface() { enroll_options(); }

void SystemInterface::enroll_options()
{
  options_.usage("[options] input_database output_matlab_script_file\n"
                 "\t\tIf output name not specified, then output file will be the\n"
                 "\t\tbasename of the input file with suffix '.m'");

  options_.enroll("help", GetLongOption::NoValue, "Print this summary and exit", nullptr);

  options_.enroll("version", GetLongOption::NoValue, "Print version and exit", nullptr);

  options_.enroll(
      "field_suffix", GetLongOption::MandatoryValue,
      "Character used to separate a field component suffix from the field name.\n"
      "\t\tEnter '_' to treat field_x, field_y, field_z as a 3-component vector 'field'.\n"
      "\t\tDefault = none (field_x, field_y, field_z are different fields)",
      "none");

  options_.enroll("minimum_time", GetLongOption::MandatoryValue,
                  "Minimum timestep for which to transfer data to matlab file.", nullptr);

  options_.enroll("maximum_time", GetLongOption::MandatoryValue,
                  "Maximum timestep for which to transfer data to matlab file.", nullptr);

  options_.enroll("list", GetLongOption::MandatoryValue,
                  "List global, nodal, element, nodeset, or sideset variables.\n\t\tEnter 'all' to "
                  "list all types.\n"
                  "\t\tCode exits after listing variable names.",
                  nullptr);

  options_.enroll("gvar", GetLongOption::MandatoryValue,
                  "Comma-separated list of global variables to be output or ALL or NONE.", "ALL");

#if 0
  options_.enroll("evar", GetLongOption::MandatoryValue,
                  "(Not Yet Implemented) Comma-separated list of element variables to be output or ALL or NONE.\n"
                  "\t\tVariables can be limited to certain blocks by appending a\n"
                  "\t\tcolon followed by the block id.  E.g. -evar sigxx:10:20",
                  nullptr);

  options_.enroll("nvar", GetLongOption::MandatoryValue,
                  "(Not Yet Implemented) Comma-separated list of nodal variables to be output or ALL or NONE.\n"
                  "\t\tVariables can be limited to certain nodes by appending a\n"
                  "\t\tcolon followed by the node id.  E.g. -nvar disp:10:20",
                  nullptr);

  options_.enroll("nsetvar", GetLongOption::MandatoryValue,
                  "(Not Yet Implemented) Comma-separated list of nodeset variables to be output or ALL or NONE.",
                  nullptr);

  options_.enroll("ssetvar", GetLongOption::MandatoryValue,
                  "(Not Yet Implemented) Comma-separated list of sideset variables to be output or ALL or NONE.",
                  nullptr);
#endif

  options_.enroll("copyright", GetLongOption::NoValue, "Show copyright and license data.", nullptr);
}

bool SystemInterface::parse_options(int argc, char **argv)
{
  int option_index = options_.parse(argc, argv);
  if (option_index < 1) {
    return false;
  }

  // Get options from environment variable also...
  char *options = getenv("exomatlab");
  if (options != nullptr) {
    fmt::print(stderr,
               "\nThe following options were specified via the EXOMATLAB_OPTIONS environment "
               "variable:\n\t{}\n\n",
               options);
    options_.parse(options, GetLongOption::basename(*argv));
  }

  if (options_.retrieve("help") != nullptr) {
    options_.usage();
    fmt::print(stderr,
               "\n\tCan also set options via EXOMATLAB_OPTIONS environment variable.\n"
               "\n\tDocumentation: "
               "https://sandialabs.github.io/seacas-docs/sphinx/html/index.html#exomatlab\n"
               "\n\t->->-> Send email to gdsjaar@sandia.gov for exomatlab support.<-<-<-\n");

    exit(EXIT_SUCCESS);
  }

  if (options_.retrieve("version") != nullptr) {
    // Version is printed up front, just exit...
    exit(0);
  }

  {
    const char *temp = options_.retrieve("field_suffix");
    if (temp != nullptr) {
      if (strcmp("none", temp) == 0) {
        // This is ASCII 1 which means it won't be found so
        // vector/tensor won't be recognized by default.
        fieldSuffix_ = 1;
      }
      else {
        fieldSuffix_ = temp[0];
      }
    }
  }

  {
    const char *temp = options_.retrieve("list");
    if (temp != nullptr) {
      listVars_ = true;
      parse_variable_names(temp, &varsToList_);
    }
  }

  {
    const char *temp = options_.retrieve("gvar");
    parse_variable_names(temp, &globalVarNames_);
  }

#if 0
  {
    const char *temp = options_.retrieve("nvar");
    parse_variable_names(temp, &nodeVarNames_);
  }

  {
    const char *temp = options_.retrieve("evar");
    parse_variable_names(temp, &elemVarNames_);
  }

  {
    const char *temp = options_.retrieve("nsetvar");
    parse_variable_names(temp, &nsetVarNames_);
  }

  {
    const char *temp = options_.retrieve("ssetvar");
    parse_variable_names(temp, &ssetVarNames_);
  }
#endif

  {
    const char *temp = options_.retrieve("minimum_time");
    if (temp != nullptr) {
      minimumTime_ = strtod(temp, nullptr);
    }
  }

  {
    const char *temp = options_.retrieve("maximum_time");
    if (temp != nullptr) {
      maximumTime_ = strtod(temp, nullptr);
    }
  }

  if (options_.retrieve("copyright") != nullptr) {
    fmt::print("{}", copyright("2011-2021"));
    exit(EXIT_SUCCESS);
  }

  // Parse remaining options as input file.
  if (option_index < argc) {
    inputFile_ = argv[option_index++];
    if (option_index < argc) {
      outputFile_ = argv[option_index++];
    }
    else {
      outputFile_ = Ioss::FileInfo(inputFile_).basename() + ".m";
    }
  }
  else {
    options_.usage();
    fmt::print(stderr, "\nERROR: no files specified\n\n");
    return false;
  }
  return true;
}

void SystemInterface::show_version()
{
  fmt::print("{}\n"
             "\t(A code for outputting exodusII global variable data for use in matlab.)\n"
             "\t(Version: {}) Modified: {}\n",
             qainfo[0], qainfo[2], qainfo[1]);
}

namespace {
  std::string LowerCase(const std::string &name)
  {
    std::string s = name;
    std::transform(s.begin(), s.end(), // source
                   s.begin(),          // destination
                   ::tolower);         // operation
    return s;
  }

  using StringVector = std::vector<std::string>;
  bool string_id_sort(const std::pair<std::string, int> &t1, const std::pair<std::string, int> &t2)
  {
    return t1.first < t2.first || (!(t2.first < t1.first) && t1.second < t2.second);
  }

  void parse_variable_names(const char *tokens, StringIdVector *variable_list)
  {
    // Break into tokens separated by ","
    if (tokens != nullptr) {
      std::string  token_string(tokens);
      StringVector var_list = SLIB::tokenize(token_string, ",");

      // At this point, var_list is either a single string, or a
      // string separated from 1 or more ids with ":" delimiter.  For
      // example, sigxx:1:10:100 would indicate that the variable
      // "sigxx" should be written only for elements with id 1, 10,
      // and 100.  "sigxx" would indicate that the variable should be
      // written for all elements.
      auto I = var_list.begin();
      while (I != var_list.end()) {
        StringVector name_id  = SLIB::tokenize(*I, ":");
        std::string  var_name = LowerCase(name_id[0]);
        if (name_id.size() == 1) {
          (*variable_list).emplace_back(var_name, 0);
        }
        else {
          for (size_t i = 1; i < name_id.size(); i++) {
            // Convert string to integer...
            int id = std::stoi(name_id[i]);
            (*variable_list).emplace_back(var_name, id);
          }
        }
        ++I;
      }
      // Sort the list...
      std::sort(variable_list->begin(), variable_list->end(), string_id_sort);
    }
  }
} // namespace
