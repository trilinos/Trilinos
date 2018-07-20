/*
 * Copyright(C) 2010-2017 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *     * Neither the name of NTESS nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */
#include "EP_SystemInterface.h"
#include "EP_Version.h"  // for qainfo
#include "GetLongOpt.h"  // for GetLongOption, etc
#include "SL_tokenize.h" // for tokenize
#include <algorithm>     // for sort, transform
#include <cctype>        // for tolower
#include <climits>       // for INT_MAX
#include <cstddef>       // for size_t
#include <cstdlib>       // for strtol, abs, exit, strtoul, etc
#include <cstring>       // for strchr, strlen
#include <iostream>      // for operator<<, basic_ostream, etc
#include <string>        // for string, char_traits, etc
#include <utility>       // for pair, make_pair
#include <vector>        // for vector

namespace {
  int case_strcmp(const std::string &s1, const std::string &s2)
  {
    const char *c1 = s1.c_str();
    const char *c2 = s2.c_str();
    for (;; c1++, c2++) {
      if (std::tolower(*c1) != std::tolower(*c2)) {
        return (std::tolower(*c1) - std::tolower(*c2));
      }
      if (*c1 == '\0') {
        return 0;
      }
    }
  }
  void parse_variable_names(const char *tokens, Excn::StringIdVector *variable_list);
} // namespace

Excn::SystemInterface::SystemInterface()
    : inExtension_(""), outExtension_(""), cwd_(""), rootDirectory_(), subDirectory_(""),
      basename_(""), raidOffset_(0), raidCount_(0), processorCount_(1), startPart_(0),
      partCount_(-1), debugLevel_(0), screenWidth_(0), stepMin_(1), stepMax_(INT_MAX),
      stepInterval_(1), subcycle_(-1), cycle_(-1), compressData_(0), sumSharedNodes_(false),
      addProcessorId_(false), mapIds_(true), omitNodesets_(false), omitSidesets_(false),
      useNetcdf4_(false), append_(false), intIs64Bit_(false), subcycleJoin_(false),
      outputSharedNodes_(false)
{
  enroll_options();
}

Excn::SystemInterface::~SystemInterface() = default;

void Excn::SystemInterface::enroll_options()
{
  options_.usage("[options] basename");

  options_.enroll("help", GetLongOption::NoValue, "Print this summary and exit", nullptr);

  options_.enroll("version", GetLongOption::NoValue, "Print version and exit", nullptr);

  options_.enroll("auto", GetLongOption::NoValue,
                  "Automatically set Root, Proc, Ext from filename 'Root/basename.ext.#p.00'.",
                  nullptr);
  options_.enroll("map", GetLongOption::NoValue,
                  "Map element ids to original order if possible [default]", nullptr);

  options_.enroll("nomap", GetLongOption::NoValue, "Do not map element ids to original order",
                  nullptr);

  options_.enroll("extension", GetLongOption::MandatoryValue,
                  "Exodus database extension for the input files", "e");

  options_.enroll("output_extension", GetLongOption::MandatoryValue,
                  "Exodus database extension for the output file", nullptr);

  options_.enroll("offset", GetLongOption::MandatoryValue, "Raid Offset", nullptr);

  options_.enroll("raid_count", GetLongOption::MandatoryValue, "Number of raids", "0");

  options_.enroll("processor_count", GetLongOption::MandatoryValue, "Number of processors", "1");

  options_.enroll("current_directory", GetLongOption::MandatoryValue, "Current Directory", ".");

  options_.enroll("Root_directory", GetLongOption::MandatoryValue, "Root directory", nullptr);

  options_.enroll("Subdirectory", GetLongOption::MandatoryValue,
                  "subdirectory containing input exodusII files", nullptr);

  options_.enroll("width", GetLongOption::MandatoryValue, "Width of output screen, default = 80",
                  "80");

  options_.enroll("add_processor_id", GetLongOption::NoValue,
                  "Add 'processor_id' element variable to the output file", nullptr);

  options_.enroll("netcdf4", GetLongOption::NoValue, "Create output database using the HDF5-based "
                                                     "netcdf which allows for up to 2.1 GB "
                                                     "nodes/elements",
                  nullptr);

  options_.enroll("large_model", GetLongOption::NoValue,
                  "Create output database using the HDF5-based netcdf which allows for up to 2.1 "
                  "GB nodes/elements (deprecated; use netcdf4 instead)",
                  nullptr);

  options_.enroll("append", GetLongOption::NoValue,
                  "Append to database instead of opening a new database.\n"
                  "\t\tTimestep transfer will start after last timestep on database",
                  nullptr);

  options_.enroll("64", GetLongOption::NoValue,
                  "The output database will be written in the 64-bit integer mode", nullptr);

  options_.enroll("compress_data", GetLongOption::MandatoryValue,
                  "The output database will be written using compression (netcdf-4 mode only).\n"
                  "\t\tValue ranges from 0 (no compression) to 9 (max compression) inclusive.",
                  nullptr);

  options_.enroll("steps", GetLongOption::MandatoryValue,
                  "Specify subset of timesteps to transfer to output file.\n"
                  "\t\tFormat is beg:end:step. 1:10:2 --> 1,3,5,7,9\n"
                  "\t\tEnter LAST for last step",
                  "1:");

  options_.enroll("Part_count", GetLongOption::MandatoryValue,
                  "How many pieces (files) of the model should be joined.", "0");

  options_.enroll("start_part", GetLongOption::MandatoryValue, "Start with piece {n} (file)", "0");

  options_.enroll("subcycle", GetLongOption::OptionalValue,
                  "Subcycle. Create $val subparts if $val is specified.\n"
                  "\t\tOtherwise, create multiple parts each of size 'Part_count'.\n"
                  "\t\tThe subparts can then be joined by a subsequent run of epu.\n"
                  "\t\tUseful if the maximum number of open files is less\n"
                  "\t\tthan the processor count.",
                  nullptr, "0");

  options_.enroll("cycle", GetLongOption::MandatoryValue,
                  "Cycle number. If subcycle # is specified, then only execute\n"
                  "\t\tcycle $val ($val < #).  The cycle number is 0-based.",
                  "-1");

  options_.enroll("join_subcycles", GetLongOption::NoValue,
                  "If -subcycle is specified, then after the subcycle files are processed,\n"
                  "\t\trun epu one more time and join the subcycle files into a single file.",
                  nullptr);

  options_.enroll("sum_shared_nodes", GetLongOption::NoValue,
                  "The nodal results data on all shared nodes (nodes on processor boundaries)\n"
                  "\t\twill be the sum of the individual nodal results data on each shared node.\n"
                  "\t\tThe default behavior assumes that the values are equal.",
                  nullptr);

  options_.enroll("gvar", GetLongOption::MandatoryValue,
                  "Comma-separated list of global variables to be joined or ALL or NONE.", nullptr);

  options_.enroll("evar", GetLongOption::MandatoryValue,
                  "Comma-separated list of element variables to be joined or ALL or NONE.\n"
                  "\t\tVariables can be limited to certain blocks by appending a\n"
                  "\t\tcolon followed by the block id.  E.g. -evar sigxx:10:20",
                  nullptr);

  options_.enroll("nvar", GetLongOption::MandatoryValue,
                  "Comma-separated list of nodal variables to be joined or ALL or NONE.", nullptr);

  options_.enroll("nsetvar", GetLongOption::MandatoryValue,
                  "Comma-separated list of nodeset variables to be joined or ALL or NONE.",
                  nullptr);

  options_.enroll("ssetvar", GetLongOption::MandatoryValue,
                  "Comma-separated list of sideset variables to be joined or ALL or NONE.",
                  nullptr);

  options_.enroll("omit_nodesets", GetLongOption::NoValue,
                  "Don't transfer nodesets to output file.", nullptr);

  options_.enroll("omit_sidesets", GetLongOption::NoValue,
                  "Don't transfer sidesets to output file.", nullptr);

  options_.enroll("output_shared_nodes", GetLongOption::NoValue,
                  "Output list of shared nodes and the processors they are shared with.", nullptr);

  options_.enroll("debug", GetLongOption::MandatoryValue,
                  "debug level (values are or'd)\n"
                  "\t\t  1 = timing information.\n"
                  "\t\t  2 = Check consistent nodal field values between processors.\n"
                  "\t\t  4 = Verbose Element block information.\n"
                  "\t\t  8 = Check consistent nodal coordinates between processors.\n"
                  "\t\t 16 = Verbose Sideset information.\n"
                  "\t\t 32 = Verbose Nodeset information.\n"
                  "\t\t 64 = put exodus library into verbose mode.\n"
                  "\t\t128 = Check consistent global field values between processors.",
                  "0");

  options_.enroll("copyright", GetLongOption::NoValue, "Show copyright and license data.", nullptr);
}

bool Excn::SystemInterface::parse_options(int argc, char **argv)
{
  // Get options from environment variable also...
  char *options = getenv("EPU_OPTIONS");
  if (options != nullptr) {
    std::cout
        << "\nThe following options were specified via the EPU_OPTIONS environment variable:\n"
        << "\t" << options << "\n\n";
    options_.parse(options, options_.basename(*argv));
  }

  int option_index = options_.parse(argc, argv);
  if (option_index < 1) {
    return false;
  }

  if (options_.retrieve("help") != nullptr) {
    options_.usage();
    std::cout << "\n\tCan also set options via EPU_OPTIONS environment variable.\n\n"
              << "\tWrites: current_directory/basename.suf\n"
              << "\tReads:  root#o/sub/basename.suf.#p.0 to\n"
              << "\t\troot(#o+#p)%#r/sub/basename.suf.#p.#p\n";
    std::cout << "\n\t->->-> Send email to gdsjaar@sandia.gov for epu support.<-<-<-\n";
    exit(EXIT_SUCCESS);
  }

  if (options_.retrieve("version") != nullptr) {
    // Version is printed up front, just exit...
    exit(0);
  }

  {
    const char *temp = options_.retrieve("extension");
    if (temp != nullptr) {
      inExtension_ = temp;
    }
  }

  {
    const char *temp = options_.retrieve("output_extension");
    if (temp != nullptr) {
      outExtension_ = temp;
    }
  }

  {
    const char *temp = options_.retrieve("offset");
    if (temp != nullptr) {
      raidOffset_ = strtol(temp, nullptr, 10);
    }
  }

  {
    const char *temp = options_.retrieve("raid_count");
    if (temp != nullptr) {
      raidCount_ = strtol(temp, nullptr, 10);
    }
  }

  {
    const char *temp = options_.retrieve("processor_count");
    if (temp != nullptr) {
      processorCount_ = strtol(temp, nullptr, 10);
    }
  }

  {
    const char *temp = options_.retrieve("Part_count");
    if (temp != nullptr) {
      partCount_ = strtol(temp, nullptr, 10);
    }
  }

  {
    const char *temp = options_.retrieve("start_part");
    if (temp != nullptr) {
      startPart_ = strtol(temp, nullptr, 10);
    }
  }

  {
    const char *temp = options_.retrieve("current_directory");
    if (temp != nullptr) {
      cwd_ = temp;
    }
  }

  {
    const char *temp = options_.retrieve("Root_directory");
    if (temp != nullptr) {
      rootDirectory_ = temp;
    }
  }

  {
    const char *temp = options_.retrieve("Subdirectory");
    if (temp != nullptr) {
      subDirectory_ = temp;
    }
  }

  {
    const char *temp = options_.retrieve("debug");
    if (temp != nullptr) {
      debugLevel_ = strtol(temp, nullptr, 10);
    }
  }

  {
    const char *temp = options_.retrieve("width");
    if (temp != nullptr) {
      screenWidth_ = strtol(temp, nullptr, 10);
    }
  }

  {
    const char *temp = options_.retrieve("steps");
    if (temp != nullptr) {
      parse_step_option(temp);
    }
  }

  {
    const char *temp = options_.retrieve("gvar");
    parse_variable_names(temp, &globalVarNames_);
  }

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

  if (options_.retrieve("add_processor_id") != nullptr) {
    addProcessorId_ = true;
  }
  else {
    addProcessorId_ = false;
  }

  if (options_.retrieve("large_model") != nullptr) {
    useNetcdf4_ = true;
    std::cerr << "\nWARNING: the -large_model option is deprecated; please use -netcdf4 instead.\n";
  }

  if (options_.retrieve("netcdf4") != nullptr) {
    useNetcdf4_ = true;
  }

  if (options_.retrieve("append") != nullptr) {
    append_ = true;
  }

  if (options_.retrieve("64") != nullptr) {
    intIs64Bit_ = true;
  }

  {
    const char *temp = options_.retrieve("compress_data");
    if (temp != nullptr) {
      compressData_ = strtol(temp, nullptr, 10);
    }
  }

  if (options_.retrieve("sum_shared_nodes") != nullptr) {
    sumSharedNodes_ = true;
  }

  if (options_.retrieve("append") != nullptr) {
    append_ = true;
  }

  {
    const char *temp = options_.retrieve("subcycle");
    if (temp != nullptr) {
      subcycle_ = strtol(temp, nullptr, 10);
    }
  }

  {
    const char *temp = options_.retrieve("cycle");
    if (temp != nullptr) {
      cycle_ = strtol(temp, nullptr, 10);
    }
  }

  if (options_.retrieve("join_subcycles") != nullptr) {
    subcycleJoin_ = true;
  }

  if (options_.retrieve("map") != nullptr) {
    mapIds_ = true;
  }

  if (options_.retrieve("nomap") != nullptr) {
    mapIds_ = false;
  }

  if (options_.retrieve("omit_nodesets") != nullptr) {
    omitNodesets_ = true;
  }
  else {
    omitNodesets_ = false;
  }

  if (options_.retrieve("omit_sidesets") != nullptr) {
    omitSidesets_ = true;
  }
  else {
    omitSidesets_ = false;
  }

  if (options_.retrieve("output_shared_nodes") != nullptr) {
    outputSharedNodes_ = true;
  }

  if (options_.retrieve("copyright") != nullptr) {
    std::cout << "\n"
	      << "Copyright(C) 2010-2017 National Technology & Engineering Solutions of Sandia, LLC (NTESS).\n"
	      << "Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government\n"
	      << "retains certain rights in this software\n."
	      << "\n"
	      << "Redistribution and use in source and binary forms, with or without\n"
	      << "modification, are permitted provided that the following conditions are\n"
	      << "met:\n"
	      << "\n"
	      << "    * Redistributions of source code must retain the above copyright\n"
	      << "      notice, this list of conditions and the following disclaimer.\n"
	      << "\n"
	      << "    * Redistributions in binary form must reproduce the above\n"
	      << "      copyright notice, this list of conditions and the following\n"
	      << "      disclaimer in the documentation and/or other materials provided\n"
	      << "      with the distribution.\n"
	      << "\n"
	      << "    * Neither the name of NTESS nor the names of its\n"
	      << "      contributors may be used to endorse or promote products derived\n"
	      << "      from this software without specific prior written permission.\n"
	      << "\n"
	      << "THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS\n"
	      << "'AS IS' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT\n"
	      << "LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR\n"
	      << "A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT\n"
	      << "OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,\n"
	      << "SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT\n"
	      << "LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,\n"
	      << "DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY\n"
	      << "THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT\n"
	      << "(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE\n"
	      << "OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.\n\n";
    exit(EXIT_SUCCESS);
  }

  // Parse remaining options as directory paths.
  if (option_index < argc) {
    basename_ = argv[option_index];

    if (options_.retrieve("auto") != nullptr) {
      // Determine Root, Proc, Extension, and Basename automatically
      // by parsing the basename_ entered by the user.  Assumed to be
      // in the form: "/directory/sub/basename.ext.#proc.34"
      bool success = decompose_filename(basename_);
      if (!success) {
        std::cerr
            << "\nERROR: (EPU) If the '-auto' option is specified, the basename must specify an "
               "existing filename.\n"
            << "       The entered basename does not contain an extension or processor count.\n";
        return false;
      }
    }
  }
  else {
    std::cerr << "\nERROR: (EPU) basename not specified\n\n";
    return false;
  }
  return true;
}

void Excn::SystemInterface::dump(std::ostream & /*unused*/) const {}

std::string Excn::SystemInterface::output_suffix() const
{
  if (outExtension_ == "") {
    return inExtension_;
  }

  return outExtension_;
}

void Excn::SystemInterface::show_version()
{
  std::cout << qainfo[0] << "\n"
            << "\t(Out of Many One -- see http://www.greatseal.com/mottoes/unum.html)\n"
            << "\tExodusII Parallel Unification Program\n"
            << "\t(Version: " << qainfo[2] << ") Modified: " << qainfo[1] << '\n';
}

void Excn::SystemInterface::parse_step_option(const char *tokens)
{
  //: The defined formats for the count attribute are:<br>
  //:  <ul>
  //:    <li><missing> -- default -- 1 <= count <= oo  (all steps)</li>
  //:    <li>"X"                  -- X <= count <= X  (just step X) LAST for last step.</li>
  //:    <li>"X:Y"                -- X to Y by 1</li>
  //:    <li>"X:"                 -- X to oo by 1</li>
  //:    <li>":Y"                 -- 1 to Y by 1</li>
  //:    <li>"::Z"                -- 1 to oo by Z</li>
  //:  </ul>
  //: The count and step must always be >= 0

  // Break into tokens separated by ":"

  // Default is given in constructor above...

  if (tokens != nullptr) {
    if (strchr(tokens, ':') != nullptr) {
      // The string contains a separator

      int vals[3];
      vals[0] = stepMin_;
      vals[1] = stepMax_;
      vals[2] = stepInterval_;

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
          val = strtoul(tmp_str, nullptr, 0);
        }

        if (tokens[j++] == '\0') {
          break; // Reached end of string
        }
      }
      stepMin_      = abs(vals[0]);
      stepMax_      = abs(vals[1]);
      stepInterval_ = abs(vals[2]);
    }
    else if (case_strcmp("LAST", tokens) == 0) {
      stepMin_ = stepMax_ = -1;
    }
    else {
      // Does not contain a separator, min == max
      stepMin_ = stepMax_ = strtol(tokens, nullptr, 0);
    }
  }
}

bool Excn::SystemInterface::decompose_filename(const std::string &cs)
{
  std::string s(cs);
  // Decompose the input string 's' into processor count, extension, basename, and root.
  // Input string should be of the form:
  // "root/basename.ext.proc.nn"
  // 'root/' is optional and is all characters preceding the last (if any) '/';
  // 'basename' can contain multiple '.'

  // NOTE: This used to use the tokenize function, but that didn't work since we need
  // to handle leading and embedded '..' which tokenize threw away...

  // Get rid of the 'nn' which is not used at this time...
  size_t ind = s.find_last_of('.', std::string::npos); // last '.'
  if (ind == std::string::npos) {
    return false;
  }
  s.erase(ind);

  // Now find the processor count...
  ind = s.find_last_of('.', std::string::npos);
  if (ind == std::string::npos) {
    return false;
  }

  std::string tmp = s.substr(ind + 1); // Skip the '.'
  processorCount_ = strtol(tmp.c_str(), nullptr, 10);
  if (processorCount_ <= 0) {
    std::cerr << "\nERROR: (EPU) Invalid processor count specified: '" << processorCount_
              << "'. Must be greater than zero.\n";
    return false;
  }
  s.erase(ind);

  // Should now be an extension...
  ind = s.find_last_of('.', std::string::npos);
  if (ind == std::string::npos) {
    inExtension_ = "";
  }
  else {
    inExtension_ = s.substr(ind + 1);
    s.erase(ind);
  }

  // Remainder of 's' consists of the basename_ and the rootDirectory_
  // If there is no '/', then it is all basename_; otherwise the
  // basename_ is the portion following the '/' and the rootDirectory_
  // is the portion preceding the '/'
  ind = s.find_last_of('/', std::string::npos);
  if (ind != std::string::npos) {
    basename_      = s.substr(ind + 1, std::string::npos);
    rootDirectory_ = s.substr(0, ind);
  }
  else {
    basename_ = s;
  }

  std::cout << "\nThe following options were determined automatically:\n"
            << "\t basename = '" << basename_ << "'\n"
            << "\t-processor_count " << processorCount_ << "\n"
            << "\t-extension " << inExtension_ << "\n"
            << "\t-Root_directory " << rootDirectory_ << "\n\n";

  return true;
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

  void parse_variable_names(const char *tokens, Excn::StringIdVector *variable_list)
  {
    // Break into tokens separated by ","

    // Value of num_vars includes optional add_processor_id

    if (tokens != nullptr) {
      std::string  token_string(tokens);
      StringVector var_list = SLIB::tokenize(token_string, ",");

      // At this point, var_list is either a single string, or a string
      // separated from 1 or more block ids with ":" delimiter.
      // For example, sigxx:1:10:100 would indicate that the variable
      // "sigxx" should be written only for blocks with id 1, 10, and
      // 100.  "sigxx" would indicate that the variable should be
      // written for all blocks.
      auto I = var_list.begin();
      while (I != var_list.end()) {
        StringVector name_id  = SLIB::tokenize(*I, ":");
        std::string  var_name = LowerCase(name_id[0]);
        if (name_id.size() == 1) {
          (*variable_list).push_back(std::make_pair(var_name, 0));
        }
        else {
          for (size_t i = 1; i < name_id.size(); i++) {
            // Convert string to integer...
            int id = strtoul(name_id[i].c_str(), nullptr, 0);
            (*variable_list).push_back(std::make_pair(var_name, id));
          }
        }
        ++I;
      }
      // Sort the list...
      std::sort(variable_list->begin(), variable_list->end(), string_id_sort);
    }
  }
} // namespace
