// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details
#include "SL_SystemInterface.h"

#include <algorithm>
#include <functional>
#include <iostream>
#include <vector>

#include <climits>
#include <cstdlib>
#include <cstring>

#include "SL_Version.h"
#include <SL_tokenize.h>
#include <copyright.h>
#include <fmt/format.h>

namespace {
  int  get_free_descriptor_count();
  bool str_equal(const std::string &s1, const std::string &s2)
  {
    return (s1.size() == s2.size()) &&
           std::equal(s1.begin(), s1.end(), s2.begin(),
                      [](char a, char b) { return std::tolower(a) == std::tolower(b); });
  }

#if 0
  void parse_variable_names(const char *tokens, StringIdVector *variable_list);
  void parse_integer_list(const char *tokens, std::vector<int> *list);
  void parse_omissions(const char *tokens, Omissions *omissions,
                       const std::string &basename, bool require_ids);
#endif
} // namespace

SystemInterface::SystemInterface() { enroll_options(); }

SystemInterface::~SystemInterface() = default;

void SystemInterface::enroll_options()
{
  options_.usage("[options] list_of_files_to_join");

  options_.enroll("help", GetLongOption::NoValue, "Print this summary and exit", nullptr);

  options_.enroll("version", GetLongOption::NoValue, "Print version and exit", nullptr);

  options_.enroll("processors", GetLongOption::MandatoryValue,
                  "Number of processors to decompose the mesh for", "1");

  options_.enroll("debug", GetLongOption::MandatoryValue, "Debug level: 0, 1, 2, 4 or'd", "0");

  options_.enroll("in_type", GetLongOption::MandatoryValue,
                  "File format for input mesh file (default = exodus)", "exodusii");

  options_.enroll("method", GetLongOption::MandatoryValue,
                  "Decomposition method\n"
                  "\t\t'linear'   : #elem/#proc to each processor\n"
                  "\t\t'scattered': Shuffle elements to each processor (cyclic)\n"
                  "\t\t'random'   : Random distribution of elements, maintains balance\n"
                  "\t\t'rb'       : Metis multilevel recursive bisection\n"
                  "\t\t'kway'     : Metis multilevel k-way graph partitioning\n"
                  "\t\t'file'     : Read element-processor assignment from file",
                  "linear");

  options_.enroll("decomposition_file", GetLongOption::MandatoryValue,
                  "File containing element to processor mapping\n"
                  "\t\twhen decomposition method 'file' specified\n"
                  "\t\tThe file contains multiple lines, each line has 1 or 2 integers.\n"
                  "\t\tIf a single integer, it is the processor for the current element\n"
                  "\t\tIf two integers (count proc), they specify that the next\n"
                  "\t\t\t'count' elements are on processor 'proc'",
                  nullptr);

  options_.enroll("output_path", GetLongOption::MandatoryValue,
                  "Path to where decomposed files will be written.\n"
                  "\t\tThe string %P will be replaced with the processor count\n"
                  "\t\tThe string %M will be replaced with the decomposition method.\n"
                  "\t\tDefault is the location of the input mesh",
                  nullptr);

  options_.enroll("Partial_read_count", GetLongOption::MandatoryValue,
                  "Split the coordinate and connetivity reads into a\n"
                  "\t\tmaximum of this many nodes or elements at a time to reduce memory.",
                  "1000000000");

  options_.enroll("max-files", GetLongOption::MandatoryValue,
                  "Specify maximum number of processor files to write at one time.\n"
                  "\t\tUsually use default value; this is typically used for debugging.",
                  nullptr);

  options_.enroll("netcdf4", GetLongOption::NoValue,
                  "Output database will be a netcdf4 "
                  "hdf5-based file instead of the "
                  "classical netcdf file format",
                  nullptr);

  options_.enroll("64-bit", GetLongOption::NoValue, "Use 64-bit integers on output database",
                  nullptr);

  options_.enroll("netcdf5", GetLongOption::NoValue,
                  "Output database will be a netcdf5 (CDF5) "
                  "file instead of the classical netcdf file format",
                  nullptr);

  options_.enroll("shuffle", GetLongOption::NoValue,
                  "Use a netcdf4 hdf5-based file and use hdf5s shuffle mode with compression.",
                  nullptr);

  options_.enroll("compress", GetLongOption::MandatoryValue,
                  "Specify the hdf5 compression level [0..9] to be used on the output file.",
                  nullptr);

#if 0
  options_.enroll("omit_blocks", GetLongOption::MandatoryValue,
                  "Omit the specified part/block pairs. The specification is\n"
                  "\t\tp#:block_id1:block_id2,p#:block_id1. For example, to\n"
                  "\t\tOmit block id 1,3,4 from part 1; block 2 3 4 from part 2;\n"
                  "\t\tand block 8 from part5, specify\n"
                  "\t\t\t '-omit_blocks p1:1:3:4,p2:2:3:4,p5:8'",
                  nullptr);

  options_.enroll("omit_nodesets", GetLongOption::OptionalValue,
                  "If no value, then don't transfer any nodesets to output file.\n"
                  "\t\tIf just p#,p#,... specified, then omit sets on specified parts\n"
                  "\t\tIf p#:id1:id2,p#:id2,id4... then omit the sets with the specified\n"
                  "\t\tid in the specified parts.",
                  0, "ALL");

  options_.enroll("omit_sidesets", GetLongOption::OptionalValue,
                  "If no value, then don't transfer any sidesets to output file.\n"
                  "\t\tIf just p#,p#,... specified, then omit sets on specified parts\n"
                  "\t\tIf p#:id1:id2,p#:id2,id4... then omit the sets with the specified\n"
                  "\t\tid in the specified parts.",
                  0, "ALL");

  options_.enroll("steps", GetLongOption::MandatoryValue,
                  "Specify subset of timesteps to transfer to output file.\n"
                  "\t\tFormat is beg:end:step. 1:10:2 --> 1,3,5,7,9\n"
                  "\t\tTo only transfer last step, use '-steps LAST'",
                  "1:");

  options_.enroll("disable_field_recognition", GetLongOption::NoValue,
                  "Do not try to combine scalar fields into higher-order fields such as\n"
                  "\t\tvectors or tensors based on the field suffix",
                  nullptr);

#endif
  options_.enroll("contiguous_decomposition", GetLongOption::NoValue,
                  "If the input mesh is contiguous, create contiguous decompositions", nullptr);

  options_.enroll("copyright", GetLongOption::NoValue, "Show copyright and license data.", nullptr);
}

bool SystemInterface::parse_options(int argc, char **argv)
{
#if (__SUNPRO_CC == 0x500)
  using namespace std;
#endif

  int option_index = options_.parse(argc, argv);
  if (option_index < 1) {
    return false;
  }

  if (options_.retrieve("help") != nullptr) {
    options_.usage();
    fmt::print(stderr, "\n\t   Can also set options via SLICE_OPTIONS environment variable.\n");
    fmt::print(stderr, "\n\t->->-> Send email to gsjaardema@gmail.com for slice support.<-<-<-\n");
    exit(EXIT_SUCCESS);
  }

  if (options_.retrieve("version") != nullptr) {
    // Version is printed up front, just exit...
    exit(0);
  }

  if (options_.retrieve("copyright") != nullptr) {
    fmt::print("{}", copyright("2016-2019"));
    exit(EXIT_SUCCESS);
  }

  if (option_index < argc) {
    inputFile_ = argv[option_index++];
  }
  else {
    fmt::print(stderr, "\nERROR: no input mesh file specified\n\n");
    return false;
  }

  if (option_index < argc) {
    nemesisFile_ = argv[option_index++];
  }
  else {
    nemesisFile_ = inputFile_;
  }

  // Get options from environment variable also...
  char *options = getenv("SLICE_OPTIONS");
  if (options != nullptr) {
    fmt::print(
        stderr,
        "\nThe following options were specified via the SLICE_OPTIONS environment variable:\n"
        "\t{}\n\n",
        options);
    options_.parse(options, options_.basename(*argv));
  }

  {
    const char *temp = options_.retrieve("processors");
    processorCount_  = strtoul(temp, nullptr, 0);
  }

  {
    const char *temp  = options_.retrieve("Partial_read_count");
    partialReadCount_ = strtoul(temp, nullptr, 0);
  }

  {
    const char *temp = options_.retrieve("max-files");
    if (temp != nullptr) {
      maxFiles_ = strtoul(temp, nullptr, 0);
    }
    else {
      maxFiles_ = get_free_descriptor_count();
    }
  }

  {
    const char *temp = options_.retrieve("debug");
    debugLevel_      = strtoul(temp, nullptr, 0);
  }

  {
    const char *temp = options_.retrieve("in_type");
    inputFormat_     = temp;
  }

  {
    const char *temp = options_.retrieve("method");
    decompMethod_    = temp;
  }

  {
    if (decompMethod_ == "file") {
      const char *temp = options_.retrieve("decomposition_file");
      if (temp != nullptr) {
        decompFile_ = temp;
      }
      else {
        fmt::print(stderr,
                   "\nThe 'file' decompositon method was specified, but no element "
                   "to processor mapping file was specified via the -decomposition_file option\n");
        return false;
      }
    }
  }

  {
    const char *temp = options_.retrieve("output_path");
    if (temp != nullptr) {
      outputPath_ = temp;
    }
  }

  ints64Bit_ = (options_.retrieve("64-bit") != nullptr);

  if (options_.retrieve("netcdf4") != nullptr) {
    netcdf4_ = true;
    netcdf5_ = false;
  }

  if (options_.retrieve("netcdf5") != nullptr) {
    netcdf5_ = true;
    netcdf4_ = false;
  }

  shuffle_ = (options_.retrieve("shuffle") != nullptr);

  {
    const char *temp = options_.retrieve("compress");
    if (temp != nullptr) {
      compressionLevel_ = std::strtol(temp, nullptr, 10);
    }
  }

#if 0
 {
    const char *temp = options_.retrieve("steps");
    if (temp != nullptr) {
      parse_step_option(temp);
    }
  }

  {
    const char *temp = options_.retrieve("omit_blocks");
    parse_omissions(temp, &blockOmissions_, "block", true);
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

  if (options_.retrieve("disable_field_recognition")) {
    disableFieldRecognition_ = true;
  } else {
    disableFieldRecognition_ = false;
  }
#endif

  if (options_.retrieve("contiguous_decomposition") != nullptr) {
    contig_ = true;
  }
  else {
    contig_ = false;
  }
  return true;
}

void SystemInterface::parse_step_option(const char *tokens)
{
  //: The defined formats for the count attribute are:<br>
  //:  <ul>
  //:    <li><missing> -- default -- 1 <= count <= oo  (all steps)</li>
  //:    <li>"X"                  -- X <= count <= X  (just step X). If X == LAST, last step
  // only</li>
  //:    <li>"X:Y"                -- X to Y by 1</li>
  //:    <li>"X:"                 -- X to oo by 1</li>
  //:    <li>":Y"                 -- 1 to Y by 1</li>
  //:    <li>"::Z"                -- 1 to oo by Z</li>
  //:    <li>"LAST"               -- last step only</li>
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
      for (int &val : vals) {
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
    else if (str_equal("LAST", tokens)) {
      stepMin_ = stepMax_ = -1;
    }
    else {
      // Does not contain a separator, min == max
      stepMin_ = stepMax_ = strtol(tokens, nullptr, 0);
    }
  }
}
void SystemInterface::dump(std::ostream & /*unused*/) const {}

void SystemInterface::show_version()
{
  fmt::print("Slice\n"
             "\t(A code for decomposing finite element meshes for running parallel analyses.)\n"
             "\t(Version: {}) Modified: {}\n",
             qainfo[2], qainfo[1]);
}

namespace {
#if 0
  std::string LowerCase(const std::string &name)
  {
    std::string s = name;
    std::transform(s.begin(), s.end(), // source
                   s.begin(),          // destination
                   ::tolower);         // operation
    return s;
  }

  typedef std::vector<std::string> StringVector;
  bool string_id_sort(const std::pair<std::string, int> &t1, const std::pair<std::string, int> &t2)
  {
    return t1.first < t2.first || (!(t2.first < t1.first) && t1.second < t2.second);
  }

  void parse_variable_names(const char *tokens, StringIdVector *variable_list)
  {
    // Break into tokens separated by ","
    if (tokens != nullptr) {
      std::string token_string(tokens);
      StringVector var_list;
      var_list = SLIB::tokenize(token_string, ",");

      // At this point, var_list is either a single string, or a string
      // separated from 1 or more block ids with ":" delimiter.
      // For example, sigxx:1:10:100 would indicate that the variable
      // "sigxx" should be written only for blocks with id 1, 10, and
      // 100.  "sigxx" would indicate that the variable should be
      // written for all blocks.
      std::vector<std::string>::iterator I = var_list.begin();
      while (I != var_list.end()) {
        StringVector name_id;
        name_id = SLIB::tokenize(*I, ":");
        std::string var_name = LowerCase(name_id[0]);
        if (name_id.size() == 1) {
          (*variable_list).push_back(std::make_pair(var_name,0));
        } else {
          for (size_t i=1; i < name_id.size(); i++) {
            // Convert string to integer...
            int id = std::stoi(name_id[i]);
            (*variable_list).push_back(std::make_pair(var_name,id));
          }
        }
        ++I;
      }
      // Sort the list...
      std::sort(variable_list->begin(), variable_list->end(), string_id_sort);
    }
  }

  void parse_integer_list(const char *tokens, std::vector<int> *list)
  {
    // Break into tokens separated by ","
    if (tokens != nullptr) {
      if (LowerCase(tokens) == "all") {
        (*list).push_back(0);
        return;
      }

      std::string token_string(tokens);
      StringVector part_list;
      part_list = SLIB::tokenize(token_string, ",");

      std::vector<std::string>::iterator I = part_list.begin();
      while (I != part_list.end()) {
        int id = std::stoi(*I);
        (*list).push_back(id);
        ++I;
      }
    }
  }

  void parse_omissions(const char *tokens, Omissions *omissions,
                       const std::string &basename, bool require_ids)
  {
    //  to Omit block id 1,3,4 from part 1; block 2 3 4 from part 2;
    //  and block 8 from part5, specify
    // '-omit_blocks p1:1:3:4,p2:2:3:4,p5:8'

    // Break into tokens separated by "," Each token will then be a
    // ":" separated list of blocks to be omitted for the specified
    // part.

    // If "require_ids" is true, then there must be at least one id
    // following the part specification.  If false, then it is OK to
    // just specify a part number and all entities (typically nset or
    // sset) will be omitted on that part.

    if (tokens == nullptr)
      return;

    std::string token_string(tokens);
    StringVector part_block_list;
    part_block_list = SLIB::tokenize(token_string, ",");

    // Now, for each token in 'part_block_list', split by ":"
    // The result should be a string starting with 'p' followed by an
    // integer indicating the part followed by 1 or more strings which
    // are the block ids of the blocks to be omitted for that part.
    //
    // Parts are 1-based.  Store the results in the 'omissions' set as
    // the pair (part#, block_id).
    std::vector<std::string>::iterator I = part_block_list.begin();
    while (I != part_block_list.end()) {
      StringVector part_block;
      part_block = SLIB::tokenize(*I, ":");
      if (part_block.empty() || (part_block[0][0] != 'p' && part_block[0][0] != 'P')) {
        fmt::print(stderr, "ERROR: Bad syntax specifying the part number.  Use 'p' + part number\n"
                   "       For example -omit_blocks p1:1:2:3,p2:2:3:4\n");
        exit(EXIT_FAILURE);
      }
      if (require_ids && part_block.size() == 1) {
        fmt::print(stderr, "ERROR: No block ids were found following the part specification.\n"
                   "       for part {}\n", part_block[0]);
        exit(EXIT_FAILURE);
      }

      // Extract the part number...
      std::string part(part_block[0],1);
      int part_num = std::stoi(part) - 1;

      // If no blocks specified for a part, then omit all entities for
      // this part.  Since don't know how many entities there are,
      // store the id as '0' to signify all.
      if (part_block.size() == 1) {
        (*omissions)[part_num].push_back("ALL");
      } else {
        // Get the list of blocks to omit for this part...
        std::vector<std::string>::iterator J = part_block.begin(); ++J; // Skip p#
        while (J != part_block.end()) {
          std::string block = *J;
          std::string name = basename + '_'+ block;
          (*omissions)[part_num].push_back(name);
          ++J;
        }
      }
      ++I;
    }
  }
#endif
#include <climits>
#include <unistd.h>

  int get_free_descriptor_count()
  {
// Returns maximum number of files that one process can have open
// at one time. (POSIX)
#ifndef _MSC_VER
    int fdmax = sysconf(_SC_OPEN_MAX);
    if (fdmax == -1) {
      /* POSIX indication that there is no limit on open files... */
      fdmax = INT_MAX;
    }
#else
    int fdmax = _getmaxstdio();
#endif
    // File descriptors are assigned in order (0,1,2,3,...) on a per-process
    // basis.

    // Assume that we have stdin, stdout, stderr, and output exodus
    // file (4 total).

    return fdmax - 4;

    // Could iterate from 0..fdmax and check for the first
    // EBADF (bad file descriptor) error return from fcntl, but that takes
    // too long and may cause other problems.  There is a isastream(filedes)
    // call on Solaris that can be used for system-dependent code.
    //
    // Another possibility is to do an open and see which descriptor is
    // returned -- take that as 1 more than the current count of open files.
    //
  }

} // namespace
