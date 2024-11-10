// Copyright(C) 1999-2024 National Technology & Engineering Solutions
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
#include <open_file_limit.h>

namespace {
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

void SystemInterface::enroll_options()
{
  options_.usage("[options] file_to_split [output_file]");

  options_.enroll("help", GetLongOption::NoValue, "Print this summary and exit", nullptr);

  options_.enroll("in_type", GetLongOption::MandatoryValue,
                  "File format for input mesh file (default = exodus)", "exodusii", nullptr, true);

  options_.enroll("processors", GetLongOption::MandatoryValue,
                  "Number of processors to decompose the mesh for", "1");

  options_.enroll(
      "method", GetLongOption::MandatoryValue,
      "Decomposition method\n"
      "\t\t'linear'   : #elem/#proc to each processor\n"
      "\t\t'scattered': Shuffle elements to each processor (cyclic)\n"
      "\t\t'random'   : Random distribution of elements, maintains balance\n"
#if USE_METIS
      "\t\t'rb'       : Metis multilevel recursive bisection\n"
      "\t\t'kway'     : Metis multilevel k-way graph partitioning\n"
#endif
#if USE_ZOLTAN
      "\t\t'rib'      : Zoltan recursive-inertial-bisection\n"
      "\t\t'rcb'      : Zoltan recursive-coordinate-bisection\n"
      "\t\t'hsfc'     : Zoltan hilbert-space-filling curve\n"
#endif
      "\t\t'variable' : Read element-processor assignment from an element variable\n"
      "\t\t'map'      : Read element-processor assignment from an element map [processor_id]\n"
      "\t\t'file'     : Read element-processor assignment from file",
      "linear");

  options_.enroll(
      "decomposition_name", GetLongOption::MandatoryValue,
      "The name of the element variable (method = `variable`)\n"
      "\t\tor element map (method = `map`) containing the element to processor mapping.\n"
      "\t\tIf no name is specified, then `processor_id` will be used.\n"
      "\t\tIf the name is followed by a ',' and an integer or 'auto', then\n"
      "\t\tthe entries in the map will be divided by the integer value or\n"
      "\t\t(if auto) by `int((max_entry+1)/proc_count)`.",
      nullptr);

  options_.enroll("decomposition_file", GetLongOption::MandatoryValue,
                  "File containing element to processor mapping\n"
                  "\t\twhen decomposition method 'file' specified\n"
                  "\t\tThe file contains multiple lines, each line has 1 or 2 integers.\n"
                  "\t\tIf a single integer, it is the processor for the current element\n"
                  "\t\tIf two integers (count proc), they specify that the next\n"
                  "\t\t\t'count' elements are on processor 'proc'",
                  nullptr);

#if USE_ZOLTAN
  options_.enroll(
      "ignore_x", GetLongOption::NoValue,
      "If using `rcb`, `rib`, or `hsfc`, decompose as if mesh in yz plane. Ignore x dimension.",
      nullptr);
  options_.enroll(
      "ignore_y", GetLongOption::NoValue,
      "If using `rcb`, `rib`, or `hsfc`, decompose as if mesh in xz plane. Ignore y dimension.",
      nullptr);
  options_.enroll(
      "ignore_z", GetLongOption::NoValue,
      "If using `rcb`, `rib`, or `hsfc`, decompose as if mesh in xy plane. Ignore z dimension.",
      nullptr);
#endif
  options_.enroll("contiguous_decomposition", GetLongOption::NoValue,
                  "If the input mesh is contiguous, create contiguous decompositions", nullptr,
                  nullptr);

  options_.enroll("line_decomp", GetLongOption::OptionalValue,
                  "Generate the `lines` or `columns` of elements from the specified surface(s).\n"
                  "\t\tSpecify a comma-separated list of surface/sideset names from which the "
                  "lines will grow.\n"
                  "\t\tOmit or enter 'ALL' for all surfaces in model\n"
                  "\t\tDo not split a line/column across processors.",
                  nullptr, "ALL", true);
  options_.enroll("output_decomp_map", GetLongOption::NoValue,
                  "Do not output the split files; instead write the decomposition information to "
                  "an element map.\n"
                  "\t\tThe name of the map is specified by `-decomposition_name`",
                  nullptr);

  options_.enroll("output_decomp_field", GetLongOption::NoValue,
                  "Do not output the split files; instead write the decomposition information to "
                  "an element map.\n"
                  "\t\tThe name of the field is specified by `-decomposition_name`",
                  nullptr);

  options_.enroll("output_path", GetLongOption::MandatoryValue,
                  "Path to where decomposed files will be written.\n"
                  "\t\tThe string %P will be replaced with the processor count\n"
                  "\t\tThe string %M will be replaced with the decomposition method.\n"
                  "\t\tDefault is the location of the input mesh",
                  nullptr, nullptr, true);

  options_.enroll("Partial_read_count", GetLongOption::MandatoryValue,
                  "Split the coordinate and connectivity reads into a\n"
                  "\t\tmaximum of this many nodes or elements at a time to reduce memory.",
                  "1000000000");

  options_.enroll("max-files", GetLongOption::MandatoryValue,
                  "Specify maximum number of processor files to write at one time.\n"
                  "\t\tUsually use default value; this is typically used for debugging.",
                  nullptr, nullptr, true);

  options_.enroll("netcdf4", GetLongOption::NoValue,
                  "Output database will be a netcdf4 "
                  "hdf5-based file instead of the "
                  "classical netcdf file format",
                  nullptr);

  options_.enroll("netcdf5", GetLongOption::NoValue,
                  "Output database will be a netcdf5 (CDF5) "
                  "file instead of the classical netcdf file format",
                  nullptr);

  options_.enroll("64-bit", GetLongOption::NoValue, "Use 64-bit integers on output database",
                  nullptr, nullptr, true);

  options_.enroll("shuffle", GetLongOption::NoValue,
                  "Use a netcdf4 hdf5-based file and use hdf5s shuffle mode with compression.",
                  nullptr);

  options_.enroll(
      "zlib", GetLongOption::NoValue,
      "Use the Zlib / libz compression method if compression is enabled (default) [exodus only].",
      nullptr);

  options_.enroll("szip", GetLongOption::NoValue,
                  "Use SZip compression. [exodus only, enables netcdf-4]", nullptr);

  options_.enroll("compress", GetLongOption::MandatoryValue,
                  "Specify the hdf5 compression level [0..9] to be used on the output file.",
                  nullptr, nullptr, true);

  options_.enroll("debug", GetLongOption::MandatoryValue,
                  "debug level (values are or'd)\n"
                  "\t\t  1 = timing information.\n"
                  "\t\t  2 = Communication, NodeSet, Sideset information.\n"
                  "\t\t  4 = Progress information in File/Rank.\n"
                  "\t\t  8 = File/Rank Decomposition information.\n"
                  "\t\t 16 = Chain/Line generation/decomp information.\n"
                  "\t\t 32 = Show decomposition histogram (elements / rank).",
                  "0");

  options_.enroll("version", GetLongOption::NoValue, "Print version and exit", nullptr);

  options_.enroll("copyright", GetLongOption::NoValue, "Show copyright and license data.", nullptr);

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
    fmt::print(stderr, "\n\tDocumentation: "
                       "https://sandialabs.github.io/seacas-docs/sphinx/html/index.html#slice\n");
    fmt::print(stderr, "\n\t->->-> Send email to gsjaardema@gmail.com for slice support.<-<-<-\n");
    exit(EXIT_SUCCESS);
  }

  if (options_.retrieve("version") != nullptr) {
    // Version is printed up front, just exit...
    exit(0);
  }

  if (options_.retrieve("copyright") != nullptr) {
    fmt::print("{}", copyright("2016-2021"));
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
    options_.parse(options, GetLongOption::basename(*argv));
  }

  processorCount_   = options_.get_option_value("processors", processorCount_);
  partialReadCount_ = options_.get_option_value("Partial_read_count", partialReadCount_);
  maxFiles_ =
      options_.get_option_value("max-files", open_file_limit() - 1); // -1 for output exodus file.
  debugLevel_   = options_.get_option_value("debug", debugLevel_);
  inputFormat_  = options_.get_option_value("in_type", inputFormat_);
  decompMethod_ = options_.get_option_value("method", decompMethod_);

  {
    if (decompMethod_ == "file") {
      const char *temp = options_.retrieve("decomposition_file");
      if (temp != nullptr) {
        decompFile_ = temp;
      }
      else {
        fmt::print(stderr,
                   "\nThe 'file' decomposition method was specified, but no element "
                   "to processor mapping file was specified via the -decomposition_file option\n");
        return false;
      }
    }
  }

  // Only used in a few methods, but see if set anyway...
  decompVariable_ = options_.get_option_value("decomposition_name", decompVariable_);

  {
    const char *temp = options_.retrieve("line_decomp");
    if (temp != nullptr) {
      lineSurfaceList_ = temp;
      lineDecomp_      = true;
    }
  }

  outputPath_ = options_.get_option_value("output_path", outputPath_);
  ints64Bit_  = (options_.retrieve("64-bit") != nullptr);

  if (options_.retrieve("netcdf4") != nullptr) {
    netcdf4_ = true;
    netcdf5_ = false;
  }

  if (options_.retrieve("netcdf5") != nullptr) {
    netcdf5_ = true;
    netcdf4_ = false;
  }

  shuffle_ = (options_.retrieve("shuffle") != nullptr);

  if (options_.retrieve("szip") != nullptr) {
    szip_ = true;
    zlib_ = false;
  }
  zlib_ = (options_.retrieve("zlib") != nullptr);

  if (szip_ && zlib_) {
    fmt::print(stderr, "ERROR: Only one of 'szip' or 'zlib' can be specified.\n");
  }

  compressionLevel_  = options_.get_option_value("compress", compressionLevel_);
  contig_            = options_.retrieve("contiguous_decomposition") != nullptr;
  outputDecompMap_   = options_.retrieve("output_decomp_map") != nullptr;
  outputDecompField_ = options_.retrieve("output_decomp_field") != nullptr;
#if USE_ZOLTAN
  ignore_x_ = options_.retrieve("ignore_x") != nullptr;
  ignore_y_ = options_.retrieve("ignore_y") != nullptr;
  ignore_z_ = options_.retrieve("ignore_z") != nullptr;
#endif

  if (outputDecompMap_ && outputDecompField_) {
    fmt::print(
        stderr,
        "\nERROR: Cannot specify BOTH `output_decomp_map` and `output_decomp_field` options.\n"
        "       Can only specify one of the two options.\n\n");
    exit(EXIT_FAILURE);
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

  disableFieldRecognition_ = options_.retrieve("disable_field_recognition") != nullptr;
#endif

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

      std::array<int, 3> vals = {stepMin_, stepMax_, stepInterval_};

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

  using StringVector = std::vector<std::string>;
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
} // namespace
