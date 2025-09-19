// Copyright(C) 1999-2025 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "EJ_SystemInterface.h"
#include "EJ_Version.h"
#include "EJ_vector3d.h"
#include <SL_tokenize.h>
#include <algorithm>
#include <cctype>
#include <copyright.h>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <fmt/format.h>
#include <iosfwd>
#include <utility>
#include <vector>

namespace {
  bool str_equal(const std::string &s1, const std::string &s2)
  {
    return (s1.size() == s2.size()) &&
           std::equal(s1.begin(), s1.end(), s2.begin(),
                      [](char a, char b) { return std::tolower(a) == std::tolower(b); });
  }

  void parse_variable_names(const char *tokens, StringIdVector *variable_list);
  void parse_offset(const char *tokens, std::vector<vector3d> &offset, bool is_offset);
  void parse_integer_list(const char *tokens, std::vector<int> *list);
  void parse_part_list(const char *tokens, std::vector<int> *list);
  void parse_omissions(const char *tokens, Omissions *omissions, const std::string &basename,
                       bool require_ids);
} // namespace

SystemInterface::SystemInterface() { enroll_options(); }

void SystemInterface::enroll_options()
{
  options_.usage("[options] list_of_files_to_join");

  options_.enroll("help", GetLongOption::NoValue, "Print this summary and exit", nullptr);

  options_.enroll("version", GetLongOption::NoValue, "Print version and exit", nullptr);

  options_.enroll("output", GetLongOption::MandatoryValue, "Name of output file to create",
                  "ejoin-out.e", nullptr, true);

  options_.enroll(
      "extract_blocks", GetLongOption::MandatoryValue,
      "Use only the specified part/block pairs. The specification is\n"
      "\t\tp#:block_id1:block_id2,p#:block_id1. For example, to\n"
      "\t\tExtract block ids 1,3,4 from part 1; blocks 2 3 4 from part 2;\n"
      "\t\tand block 8 from part5, specify\n"
      "\t\t\t '-extract_blocks p1:1:3:4,p2:2:3:4,p5:8'\n"
      "\t\tIf an extract is specified, then only that id(s) will be used for that part.",
      nullptr);

  options_.enroll("omit_blocks", GetLongOption::MandatoryValue,
                  "Omit the specified part/block pairs. The specification is\n"
                  "\t\tp#:block_id1:block_id2,p#:block_id1. For example, to\n"
                  "\t\tOmit block id 1,3,4 from part 1; block 2 3 4 from part 2;\n"
                  "\t\tand block 8 from part5, specify\n"
                  "\t\t\t '-omit_blocks p1:1:3:4,p2:2:3:4,p5:8'",
                  nullptr);

  options_.enroll("omit_assemblies", GetLongOption::OptionalValue,
                  "If no value, then don't transfer any assemblies to output file.\n"
                  "\t\tIf just p#,p#,... specified, then omit assemblies on specified parts\n"
                  "\t\tIf p#:id1:id2,p#:id2,id4... then omit the assemblies with the specified\n"
                  "\t\tid in the specified parts.",
                  nullptr, "ALL");

  options_.enroll(
      "omit_part_assemblies", GetLongOption::NoValue,
      "Do not create an assembly for each input part containing the blocks in that part.\n"
      "\t\tDefault is to create the part assemblies.",
      nullptr);

  options_.enroll("omit_nodesets", GetLongOption::OptionalValue,
                  "If no value, then don't transfer any nodesets to output file.\n"
                  "\t\tIf just p#,p#,... specified, then omit sets on specified parts\n"
                  "\t\tIf p#:id1:id2,p#:id2,id4... then omit the sets with the specified\n"
                  "\t\tid in the specified parts.",
                  nullptr, "ALL");

  options_.enroll("omit_sidesets", GetLongOption::OptionalValue,
                  "If no value, then don't transfer any sidesets to output file.\n"
                  "\t\tIf just p#,p#,... specified, then omit sets on specified parts\n"
                  "\t\tIf p#:id1:id2,p#:id2,id4... then omit the sets with the specified\n"
                  "\t\tid in the specified parts.",
                  nullptr, "ALL");

  options_.enroll("convert_nodal_to_nodesets", GetLongOption::MandatoryValue,
                  "For each part listed (or ALL),\n"
                  "\t\tcreate a nodeset containing the nodes of that part\n"
                  "\t\tand output the nodal fields as nodeset fields instead of nodal fields.\n"
                  "\t\tFormat is comma-separated list of parts (1-based), or ALL",
                  nullptr, nullptr, true);

  options_.enroll("match_node_ids", GetLongOption::NoValue,
                  "Combine nodes if their global ids match.", nullptr);

  options_.enroll("match_node_coordinates", GetLongOption::NoValue,
                  "Combine nodes if they are within tolerance distance of each other.", nullptr);

  options_.enroll("match_nodeset_nodes", GetLongOption::MandatoryValue,
                  "Combine nodes in the specified nodeset(s) if they are within\n"
                  "\t\t`tolerance` distance of each other.\n"
                  "\t\tSpecify nodesets in each part as p#:id1:id2,p#:id2,id4...",
                  nullptr);

  options_.enroll("tolerance", GetLongOption::MandatoryValue,
                  "Maximum distance between two nodes to be considered colocated.", nullptr,
                  nullptr, true);

  options_.enroll(
      "combine_nodesets", GetLongOption::NoValue,
      "Input nodesets with the same name will be combined into a single nodeset on output.",
      nullptr);
  options_.enroll("combine_sidesets", GetLongOption::NoValue,
                  "Input sidesets with the same name will be combined into a "
                  "single sideset on output.",
                  nullptr);
  options_.enroll("combine_element_blocks", GetLongOption::NoValue,
                  "Element blocks with the same name and topology will be "
                  "combined into a\n"
                  "\t\tsingle element block on output.",
                  nullptr);

  options_.enroll(
      "nodeset_combines", GetLongOption::MandatoryValue,
      "List of names of output nodesets and the input nodesets which will be combined into that "
      "output.\n"
      "\t\tSyntax: out1:in1,in2,..,inX;out2:inA,inB,...,inZ\n"
      "\t\t        Nodeset 'out1' will contain input nodesets 'in1', 'in2', ..., 'inX'\n"
      "\t\t       Out name separated by ':' from comma-separated list of input.  Multiple outs "
      "separated by ';'",
      nullptr);
  options_.enroll("sideset_combines", GetLongOption::MandatoryValue,
                  "List of names of output sidesets and the input sidesets which will be combined "
                  "into that output.\n"
                  "\t\t See `nodeset_combines` for syntax.",
                  nullptr);
  options_.enroll(
      "element_block_combines", GetLongOption::MandatoryValue,
      "List of names of output element blocks and the input element blocks which will be\n"
      "\t\tcombined into that output. See `nodeset_combines` for syntax.",
      nullptr, nullptr, true);

#if 0
  options_.enroll("match_elem_ids", GetLongOption::NoValue,
                  "Combine elements if their global ids match and they are compatible.\n"
                  "\t\tCompatible = same element type, nodes of the two elements match",
                  nullptr);

  options_.enroll("match_element_coordinates", GetLongOption::NoValue,
                  "Combine elements if their centroids are within tolerance distance\n"
                  "\t\tand they are compatible (same element type, nodes match).",
                  nullptr);
#endif

  options_.enroll(
      "block_prefix", GetLongOption::MandatoryValue,
      "Prefix used on the input block names of second and subsequent part meshes to make them\n"
      "\t\tunique.  Default is 'p'.  Example: block1, p2_block1, p3_block1.",
      "p");

  options_.enroll("offset", GetLongOption::MandatoryValue,
                  "Comma-separated x,y,z offset for coordinates of second and subsequent parts.\n"
                  "\t\tIf there are only 3 values specified, then The offset will be multiplied by "
                  "the part number-1 so:\n"
                  "\t\tP1: no offset; P2: 1x, 1y, 1z; P3: 2x, 2y, 2z; P(n+1): nx, ny, nz\n"
                  "\t\tYou can also specify the offset of specific parts using the syntax:\n"
                  "\t\tpn:xn,yn,zn:pm:xm,ym,zm:pk:xk,yk,zk. (note ':', ',')  For example: `-offset "
                  "p1:1.1,2.2,3.3:p3:2.2,1.0,3.0`\n"
                  "\t\tThe final coordinates are `scale * orig + offset`",
                  nullptr, nullptr, true);

  options_.enroll("scale", GetLongOption::MandatoryValue,
                  "Comma-separated x,y,z scale for coordinates of input parts.\n"
                  "\t\tIf there are only 3 values specified, then The same scale will be used by "
                  "all parts (including the first)\n"
                  "\t\tYou can also specify the scale of specific parts using the syntax:\n"
                  "\t\tpn:xn,yn,zn:pm:xm,ym,zm:pk:xk,yk,zk. (note ':', ',')  For example: `-scale "
                  "p1:1.1,2.2,3.3:p3:2.2,1.0,3.0`\n"
                  "\t\tThe final coordinates are `scale * orig + offset`",
                  nullptr, nullptr, true);

  options_.enroll("steps", GetLongOption::MandatoryValue,
                  "Specify subset of timesteps to transfer to output file.\n"
                  "\t\tFormat is beg:end:step. 1:10:2 --> 1,3,5,7,9\n"
                  "\t\tIf the 'beg' or 'end' is < 0, then it is the \"-Nth\" step...\n"
                  "\t\t-1 is \"first last\" or last, -3 is \"third last\"\n"
                  "\t\tTo copy just the last 3 steps, do: `-steps -3:-1`\n"
                  "\t\tEnter LAST for last step",
                  "1:", nullptr, true);

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
                  "Comma-separated list of nodeset variables to be joined or ALL or NONE.\n"
                  "\t\tVariables can be limited to certain sets by appending a\n"
                  "\t\tcolon followed by the nodeset id.  E.g. -nsetvar sigxx:10:20",
                  nullptr);

  options_.enroll("ssetvar", GetLongOption::MandatoryValue,
                  "Comma-separated list of sideset variables to be joined or ALL or NONE.\n"
                  "\t\tVariables can be limited to certain sidesets by appending a\n"
                  "\t\tcolon followed by the sideset id.  E.g. -ssetvar sigxx:10:20",
                  nullptr);

  options_.enroll(
      "info_records", GetLongOption::OptionalValue,
      "If no value specified or not present,\n"
      "\t\tthen don't transfer any information records to output file.\n"
      "\t\tIf 'p#,p#,...' specified, then only transfer information records on specified parts\n"
      "\t\tIf 'all' specified, then transfer all information records",
      nullptr, "NONE");

  options_.enroll("ignore_element_ids", GetLongOption::NoValue,
                  "Ignore the element id maps on the input database and just use a 1..#elements id "
                  "map on output.\n"
                  "\t\tMuch faster for large models if do not need specific element global ids",
                  nullptr, nullptr, true);

  options_.enroll("netcdf4", GetLongOption::NoValue,
                  "Create output database using the HDF5-based "
                  "netcdf which allows for up to 2.1 GB "
                  "nodes and elements",
                  nullptr);

  options_.enroll("64-bit", GetLongOption::NoValue,
                  "True if forcing the use of 64-bit integers for the output file", nullptr);

  options_.enroll("zlib", GetLongOption::NoValue,
                  "Use the Zlib / libz compression method if compression is enabled (default) "
                  "[exodus only, enables netcdf-4].",
                  nullptr);

  options_.enroll("szip", GetLongOption::NoValue,
                  "Use SZip compression. [exodus only, enables netcdf-4]", nullptr);
  options_.enroll("zstd", GetLongOption::NoValue,
                  "Use Zstd compression. [exodus only, enables netcdf-4, experimental]", nullptr);
  options_.enroll("bzip2", GetLongOption::NoValue,
                  "Use Bzip2 compression. [exodus only, enables netcdf-4, experimental]", nullptr);

  options_.enroll("compress", GetLongOption::MandatoryValue,
                  "Specify the compression level to be used.  Values depend on algorithm:\n"
                  "\t\tzlib/bzip2:  0..9\t\tszip:  even, 4..32\t\tzstd:  -131072..22",
                  nullptr);

  options_.enroll("quantize_nsd", GetLongOption::MandatoryValue,
                  "Use the lossy quantize compression method.\n"
                  "\t\tValue specifies number of digits to "
                  "retain (1..15) [exodus only]",
                  nullptr, nullptr, true);

  options_.enroll("disable_field_recognition", GetLongOption::NoValue,
                  "Do not try to combine scalar fields into higher-order fields such as\n"
                  "\t\tvectors or tensors based on the field suffix",
                  nullptr);

  options_.enroll("copyright", GetLongOption::NoValue, "Show copyright and license data.", nullptr);
}

bool SystemInterface::parse_options(int argc, char **argv)
{
  int option_index = options_.parse(argc, argv);
  if (option_index < 1) {
    return false;
  }

  if (options_.retrieve("help") != nullptr) {
    options_.usage();
    fmt::print(
        stderr,
        "\n\tCan also set options via EJOIN_OPTIONS environment variable.\n"
        "\n\tDocumentation: https://sandialabs.github.io/seacas-docs/sphinx/html/index.html#ejoin\n"
        "\n\t->->-> Send email to sierra-help@sandia.gov for ejoin support.<-<-<-\n");
    exit(EXIT_SUCCESS);
  }

  if (options_.retrieve("version") != nullptr) {
    // Version is printed up front, just exit...
    exit(0);
  }

  if (options_.retrieve("copyright") != nullptr) {
    fmt::print("{}", copyright("2010-2021"));
    exit(EXIT_SUCCESS);
  }

  // Parse remaining options as directory paths.
  // Note that inputFiles_.size() is the number of parts which can be used in
  // parsing of other options...

  if (option_index < argc) {
    while (option_index < argc) {
      inputFiles_.emplace_back(argv[option_index++]);
    }
  }
  else {
    fmt::print(stderr, "\nERROR: no files specified\n\n");
    return false;
  }

  size_t part_count = inputFiles_.size();
  blockOmissions_.resize(part_count);
  blockInclusions_.resize(part_count);
  nodesetOmissions_.resize(part_count);
  sidesetOmissions_.resize(part_count);
  assemblyOmissions_.resize(part_count);
  nodesetMatch_.resize(part_count);
  offset_.resize(part_count);
  scale_.resize(part_count, {1.0, 1.0, 1.0});

  // Get options from environment variable also...
  char *options = getenv("EJOIN_OPTIONS");
  if (options != nullptr) {
    fmt::print(
        stderr,
        "\nThe following options were specified via the EJOIN_OPTIONS environment variable:\n"
        "\t{}\n\n",
        options);
    options_.parse(options, GetLongOption::basename(*argv));
  }

  outputName_  = options_.get_option_value("output", outputName_);
  blockPrefix_ = options_.get_option_value("block_prefix", blockPrefix_);

  {
    const char *temp = options_.retrieve("offset");
    if (temp != nullptr) {
      parse_offset(temp, offset_, true);
    }
  }

  {
    const char *temp = options_.retrieve("scale");
    if (temp != nullptr) {
      parse_offset(temp, scale_, false);
    }
  }

  tolerance_ = options_.get_option_value("tolerance", tolerance_);

  {
    const char *temp = options_.retrieve("steps");
    if (temp != nullptr) {
      parse_step_option(temp);
    }
  }

  {
    const char *temp = options_.retrieve("convert_nodal_to_nodesets");
    if (temp != nullptr) {
      parse_integer_list(temp, &nodesetConvertParts_);
      std::sort(nodesetConvertParts_.begin(), nodesetConvertParts_.end());
    }
  }

  {
    const char *temp = options_.retrieve("info_records");
    if (temp != nullptr) {
      parse_part_list(temp, &infoRecordParts_);
      std::sort(infoRecordParts_.begin(), infoRecordParts_.end());
    }
  }

  {
    const char *temp = options_.retrieve("omit_blocks");
    if (temp != nullptr) {
      parse_omissions(temp, &blockOmissions_, "block", true);
    }
  }

  {
    const char *temp = options_.retrieve("omit_assemblies");
    if (temp != nullptr) {
      if (str_equal("ALL", temp)) {
        omitAssemblies_ = true;
      }
      else {
        parse_omissions(temp, &assemblyOmissions_, "assembly", true);
      }
    }
    else {
      omitAssemblies_ = false;
    }
  }

  {
    const char *temp = options_.retrieve("extract_blocks");
    if (temp != nullptr) {
      parse_omissions(temp, &blockInclusions_, "block", true);
    }
  }

  {
    const char *temp = options_.retrieve("match_nodeset_nodes");
    if (temp != nullptr) {
      parse_omissions(temp, &nodesetMatch_, "nodelist", true);
    }
  }

  {
    const char *temp = options_.retrieve("omit_nodesets");
    if (temp != nullptr) {
      if (str_equal("ALL", temp)) {
        omitNodesets_ = true;
      }
      else {
        parse_omissions(temp, &nodesetOmissions_, "nodelist", false);
      }
    }
    else {
      omitNodesets_ = false;
    }
  }

  {
    const char *temp = options_.retrieve("omit_sidesets");
    if (temp != nullptr) {
      if (str_equal("ALL", temp)) {
        omitSidesets_ = true;
      }
      else {
        parse_omissions(temp, &sidesetOmissions_, "surface", false);
      }
    }
    else {
      omitSidesets_ = false;
    }
  }

  {
    const char *temp = options_.retrieve("gvar");
    if (temp != nullptr) {
      parse_variable_names(temp, &globalVarNames_);
    }
  }

  {
    const char *temp = options_.retrieve("nvar");
    if (temp != nullptr) {
      parse_variable_names(temp, &nodeVarNames_);
    }
  }

  {
    const char *temp = options_.retrieve("evar");
    if (temp != nullptr) {
      parse_variable_names(temp, &elemVarNames_);
    }
  }

  {
    const char *temp = options_.retrieve("nsetvar");
    if (temp != nullptr) {
      parse_variable_names(temp, &nodesetVarNames_);
    }
  }

  {
    const char *temp = options_.retrieve("ssetvar");
    if (temp != nullptr) {
      parse_variable_names(temp, &sidesetVarNames_);
    }
  }

  createAssemblies_        = options_.retrieve("omit_part_assemblies") != nullptr;
  disableFieldRecognition_ = options_.retrieve("disable_field_recognition") != nullptr;
  useNetcdf4_              = options_.retrieve("netcdf4") != nullptr;
  ignoreElementIds_        = options_.retrieve("ignore_element_ids") != nullptr;
  combineNodesets_         = options_.retrieve("combine_nodesets") != nullptr;
  combineSidesets_         = options_.retrieve("combine_sidesets") != nullptr;
  combineElementBlocks_    = options_.retrieve("combine_element_blocks") != nullptr;
  ints64bit_               = options_.retrieve("64-bit") != nullptr;

  {
    const char *temp = options_.retrieve("nodeset_combines");
    if (temp != nullptr) {
      nodesetCombines_ = temp;
    }
  }

  {
    const char *temp = options_.retrieve("sideset_combines");
    if (temp != nullptr) {
      sidesetCombines_ = temp;
    }
  }

  {
    const char *temp = options_.retrieve("element_block_combines");
    if (temp != nullptr) {
      elementBlockCombines_ = temp;
    }
  }

  zlib_ = (options_.retrieve("zlib") != nullptr);
  szip_ = (options_.retrieve("szip") != nullptr);
  zstd_ = (options_.retrieve("zstd") != nullptr);
  bz2_  = (options_.retrieve("bzip2") != nullptr);

  if (szip_ + zlib_ + zstd_ + bz2_ > 1) {
    fmt::print(stderr,
               "ERROR: Only one of 'szip' or 'zlib' or 'zstd' or 'bzip2' can be specified.\n");
  }

  {
    const char *temp = options_.retrieve("compress");
    if (temp != nullptr) {
      compressionLevel_ = std::strtol(temp, nullptr, 10);
      if (!szip_ && !zlib_ && !zstd_ && !bz2_) {
        zlib_ = true;
      }

      if (zlib_ || bz2_) {
        if (compressionLevel_ < 0 || compressionLevel_ > 9) {
          fmt::print(stderr,
                     "ERROR: Bad compression level {}, valid value is between 0 and 9 inclusive "
                     "for gzip/zlib/bzip2 compression.\n",
                     compressionLevel_);
          return false;
        }
      }
      else if (szip_) {
        if (compressionLevel_ % 2 != 0) {
          fmt::print(
              stderr,
              "ERROR: Bad compression level {}. Must be an even value for szip compression.\n",
              compressionLevel_);
          return false;
        }
        if (compressionLevel_ < 4 || compressionLevel_ > 32) {
          fmt::print(stderr,
                     "ERROR: Bad compression level {}, valid value is between 4 and 32 inclusive "
                     "for szip compression.\n",
                     compressionLevel_);
          return false;
        }
      }
    }
  }

  {
    const char *temp = options_.retrieve("quantize_nsd");
    if (temp != nullptr) {
      quantizeNSD_ = std::strtol(temp, nullptr, 10);
      if (!szip_ && !zlib_ && !zstd_ && !bz2_) {
        zlib_ = true;
      }
    }
  }

  if (options_.retrieve("match_node_ids") != nullptr) {
    matchNodeIds_ = true;
    matchNodeXYZ_ = false;
  }
  else {
    matchNodeIds_ = false;
  }

  if (options_.retrieve("match_node_coordinates") != nullptr) {
    matchNodeXYZ_ = true;
    matchNodeIds_ = false;
  }
  else {
    matchNodeXYZ_ = false;
  }

#if 0
  if (options_.retrieve("match_elem_ids")) {
    matchElemIds_ = true;
    matchElemXYZ_ = false;
  } else {
    matchElemIds_ = false;
  }

  if (options_.retrieve("match_elem_coordinates")) {
    matchElemXYZ_ = true;
    matchElemIds_ = false;
  } else {
    matchElemXYZ_ = false;
  }
#endif

  return true;
}

bool SystemInterface::convert_nodes_to_nodesets(int part_number) const
{
  if (nodesetConvertParts_.empty()) {
    return false;
  }
  if (nodesetConvertParts_[0] == 0) {
    return true;
  }
  return std::find(nodesetConvertParts_.cbegin(), nodesetConvertParts_.cend(), part_number) !=
         nodesetConvertParts_.cend();
}

void SystemInterface::parse_step_option(const char *tokens)
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
  //: The step must always be > 0
  //: If the 'from' or 'to' is < 0, then it is the "-Nth" step...
  //: -1 is "first last" or last
  //: -4 is "fourth last step"
  //: To copy just the last 3 steps, do: `-steps -3:-1`

  // Break into tokens separated by ":"

  // Default is given in constructor above...

  if (tokens != nullptr) {
    if (strchr(tokens, ':') != nullptr) {
      // The string contains a separator

      std::array<int, 3> vals{stepMin_, stepMax_, stepInterval_};

      int j = 0;
      for (auto &val : vals) {
        // Parse 'i'th field
        char tmp_str[128];

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
      stepMin_      = vals[0];
      stepMax_      = vals[1];
      stepInterval_ = abs(vals[2]); // step is always positive...
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

void SystemInterface::show_version()
{
  fmt::print("EJoin\n"
             "\t(A code for merging Exodus databases; with or without results data.)\n"
             "\t(Version: {}) Modified: {}\n\n",
             qainfo[2], qainfo[1]);
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
  void parse_variable_names(const char *tokens, StringIdVector *variable_list)
  {
    // Break into tokens separated by ","
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
      std::sort(
          variable_list->begin(), variable_list->end(),
          [](const std::pair<std::string, size_t> &t1, const std::pair<std::string, size_t> &t2) {
            return t1.first < t2.first || (!(t2.first < t1.first) && t1.second < t2.second);
          });
    }
  }

  void parse_offset(const char *tokens, std::vector<vector3d> &offset, bool is_offset)
  {
    // Sets the `offset` or `scale`
    // Break into tokens separated by ","
    if (tokens != nullptr) {
      std::string token_string(tokens);
      if (token_string.find(':') == std::string::npos) {
        // This is specifying just 3 values which are applied to all parts
        StringVector var_list = SLIB::tokenize(token_string, ",");

        // At this point, var_list should contain 3 strings
        // corresponding to the x, y, and z coordinate offsets/scales.
        if (var_list.size() != 3) {
          fmt::print(stderr,
                     "ERROR: Incorrect number of offset components specified--3 required.\n\n");
          return;
        }

        std::string offx = var_list[0];
        std::string offy = var_list[1];
        std::string offz = var_list[2];
        double      x    = std::stod(offx);
        double      y    = std::stod(offy);
        double      z    = std::stod(offz);

        for (size_t i = 0; i < offset.size(); i++) {
          double di = is_offset ? (double)i : 1.0;
          offset[i] = {di * x, di * y, di * z};
        }
      }
      else {
        // Tokens specify explicit offset/scale(s) for 1 or more parts...
        // Form is:  `p1:x1,y1,z1:p3:x3,y3,z3:pN:xN,yN,zN`
        // colon separates part from comma-separated. x,y,z
        // colon also separates the part groups.
        auto groups = SLIB::tokenize(token_string, ":");
        if (groups.size() % 2 != 0) {
          fmt::print(
              stderr,
              "ERROR: Invalid syntax for offset/scale.  Make sure parts are surrounded by ':'\n");
          exit(EXIT_FAILURE);
        }
        for (size_t i = 0; i < groups.size(); i += 2) {
          auto &part_string = groups[i];
          auto &off_string  = groups[i + 1];

          int part_num = -1;
          if (part_string[0] == 'p' || part_string[0] == 'P') {
            part_num = std::stoi(part_string.substr(1));
            if ((size_t)part_num > offset.size()) {
              fmt::print(
                  stderr,
                  "ERROR: Part number too large in offset/scale command ({} must be less or equal "
                  "to {})\n",
                  part_num, offset.size());
              exit(EXIT_FAILURE);
            }
          }
          else {
            fmt::print(stderr,
                       "ERROR: Bad syntax ({}) specifying part number. Use 'p'+ part_number\n"
                       "       For example -offset p1:0,1,2\n",
                       part_string);
            exit(EXIT_FAILURE);
          }

          auto        soff     = SLIB::tokenize(off_string, ",");
          std::string offx     = soff[0];
          std::string offy     = soff[1];
          std::string offz     = soff[2];
          double      x        = std::stod(offx);
          double      y        = std::stod(offy);
          double      z        = std::stod(offz);
          offset[part_num - 1] = {x, y, z};
        }
      }
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

      std::string  token_string(tokens);
      StringVector part_list = SLIB::tokenize(token_string, ",");

      auto I = part_list.begin();
      while (I != part_list.end()) {
        int id = std::stoi(*I);
        (*list).push_back(id);
        ++I;
      }
    }
  }

  void parse_part_list(const char *tokens, std::vector<int> *list)
  {
    // Break into tokens separated by ","
    // Tokens will be of the form "p$" or "all"
    if (tokens != nullptr) {
      if (LowerCase(tokens) == "all") {
        (*list).push_back(0);
        return;
      }

      std::string  token_string(tokens);
      StringVector part_list = SLIB::tokenize(token_string, ",");

      auto I = part_list.begin();
      while (I != part_list.end()) {
        std::string part = *I;
        if (part[0] == 'p' || part[0] == 'P') {
          std::string part_id(part, 1);
          int         part_num = std::stoi(part_id);
          list->push_back(part_num);
        }
        else {
          fmt::print(stderr,
                     "ERROR: Bad syntax ({}) specifying part number. Use 'p'+ part_number\n"
                     "       For example -info_records p1,p2,p7\n",
                     part);
          exit(EXIT_FAILURE);
        }
        ++I;
      }
    }
  }

  void parse_omissions(const char *tokens, Omissions *omissions, const std::string &basename,
                       bool require_ids)
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

    if (tokens == nullptr) {
      return;
    }

    std::string  token_string(tokens);
    StringVector part_block_list = SLIB::tokenize(token_string, ",");

    // Now, for each token in 'part_block_list', split by ":"
    // The result should be a string starting with 'p' followed by an
    // integer indicating the part followed by 1 or more strings which
    // are the block ids of the blocks to be omitted for that part.
    //
    // Parts are 1-based.  Store the results in the 'omissions' set as
    // the pair (part#, block_id).
    auto I = part_block_list.begin();
    while (I != part_block_list.end()) {
      StringVector part_block = SLIB::tokenize(*I, ":");
      if (part_block.empty() || (part_block[0][0] != 'p' && part_block[0][0] != 'P')) {
        fmt::print(stderr, "ERROR: Bad syntax specifying the part number.  Use 'p' + part number\n"
                           "       For example -omit_blocks p1:1:2:3,p2:2:3:4\n");
        exit(EXIT_FAILURE);
      }
      if (require_ids && part_block.size() == 1) {
        fmt::print(stderr,
                   "ERROR: No block ids were found following the part specification.\n"
                   "       for part {}\n",
                   part_block[0]);
        exit(EXIT_FAILURE);
      }

      // Extract the part number...
      std::string part(part_block[0], 1);
      int         part_num = std::stoi(part) - 1;

      // If no blocks specified for a part, then omit all entities for
      // this part.  Since don't know how many entities there are,
      // store the id as '0' to signify all.
      if (part_block.size() == 1) {
        (*omissions)[part_num].emplace_back("ALL");
      }
      else {
        // Get the list of blocks to omit for this part...
        auto J = part_block.begin();
        ++J; // Skip p#
        while (J != part_block.end()) {
          std::string block = *J;
          std::string name  = basename + '_' + block;
          (*omissions)[part_num].push_back(name);
          ++J;
        }
      }
      ++I;
    }
  }
} // namespace
