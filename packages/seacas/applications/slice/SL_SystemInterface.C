// Copyright(C) 2016-2017 National Technology & Engineering Solutions of
// Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above
//   copyright notice, this list of conditions and the following
//   disclaimer in the documentation and/or other materials provided
//   with the distribution.
//
// * Neither the name of NTESS nor the names of its
//   contributors may be used to endorse or promote products derived
//   from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
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

#if defined(__PUMAGON__)
#define NPOS (size_t) - 1
#else
#define NPOS std::string::npos
#endif

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
#if 0
  void parse_variable_names(const char *tokens, StringIdVector *variable_list);
  void parse_integer_list(const char *tokens, std::vector<int> *list);
  void parse_omissions(const char *tokens, Omissions *omissions,
		       const std::string &basename, bool require_ids);
#endif
} // namespace

SystemInterface::SystemInterface()
    : decompMethod_("linear"), partialReadCount_(1000000000), processorCount_(1), debugLevel_(0),
      screenWidth_(0), stepMin_(1), stepMax_(INT_MAX), stepInterval_(1), omitNodesets_(false),
      omitSidesets_(false), disableFieldRecognition_(false), contig_(false)
{
  enroll_options();
}

SystemInterface::~SystemInterface() = default;

void SystemInterface::enroll_options()
{
  options_.usage("[options] list_of_files_to_join");

  options_.enroll("help", GetLongOption::NoValue, "Print this summary and exit", nullptr);

  options_.enroll("version", GetLongOption::NoValue, "Print version and exit", nullptr);

  options_.enroll("processors", GetLongOption::MandatoryValue,
                  "Number of processors to decompose the mesh for", "1");

  options_.enroll("debug", GetLongOption::MandatoryValue, "Debug level: 0, 1, 2", "0");

  options_.enroll("input_type", GetLongOption::MandatoryValue,
                  "File format for input mesh file (default = exodus)", "exodusii");

  options_.enroll("method", GetLongOption::MandatoryValue,
                  "Decomposition method\n"
                  "\t\t'linear'   : #elem/#proc to each processor\n"
                  "\t\t'scattered': Shuffle elements to each processor\n"
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
                  "Split the coordinate and connetivity reads into a maximum of this many"
                  " nodes or elements at a time to reduce memory.",
                  "1000000000");

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
    std::cerr << "\n\tCan also set options via SLICE_OPTIONS environment variable.\n";
    std::cerr << "\n\t->->-> Send email to gdsjaar@sandia.gov for slice support.<-<-<-\n";
    exit(EXIT_SUCCESS);
  }

  if (options_.retrieve("version") != nullptr) {
    // Version is printed up front, just exit...
    exit(0);
  }

  if (options_.retrieve("copyright") != nullptr) {
    std::cerr << "\n"
              << "Copyright(C) 2016-2017 National Technology & Engineering Solutions of\n"
              << "Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with\n"
              << "NTESS, the U.S. Government retains certain rights in this software.\n"
              << "\n"
              << "Redistribution and use in source and binary forms, with or without\n"
              << "modification, are permitted provided that the following conditions are\n"
              << "met:\n"
              << "\n"
              << "* Redistributions of source code must retain the above copyright\n"
              << "   notice, this list of conditions and the following disclaimer.\n"
              << "\n"
              << "* Redistributions in binary form must reproduce the above\n"
              << "  copyright notice, this list of conditions and the following\n"
              << "  disclaimer in the documentation and/or other materials provided\n"
              << "  with the distribution.\n"
              << "\n"
              << "* Neither the name of NTESS nor the names of its\n"
              << "  contributors may be used to endorse or promote products derived\n"
              << "  from this software without specific prior written permission.\n"
              << "\n"
              << "THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS\n"
              << "\" AS IS \" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT\n"
              << "LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR\n"
              << "A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT\n"
              << "OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,\n"
              << "SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT\n"
              << "LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,\n"
              << "DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY\n"
              << "THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT\n"
              << "(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE\n"
              << "OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.\n"
              << "\n";
    exit(EXIT_SUCCESS);
  }

  if (option_index < argc) {
    inputFile_ = argv[option_index++];
  }
  else {
    std::cerr << "\nERROR: no input mesh file specified\n\n";
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
    std::cerr
        << "\nThe following options were specified via the SLICE_OPTIONS environment variable:\n"
        << "\t" << options << "\n\n";
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
    const char *temp = options_.retrieve("debug");
    debugLevel_      = strtoul(temp, nullptr, 0);
  }

  {
    const char *temp = options_.retrieve("input_type");
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
        std::cerr << "\nThe 'file' decompositon method was specified, but no element "
                     "to processor mapping file was specified via the -decomposition_file option\n";
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
    else if (case_strcmp("LAST", tokens) == 0) {
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
  std::cout << "Slice"
            << "\n"
            << "\t(A code for decomposing finite element meshes for running parallel analyses.)\n"
            << "\t(Version: " << qainfo[2] << ") Modified: " << qainfo[1] << '\n';
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
	    int id = strtoul(name_id[i].c_str(), nullptr, 0);
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
	int id = strtol((*I).c_str(), nullptr, 0);
	(*list).push_back(id);
	++I;
      }
    }
  }

  void parse_omissions(const char *tokens, Omissions *omissions,
		       const std::string &basename, bool require_ids)
  {
    //	to Omit block id 1,3,4 from part 1; block 2 3 4 from part 2;
    //	and block 8 from part5, specify
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
	std::cerr << "ERROR: Bad syntax specifying the part number.  Use 'p' + part number\n"
		  << "       For example -omit_blocks p1:1:2:3,p2:2:3:4\n";
	exit(EXIT_FAILURE);
      }
      if (require_ids && part_block.size() == 1) {
	std::cerr << "ERROR: No block ids were found following the part specification.\n"
		  << "       for part " << part_block[0] << "\n";
	exit(EXIT_FAILURE);
      }
      
      // Extract the part number...
      std::string part(part_block[0],1);
      int part_num = strtoul(part.c_str(), nullptr, 0) - 1;

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
