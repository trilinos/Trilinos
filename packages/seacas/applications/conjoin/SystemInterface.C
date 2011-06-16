#include <SystemInterface.h>

#include <iostream>
#include <algorithm>
#include <functional>
#include <vector>

#include <stdlib.h>
#include <limits.h>
#include <string.h>

#include "Version.h"
#include "tokenize.h"

#if defined(__PUMAGON__)
#define NPOS (size_t)-1
#else
#define NPOS std::string::npos
#endif

namespace {
  void parse_variable_names(const char *tokens, StringIdVector *variable_list);
}

Excn::SystemInterface::SystemInterface()
  : outputName_(),
    debugLevel_(0), screenWidth_(0),
    omitNodesets_(false), omitSidesets_(false),
    aliveValue_(-1.0)
{
  enroll_options();
}

Excn::SystemInterface::~SystemInterface() {}

void Excn::SystemInterface::enroll_options()
{
  options_.usage("[options] list_of_files_to_join");

  options_.enroll("help", GetLongOpt::NoValue,
		  "Print this summary and exit", 0);

  options_.enroll("version", GetLongOpt::NoValue,
		  "Print version and exit", NULL);

  options_.enroll("output", GetLongOpt::MandatoryValue,
		  "Name of output file to create",
		  "conjoin-out.e");
  
  options_.enroll("alive_value", GetLongOpt::MandatoryValue,
		  "Value (1 or 0) to indicate that an element is alive, default = 0",
		  "0");
  
  options_.enroll("combine_status_variables", GetLongOpt::MandatoryValue,
		  "The conjoin elem_status variable will be combined\n"
		  "\t\twith the specified status variable ($val) existing on the mesh.\n"
		  "\t\tBoth variables must have the same value (1 or 0) for 'alive'.\n"
		  "\t\t\tIf 1 is alive, then the combined variable is the min of the two values.\n"
		  "\t\t\tIf 0 is alive, then the combined variable is the max of the two values.\n"
		  "\t\tUse the 'alive_value' option to set conjoin's alive value",
		  "");
  options_.enroll("element_status_variable", GetLongOpt::MandatoryValue,
		  "Name to use as element existance status variable;\n"
                  "\t\tmust not exist on input files. If NONE, then not created.\n"
		  "\t\tDefault = elem_status",
		  "elem_status");
  
  options_.enroll("nodal_status_variable", GetLongOpt::MandatoryValue,
		  "Name to use as nodal status variable;\n\t\tmust not exist on input files.\n"
		  "\t\tIf NONE, then not created. Default = node_status",
		  "node_status");
  
  options_.enroll("omit_nodesets", GetLongOpt::NoValue,
		  "Don't transfer nodesets to output file.",
		  NULL);

  options_.enroll("omit_sidesets", GetLongOpt::NoValue,
		  "Don't transfer sidesets to output file.",
		  NULL);

  options_.enroll("gvar", GetLongOpt::MandatoryValue,
		  "Comma-separated list of global variables to be joined or ALL or NONE.",
		  0);

  options_.enroll("evar", GetLongOpt::MandatoryValue,
		  "Comma-separated list of element variables to be joined or ALL or NONE.\n"
		  "\t\tVariables can be limited to certain blocks by appending a\n"
		  "\t\tcolon followed by the block id.  E.g. -evar sigxx:10:20",
		  0);

  options_.enroll("nvar", GetLongOpt::MandatoryValue,
		  "Comma-separated list of nodal variables to be joined or ALL or NONE.",
		  0);

  options_.enroll("nsetvar", GetLongOpt::MandatoryValue,
		  "Comma-separated list of nodeset variables to be joined or ALL or NONE.",
		  0);

  options_.enroll("ssetvar", GetLongOpt::MandatoryValue,
		  "Comma-separated list of sideset variables to be joined or ALL or NONE.",
		  0);

  options_.enroll("debug", GetLongOpt::MandatoryValue,
		  "debug level (values are or'd)\n"
		  "\t\t  1 = timing information.\n"
		  "\t\t  4 = Verbose Element block information.\n"
		  "\t\t  8 = Check consistent nodal coordinates between parts.\n"
		  "\t\t 16 = Verbose Sideset information.\n"
		  "\t\t 32 = Verbose Nodeset information.\n"
		  "\t\t 64 = put exodus library into verbose mode.\n"
		  "\t\t128 = Check consistent global field values between parts.",
		  "0");

  options_.enroll("width", GetLongOpt::MandatoryValue,
		  "Width of output screen, default = 80",
		  "80");
  
  options_.enroll("copyright", GetLongOpt::NoValue,
		  "Show copyright and license data.",
		  NULL);
}

bool Excn::SystemInterface::parse_options(int argc, char **argv)
{
#if (__SUNPRO_CC == 0x500)
  using namespace std;
#endif

  int option_index = options_.parse(argc, argv);
  if ( option_index < 1 )
    return false;

  // Get options from environment variable also...
  char *options = getenv("CONJOIN_OPTIONS");
  if (options != NULL) {
    std::cerr << "\nThe following options were specified via the CONJOIN_OPTIONS environment variable:\n"
	      << "\t" << options << "\n\n";
    options_.parse(options, options_.basename(*argv));
  }

  if (options_.retrieve("help")) {
    options_.usage();
    std::cerr << "\n\t->->-> Send email to gdsjaar@sandia.gov for conjoin support.<-<-<-\n";
    exit(EXIT_SUCCESS);
  }

  if (options_.retrieve("version")) {
    // Version is printed up front, just exit...
    exit(0);
  }
  
  {
    const char *temp = options_.retrieve("debug");
    debugLevel_ = strtol(temp, NULL, 10);
  }

  {
    const char *temp = options_.retrieve("alive_value");
    int value = strtol(temp, NULL, 10);
    if (value == 1 || value == 0) {
      aliveValue_ = value;
    } else {
      std::cerr << "\nInvalid value specified for node and element status."
		<< "\nValid values are '1' or '0'.  Found '" << value << "'\n";
      exit(EXIT_FAILURE);
    }
  }

  {
    const char *temp = options_.retrieve("element_status_variable");
    elementStatusVariable_ = temp;
  }

  {
    const char *temp = options_.retrieve("nodal_status_variable");
    nodalStatusVariable_ = temp;
  }

  {
    const char *temp = options_.retrieve("combine_status_variables");
    meshCombineStatusVariable_ = temp;
  }

  {
    const char *temp = options_.retrieve("width");
    screenWidth_ = strtol(temp, NULL, 10);
  }

  {
    const char *temp = options_.retrieve("output");
    outputName_ = temp;
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

  if (options_.retrieve("omit_nodesets")) {
    omitNodesets_ = true;
  } else {
    omitNodesets_ = false;
  }
  
  if (options_.retrieve("omit_sidesets")) {
    omitSidesets_ = true;
  } else {
    omitSidesets_ = false;
  }
  
  if (options_.retrieve("copyright")) {
    std::cerr << "\n"
	      << "Copyright(C) 2009-2010 Sandia Corporation.\n"
	      << "\n"
	      << "Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,\n"
	      << "the U.S. Government retains certain rights in this software.\n"
	      << "        \n"
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
	      << "    * Neither the name of Sandia Corporation nor the names of its\n"
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
    while (option_index < argc) {
      inputFiles_.push_back(argv[option_index++]);
    }
  } else {
    std::cerr << "\nERROR: no files specified\n\n";
    return false;
  }
  return true;
}

void Excn::SystemInterface::dump(std::ostream &) const
{
}

void Excn::SystemInterface::show_version()
{
  std::cout << qainfo[0] << "\n"
	    << "\t(A code for sequentially appending Exodus II databases. Supercedes conex and conex2.)\n"
	    << "\t(Version: " << qainfo[2] << ") Modified: " << qainfo[1] << '\n';
}

namespace {
  std::string LowerCase(const std::string &name)
  {
    std::string s = name;
    std::transform (s.begin(), s.end(),    // source
		    s.begin(),             // destination
		    ::tolower);            // operation
    return s;
  }

  typedef std::vector<std::string> StringVector;
  bool string_id_sort(const std::pair<std::string,int> &t1,
		      const std::pair<std::string,int> &t2)
  {
    return t1.first < t2.first || (!(t2.first < t1.first) &&
				   t1.second < t2.second);
  }

  void parse_variable_names(const char *tokens, StringIdVector *variable_list)
  {
    // Break into tokens separated by ","
    if (tokens != NULL) {
      std::string token_string(tokens);
      StringVector var_list;
      tokenize(token_string, ",", var_list);
    
      // At this point, var_list is either a single string, or a string
      // separated from 1 or more block ids with ":" delimiter.
      // For example, sigxx:1:10:100 would indicate that the variable
      // "sigxx" should be written only for blocks with id 1, 10, and
      // 100.  "sigxx" would indicate that the variable should be
      // written for all blocks.
      std::vector<std::string>::iterator I = var_list.begin();
      while (I != var_list.end()) {
	StringVector name_id;
	tokenize(*I, ":", name_id);
	std::string var_name = LowerCase(name_id[0]);
	if (name_id.size() == 1) {
	  (*variable_list).push_back(std::make_pair(var_name,0));
	} else {
	  for (size_t i=1; i < name_id.size(); i++) {
	    // Convert string to integer...
	    int id = strtoul(name_id[i].c_str(), NULL, 0);
	    (*variable_list).push_back(std::make_pair(var_name,id));
	  }
	}
	++I;
      }
      // Sort the list...
      std::sort(variable_list->begin(), variable_list->end(), string_id_sort);
    }
  }
}
