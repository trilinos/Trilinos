#include "SystemInterface.h"

#include <iostream>
#include <algorithm>
#include <functional>
#include <vector>

#include <stdlib.h>
#include <limits.h>
#include <string.h>

#include "Version.h"
#include "tokenize.h"

#include <Ioss_FileInfo.h>

#if defined(__PUMAGON__)
#define NPOS (size_t)-1
#else
#define NPOS std::string::npos
#endif

namespace {
  void parse_variable_names(const char *tokens, StringIdVector *variable_list);
}

SystemInterface::SystemInterface()
  : minimumTime_(0.0), maximumTime_(-1.0), inputFile_(), outputFile_(),
    listVars_(false), fieldSuffix_('_')
{
  enroll_options();
}

SystemInterface::~SystemInterface() {}

void SystemInterface::enroll_options()
{
  options_.usage("[options] input_database output_matlab_script_file\n"
		 "\t\tIf output name not specified, then output file will be the\n"
		 "\t\tbasename of the input file with suffix '.m'");

  options_.enroll("help", GetLongOpt::NoValue,
		  "Print this summary and exit", 0);

  options_.enroll("version", GetLongOpt::NoValue,
		  "Print version and exit", NULL);

  options_.enroll("field_suffix", GetLongOpt::MandatoryValue,
		  "Character used to separate a field component suffix from the field name.\n"
		  "\t\tEnter 'none' for no separator (fieldx, fieldy fieldz).\n"
		  "\t\tDefault = '_' (field_x, field_y, field_z)", "_");
  
  options_.enroll("minimum_time", GetLongOpt::MandatoryValue,
		  "Minimum timestep for which to transfer data to matlab file.", 0);
  
  options_.enroll("maximum_time", GetLongOpt::MandatoryValue,
		  "Maximum timestep for which to transfer data to matlab file.", 0);
  
  options_.enroll("list", GetLongOpt::MandatoryValue,
		  "List global, nodal, element, nodeset, or sideset variables.\n\t\tEnter 'all' to list all types.\n"
		  "\t\tCode exits after listing variable names.", NULL);
  
  options_.enroll("gvar", GetLongOpt::MandatoryValue,
		  "Comma-separated list of global variables to be output or ALL or NONE.",
		  0);

  options_.enroll("evar", GetLongOpt::MandatoryValue,
		  "(NI) Comma-separated list of element variables to be output or ALL or NONE.\n"
		  "\t\tVariables can be limited to certain blocks by appending a\n"
		  "\t\tcolon followed by the block id.  E.g. -evar sigxx:10:20",
		  0);

  options_.enroll("nvar", GetLongOpt::MandatoryValue,
		  "(NI) Comma-separated list of nodal variables to be output or ALL or NONE."
		  "\t\tVariables can be limited to certain nodes by appending a\n"
		  "\t\tcolon followed by the node id.  E.g. -nvar disp:10:20",
		  0);

  options_.enroll("nsetvar", GetLongOpt::MandatoryValue,
		  "(NI) Comma-separated list of nodeset variables to be output or ALL or NONE.",
		  0);

  options_.enroll("ssetvar", GetLongOpt::MandatoryValue,
		  "(NI) Comma-separated list of sideset variables to be output or ALL or NONE.",
		  0);

  options_.enroll("copyright", GetLongOpt::NoValue,
		  "Show copyright and license data.",
		  NULL);
}

bool SystemInterface::parse_options(int argc, char **argv)
{
  int option_index = options_.parse(argc, argv);
  if ( option_index < 1 )
    return false;

  // Get options from environment variable also...
  char *options = getenv("exomatlab");
  if (options != NULL) {
    std::cerr << "\nThe following options were specified via the EXOMATLAB_OPTIONS environment variable:\n"
	      << "\t" << options << "\n\n";
    options_.parse(options, options_.basename(*argv));
  }

  if (options_.retrieve("help")) {
    options_.usage();
    exit(EXIT_SUCCESS);
  }

  if (options_.retrieve("version")) {
    // Version is printed up front, just exit...
    exit(0);
  }
  
  {
    const char *temp = options_.retrieve("field_suffix");
    if (strcmp("none", temp) == 0) {
      fieldSuffix_ = 0;
    } else {
      fieldSuffix_ = temp[0];
    }
  }

  {
    const char *temp = options_.retrieve("list");
    if (temp != NULL) {
      listVars_ = true;
      parse_variable_names(temp, &varsToList_);
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

  {
    const char *temp = options_.retrieve("minimum_time");
    if (temp != NULL)
      minimumTime_  = strtod(temp, NULL);
  }
  
  {
    const char *temp = options_.retrieve("maximum_time");
    if (temp != NULL)
      maximumTime_  = strtod(temp, NULL);
  }
  
  if (options_.retrieve("copyright")) {
    std::cerr << "\n"
	      << "Copyright(C) 2010 Sandia Corporation.\n"
	      << "***not yet submitted through the copyright process***\n\n";
    exit(EXIT_SUCCESS);
  }  
  
  // Parse remaining options as input file.
  if (option_index < argc) {
    inputFile_ = argv[option_index++];
    if (option_index < argc) {
      outputFile_ = argv[option_index++];
    } else {
      outputFile_ = Ioss::FileInfo(inputFile_).basename() + ".m";
    }
  } else {
    options_.usage();
    std::cerr << "\nERROR: no files specified\n\n";
    return false;
  }
  return true;
}

void SystemInterface::show_version()
{
  std::cout << qainfo[0] << "\n"
	    << "\t(A code for outputting exodusII global variable data for use in matlab.)\n"
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
    
      // At this point, var_list is either a single string, or a
      // string separated from 1 or more ids with ":" delimiter.  For
      // example, sigxx:1:10:100 would indicate that the variable
      // "sigxx" should be written only for elements with id 1, 10,
      // and 100.  "sigxx" would indicate that the variable should be
      // written for all elements.
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
