#include <stk_util/diag/FileUtils.hpp>
#include <stk_util/environment/ProgramOptions.hpp>
#include <stk_util/parallel/ParallelComm.hpp>

#include <stk_util/diag/Env.hpp>
#include <stk_util/util/tokenize.hpp>

#include <Ioss_Utils.h>
#include <Ioss_FileInfo.h>
#include <algorithm>

typedef std::vector<std::string> StringVector;

namespace {
  std::string  get_input_file_basename()
  {
    std::string filename;
    
    std::string input_file_name = "stdin";
    if (stk::get_variables_map().count("input-deck")) {
      input_file_name = sierra::Env::getInputFileName();
    }

    Ioss::FileInfo input_file(input_file_name);
    std::string basename = input_file.basename();

    // See if name contains ".aprepro"
    size_t pos = basename.find(".aprepro");
    if (pos != std::string::npos) {
      // Strip if off...
      filename = basename.substr(0,pos) + basename.substr(pos+8);
    } else {
      filename = basename;
    }
    return filename;
  }
}

namespace stk {
  namespace util {
    void filename_substitution(std::string &filename)
    {
      // See if filename contains "%P" which is replaced by the number of processors...
      // Assumes that %P only occurs once...
      // filename is changed.
      size_t pos = filename.find("%P");
      if (pos != std::string::npos) {
	// Found the characters...  Replace with the processor count...
	int num_proc = std::max(1, sierra::Env::parallel_size());
	std::string tmp(filename, 0, pos);
	tmp += Ioss::Utils::to_string(num_proc);
	tmp += filename.substr(pos+2);
	filename = tmp;
      }

      pos = filename.find("%B");
      if (pos != std::string::npos) {
	// Found the characters...  Replace with the input file basename...
	std::string basename = get_input_file_basename();
	if (!basename.empty()) {
	  std::string tmp(filename, 0, pos);
	  tmp += basename;
	  tmp += filename.substr(pos+2);
	  filename = tmp;
	}
      }
    }
  }
}
