
#include <stk_util/environment/FileUtils.hpp>
#include <Ioss_FileInfo.h>              // for FileInfo
#include <Ioss_Utils.h>                 // for Utils
#include <stddef.h>                     // for size_t
#include <algorithm>                    // for max
#include <stk_util/environment/EnvData.hpp>  // for EnvData
#include <stk_util/environment/ProgramOptions.hpp>
#include <vector>                       // for vector
#include "boost/program_options/variables_map.hpp"  // for variables_map



typedef std::vector<std::string> StringVector;

namespace {
  std::string  get_input_file_basename()
  {
    std::string filename;
    
    std::string input_file_name = "stdin";
    if (stk::get_variables_map().count("input-deck")) {
      input_file_name = stk::EnvData::instance().m_inputFile;
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
	int num_proc = std::max(1, EnvData::instance().m_parallelSize);
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
