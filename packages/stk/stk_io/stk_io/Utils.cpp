#include <stk_io/Utils.hpp>
#include <stk_io/DatabasePurpose.hpp>

#include <stk_util/diag/Env.hpp>
#include <stk_util/util/tokenize.hpp>

#include <Ioss_Utils.h>
#include <Ioss_FileInfo.h>

#include <fstream>
#include <algorithm>

typedef std::vector<std::string> StringVector;

namespace {
  inline int to_lower(int c) { return std::tolower(c); }
  
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

  StringVector &get_control_shutdown_filenames()
  {
    static StringVector executable_names;
    static bool first = true;
    
    if (sierra::Env::parallel_rank() == 0) {
      if (first) {
	std::string work_dir = "./";
	if (stk::get_variables_map().count("directory")) {
	  work_dir = stk::get_variables_map()["directory"].as<std::string>();
	}
	executable_names.push_back(work_dir);
	executable_names.push_back("sierra");
	
	// product name may have slashes or capitalization. E.g. "Adagio/Presto" 
	// Try to break into individual names and lowercase...
	std::string product_name = sierra::Env::product_name();
	StringVector tokens;
	stk::util::tokenize(product_name,"/-",tokens);
	for (size_t i=0; i < tokens.size(); i++) {
	  std::string product = tokens[i];
	  std::transform(product.begin(), product.end(), product.begin(), to_lower);
	  executable_names.push_back(product);
	}
	std::string filename = get_input_file_basename();
	if (!filename.empty()) {
	  executable_names.push_back(filename);
	}
	first = false;
      }
    }
    return executable_names;
  }
}

namespace stk {
  namespace io {
    void Utils::filename_substitution(std::string &filename)
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

    void Utils::check_shutdown_request()
    {
      StringVector &executable_names = get_control_shutdown_filenames();

      if (sierra::Env::parallel_rank() == 0) {
	std::string work_dir = executable_names[0];
	for (size_t i=1; i<executable_names.size(); i++) {
	  // Try each file and see if it exists...
	  std::string filename = work_dir + executable_names[i] + ".shutdown";
	  Ioss::FileInfo file_info(filename);
	  if (file_info.exists()) {
	    sierra::Env::outputP0() << "Shutdown requested due to presence of '"
				    << executable_names[i]
				    << ".shutdown' file. Output files will be dumped and execution halted."
				    << std::endl;
	    // OK to call only on 1 processor.
	    sierra::Env::request_shutdown(true);
	  
	    // Delete the file...
	    file_info.remove_file();
	    break;
	  }
	}
      }
    }

    unsigned int Utils::check_control_request(double time, int step)
    {
      StringVector &executable_names = get_control_shutdown_filenames();
      
      unsigned int result = 0;
      if (sierra::Env::parallel_rank() == 0) {
	std::string work_dir = executable_names[0];
	for (size_t i=1; i<executable_names.size(); i++) {
	  // Try each file and see if it exists...
	  std::string filename = work_dir + executable_names[i] + ".control";

	  Ioss::FileInfo file_info(filename);
	  if (file_info.exists()) {
	  
	    // Open the file and read the first line...
	    std::string input_line;
	    {
	      std::ifstream control_file(filename.c_str());
	      std::getline(control_file, input_line);
	      std::transform(input_line.begin(), input_line.end(), input_line.begin(), to_lower);
	    }
	    file_info.remove_file();

	    // Act on the data.  Very simple parsing.  Support dyna and sierra syntax
	    // DUMP RESTART|RESULTS|HEARTBEAT|HISTORY ABORT|STOP|SHUTDOWN|CONTINUE
	    // SW1 - (dump restart shutdown)
	    // SW2 - (dump step    continue)
	    // SW3 - (dump restart continue)
	    // Sw4 - (dump results continue)
	    // If just "dump" -- dump all and continue
	    // DUMP required.  Can combine outputs "restart results history".
	    // Assumes CONTINUE.
	    
	    // Tokenize input line:
	    bool shutdown = false;
	    StringVector tokens;
	    stk::util::tokenize(input_line, " \t", tokens);
	    if (tokens[0] == "dump") {
	      for (size_t j = 1; j < tokens.size(); j++) {
		if (tokens[j] == "restart")
		  result |= stk::io::WRITE_RESTART;
		else if (tokens[j] == "results")
		  result |= stk::io::WRITE_RESULTS;
#if 0
		else if (tokens[j] == "heartbeat")
		  result |= stk::io::WRITE_HEARTBEAT;
		else if (tokens[j] == "history")
		  result |= stk::io::WRITE_HISTORY;
#endif
		else if (tokens[j] == "step" || tokens[j] == "time")
		  ; // Step and time echoed below.
		else if (tokens[j] == "abort" ||
			 tokens[j] == "stop" ||
			 tokens[j] == "shutdown")
		  shutdown = true;
	      }
	    }
	    else if (tokens[0] == "sw1") {
	      result |= stk::io::WRITE_RESTART;
	      shutdown = true;
	    }	      
	    else if (tokens[0] == "sw2") {
	      ; // Step and time echoed below.
	    }	      
	    else if (tokens[0] == "sw3") {
	      result |= stk::io::WRITE_RESTART;
	    }	      
	    else if (tokens[0] == "sw4") {
	      result |= stk::io::WRITE_RESULTS;
	    }
	    // OK to call this on only a single processor.
	    if (shutdown) sierra::Env::request_shutdown(true);

	    sierra::Env::outputP0() << step << "\t  " << time
				    << "\tcontrol file "
				    << executable_names[i]
				    << ".control: '"
				    << input_line << "'"
				    << std::endl;
	  }
	}
      }
      MPI_Bcast(&result, 1, MPI_UNSIGNED, 0, sierra::Env::parallel_comm());
      return result;
    }
  }
}
