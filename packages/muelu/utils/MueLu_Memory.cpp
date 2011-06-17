#include <sstream>
#include <fstream>
#include "MueLu_Memory.hpp"

namespace MueLu {
  
  namespace MemUtils {
    
    std::string PrintMemoryUsage() {
      std::ostringstream mem;
      std::ifstream proc("/proc/self/status");
      std::string s;
      while(getline(proc, s), !proc.fail()) {
	if(s.substr(0, 6) == "VmSize") {
	  mem << s;
	  return mem.str();
	}
      }
      return mem.str();
    }
    
  } //namespace MemUtils
  
} //namespace MueLu
