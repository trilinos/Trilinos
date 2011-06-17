#include <sstream>
#include <fstream>
#include "MueLu_Memory.hpp"

#include <iostream> // TODO: remove

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
      
      std::cout << "Free memory of the node:" << PrintMemoryInfo() << std::endl; //TODO: remove

      return mem.str();
    }

    std::string PrintMemoryInfo() {
      std::ostringstream mem;
      std::ifstream proc("/proc/meminfo");
      std::string s;
      while(getline(proc, s), !proc.fail()) {
	if(s.substr(0, 7) == "MemFree") {
	  mem << s;
	  return mem.str();
	}
      }
      return mem.str();
    }
    
  } //namespace MemUtils
  
} //namespace MueLu
