#ifndef PANZER_STRING_UTILITIES_HPP
#define PANZER_STRING_UTILITIES_HPP

#include <vector>
#include <string>

namespace panzer {

  void StringTokenizer(std::vector<std::string>& tokens,
		       const std::string& str,
		       const std::string delimiter = ",");
  
}

#endif
