#ifndef PANZER_STRING_UTILITIES_HPP
#define PANZER_STRING_UTILITIES_HPP

#include <vector>
#include <string>

namespace panzer {

  //! Tokenize a string, put tokens in a vector
  void StringTokenizer(std::vector<std::string>& tokens,
		       const std::string& str,
		       const std::string delimiter = ",",bool trim=false);

  //! Turn a vector of tokens into a vector of doubles
  void TokensToDoubles(std::vector<double> & values,const std::vector<std::string> & tokens);

}

#endif
