#include "tokenize.hpp"

void stk_classic::util::tokenize(const std::string& str, const std::string& separators,
			 std::vector<std::string>& tokens)
{
  std::string curr_token = "";
  for (size_t i = 0; i < str.length(); ++i) {
    char curr_char = str[i];

    // determine if current character is a separator
    bool is_separator = false;
    for (size_t j = 0; j < separators.length(); ++j) {
      if (curr_char == separators[j]) {
        is_separator = true;
        break;
      }
    }

    if (is_separator && curr_token != "") {
      // we just completed a token
      tokens.push_back(curr_token);
      curr_token.clear();
    }
    else if (!is_separator) {
      curr_token += curr_char;
    }
  }
  if (curr_token != "") {
    tokens.push_back(curr_token);
  }
}
