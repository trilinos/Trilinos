#ifndef STK_UTIL_UTIL_CONCAT_VARIABLE_NAME_H
#define STK_UTIL_UTIL_CONCAT_VARIABLE_NAME_H

#include <string>
namespace stk {
namespace util {
  /*!
   *  Parser reads a string like: 'stress(1,2)' as two strings,
   *  'stress(1' and '2)'.  This routines checks if two strings
   *  form a single variable name and if so concats them back
   *  together.  Returns true if the strings were concated,
   *  false otherwise.
   */
  bool concat_variable_name(const std::string& string1, const std::string& string2,
			    std::string& concat_string);
}
}

#endif
