#include <stk_util/util/concat_variable_name.hpp>

namespace stk {
  namespace util {
    bool concat_variable_name(const std::string& first_string,
			      const std::string& second_string,
			      std::string& concat_string) {
      int num_left_paren_first_string = 0;
      int num_right_paren_first_string = 0;
      int num_left_paren_second_string = 0;
      int num_right_paren_second_string = 0;
      for(size_t ifirst=0; ifirst<first_string.length(); ++ifirst) {
	if(first_string[ifirst] == '(') num_left_paren_first_string++;
	if(first_string[ifirst] == ')') num_right_paren_first_string++;
      }
      for(size_t isecond=0; isecond<second_string.length(); ++isecond) {
	if(second_string[isecond] == '(') num_left_paren_second_string++;
	if(second_string[isecond] == ')') num_right_paren_second_string++;
      }
      if(num_left_paren_first_string==1 && num_right_paren_first_string==0 &&
	 num_left_paren_second_string==0 && num_left_paren_first_string==1) {
	concat_string = first_string + "," + second_string;
	return true;
      } else {
	return false;
      }
    }
  }
}
