// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include "stk_util/util/concat_variable_name.hpp"
#include <cstddef>  // for size_t

namespace stk {
namespace util {
  bool concat_variable_name(const std::string& first_string,
                            const std::string& second_string,
                            std::string& concat_string)
  {
    int num_left_paren_first_string = 0;
    int num_right_paren_first_string = 0;
    int num_left_paren_second_string = 0;
    for(size_t ifirst=0; ifirst<first_string.length(); ++ifirst) {
      if(first_string[ifirst] == '(') num_left_paren_first_string++;
      if(first_string[ifirst] == ')') num_right_paren_first_string++;
    }
    for(size_t isecond=0; isecond<second_string.length(); ++isecond) {
      if(second_string[isecond] == '(') num_left_paren_second_string++;
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
