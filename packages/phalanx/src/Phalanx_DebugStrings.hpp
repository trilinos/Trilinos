#ifndef PHX_DEBUG_STRINGS_HPP
#define PHX_DEBUG_STRINGS_HPP

#include <string>

namespace PHX {
  
  template<typename ScalarT, typename Traits>
  std::string getTypeString() 
  { 
    return Traits::template TypeString<ScalarT>::value;
  }

}
   
#endif
