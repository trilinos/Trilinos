#ifndef PHX_DEBUG_STRINGS_HPP
#define PHX_DEBUG_STRINGS_HPP

#include <string>

namespace PHX {
  
  template<typename ObjectT, typename Traits>
  std::string getTypeString() 
  { 
    return Traits::template TypeString<ObjectT>::value;
  }

}
   
#endif
