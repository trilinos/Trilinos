#ifndef PHALANX_PRINT_HPP
#define PHALANX_PRINT_HPP

#include "Teuchos_TypeNameTraits.hpp"
#include <string>

namespace PHX {

  /// Default print function for a user defined object, typically used for user defined MDField extents. Users can specialize this.
  template<typename T> std::string print()
  {return Teuchos::demangleName(typeid(T).name());}

}

#endif
