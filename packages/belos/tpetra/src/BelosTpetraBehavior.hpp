#ifndef BELOSTPETRABEHAVIOR_HPP
#define BELOSTPETRABEHAVIOR_HPP

#include "Teuchos_EnvVariables.hpp"

namespace Belos {

  class Behavior {
  public:

    static std::string denseMatrixAbstraction() {
      constexpr const std::string_view optionName = "BELOS_DENSE_MATRIX_ABSTRACTION";
      const std::string defaultValue("Teuchos");

      static std::string value_ = defaultValue;
      static bool initialized_  = false;
      return Teuchos::idempotentlyGetEnvironmentVariable(value_, initialized_, optionName, defaultValue);
    }
  };

}
#endif
