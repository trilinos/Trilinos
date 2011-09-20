#include "MueLu_VerbosityLevel.hpp"
#include "MueLu_Exceptions.hpp"

namespace MueLu {
  
  VerbLevel toMueLuVerbLevel(const Teuchos::EVerbosityLevel verbLevel) {
      switch(verbLevel)
        {  
        case Teuchos::VERB_NONE:
          return None;
        case Teuchos::VERB_DEFAULT:
          return Default;
        case Teuchos::VERB_LOW:
          return Low;
        case Teuchos::VERB_MEDIUM:
          return Medium;
        case Teuchos::VERB_HIGH:
          return High;
        case Teuchos::VERB_EXTREME:
          return Extreme;
        default:  
          TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "Unknown enum value found.");  
        }  
  }

} // namespace MueLu
