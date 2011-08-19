#include "MueLu_Level.hpp"

namespace MueLu {

  //Print function.  Not a friend b/c it uses only public interfaces for data access.
  std::ostream& operator<<(std::ostream& os, Level const &level) {
  
    os << "Printing Level object " << level.GetLevelID() << std::endl;
    //     if (level.template Get< Teuchos::RCP<Operator> >("A") != Teuchos::null) os << *level.template Get< Teuchos::RCP<Operator> >("A") << std::endl;
    //     if (level.template Get< Teuchos::RCP<Operator> >("R") != Teuchos::null) os << *level.template Get< Teuchos::RCP<Operator> >("R") << std::endl;
    //     if (level.template Get< Teuchos::RCP<Operator> >("P") != Teuchos::null) os << *level.template Get< Teuchos::RCP<Operator> >("P") << std::endl;
    return os;
  }

}
