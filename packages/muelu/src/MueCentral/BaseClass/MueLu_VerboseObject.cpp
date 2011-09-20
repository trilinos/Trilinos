#include "MueLu_VerboseObject.hpp"

namespace MueLu {
 
  RCP<Teuchos::FancyOStream> VerboseObject::blackHole_ = Teuchos::getFancyOStream(rcp(new Teuchos::oblackholestream()));

}
