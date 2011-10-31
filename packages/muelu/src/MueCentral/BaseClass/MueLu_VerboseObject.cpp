#include "MueLu_VerboseObject_def.hpp"

namespace MueLu {
 
  RCP<Teuchos::FancyOStream> VerboseObject::blackHole_ = Teuchos::getFancyOStream(rcp(new Teuchos::oblackholestream()));

}
