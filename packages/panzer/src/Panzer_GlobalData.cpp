
#include "Panzer_GlobalData.hpp"
#include <ostream>

namespace panzer {
  
  Teuchos::RCP<panzer::GlobalData> 
  createGlobalData(bool build_default_os, int print_process)
  {
    Teuchos::RCP<panzer::GlobalData> gd = Teuchos::rcp(new panzer::GlobalData);
    gd->pl = Teuchos::rcp(new panzer::ParamLib);

    if (build_default_os) {
      gd->os = Teuchos::rcp(new Teuchos::FancyOStream(Teuchos::rcpFromRef(std::cout)));
      gd->os->setOutputToRootOnly(print_process);
    }

    return gd;
  }

}
