
#include "Panzer_GlobalData.hpp"

namespace panzer {
  
  Teuchos::RCP<panzer::GlobalData> createGlobalData()
  {
    Teuchos::RCP<panzer::GlobalData> gd = Teuchos::rcp(new panzer::GlobalData);
    gd->pl = Teuchos::rcp(new panzer::ParamLib);
    gd->pv = Teuchos::rcp(new panzer::ParamVec);
    return gd;
  }

}
