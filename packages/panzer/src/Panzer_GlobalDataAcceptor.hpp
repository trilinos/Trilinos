#ifndef PANZER_GLOBAL_DATA_ACCEPTOR_HPP
#define PANZER_GLOBAL_DATA_ACCEPTOR_HPP

#include "Teuchos_RCP.hpp"

namespace panzer {

  class GlobalData;
  
  /** \brief Interface for accessing the GlobalData object.  */
  class GlobalDataAcceptor {

  public:
    
    virtual ~GlobalDataAcceptor() {}

    virtual void setGlobalData(const Teuchos::RCP<panzer::GlobalData>& gd) = 0;

    virtual Teuchos::RCP<panzer::GlobalData> getGlobalData() const = 0;

  };

}

#endif
