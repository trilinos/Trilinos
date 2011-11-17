#ifndef PANZER_GLOBAL_DATA_ACCEPTOR_DEFAULT_IMPL_HPP
#define PANZER_GLOBAL_DATA_ACCEPTOR_DEFAULT_IMPL_HPP

#include "Teuchos_RCP.hpp"
#include "Panzer_GlobalDataAcceptor.hpp"

namespace panzer {

  /** \brief Default implementation for accessing the GlobalData object.  */
  class GlobalDataAcceptorDefaultImpl : public GlobalDataAcceptor {

  public:

    GlobalDataAcceptorDefaultImpl();

    GlobalDataAcceptorDefaultImpl(const Teuchos::RCP<const panzer::GlobalData>& gd);

    ~GlobalDataAcceptorDefaultImpl();

    void setGlobalData(const Teuchos::RCP<const panzer::GlobalData>& gd);

    Teuchos::RCP<const panzer::GlobalData> getGlobalData() const;

  private:

    Teuchos::RCP<const panzer::GlobalData> m_gd;

  };

}

#endif
