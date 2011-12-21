
#include "Panzer_GlobalDataAcceptor_DefaultImpl.hpp"
#include "Teuchos_Assert.hpp"

namespace panzer {

  GlobalDataAcceptorDefaultImpl::GlobalDataAcceptorDefaultImpl()
  { }  
  
  
  GlobalDataAcceptorDefaultImpl::GlobalDataAcceptorDefaultImpl(const Teuchos::RCP<panzer::GlobalData>& gd) : m_gd(gd)
  { }
  
  GlobalDataAcceptorDefaultImpl::~GlobalDataAcceptorDefaultImpl()
  { }
  
  void GlobalDataAcceptorDefaultImpl::setGlobalData(const Teuchos::RCP<panzer::GlobalData>& gd)
  {
    m_gd = gd;
  }
  
  Teuchos::RCP<panzer::GlobalData> GlobalDataAcceptorDefaultImpl::getGlobalData() const
  {
    TEUCHOS_ASSERT(nonnull(m_gd));
    return m_gd;
  }

}
