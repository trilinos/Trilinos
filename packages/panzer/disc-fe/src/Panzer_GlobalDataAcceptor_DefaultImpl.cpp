// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
