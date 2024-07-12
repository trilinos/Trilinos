// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_GLOBAL_DATA_ACCEPTOR_DEFAULT_IMPL_HPP
#define PANZER_GLOBAL_DATA_ACCEPTOR_DEFAULT_IMPL_HPP

#include "Teuchos_RCP.hpp"
#include "Panzer_GlobalDataAcceptor.hpp"

namespace panzer {

  /** \brief Default implementation for accessing the GlobalData object.  */
  class GlobalDataAcceptorDefaultImpl : public GlobalDataAcceptor {

  public:

    GlobalDataAcceptorDefaultImpl();

    GlobalDataAcceptorDefaultImpl(const Teuchos::RCP<panzer::GlobalData>& gd);

    ~GlobalDataAcceptorDefaultImpl();

    void setGlobalData(const Teuchos::RCP<panzer::GlobalData>& gd);

    Teuchos::RCP<panzer::GlobalData> getGlobalData() const;

  private:

    Teuchos::RCP<panzer::GlobalData> m_gd;

  };

}

#endif
