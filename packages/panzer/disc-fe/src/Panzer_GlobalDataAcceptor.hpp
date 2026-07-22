// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_GLOBAL_DATA_ACCEPTOR_HPP
#define PANZER_GLOBAL_DATA_ACCEPTOR_HPP

#include "Teuchos_RCP.hpp"

namespace panzer {

  struct GlobalData;
  
  /** \brief Interface for accessing the GlobalData object.  */
  class GlobalDataAcceptor {

  public:
    
    virtual ~GlobalDataAcceptor() {}

    virtual void setGlobalData(const Teuchos::RCP<panzer::GlobalData>& gd) = 0;

    virtual Teuchos::RCP<panzer::GlobalData> getGlobalData() const = 0;

  };

}

#endif
