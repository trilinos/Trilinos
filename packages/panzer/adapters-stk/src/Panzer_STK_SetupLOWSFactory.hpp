// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_STK_SetupLOWSFactory_hpp__
#define __Panzer_STK_SetupLOWSFactory_hpp__

#include <ostream>
#include <string>
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_DefaultMpiComm.hpp"

#include "PanzerAdaptersSTK_config.hpp"

#include "Panzer_GlobalIndexer.hpp"
#include "Panzer_ConnManager.hpp"
#include "Panzer_BlockedDOFManager.hpp"

#include "Panzer_STKConnManager.hpp"

#include "Thyra_LinearOpWithSolveFactoryBase.hpp"

#ifdef PANZER_HAVE_TEKO 
#include "Teko_RequestHandler.hpp"
#endif

namespace panzer_stk {

/** Build LOWS factory. */
Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> >
buildLOWSFactory(bool blockedAssembly,
                 const Teuchos::RCP<const panzer::GlobalIndexer> & globalIndexer,
                 const Teuchos::RCP<panzer::ConnManager> & conn_manager,
                 int spatialDim,
                 const Teuchos::RCP<const Teuchos::MpiComm<int> > & mpi_comm,
                 const Teuchos::RCP<Teuchos::ParameterList> & strat_params,
                 #ifdef PANZER_HAVE_TEKO 
                 const Teuchos::RCP<Teko::RequestHandler> & req_handler=Teuchos::null,
                 #endif 
                 bool writeCoordinates=false,
                 bool writeTopo=false,
                 const Teuchos::RCP<const panzer::GlobalIndexer> & auxGlobalIndexer=Teuchos::null,
                 bool useCoordinates=true
                 );

/** Build LOWS factory. */
Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > 
buildLOWSFactory(bool blockedAssembly,
                 const Teuchos::RCP<const panzer::GlobalIndexer> & globalIndexer,
                 const Teuchos::RCP<panzer_stk::STKConnManager> & stkConn_manager,
                 int spatialDim,
                 const Teuchos::RCP<const Teuchos::MpiComm<int> > & mpi_comm,
                 const Teuchos::RCP<Teuchos::ParameterList> & strat_params,
                 #ifdef PANZER_HAVE_TEKO 
                 const Teuchos::RCP<Teko::RequestHandler> & req_handler,
                 #endif 
                 bool writeCoordinates=false,
                 bool writeTopo=false,
                 const Teuchos::RCP<const panzer::GlobalIndexer> & auxGlobalIndexer=Teuchos::null,
                 bool useCoordinates=true
                 );

}

#endif
