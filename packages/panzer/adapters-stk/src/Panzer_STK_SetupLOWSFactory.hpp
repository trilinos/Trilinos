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

/** \brief Convenience overload accepting a generic panzer::ConnManager,
  * which must actually be a panzer_stk::STKConnManager; dynamic-casts
  * conn_manager and delegates to the panzer_stk::STKConnManager
  * overload of buildLOWSFactory() below.
  */
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

/** \brief Builds a Thyra LOWS (linear-op-with-solve) factory from a
  * Stratimikos parameter list, configured for the given global indexer
  * and STK connection manager.
  *
  * \param[in] blockedAssembly true if the linear system is assembled in blocked form (one block per physics/DOF group).
  * \param[in] globalIndexer the DOF manager describing the linear system's unknowns.
  * \param[in] stkConn_manager the mesh connection manager.
  * \param[in] spatialDim the spatial dimension of the mesh.
  * \param[in] mpi_comm the MPI communicator to build the solver on.
  * \param[in] strat_params the Stratimikos parameter list describing the requested solver/preconditioner.
  * \param[in] req_handler (only if PANZER_HAVE_TEKO) an optional Teko request handler; a default one is created if null.
  * \param[in] writeCoordinates if true, writes mesh coordinates to file (e.g. for preconditioner setup diagnostics).
  * \param[in] writeTopo if true, writes the mesh topology to file.
  * \param[in] auxGlobalIndexer an optional auxiliary DOF manager, e.g. for preconditioners needing extra unknowns.
  * \param[in] useCoordinates if true, supplies mesh coordinates to the solver/preconditioner (e.g. for algebraic multigrid).
  * \return the constructed LOWS factory.
  */
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
