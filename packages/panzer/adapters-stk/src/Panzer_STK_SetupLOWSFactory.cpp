// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#include "PanzerAdaptersSTK_config.hpp"

#include "Panzer_STK_SetupLOWSFactory.hpp"
#include "Panzer_STK_SetupLOWSFactory_impl.hpp"

namespace panzer_stk {

  template
  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> >
  buildLOWSFactory<int>(bool blockedAssembly,
                   const Teuchos::RCP<const panzer::UniqueGlobalIndexerBase> & globalIndexer,
                   const Teuchos::RCP<panzer_stk::STKConnManager> & stkConn_manager,
                   int spatialDim,
                   const Teuchos::RCP<const Teuchos::MpiComm<int> > & mpi_comm,
                   const Teuchos::RCP<Teuchos::ParameterList> & strat_params,
                   #ifdef PANZER_HAVE_TEKO
                   const Teuchos::RCP<Teko::RequestHandler> & reqHandler,
                   #endif
                   bool writeCoordinates,
                   bool writeTopo,
                   const Teuchos::RCP<const panzer::UniqueGlobalIndexerBase> & auxGlobalIndexer,
                   bool useCoordinates
                   );

#ifndef PANZER_ORDINAL64_IS_INT
  template
  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> >
  buildLOWSFactory<panzer::Ordinal64>(bool blockedAssembly,
                   const Teuchos::RCP<const panzer::UniqueGlobalIndexerBase> & globalIndexer,
                   const Teuchos::RCP<panzer_stk::STKConnManager> & stkConn_manager,
                   int spatialDim,
                   const Teuchos::RCP<const Teuchos::MpiComm<int> > & mpi_comm,
                   const Teuchos::RCP<Teuchos::ParameterList> & strat_params,
                   #ifdef PANZER_HAVE_TEKO
                   const Teuchos::RCP<Teko::RequestHandler> & reqHandler,
                   #endif
                   bool writeCoordinates,
                   bool writeTopo,
                   const Teuchos::RCP<const panzer::UniqueGlobalIndexerBase> & auxGlobalIndexer,
                   bool useCoordinates
                   );
#endif

  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> >
  buildLOWSFactory(bool blockedAssembly,
                   const Teuchos::RCP<const panzer::UniqueGlobalIndexerBase> & globalIndexer,
                   const Teuchos::RCP<panzer::ConnManager> & conn_manager,
                   int spatialDim,
                   const Teuchos::RCP<const Teuchos::MpiComm<int> > & mpi_comm,
                   const Teuchos::RCP<Teuchos::ParameterList> & strat_params,
                   #ifdef PANZER_HAVE_TEKO
                   const Teuchos::RCP<Teko::RequestHandler> & reqHandler,
                   #endif
                   bool writeCoordinates,
                   bool writeTopo,
                   const Teuchos::RCP<const panzer::UniqueGlobalIndexerBase> & auxGlobalIndexer,
                   bool useCoordinates
                   )
  {
    #ifdef PANZER_HAVE_TEKO
    Teuchos::RCP<Teko::RequestHandler> reqHandler_local = reqHandler;
    if(reqHandler_local==Teuchos::null)
      reqHandler_local = Teuchos::rcp(new Teko::RequestHandler);
    #endif

    auto stk_conn_manager = Teuchos::rcp_dynamic_cast<panzer_stk::STKConnManager>(conn_manager,true);

#ifndef PANZER_ORDINAL64_IS_INT
    auto globalOrdinalIsOrdinal64 = Teuchos::rcp_dynamic_cast<const panzer::UniqueGlobalIndexer<int,panzer::Ordinal64>>(globalIndexer);
    auto globalOrdinalIsOrdinal64_BLOCKED = Teuchos::rcp_dynamic_cast<const panzer::UniqueGlobalIndexer<int,std::pair<int,panzer::Ordinal64>>>(globalIndexer);
    if(nonnull(globalOrdinalIsOrdinal64) or nonnull(globalOrdinalIsOrdinal64_BLOCKED))
      return buildLOWSFactory<panzer::Ordinal64>(blockedAssembly,globalIndexer,stk_conn_manager,spatialDim,mpi_comm,strat_params,
#ifdef PANZER_HAVE_TEKO
                                                 reqHandler_local,
#endif
                                                 writeCoordinates,
                                                 writeTopo,
                                                 auxGlobalIndexer,
                                                 useCoordinates
                                                 );
#endif

    auto globalOrdinalIsInt = Teuchos::rcp_dynamic_cast<const panzer::UniqueGlobalIndexer<int,int>>(globalIndexer);
    auto globalOrdinalIsInt_BLOCKED = Teuchos::rcp_dynamic_cast<const panzer::UniqueGlobalIndexer<int,std::pair<int,int>>>(globalIndexer);
    if(nonnull(globalOrdinalIsInt) or nonnull(globalOrdinalIsInt_BLOCKED))
      return buildLOWSFactory<int>(blockedAssembly,globalIndexer,stk_conn_manager,spatialDim,mpi_comm,strat_params,
#ifdef PANZER_HAVE_TEKO
                                   reqHandler_local,
#endif
                                   writeCoordinates,
                                   writeTopo,
                                   auxGlobalIndexer,
                                   useCoordinates
                                   );

    // should never reach this
    TEUCHOS_ASSERT(false);
    return Teuchos::null;
  }
}
