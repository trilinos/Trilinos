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
                   const Teuchos::RCP<panzer_stk::STKConnManager<int> > & stkConn_manager,
                   int spatialDim,
                   const Teuchos::RCP<const Teuchos::MpiComm<int> > & mpi_comm,
                   const Teuchos::RCP<Teuchos::ParameterList> & strat_params,
                   #ifdef PANZER_HAVE_TEKO
                   const Teuchos::RCP<Teko::RequestHandler> & reqHandler,
                   #endif
                   bool writeCoordinates,
                   bool writeTopo
                   );

#ifndef PANZER_ORDINAL64_IS_INT
  template
  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> >
  buildLOWSFactory<panzer::Ordinal64>(bool blockedAssembly,
                   const Teuchos::RCP<const panzer::UniqueGlobalIndexerBase> & globalIndexer,
                   const Teuchos::RCP<panzer_stk::STKConnManager<panzer::Ordinal64> > & stkConn_manager,
                   int spatialDim,
                   const Teuchos::RCP<const Teuchos::MpiComm<int> > & mpi_comm,
                   const Teuchos::RCP<Teuchos::ParameterList> & strat_params,
                   #ifdef PANZER_HAVE_TEKO
                   const Teuchos::RCP<Teko::RequestHandler> & reqHandler,
                   #endif
                   bool writeCoordinates,
                   bool writeTopo
                   );
#endif

  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> >
  buildLOWSFactory(bool blockedAssembly,
                   const Teuchos::RCP<const panzer::UniqueGlobalIndexerBase> & globalIndexer,
                   const Teuchos::RCP<panzer::ConnManagerBase<int> > & conn_manager,
                   int spatialDim,
                   const Teuchos::RCP<const Teuchos::MpiComm<int> > & mpi_comm,
                   const Teuchos::RCP<Teuchos::ParameterList> & strat_params,
                   #ifdef PANZER_HAVE_TEKO
                   const Teuchos::RCP<Teko::RequestHandler> & reqHandler,
                   #endif
                   bool writeCoordinates,
                   bool writeTopo
                   )
  {
    #ifdef PANZER_HAVE_TEKO
    Teuchos::RCP<Teko::RequestHandler> reqHandler_local = reqHandler;
    if(reqHandler_local==Teuchos::null)
      reqHandler_local = Teuchos::rcp(new Teko::RequestHandler);
    #endif

#ifndef PANZER_ORDINAL64_IS_INT
    Teuchos::RCP<panzer_stk::STKConnManager<panzer::Ordinal64> > long_conn = Teuchos::rcp_dynamic_cast<panzer_stk::STKConnManager<panzer::Ordinal64> >(conn_manager);
    if(long_conn!=Teuchos::null)
      return buildLOWSFactory(blockedAssembly,globalIndexer,long_conn,spatialDim,mpi_comm,strat_params,
                              #ifdef PANZER_HAVE_TEKO
                              reqHandler_local,
                              #endif
                              writeCoordinates,
                              writeTopo
                              );
#endif

    Teuchos::RCP<panzer_stk::STKConnManager<int> > int_conn = Teuchos::rcp_dynamic_cast<panzer_stk::STKConnManager<int> >(conn_manager);
    if(int_conn!=Teuchos::null)
      return buildLOWSFactory(blockedAssembly,globalIndexer,int_conn,spatialDim,mpi_comm,strat_params,
                              #ifdef PANZER_HAVE_TEKO
                              reqHandler_local,
                              #endif
                              writeCoordinates,
                              writeTopo
                              );

    // should never reach this
    TEUCHOS_ASSERT(false);
    return Teuchos::null;
  }
}
