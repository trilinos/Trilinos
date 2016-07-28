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

#ifndef __Panzer_STK_SetupLOWSFactory_hpp__
#define __Panzer_STK_SetupLOWSFactory_hpp__

#include <ostream>
#include <string>
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_DefaultMpiComm.hpp"

#include "PanzerAdaptersSTK_config.hpp"

#include "Panzer_UniqueGlobalIndexer.hpp"
#include "Panzer_ConnManager.hpp"
#include "Panzer_BlockedDOFManager.hpp"
#ifdef PANZER_HAVE_FEI
#include "Panzer_DOFManagerFEI.hpp"
#endif

#include "Panzer_STKConnManager.hpp"

#include "Thyra_LinearOpWithSolveFactoryBase.hpp"

#ifdef PANZER_HAVE_TEKO 
#include "Teko_RequestHandler.hpp"
#endif

namespace panzer_stk {


Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> >
buildLOWSFactory(bool blockedAssembly,
                 const Teuchos::RCP<const panzer::UniqueGlobalIndexerBase> & globalIndexer,
                 const Teuchos::RCP<panzer::ConnManagerBase<int> > & conn_manager,
                 int spatialDim,
                 const Teuchos::RCP<const Teuchos::MpiComm<int> > & mpi_comm,
                 const Teuchos::RCP<Teuchos::ParameterList> & strat_params,
                 #ifdef PANZER_HAVE_TEKO 
                 const Teuchos::RCP<Teko::RequestHandler> & req_handler=Teuchos::null,
                 #endif 
                 bool writeCoordinates=false,
                 bool writeTopo=false
                 );

/** Build LOWS factory.
  */
template <typename GO> 
Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > 
buildLOWSFactory(bool blockedAssembly,
                 const Teuchos::RCP<const panzer::UniqueGlobalIndexerBase> & globalIndexer,
                 const Teuchos::RCP<panzer_stk::STKConnManager<GO> > & stkConn_manager,
                 int spatialDim,
                 const Teuchos::RCP<const Teuchos::MpiComm<int> > & mpi_comm,
                 const Teuchos::RCP<Teuchos::ParameterList> & strat_params,
                 #ifdef PANZER_HAVE_TEKO 
                 const Teuchos::RCP<Teko::RequestHandler> & req_handler,
                 #endif 
                 bool writeCoordinates=false,
                 bool writeTopo=false
                 );

template <typename GO>
void writeTopology(const panzer::BlockedDOFManager<int,GO> & blkDofs);

#ifdef PANZER_HAVE_FEI
template <typename GO> 
void writeTopology(const panzer::DOFManagerFEI<int,GO> & dofs,const std::string & block,std::ostream & os);
#endif

}

#endif
