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

#include "PanzerDiscFE_config.hpp"
#include "Panzer_UniqueGlobalIndexer.hpp"
#include "Panzer_IntrepidOrientation.hpp"

namespace panzer {

  Teuchos::RCP<std::vector<Intrepid2::Orientation> >
  buildIntrepidOrientation(const Teuchos::RCP<const UniqueGlobalIndexerBase> globalIndexer)
  {
    using Teuchos::rcp_dynamic_cast;
    using Teuchos::RCP;
    using Teuchos::rcp;
  
    auto orientation = rcp(new std::vector<Intrepid2::Orientation>);

    // int, int 
    {
      typedef int LO;
      typedef int GO;

      RCP<const UniqueGlobalIndexer<LO,GO> > ugi 
        = rcp_dynamic_cast<const UniqueGlobalIndexer<LO,GO> >(globalIndexer);

      if (ugi!=Teuchos::null) {
        const auto connMgrBase = ugi->getConnManagerBase();
        const auto connMgr = rcp_dynamic_cast<ConnManager<LO,GO> >(connMgrBase->noConnectivityClone());
        
        TEUCHOS_TEST_FOR_EXCEPTION(connMgr == Teuchos::null,std::logic_error,
                                   "panzer::buildIntrepidOrientation: Could not cast ConnManagerBase");
        
        buildIntrepidOrientation(*orientation, *connMgr);
        return orientation;
      }
    }

    // int, Ordinal64
    {
      typedef int LO;
      typedef Ordinal64 GO;
      
      RCP<const UniqueGlobalIndexer<LO,GO> > ugi 
        = rcp_dynamic_cast<const UniqueGlobalIndexer<LO,GO> >(globalIndexer);
      if (ugi!=Teuchos::null) {
        const auto connMgrBase = ugi->getConnManagerBase();
        const auto connMgr = rcp_dynamic_cast<ConnManager<LO,GO> >(connMgrBase->noConnectivityClone());
        
        TEUCHOS_TEST_FOR_EXCEPTION(connMgr == Teuchos::null,std::logic_error,
                                   "panzer::buildIntrepidOrientation: Could not cast ConnManagerBase");

        buildIntrepidOrientation(*orientation, *connMgr);
        return orientation;
      }
    }

    // int, pair<int,int>
    {
      typedef int LO;
      typedef std::pair<int,int> GO;

      RCP<const UniqueGlobalIndexer<LO,GO> > ugi 
        = rcp_dynamic_cast<const UniqueGlobalIndexer<LO,GO> >(globalIndexer);
      if(ugi!=Teuchos::null) {
        const auto connMgrBase = ugi->getConnManagerBase();
        const auto connMgr = rcp_dynamic_cast<ConnManager<LO,int> >(connMgrBase->noConnectivityClone());
        
        TEUCHOS_TEST_FOR_EXCEPTION(connMgr == Teuchos::null,std::logic_error,
                                   "panzer::buildIntrepidOrientation: Could not cast ConnManagerBase");

        buildIntrepidOrientation(*orientation, *connMgr);
        return orientation;
      }
    }

    // int, pair<int,Ordinal64>
    {
      typedef int LO;
      typedef std::pair<int,Ordinal64> GO;

      RCP<const UniqueGlobalIndexer<LO,GO> > ugi 
        = rcp_dynamic_cast<const UniqueGlobalIndexer<LO,GO> >(globalIndexer);
      if(ugi!=Teuchos::null) {
        const auto connMgrBase = ugi->getConnManagerBase();
        const auto connMgr = rcp_dynamic_cast<ConnManager<LO,Ordinal64> >(connMgrBase->noConnectivityClone());
        
        TEUCHOS_TEST_FOR_EXCEPTION(connMgr == Teuchos::null,std::logic_error,
                                   "panzer::buildIntrepidOrientation: Could not cast ConnManagerBase");

        buildIntrepidOrientation(*orientation, *connMgr);
        return orientation;
      }
    }
    
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                               "panzer::buildIntrepidOrientation: Could not cast UniqueGlobalIndexerBase");    
  }

} // end namespace panzer
