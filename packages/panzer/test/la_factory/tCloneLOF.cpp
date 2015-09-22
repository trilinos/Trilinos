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

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <string>
#include <iostream>

#include "Phalanx_KokkosUtilities.hpp"

#include "Panzer_config.hpp"
#include "Panzer_EpetraLinearObjFactory.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_LinearObjFactory_Utilities.hpp"

#include "UnitTest_ConnManager.hpp"

#ifdef HAVE_MPI
   #include "Epetra_MpiComm.h"
#else
   #include "Epetra_SerialComm.h"
#endif

#include "Epetra_Map.h"

#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;
using Teuchos::rcpFromRef;

typedef Intrepid::FieldContainer<double> FieldContainer;

namespace panzer {

template <typename IntrepidType>
RCP<const panzer::FieldPattern> buildFieldPattern()
{
   // build a geometric pattern from a single basis
   RCP<Intrepid::Basis<double,FieldContainer> > basis = rcp(new IntrepidType);
   RCP<const panzer::FieldPattern> pattern = rcp(new panzer::IntrepidFieldPattern(basis));
   return pattern;
}

TEUCHOS_UNIT_TEST(tCloneLOF, epetra)
{
   PHX::KokkosDeviceSession session;

   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
   #else
      Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_SerialComm());
   #endif

   int myRank = eComm->MyPID();
   int numProc = eComm->NumProc();

   RCP<ConnManager<int,int> > connManager = rcp(new unit_test::ConnManager(myRank,numProc));

   RCP<const FieldPattern> patternC1
         = buildFieldPattern<Intrepid::Basis_HGRAD_HEX_C1_FEM<double,FieldContainer> >();

   RCP<panzer::DOFManager<int,int> > indexer = rcp(new panzer::DOFManager<int,int>());
   indexer->setConnManager(connManager,MPI_COMM_WORLD);
   indexer->addField("U",patternC1);
   indexer->addField("V",patternC1);
   indexer->buildGlobalUnknowns();

   RCP<panzer::DOFManager<int,int> > control_indexer = rcp(new panzer::DOFManager<int,int>());
   control_indexer->setConnManager(connManager,MPI_COMM_WORLD);
   control_indexer->addField("Z",patternC1);
   control_indexer->buildGlobalUnknowns();
 
   // setup factory
   RCP<EpetraLinearObjFactory<Traits,int> > ep_lof
         = Teuchos::rcp(new EpetraLinearObjFactory<Traits,int>(eComm.getConst(),indexer));
   
   // this is the member we are testing!
   RCP<const LinearObjFactory<Traits> > control_lof = cloneWithNewDomain(*ep_lof,control_indexer);
   RCP<const EpetraLinearObjFactory<Traits,int> > ep_control_lof 
       = rcp_dynamic_cast<const EpetraLinearObjFactory<Traits,int> >(control_lof);

   std::vector<int> control_owned, owned;
   control_indexer->getOwnedIndices(control_owned);

   TEST_ASSERT(ep_control_lof->getMap()->SameAs(*ep_lof->getMap()));
   TEST_EQUALITY(ep_control_lof->getColMap()->NumMyElements(),Teuchos::as<int>(control_owned.size()));
}

}
